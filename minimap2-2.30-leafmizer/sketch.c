#include "bseq.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct
{ // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0)
		return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

/**
 * Find symmetric (w,k)-leafmizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a leafmizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      leafmizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th leafmizer,
 *               and strand indicates whether the leafmizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len / w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = {UINT64_MAX, UINT64_MAX};
		if (c < 4)
		{ // not an ambiguous base
			int z;
			if (is_hpc)
			{
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c)
				{
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k)
					kmer_span -= tq_shift(&tq);
			}
			else
				kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;			 // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1])
				continue;				   // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1] ? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256)
			{
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | z;
			}
		}
		else
			l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX)
		{ // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y)
					kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y)
					kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x)
		{ // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX)
				kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		}
		else if (buf_pos == min_pos)
		{ // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX)
				kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x)
					min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x)
					min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX)
			{									  // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y)
						kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y)
						kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w)
			buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);
}

// leafmizer_edit_20251228
void mm_sketch_query(const mm_idx_t *mi, void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};
	uint32_t buf_z[256], min_z = UINT32_MAX;
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(buf_z, 0xff, w * 4);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len / w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = {UINT64_MAX, UINT64_MAX};
		uint32_t info_z = UINT32_MAX;
		if (c < 4)
		{ // not an ambiguous base
			int z;
			if (is_hpc)
			{
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c)
				{
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k)
					kmer_span -= tq_shift(&tq);
			}
			else
				kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;			 // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1])
				continue;				   // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1] ? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256)
			{
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | z;
				info_z = mi->kmer_counts[kmer[z] & mi->kmer_counts_size]; // get kmer count from the index
			}
		}
		else
			l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		buf_z[buf_pos] = info_z;
		if (l == w + k - 1 && min_z != UINT32_MAX && min.x != UINT64_MAX)
		{ // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min_z == buf_z[j] && buf[j].y != min.y)
					kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min_z == buf_z[j] && buf[j].y != min.y)
					kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info_z <= min_z)
		{ // a new minimum; then write the old min
			if (l >= w + k && min_z != UINT32_MAX && min.x != UINT64_MAX)
				kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
			min_z = info_z;
		}
		else if (buf_pos == min_pos)
		{ // old min has moved outside the window
			if (l >= w + k - 1 && min_z != UINT32_MAX && min.x != UINT64_MAX)
				kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min_z = UINT32_MAX, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min_z >= buf_z[j])
					min_z = buf_z[j], min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min_z >= buf_z[j])
					min_z = buf_z[j], min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min_z != UINT32_MAX)
			{									  // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min_z == buf_z[j] && min.y != buf[j].y)
						kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min_z == buf_z[j] && min.y != buf[j].y)
						kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w)
			buf_pos = 0;
	}
	if (min_z != UINT32_MAX && min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);
}

// leafmizer_edit_20251228
void mm_sketch_reference(const mm_idx_t *mi, void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	/*
	When comparing three consecutive leafmizer1, leafmizer2, and leafmizer3,
	first, keep the first one and temporarily store the second one in a variable.
	Based on the position of the third one, decide whether to save the second one.
	If the difference in position between the two leafmizers is less than or equal to the window size,
	delete the middle leafmizer, update the relative position,
	and change the third one to the second one.
	Conversely, keep the middle leafmizer, update the relative positions,
	change the second one to the first one,
	and change the third one to the second one.
	*/

	unsigned long long int leafmizer1_suoyin_concisemizer = 0; // The location of leafmizer1
	unsigned long long int leafmizer2_concisemizer = 0;		   // leafmizer2
	unsigned long long int leafmizer2_suoyin_concisemizer = 0; // The location of leafmizer2
	mm128_t min2 = {UINT64_MAX, UINT64_MAX};

	uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = {0, 0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = {UINT64_MAX, UINT64_MAX};
	uint32_t buf_z[256], min_z = UINT32_MAX;
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(buf_z, 0xff, w * 4);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len / w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i)
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = {UINT64_MAX, UINT64_MAX};
		uint32_t info_z = UINT32_MAX;
		if (c < 4)
		{ // not an ambiguous base
			int z;
			if (is_hpc)
			{
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c)
				{
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k)
					kmer_span -= tq_shift(&tq);
			}
			else
				kmer_span = l + 1 < k ? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;			 // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1])
				continue;				   // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1] ? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256)
			{
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | z;
				info_z = mi->kmer_counts[kmer[z] & mi->kmer_counts_size]; // get kmer count from the index
			}
		}
		else
			l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		buf_z[buf_pos] = info_z;
		if (l == w + k - 1 && min.x != UINT64_MAX && min_z != UINT32_MAX)
		{ // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)

				if (min_z == buf_z[j] && buf[j].y != min.y)
				{

					leafmizer1_suoyin_concisemizer = buf[j].y;
					kv_push(mm128_t, km, *p, buf[j]);
				}
			for (j = 0; j < buf_pos; ++j)

				if (min_z == buf_z[j] && buf[j].y != min.y)
				{

					leafmizer1_suoyin_concisemizer = buf[j].y;
					kv_push(mm128_t, km, *p, buf[j]);
				}
		}
		if (info_z <= min_z)
		{ // a new minimum; then write the old min

			if (l >= w + k && min.x != UINT64_MAX && min_z != UINT32_MAX)
			{
				// If leafmizer3, normal comparison; If it is leafmizer2, the distance between it and leafmizer1 must be less than or equal to w.
				// So there is no need to determine which one it is, so the old version of the code is commented out.

				if (((min.y >> 1) - (leafmizer1_suoyin_concisemizer >> 1)) > w)
				{
					leafmizer1_suoyin_concisemizer = leafmizer2_suoyin_concisemizer;
					min2.x = leafmizer2_concisemizer;
					min2.y = leafmizer2_suoyin_concisemizer;
					kv_push(mm128_t, km, *p, min2);
					leafmizer2_concisemizer = min.x;
					leafmizer2_suoyin_concisemizer = min.y;
				}
				else
				{
					leafmizer2_concisemizer = min.x;
					leafmizer2_suoyin_concisemizer = min.y;
				}
			}
			min = info, min_pos = buf_pos;
			min_z = info_z;
		}
		else if (buf_pos == min_pos)
		{ // old min has moved outside the window

			if (l >= w + k - 1 && min.x != UINT64_MAX && min_z != UINT32_MAX)
			{
				if (((min.y >> 1) - (leafmizer1_suoyin_concisemizer >> 1)) > w)
				{
					leafmizer1_suoyin_concisemizer = leafmizer2_suoyin_concisemizer;
					min2.x = leafmizer2_concisemizer;
					min2.y = leafmizer2_suoyin_concisemizer;
					kv_push(mm128_t, km, *p, min2);
					leafmizer2_concisemizer = min.x;
					leafmizer2_suoyin_concisemizer = min.y;
				}
				else
				{
					leafmizer2_concisemizer = min.x;
					leafmizer2_suoyin_concisemizer = min.y;
				}
			}
			for (j = buf_pos + 1, min_z = UINT32_MAX, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min_z >= buf_z[j])
					min_z = buf_z[j], min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min_z >= buf_z[j])
					min_z = buf_z[j], min = buf[j], min_pos = j;

			if (l >= w + k - 1 && min.x != UINT64_MAX && min_z != UINT32_MAX)
			{									  // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min_z == buf_z[j] && min.y != buf[j].y)
					{
						if (((buf[j].y >> 1) - (leafmizer1_suoyin_concisemizer >> 1)) > w)
						{
							leafmizer1_suoyin_concisemizer = leafmizer2_suoyin_concisemizer;
							min2.x = leafmizer2_concisemizer;
							min2.y = leafmizer2_suoyin_concisemizer;
							kv_push(mm128_t, km, *p, min2);
							leafmizer2_concisemizer = buf[j].x;
							leafmizer2_suoyin_concisemizer = buf[j].y;
						}
						else
						{
							leafmizer2_concisemizer = buf[j].x;
							leafmizer2_suoyin_concisemizer = buf[j].y;
						}
					}
				for (j = 0; j <= buf_pos; ++j)
					if (min_z == buf_z[j] && min.y != buf[j].y)
					{
						if (((buf[j].y >> 1) - (leafmizer1_suoyin_concisemizer >> 1)) > w)
						{
							leafmizer1_suoyin_concisemizer = leafmizer2_suoyin_concisemizer;
							min2.x = leafmizer2_concisemizer;
							min2.y = leafmizer2_suoyin_concisemizer;
							kv_push(mm128_t, km, *p, min2);
							leafmizer2_concisemizer = buf[j].x;
							leafmizer2_suoyin_concisemizer = buf[j].y;
						}
						else
						{
							leafmizer2_concisemizer = buf[j].x;
							leafmizer2_suoyin_concisemizer = buf[j].y;
						}
					}
			}
		}
		if (++buf_pos == w)
			buf_pos = 0;
	}
	if (min_z != UINT32_MAX && min.x != UINT64_MAX)
	{

		kv_push(mm128_t, km, *p, min);
		if ((((min.y >> 1) - (leafmizer1_suoyin_concisemizer >> 1)) > w) && (leafmizer2_suoyin_concisemizer > 0))
		{
			min2.x = leafmizer2_concisemizer;
			min2.y = leafmizer2_suoyin_concisemizer;
			kv_push(mm128_t, km, *p, min2);
		}
	}
}

// leafmizer_edit_20251228: Independently scan the reference sequence file and generate a k-mer count table.
void mm_count_kmers_from_file(const char *fn, mm_idx_t *mi, int k)
{
	if (!fn || !mi || !mi->kmer_counts)
		return;
	mm_bseq_file_t *fp = mm_bseq_open(fn);
	if (!fp)
		return;
	mm_bseq1_t *seqs = NULL;
	int n_seq = 0;
	while ((seqs = mm_bseq_read(fp, 256, 0, &n_seq)) != NULL)
	{
		int i = 0;
		for (i = 0; i < n_seq; ++i)
		{
			const char *s = seqs[i].seq;
			int len = seqs[i].l_seq;
			uint64_t kmer_f = 0, kmer_r = 0;
			uint64_t kmer_mask = (1ULL << (2 * k)) - 1;
			uint64_t shift1 = 2 * (k - 1);
			int l = 0;
			int j = 0;
			for (j = 0; j < len; ++j)
			{
				int c = seq_nt4_table[(uint8_t)s[j]];
				if (c < 4)
				{
					kmer_f = ((kmer_f << 2) | c) & kmer_mask;
					kmer_r = (kmer_r >> 2) | ((uint64_t)(3ULL ^ c) << shift1);
					++l;
					if (l >= k)
					{
						uint64_t kmer = kmer_f < kmer_r ? kmer_f : kmer_r;
						uint64_t idx = kmer & mi->kmer_counts_size;
						mi->kmer_counts[idx]++;
					}
				}
				else
				{
					l = 0;
					kmer_f = kmer_r = 0;
				}
			}
		}
		i = 0;
		for (i = 0; i < n_seq; ++i)
		{
			free(seqs[i].seq);
			free(seqs[i].name);
		}
		free(seqs);
	}
	mm_bseq_close(fp);
	// Check mi->kmer_counts; if any value is zero, set it to the maximum value of a 32-bit unsigned integer.
	uint64_t i = 0;
	for (i = 0; i <= mi->kmer_counts_size; ++i)
	{
		if (mi->kmer_counts[i] == 0)
			mi->kmer_counts[i] = UINT32_MAX;
	}
}