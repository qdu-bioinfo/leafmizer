# Leafmizer

**Global abundance-aware seed sampling for sequence alignment with minimap2**

Leafmizer (*least frequently occurring k-mer minimizer*) is a seed sampling algorithm that incorporates reference-wide k-mer abundance into seed selection. Instead of selecting a representative k-mer solely by its hash value, Leafmizer prioritizes candidate k-mers with lower estimated abundance within each sliding window. The goal is to reduce repetitive seed occurrences and unnecessary candidate anchors while retaining the familiar minimap2 mapping workflow.

This repository provides the **Leafmizer-modified minimap2 v2.30 source code** and three supplementary workbooks associated with the study *Leafmizer: A Global Abundance-Aware Seed Sampling Algorithm for Sequence Alignment*.

> **Terminology:** `minimap2_mini` denotes the original minimap2/minimizer implementation in the study; `minimap2_leaf` denotes the Leafmizer-modified implementation. The compiled executable in this repository is named `minimap2`, not `minimap2_leaf`.

## Features

- **Abundance-aware seed ordering:** candidate k-mers are prioritized using estimated occurrence information from the reference genome.
- **Efficient approximation:** a reference-derived, bounded-size abundance table avoids maintaining exact frequencies for every possible k-mer.
- **Spatial redundancy removal:** a subsequent selection step reduces excessive local seed density.
- **minimap2 integration:** the modified seed sampling procedure is incorporated into the minimap2 indexing and mapping framework; the downstream chaining and alignment procedures are retained.
- **Long-read evaluation:** the accompanying study evaluates seed redundancy, computational cost, and mapping outcomes using reference genomes, simulated reads, and real long-read data.

## Repository contents

```text
leafmizer/
├── README.md
├── minimap2-2.30-leafmizer/                       # Leafmizer-modified minimap2 source
├── Supplementary Table S1 Evaluation_metrics.xlsx
├── Supplementary Table S2 Simulation_data_comparison_results.xlsx
└── Supplementary Table S3 Homo_sapiens_real_data_PRJEB77558.xlsx
```

## Requirements

To compile from source, you need:

- A C compiler (such as GCC or Clang)
- GNU Make
- zlib and its development headers
- A supported CPU architecture (the bundled minimap2 source uses SIMD-accelerated alignment code)

On Debian/Ubuntu, the main build prerequisites can be installed with:

```bash
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev
```

## Installation

Clone this repository and compile the **modified** source tree:

```bash
git clone https://github.com/qdu-bioinfo/leafmizer.git
cd leafmizer/minimap2-2.30-leafmizer
make -j4
./minimap2 --version
```

The resulting executable is `./minimap2` inside `minimap2-2.30-leafmizer/`.

**Important:** installing the standard `minimap2` package from an external package manager does **not** install the Leafmizer modification. Build and use the executable from this repository when evaluating `minimap2_leaf`.

## Quick start

The following commands assume you are in `leafmizer/minimap2-2.30-leafmizer/` and have a locally accessible reference FASTA (`ref.fa`) and Nanopore reads (`reads.fastq.gz`). Replace these example filenames with your own files.

### Option 1: Align directly from a reference FASTA

```bash
./minimap2 -t 8 -ax map-ont ref.fa reads.fastq.gz > alignment.sam
```

### Option 2: Build an index and reuse it

```bash
# Construct an index with the Leafmizer-modified executable
./minimap2 -t 8 -x map-ont -d ref.leaf.mmi ref.fa

# Align Nanopore reads with that index
./minimap2 -t 8 -ax map-ont ref.leaf.mmi reads.fastq.gz > alignment.sam
```

To obtain PAF output instead of SAM, omit `-a`:

```bash
./minimap2 -t 8 -x map-ont ref.leaf.mmi reads.fastq.gz > alignment.paf
```

Other minimap2 presets used in the study include `map-pb`, `map-hifi`, and `lr:hq`. Choose the preset appropriate for your sequencing technology, and use the same preset/seed settings when building and using an index.

**Index compatibility:** Leafmizer's reference index stores additional abundance-count information used during seed selection. **Build your `.mmi` index with the Leafmizer executable and use the same implementation to load it.** Do not assume that an index created by unmodified minimap2 is interchangeable with a Leafmizer index. The current reference-counting implementation also reads the reference sequence file when constructing the index; using a regular, locally readable FASTA path is recommended.

For standard minimap2 command-line flags, output formats, and mapping presets, see the [minimap2 documentation](https://github.com/lh3/minimap2) and the documentation included in `minimap2-2.30-leafmizer/`.

## Algorithm overview

The central distinction between the two implementations is the **seed ordering criterion**:

1. **Reference abundance estimation.** Leafmizer constructs an approximate k-mer abundance lookup from the reference sequence before sampling reference seeds.
2. **Abundance-aware selection.** Within each window, candidate k-mers are ranked by their estimated abundance, favoring less frequent candidates; the k-mer hash remains part of the seed representation.
3. **Spatial redundancy removal.** Locally redundant selected seeds are further reduced by a position-based filtering step.
4. **Mapping.** Query seeds are selected using the reference-derived abundance information, and the resulting anchors enter the existing minimap2 chaining and alignment stages.

By prioritizing less repetitive candidates, Leafmizer aims to lower average seed occurrence. Lower seed redundancy is an index-level property; its downstream effects on mapping sensitivity, accuracy, runtime, and memory must be assessed separately.

## Supplementary data

Three Excel workbooks are available in the repository root:

| File | Description |
| --- | --- |
| [Supplementary Table S1 — Evaluation metrics](Supplementary%20Table%20S1%20Evaluation_metrics.xlsx) | Reference-genome seed sampling and index-evaluation results, associated with the study's 1,435-genome benchmark. |
| [Supplementary Table S2 — Simulation comparison results](Supplementary%20Table%20S2%20Simulation_data_comparison_results.xlsx) | Alignment evaluation results on simulated long-read data, including the simulated-data comparisons discussed in the manuscript. |
| [Supplementary Table S3 — Human real-data results](Supplementary%20Table%20S3%20Homo_sapiens_real_data_PRJEB77558.xlsx) | Evaluation results for real human long-read samples from study accession **PRJEB77558**, including computational and mapping outcome metrics. |

The reference-genome benchmark covers 1,435 RefSeq genomes from invertebrates, plants, mammalian vertebrates, and other vertebrates. The real-data analysis in the manuscript uses Oxford Nanopore and PacBio long-read samples; the Nanopore comparison includes 161 samples.

These tables provide results accompanying the study. They are not substitutes for the underlying reference genomes or sequencing reads, which should be obtained from their original sources for independent reruns.

## Summary of the reported results

In the experiments described in the accompanying study:

- **Reference-genome seed redundancy (`map-ont`):** the per-genome median reduction in average seed occurrence was approximately **20.8%** across 1,435 reference genomes; all evaluated genomes showed a reduction relative to `minimap2_mini`.
- **Real Nanopore CPU time:** the median per-sample reduction was approximately **31.8%** across 161 samples.
- **Real Nanopore peak memory usage:** the median per-sample reduction was approximately **8.0%**; 159 of 161 samples had lower peak memory usage.
- **Real Nanopore unmapped read counts:** the median per-sample reduction was approximately **2.9%**; 156 samples had fewer unmapped reads, one was unchanged, and four had more.

For a metric where lower values are preferable, the per-sample (or per-genome) relative reduction is calculated as:

```text
Relative reduction (%) = 100 × (value_minimap2_mini − value_minimap2_leaf) / value_minimap2_mini
```

A positive value denotes a reduction with `minimap2_leaf`. The **median of paired relative reductions** is not the same as the **percentage reduction between group means**; please distinguish these statistics when reporting results. Also, a decrease in unmapped-read count does not by itself establish improved alignment correctness.

## Reproducibility notes

When comparing `minimap2_mini` with `minimap2_leaf`, we recommend recording:

- The reference FASTA and its assembly/accession version
- The source revision or Git commit of each implementation
- The selected preset and relevant `-k`, `-w`, and other mapping/indexing parameters
- CPU/thread configuration and hardware specifications
- Whether mapping used a freshly built index or a saved implementation-specific `.mmi` file
- Definitions and units for CPU time, memory usage, and mapping statistics

For fair comparisons, keep the reference, reads, presets, and hardware conditions consistent. The modified index-building stage involves additional reference abundance estimation, so index construction overhead should be considered separately from read mapping performance.

## Citation

If you use Leafmizer, please cite the accompanying manuscript when a publication record becomes available:

> Gao, H., Sun, Y., Xing, J., Yu, J., and Su, X. **Leafmizer: A Global Abundance-Aware Seed Sampling Algorithm for Sequence Alignment.** Manuscript.

Also cite the original minimap2 publication:

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

## License and acknowledgments

This implementation is derived from [minimap2](https://github.com/lh3/minimap2) by Heng Li and contributors. The bundled minimap2 source contains an MIT license in [`minimap2-2.30-leafmizer/LICENSE.txt`](minimap2-2.30-leafmizer/LICENSE.txt). Retain the applicable copyright and license notices when redistributing the source.

## Issues and questions

Please report reproducible problems or questions via this repository's [GitHub Issues](https://github.com/qdu-bioinfo/leafmizer/issues). When reporting a problem, include the exact command, input type, program revision, platform, and relevant error messages.
