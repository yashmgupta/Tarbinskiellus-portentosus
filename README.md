# *Tarbinskiellus portentosus* mitochondrial-genome reproducibility repository

This repository contains the data and code used to validate read support, nucleotide-depth coverage, primer/amplicon positions, and independent TP13 Sanger support for the completed mitochondrial genome of *Tarbinskiellus portentosus*.

**Final mitochondrial genome:** GenBank accession **PZ447390.1**  
**Genome length:** **15,495 bp**  
**Main reproducibility notebook:** `notebooks/TP_mitogenome_custom_mapping_reproducibility.ipynb`

## What this repository reproduces

The supplied notebook and Python code reproduce:

- circular mapping of the paired-end short reads to the final mitochondrial sequence;
- per-base Illumina depth for the BAN108 and BAN72 datasets and for the combined dataset;
- mapping statistics and the accepted-read mapping table;
- the diagnostic lowest-depth 325-bp window;
- remapping of the 15 overlapping PCR amplicons to the final circular genome;
- quality trimming and local alignment of the TP13 Sanger chromatogram;
- the main coverage and control-region figures used to document sequence support.

This repository is a **validation and reproducibility package**, not a byte-for-byte reconstruction of the historical primary assembly. CAP3 v10.2011 remains the primary assembly program recorded for the deposited sequence, but the original intermediate CAP3 files and any historical draft containing an N-block are not available in this package. See `REPRODUCIBILITY_SCOPE.md` for the exact scope and limitations.

## Repository structure

```text
.
├── README.md
├── REPRODUCIBILITY_SCOPE.md
├── requirements.txt
├── SHA256SUMS.txt
├── raw_data/
│   ├── sequence.gb
│   ├── TP13-TP13FR.ab1
│   ├── TP13-TP13FR.seq
│   ├── ...BAN108...R1.fq
│   ├── ...BAN108...R2.fq
│   ├── ...BAN72...R1.fasta
│   └── ...BAN72...R2.fasta
├── derived_data/
│   ├── BAN108_coverage.tsv
│   ├── BAN72_coverage.tsv
│   ├── combined_coverage.tsv
│   ├── mapping_stats.json
│   ├── read_mappings.tsv
│   ├── coverage_summary.txt
│   ├── lowest_coverage_325bp_window.tsv
│   ├── primer_amplicons_remapped_to_final_genome.csv
│   └── TP13_sanger_qtrim_alignment_summary.tsv
├── scripts/
│   └── custom_seed_mismatch_mapper.py
├── notebooks/
│   ├── TP_mitogenome_custom_mapping_reproducibility.ipynb
│   └── TP_mitogenome_custom_mapping_reproducibility_EXECUTED.ipynb
├── figures/
│   └── publication and reproducibility figures in PNG/PDF format
└── manuscript/
    ├── METHODS_custom_mapping.md
    └── FIGURE_LEGENDS.md
```

## Software requirements

The workflow was prepared with Python and the packages listed in `requirements.txt`:

- Biopython 1.86
- NumPy 2.3.5
- pandas 2.2.3
- Matplotlib 3.10.8
- JupyterLab 4+
- nbformat 5.10.4

No external BWA, Bowtie2, SAMtools, or BCFtools executable is required for the custom coverage-validation workflow.

## Installation

Create a clean Python environment from the repository root:

```bash
python -m venv tp_mito_env
```

Activate it.

Linux/macOS:

```bash
source tp_mito_env/bin/activate
```

Windows PowerShell:

```powershell
tp_mito_env\Scripts\Activate.ps1
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Then start Jupyter:

```bash
jupyter lab
```

Open:

```text
notebooks/TP_mitogenome_custom_mapping_reproducibility.ipynb
```

and choose **Run All**. The notebook detects whether Jupyter was launched from the repository root or from the `notebooks/` directory.

## Expected validation results

Running the notebook should reproduce the following values:

| Metric | Expected value |
|---|---:|
| Final genome length | 15,495 bp |
| Total paired-end reads | 62,354 |
| Accepted mappings | 54,501 |
| Mapping rate | 87.41% |
| Combined mean Illumina depth | 459.02× |
| Combined median Illumina depth | 292× |
| Bases with ≥1× Illumina depth | 15,418 / 15,495 (99.50%) |
| Lowest-depth 325-bp window | positions 14,937–15,261 |
| Mean depth in that window | 1.41× |
| Zero-depth positions in that window | 77 / 325 |
| Quality-trimmed TP13 Sanger alignment | position 14,968 → circular origin → 123 |
| Quality-trimmed TP13 length | 651 bp |
| TP13 Sanger mismatches to final sequence | 1 |

## Mapping method

The custom mapper in `scripts/custom_seed_mismatch_mapper.py` was used as an **assembly-validation and coverage-diagnostic** method. It uses:

- an exact 19-nt seed index built from the circularized reference;
- six seed positions per read;
- up to 80 candidate alignment starts;
- both forward and reverse-complement orientations;
- full-read Hamming mismatch scoring without gapped alignment;
- acceptance of the best placement when the read contains no more than 12 mismatches.

The reference is treated as circular, so placements crossing the artificial linearization origin wrap back to the beginning of the genome.

## Important interpretation of coverage

The short-read library was generated from **pooled overlapping PCR amplicons**, not from randomly fragmented genomic DNA. Therefore, large differences in read depth primarily reflect differences in PCR amplification and amplicon representation in the pool and should not be interpreted as ordinary whole-genome shotgun coverage variation.

The TP13 Sanger sequence is an independent evidence track. It is **not** converted into artificial Illumina read depth. A position with 0× Illumina depth can therefore still have independent sequence support from TP13 Sanger sequencing.

## Figures

The `figures/` directory contains both notebook-generated analytical figures and publication-ready redesigned exports. The two most useful figures for the manuscript are:

- `TP_coverage_with_amplicon_design_qtrim.png` / `.pdf` — whole-genome Illumina depth with overlapping PCR amplicons and TP13 Sanger support;
- `TP_control_region_Sanger_support_redesigned.png` / `.pdf` — focused visualization of the control region, circular origin, low Illumina support, and independent TP13 Sanger support.

Dataset-specific and additional diagnostic coverage figures are also included.

## Manuscript-support files

`manuscript/METHODS_custom_mapping.md` contains the updated reproducibility-oriented Materials and Methods text, while `manuscript/FIGURE_LEGENDS.md` contains suggested legends for the coverage figures.

## Integrity checking

`SHA256SUMS.txt` contains SHA-256 checksums for the files in this repository. On Linux/macOS, they can be checked with:

```bash
sha256sum -c SHA256SUMS.txt
```

## Citation and data availability

When using this repository, cite the associated *Tarbinskiellus portentosus* mitochondrial-genome manuscript and the deposited GenBank sequence **PZ447390.1**. After the manuscript is published, the final article citation and DOI can be added to this section.

## Reproducibility note

For the exact analytical scope, including what is **not** reconstructed by this package, read `REPRODUCIBILITY_SCOPE.md` before interpreting the coverage statistics or the 325-bp lowest-depth interval.
