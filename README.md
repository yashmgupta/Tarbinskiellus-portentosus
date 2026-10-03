# *Tarbinskiellus portentosus* mitochondrial-genome reproducibility repository

This repository contains the data, code, and notebooks used to validate read support, nucleotide-depth coverage, primer/amplicon positions, and independent TP13 Sanger support for the completed mitochondrial genome of *Tarbinskiellus portentosus*.

**Final mitochondrial genome:** GenBank accession **PZ447390.1**  
**Genome length:** **15,495 bp**  
**Main reproducibility notebook:** `notebooks/TP_mitogenome_custom_mapping_reproducibility.ipynb`

## Purpose of this repository

The files in this repository reproduce the validation analyses associated with the completed mitochondrial genome, including:

- circular mapping of paired-end short reads to the final mitochondrial sequence;
- per-base Illumina depth for the BAN108 and BAN72 datasets and for the combined dataset;
- mapping statistics and accepted-read mapping records;
- identification of the diagnostic lowest-depth 325-bp interval;
- remapping of the 15 overlapping PCR amplicons to the final circular genome;
- quality trimming and local alignment of the TP13 Sanger chromatogram; and
- numerical results used to describe sequence support across the mitochondrial genome and control-region/origin interval.

This repository is a **validation and reproducibility package**, not a byte-for-byte reconstruction of the historical primary assembly. CAP3 v10.2011 remains the primary assembly program recorded for the deposited sequence, but the original intermediate CAP3 files and any historical draft containing an N-block are not included here. See `REPRODUCIBILITY_SCOPE.md` for the exact analytical scope and limitations.

## Repository structure

```text
.
├── README.md
├── REPRODUCIBILITY_SCOPE.md
├── requirements.txt
├── .gitignore
│
├── raw_data/
│   ├── sequence.gb
│   ├── TP13-TP13FR.ab1
│   ├── TP13-TP13FR.seq
│   ├── TSN20260330-010-00157_TP-mito_20260331-BAN108_D06_D06_i29.R1.fq
│   ├── TSN20260330-010-00157_TP-mito_20260331-BAN108_D06_D06_i29.R2.fq
│   ├── TSN20260330-010-00157_TP-mito_20260403-BAN72_A05_A05_i27.R1.fasta
│   └── TSN20260330-010-00157_TP-mito_20260403-BAN72_A05_A05_i27.R2.fasta
│
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
│
├── scripts/
│   └── custom_seed_mismatch_mapper.py
│
└── notebooks/
    ├── TP_mitogenome_custom_mapping_reproducibility.ipynb
    └── TP_mitogenome_custom_mapping_reproducibility_EXECUTED.ipynb
```

## Raw data

The `raw_data/` directory contains the sequence files required for the validation workflow.

- `sequence.gb` — final annotated mitochondrial genome used as the mapping reference.
- `TP13-TP13FR.ab1` — TP13 Sanger chromatogram.
- `TP13-TP13FR.seq` — associated TP13 sequence file.
- BAN108 R1/R2 FASTQ files — one paired-end short-read dataset.
- BAN72 R1/R2 FASTA files — the second paired-end read dataset supplied for the analysis.

The completed mitochondrial sequence used in this workflow is **15,495 bp** in length and is deposited in GenBank as **PZ447390.1**.

## Derived data

The `derived_data/` directory contains outputs from the validation analysis.

- `BAN108_coverage.tsv` — per-base coverage for BAN108.
- `BAN72_coverage.tsv` — per-base coverage for BAN72.
- `combined_coverage.tsv` — combined per-base Illumina depth.
- `mapping_stats.json` — overall read-mapping statistics.
- `read_mappings.tsv` — accepted read placements, including mapping coordinates, strand, read length, mismatch count, and best-hit information.
- `coverage_summary.txt` — summary coverage statistics.
- `lowest_coverage_325bp_window.tsv` — diagnostic lowest-depth 325-bp interval.
- `primer_amplicons_remapped_to_final_genome.csv` — primer/amplicon coordinates after remapping to the final circular genome.
- `TP13_sanger_qtrim_alignment_summary.tsv` — summary of the quality-trimmed TP13 Sanger alignment to the final sequence.

## Software requirements

The workflow uses Python and the packages listed in `requirements.txt`:

- Biopython 1.86
- NumPy 2.3.5
- pandas 2.2.3
- Matplotlib 3.10.8
- JupyterLab 4+
- nbformat 5.10.4

No external BWA, Bowtie2, SAMtools, or BCFtools executable is required for this custom coverage-validation workflow.

## Installation

From the repository root, create a Python virtual environment:

```bash
python -m venv tp_mito_env
```

Activate it.

### Linux/macOS

```bash
source tp_mito_env/bin/activate
```

### Windows PowerShell

```powershell
tp_mito_env\Scripts\Activate.ps1
```

Install the required packages:

```bash
pip install -r requirements.txt
```

Start JupyterLab:

```bash
jupyter lab
```

Open:

```text
notebooks/TP_mitogenome_custom_mapping_reproducibility.ipynb
```

and select **Run All**.

An already executed version is also provided:

```text
notebooks/TP_mitogenome_custom_mapping_reproducibility_EXECUTED.ipynb
```

This allows the analysis outputs to be inspected without rerunning the notebook.

## Expected validation results

The validation workflow should reproduce the following principal results:

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

## Custom read-mapping method

The mapper in `scripts/custom_seed_mismatch_mapper.py` was used as an **assembly-validation and coverage-diagnostic tool**, not as the primary mitochondrial genome assembler.

The mapping procedure uses:

- an exact 19-nucleotide seed index constructed from the circularized reference;
- six seed positions distributed across each read;
- up to 80 highest-ranked candidate alignment starts;
- evaluation of both forward and reverse-complement orientations;
- full-read Hamming mismatch scoring without gapped alignment; and
- acceptance of the best placement when it contains no more than 12 mismatches.

The mitochondrial reference is treated as circular. Placements extending beyond the artificial linearization boundary are therefore wrapped back to the beginning of the genome.

## Interpretation of Illumina coverage

The short-read library was generated from **pooled overlapping PCR amplicons**, rather than from a randomly fragmented genomic library. Consequently, read depth is expected to vary substantially across the mitochondrial genome because of differences in amplification efficiency and representation of individual amplicons in the sequencing pool.

Coverage variation should therefore be interpreted primarily as **amplicon-representation variation**, rather than as conventional whole-genome shotgun coverage variation.

The diagnostic 325-bp interval at positions **14,937–15,261** is the empirically identified lowest-depth 325-bp window in the mapped Illumina data. It should not automatically be interpreted as the exact coordinates of any historical N-containing draft interval because the original historical draft sequence is not included in this repository.

## TP13 Sanger support

The TP13 chromatogram was quality trimmed with Biopython using the `abi-trim` reader. The quality-trimmed sequence is **651 bp** long.

Its reverse complement aligns to the final circular mitochondrial genome from position **14,968**, across the circular origin, to position **123**, with **one nucleotide mismatch** across the aligned segment.

TP13 Sanger sequencing is treated as an **independent sequence-evidence track**. It is not converted into artificial Illumina read depth. Therefore, an interval with low or zero Illumina depth may still have independent sequence support from the TP13 Sanger sequence.

## Primer and amplicon validation

The 15 primer pairs were remapped directly to the completed circular mitochondrial genome. These final-genome coordinates are reported in:

```text
derived_data/primer_amplicons_remapped_to_final_genome.csv
```

The origin-spanning amplicons include TP13, TP14, and TPJ. Their remapped coordinates can be used to verify sequence continuity across the circular genome origin and to interpret local coverage patterns.

## Figures

Separate figure files are not included in this reduced GitHub repository. Analytical plots can be regenerated from the supplied notebook using the raw and derived data.

The executed notebook also preserves the analysis outputs for inspection.

## Reproducibility scope

Please read:

```text
REPRODUCIBILITY_SCOPE.md
```

before interpreting the validation results.

In particular, this repository reproduces the **read-mapping validation, coverage calculations, primer remapping, and TP13 Sanger alignment**. It does not claim to reproduce the historical CAP3 assembly process byte-for-byte because the original intermediate assembly files are not available in this package.

## Citation and data availability

When using this repository, please cite the associated *Tarbinskiellus portentosus* mitochondrial-genome manuscript and the deposited mitochondrial sequence:

**GenBank accession: PZ447390.1**

After publication, the final article citation and DOI can be added here.

## Contact

For questions concerning this reproducibility package, please contact the corresponding authors listed in the associated manuscript.
