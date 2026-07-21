# Watzlowik et al. 2023 MNase-seq — Nextflow DSL2

This is a modular DSL2 conversion of `MNase-Seq/Snakefile` from the uploaded Watzlowik et al. repository. Every executable workflow step is defined in its own `modules/*.nf` file. All R scripts from the authors' `1_code` tree are placed in `bin/`; DANPOS3 is expected at `bin/DANPOS3/danpos.py`. Only `AdjustFDR.R` is part of the main Snakemake-derived DAG; the other R scripts remain optional downstream utilities.

## Workflow

1. FastQC on raw R1/R2 reads
2. Paired-end adapter trimming with Skewer
3. Bowtie2 reference indexing
4. Bowtie2 alignment
5. SAM-to-BAM conversion, sorting and indexing
6. Qualimap and MultiQC on mapped BAMs
7. `alignmentSieve` mononucleosome filtering (75–175 bp by default, proper pairs, MAPQ ≥20)
8. Sorting/indexing of filtered BAMs
9. Qualimap and MultiQC on filtered BAMs
10. Condition-aware staging of filtered BAMs
11. DANPOS3 `dpos` differential nucleosome analysis
12. Adjusted-FDR calculation with the authors' `AdjustFDR.R`

## Required samplesheet

CSV and TSV are both accepted. Required columns are:

| Column | Meaning |
|---|---|
| `Samples` | Unique safe sample ID |
| `R1` | R1 FASTQ/FASTQ.GZ path |
| `R2` | R2 FASTQ/FASTQ.GZ path |
| `Condition` | Group name used for DANPOS directories |
| `Replicate` | Replicate label |

Relative FASTQ paths are resolved relative to the samplesheet directory. `Samples`, `Condition`, and `Replicate` may contain letters, numbers, `.`, `_`, and `-` only.

**The comparison names must exactly match the `Condition` column.** For example, `iKO_24hr` is not the same as `iKO_24h`.

## DANPOS comparison semantics

```bash
--danpos_comparison 'iKO_24h:noKO_24h'
```

is translated to:

```text
iKO_24h/:noKO_24h/
```

and DANPOS is run as:

```bash
python3 bin/DANPOS3/danpos.py dpos iKO_24h/:noKO_24h/ -o results \
  -m 1 -c 51137214 --extend 73 -a 1 -z 20 -jd 100 -e 1
```

Reversing the parameter:

```bash
--danpos_comparison 'noKO_24h:iKO_24h'
```

creates and passes:

```text
noKO_24h/:iKO_24h/
```

The left side is treated as the DANPOS sample/treatment and the right side as the control. Additional conditions in the sheet are processed through QC/alignment but are not included in this single DANPOS contrast.

## Installation

Install Nextflow, Java 17+ and Mamba/Conda. Then install DANPOS3 into the required `bin` location:

```bash
chmod +x scripts/*.sh bin/*.R
./scripts/install_danpos3.sh
```

The source archive itself did not contain DANPOS3; its original preparation script cloned the official repository.

## Run

```bash
nextflow run main.nf \
  -profile conda \
  --input assets/samplesheet.tsv \
  --reference /path/to/PlasmoDB-52_Pfalciparum3D7_Genome.fasta \
  --annotation /path/to/PlasmoDB-52_Pfalciparum3D7.gff \
  --chrom_sizes /path/to/sizes.genome \
  --danpos_comparison 'iKO_24h:noKO_24h' \
  --outdir results
```

Resume an interrupted run with `-resume`.

## Test reads from the uploaded repository

The eight small paired-end FASTQ files are retained under `test_data/reads/`, with a matching `assets/test_samplesheet.tsv`. Reference files are not present in the uploaded archive. They can be prepared using:

```bash
./scripts/prepare_plasmodb52_reference.sh
```

Then run:

```bash
nextflow run main.nf -profile conda \
  --input assets/test_samplesheet.tsv \
  --reference reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7_Genome.fasta \
  --annotation reference/plasmoDB52/PlasmoDB-52_Pfalciparum3D7.gff \
  --chrom_sizes reference/plasmoDB52/sizes.genome \
  --danpos_comparison 'iKO_24h:noKO_24h' \
  --outdir test_results
```

## Main output

```text
results/
├── 00_index/bowtie2/
├── 00_logs/
├── 00_qc/
│   ├── preprocessing/raw_fastqc/
│   ├── mapping/{qualimap,multiqc}/
│   ├── filtered/{deeptools_alignmentSieve,qualimap,multiqc}/
│   └── nucleosome_positioning/{danpos.log,adjustFDR.log}
├── 01_trimmed_reads/
├── 02_alignment/sort/
├── 03_alignment_filtered/bam/
└── 04_differential_nucleosome_positioning/
    ├── <treatment>/
    ├── <control>/
    └── results/
        ├── pooled/
        ├── diff/
        ├── <treatment>-<control>.positions.integrative.xls
        └── adjustedFDR/
            └── <treatment>-<control>.positions.FDR_corrected.xls
```

## Static validation

Run the bundled integrity checker after unpacking:

```bash
python scripts/validate_conversion.py
```

When the original source zip is available, byte-level checks of the author R scripts and the source rule count can also be run with `--source-zip /path/to/archive.zip`. The generated `VALIDATION_REPORT.txt` records the checks performed for this delivered archive.

## Parameter overrides

The original DANPOS settings are the default. Override them only when scientifically justified:

```bash
--danpos_extra '-m 1 -c 51137214 --extend 73 -a 1 -z 20 -jd 100 -e 1'
```

Other configurable defaults are in `nextflow.config`, including adapters, fragment length, MAPQ, Qualimap memory and adjusted-FDR simulation settings.

## Optional author analysis code

The original R utilities `find_dynamic_nucleosomes.R` and `get_nucleosomeDynamics_gff.R` are in `bin/`, as requested, but are not invoked by the main workflow because they were not Snakemake rules and contain analysis-specific assumptions. The manuscript R Markdown and Docker-based nucleosomeDynamics helpers are retained under `extras/`.
