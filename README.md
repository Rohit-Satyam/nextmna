# nextmna: Nextflow workflow for differential nucleosome dynamics analysis

> Note: This pipeline is under active development.

This is a modular DSL2 pipeline made using inspiration from `Watzlowik_et_al_2023` code available [here](https://github.com/SimHolz/Watzlowik_et_al_2023) to reproduce their analysis and reuse it for in-house data analysis since snakemake pipeline needed update. Every executable workflow step is defined in its own `modules/*.nf` file. All R scripts from the authors' `1_code` tree are placed in `bin/`; DANPOS3 is expected at `bin/DANPOS3/danpos.py`. Only `AdjustFDR.R` is part of the main Snakemake-derived DAG; the other R scripts remain optional downstream utilities. This pipeline is tested using `Nextflow version 26.04.6` and might throw error related to syntax for versions `< 25`. The docker image required by the workflow can be found at https://zenodo.org/records/21476365.

## Workflow

1. FastQC on raw R1/R2 reads
2. Paired-end adapter trimming with Fastp
3. Bowtie2 reference indexing
4. Bowtie2 alignment
5. SAM-to-BAM conversion, sorting and indexing
6. Qualimap and MultiQC on mapped BAMs
7. `alignmentSieve` mononucleosome filtering (75–175 bp by default, proper pairs, MAPQ ≥20)
8. Sorting/indexing of filtered BAMs
9. Qualimap and MultiQC on filtered BAMs
10. Condition-aware staging of filtered BAMs
11. DANPOS3 `dpos` differential nucleosome analysis. We are not using the latest DANPOS which is being actively maintained [here](https://github.com/boenc28-cmyk/DANPOS) because we observed some discrepencies and therefore we continue to use the code being used by [everyone](https://github.com/sklasfeld/DANPOS3). 
12. Adjusted-FDR calculation with the authors' `AdjustFDR.R`
13. Nucleosome dynamics using `nucleR` and `NucDyn`

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

Install Nextflow, Java 17+ and Mamba/Conda. The dependencies for each modules are handled by their private `.yml` files. DANPOS3 code is included in this repository

```bash
conda install conda-forge::mamba
mamba create -n nextmna bioconda::nextflow==26.04.6
```

## Run

```bash
nextflow run main.nf \
--input di_18hpi.csv \
--outdir 18hpi_results \
--danpos_comparison='RAPA_18hpi:DMSO_18hpi' \
--nucdyn_comparison='RAPA_18hpi:DMSO_18hpi' \
--nucdyn_contrast_name=18hpi \
-with-conda \
--nucdyn_expected_replicates=3 \
-w work_18hpi \
--fastp_ext='--detect_adapter_for_pe --correction --adapter_fasta ./resources/adapters/TruSeq3-PE.fa' \
--danpos_extra='-m 1 -c 23332839 --extend 70 -a 1 -z 20 -jd 75 -jw 10 -e 1' \
--annotation="/data/di_mnaseq/Watzlowik_et_al_2023/nextmna/resources/PlasmoDB-68_Pfalciparum3D7.gff" \
--nucdyn_annotation="/data/di_mnaseq/Watzlowik_et_al_2023/nextmna/resources/pf_nucDyn_v68.gff"
```

Resume an interrupted run with `-resume`.

## Test reads from the uploaded repository

The eight small paired-end FASTQ files are retained under `data/reads/`, with a matching `assets/samplesheet.tsv`. Reference files are not present in the uploaded archive. User need to provide the reference fasta and GFF file. A separate GFF file is required for nucleodynamics step called `pf_nucDyn.gff` containing `tss` and `tts` information. This can be produced using code below. 

```r
library(tidyverse)
library(plyranges)

inGFF <- read_gff(file = "PlasmoDB-52_Pfalciparum3D7.gff")

inGFF %>% 
  filter(type %in% c("protein_coding_gene", "pseudogene", "ncRNA_gene")) %>% 
  mutate(tss = ifelse(strand == "+", start,end),
         tts = ifelse(strand == "+", end,start),
		 name = ID) %>% 
  write_gff3(file = "pf_nucDyn_v52.gff")
```

## Main output

```text
── 00_index
│   └── bowtie2
├── 00_logs
│   ├── bowtie2
│   ├── bowtie2_index
│   ├── deeptools_alignmentSieve
│   ├── fastqc
│   ├── qualimap
│   └── qualimap_filtered
├── 00_qc
│   ├── fastp
│   ├── filtered
│   ├── mapping
│   ├── nucleosome_positioning
│   └── preprocessing
├── 01_trimmed_reads
│   ├── PfSnf2L-0hpi_iKO_Rep1-trimmed-pair1.fastq.gz
│   ├── PfSnf2L-0hpi_iKO_Rep1-trimmed-pair2.fastq.gz
│   ├── PfSnf2L-0hpi_iKO_Rep2-trimmed-pair1.fastq.gz
│   ├── PfSnf2L-0hpi_iKO_Rep2-trimmed-pair2.fastq.gz
│   ├── PfSnf2L-0hpi_iKO_Rep3-trimmed-pair1.fastq.gz
│   ├── PfSnf2L-0hpi_iKO_Rep3-trimmed-pair2.fastq.gz
│   ├── PfSnf2L-0hpi_noKO_Rep1-trimmed-pair1.fastq.gz
│   ├── PfSnf2L-0hpi_noKO_Rep1-trimmed-pair2.fastq.gz
│   ├── PfSnf2L-0hpi_noKO_Rep2-trimmed-pair1.fastq.gz
│   ├── PfSnf2L-0hpi_noKO_Rep2-trimmed-pair2.fastq.gz
│   ├── PfSnf2L-0hpi_noKO_Rep3-trimmed-pair1.fastq.gz
│   └── PfSnf2L-0hpi_noKO_Rep3-trimmed-pair2.fastq.gz
├── 02_alignment
│   └── sort
├── 03_alignment_filtered
│   └── bam
├── 04_differential_nucleosome_positioning
│   ├── 0hpi_iKO
│   ├── 0hpi_noKO
│   └── results
├── NucleosomeDynamics
│   ├── bam
│   ├── NucDyn
│   ├── nucleR
│   ├── RData
│   ├── statistics
│   └── txstart
├── pipeline_dag.html
├── pipeline_report.html
├── pipeline_timeline.html
└── pipeline_trace.txt
```

## Parameter overrides

The original DANPOS settings are the default. Override them only when scientifically justified:

```bash
--danpos_extra '-m 1 -c 51137214 --extend 73 -a 1 -z 20 -jd 100 -e 1'
```

Other configurable defaults are in `nextflow.config`, including adapters, fragment length, MAPQ, Qualimap memory and adjusted-FDR simulation settings.

## Optional author analysis code

The original R utilities `find_dynamic_nucleosomes.R` has been modified and is available in `assets/filter_nucleosomedynamics_results.Rmd`. The manuscript R Markdown and Docker-based nucleosomeDynamics helpers are retained under `extras/` in case user wishes to build docker image him/herself.

## Citations

We borrow most of the code from the following article.

```
Watzlowik, M.T., Silberhorn, E., Das, S., Singhal, R., Venugopal, K., Holzinger, S., Stokes, B., Schadt, E., Sollelis, L., Bonnell, V.A. and Gow, M., 2025. Plasmodium blood stage development requires the chromatin remodeller Snf2L. Nature, 639(8056), pp.1069-1075.
```
