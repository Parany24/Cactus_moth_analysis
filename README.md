
# Cactus Moth Population Genomics

## Project Overview

This repository contains scripts for the bioinformatic analysis of RAD-seq data from cactus moth populations. The workflow covers sequencing quality control, read processing, read mapping, variant calling, population structure, and genetic diversity analyses.

## Project Structure

```text
Cactus_moth_analysis/
├── 01_quality_control/
├── 02_read_processing/
├── 03_mapping/
├── 04_variant_calling/
├── 05_population_structure/
├── 06_genetic_diversity/
└── README.md
```

### 01 — Quality Control
Scripts for assessing sequencing read quality and generating quality control reports using FastQC and MultiQC.

### 02 — Read Processing
Scripts for processing RAD-seq reads, including adapter handling and read filtering with process_radtags.

### 03 — Read Mapping
Scripts for estimating sequencing coverage and mapping reads to a reference genome using BWA.

### 04 — Variant Calling
Scripts for genotype and variant processing using Stacks, including gstacks, populations, and downstream VCF filtering and summary tables.

### 05 — Population Structure
Scripts for investigating population structure using ADMIXTURE.

### 06 — Genetic Diversity
R scripts for calculating and visualizing genetic diversity statistics, including:
- Observed heterozygosity (Ho)
- Expected heterozygosity (He)
- Nucleotide diversity (Pi)
- Individual nucleotide diversity and plots

## Software and Tools

The workflow uses a combination of command-line tools and R scripts, including:
- FastQC and MultiQC
- Stacks
- BWA
- ADMIXTURE
- R

Exact software versions and execution parameters should be documented alongside the relevant scripts or in a dedicated environment file.

## Execution

Several analyses are designed to run on a computing cluster using SLURM job scripts (`.slurm` files).

Before running a script, check its input paths, reference files, software modules, resource requirements, and output directories. Run scripts in the appropriate workflow order.

## Data and Results

Large intermediate files and analysis outputs are not tracked in this Git repository. The `.gitignore` file excludes selected data files, result files, and job logs.

Raw sequencing data, reference genomes, and other required input files must be stored separately.

## Reproducibility

This repository is intended to keep analysis scripts organized and version-controlled. Record any changes to scripts, parameters, and software versions to support reproducibility.

## Repository

GitHub: https://github.com/Parany24/Cactus_moth_analysis

