#!/bin/bash
#SBATCH -J cactus_bwa
#SBATCH -o bwa_%A_%a.out
#SBATCH -e bwa_%A_%a.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-137
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================
# BWA-MEM mapping of filtered cactus moth RADseq reads
#
# Input:
#   Clean_data/process_radtags_adapters/
#   genome_reference_cactorum/
#   metadata/samples_clean_bwa.txt
#
# Output:
#   Mapping/bam/
#   Mapping/stats/
#
# The project directory is detected automatically from the
# location of this script.
# ============================================================

# Define project directory
ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

# Define input files and directories
CLEAN_DIR="${ORIGINAL_DIR}/Clean_data/process_radtags_adapters"
REF="${ORIGINAL_DIR}/genome_reference_cactorum/Ccactorum.fcs_gx_2026_FINAL.fa"
SAMPLE_LIST="${ORIGINAL_DIR}/metadata/samples_clean_bwa.txt"

# Define output directories
BAM_DIR="${OUTPUT_DIR}/Mapping/bam"
STATS_DIR="${OUTPUT_DIR}/Mapping/stats"

# Create output directories
mkdir -p "${BAM_DIR}" "${STATS_DIR}"

# Load required software
module purge
module load bioinfo/bwa/0.7.19
module load bioinfo/samtools/1.20

# Get sample corresponding to the SLURM array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# Define paired-end input files
R1="${CLEAN_DIR}/${SAMPLE}.1.fq.gz"
R2="${CLEAN_DIR}/${SAMPLE}.2.fq.gz"

echo "======================================"
echo "SLURM array task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Reference: ${REF}"
echo "======================================"

# Check that input files exist
if [[ ! -f "${R1}" ]]; then
    echo "ERROR: R1 file not found: ${R1}"
    exit 1
fi

if [[ ! -f "${R2}" ]]; then
    echo "ERROR: R2 file not found: ${R2}"
    exit 1
fi

if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference genome not found: ${REF}"
    exit 1
fi

# Mapping with BWA-MEM and direct sorting with samtools
bwa mem \
    -t 8 \
    -M \
    "${REF}" \
    "${R1}" \
    "${R2}" \
    | samtools sort \
        -@ 8 \
        -o "${BAM_DIR}/${SAMPLE}.sorted.bam" \
        -

# Index sorted BAM
samtools index "${BAM_DIR}/${SAMPLE}.sorted.bam"

# Generate mapping statistics
samtools flagstat \
    "${BAM_DIR}/${SAMPLE}.sorted.bam" \
    > "${STATS_DIR}/${SAMPLE}.flagstat.txt"

samtools stats \
    "${BAM_DIR}/${SAMPLE}.sorted.bam" \
    > "${STATS_DIR}/${SAMPLE}.samtools_stats.txt"

echo "Mapping finished successfully for ${SAMPLE}"
