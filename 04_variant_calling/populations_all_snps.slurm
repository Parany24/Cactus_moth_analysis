#!/bin/bash
#SBATCH -J cactus_pop_allsnps
#SBATCH -o populations_all_snps_%j.out
#SBATCH -e populations_all_snps_%j.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================
# Stacks populations: all SNPs per locus
#
# Input:
#   Stacks/gstacks/
#   Stacks/popmap_biological_names.tsv
#
# Output:
#   Stacks/populations_all_snps/
#
# The project directory is detected automatically from the
# location of this script.
# ============================================================

ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_PROJECT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

GSTACKS_DIR="${OUTPUT_PROJECT_DIR}/Stacks/gstacks"
POPMAP="${ORIGINAL_DIR}/Stacks/popmap_biological_names.tsv"
POP_OUT="${OUTPUT_PROJECT_DIR}/Stacks/populations_all_snps"

mkdir -p "${POP_OUT}"

module purge
module load bioinfo/Stacks/2.68

echo "Running populations: all SNPs per locus"
echo "Input gstacks directory: ${GSTACKS_DIR}"
echo "Popmap: ${POPMAP}"
echo "Output directory: ${POP_OUT}"

populations \
    -P "${GSTACKS_DIR}" \
    -M "${POPMAP}" \
    -O "${POP_OUT}" \
    -t 16 \
    -p 2 \
    -r 0.80 \
    --min-maf 0.05 \
    --max-obs-het 0.70 \
    --vcf \
    --structure \
    --plink \
    --genepop \
    --fstats

echo "populations all-SNPs run finished successfully."
echo "Output files:"
ls -lh "${POP_OUT}"
