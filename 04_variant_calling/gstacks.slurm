#!/bin/bash
#SBATCH -J cactus_gstacks
#SBATCH -o gstacks_%j.out
#SBATCH -e gstacks_%j.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================
# Stacks gstacks for cactus moth RADseq data
#
# Input:
#   Stacks/bam_biological_names/
#   Stacks/popmap_biological_names.tsv
#
# Output:
#   Stacks/gstacks/
#
# The project directory is detected automatically from the
# location of this script.
# ============================================================

ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_PROJECT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

BAM_DIR="${ORIGINAL_DIR}/Stacks/bam_biological_names"
POPMAP="${ORIGINAL_DIR}/Stacks/popmap_biological_names.tsv"
GSTACKS_OUT="${OUTPUT_PROJECT_DIR}/Stacks/gstacks"

mkdir -p "${GSTACKS_OUT}"

module purge
module load bioinfo/Stacks/2.68

echo "Running gstacks"
echo "BAM directory: ${BAM_DIR}"
echo "BAM suffix: .sorted.bam"
echo "Popmap: ${POPMAP}"
echo "Output directory: ${GSTACKS_OUT}"

gstacks \
    -I "${BAM_DIR}" \
    -S .sorted.bam \
    -M "${POPMAP}" \
    -O "${GSTACKS_OUT}" \
    -t 16 \
    --min-mapq 10

echo "gstacks finished successfully."
