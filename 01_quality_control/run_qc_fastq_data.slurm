#!/bin/bash
#SBATCH -J cactus_fastqc
#SBATCH -o fastqc_%j.out
#SBATCH -e fastqc_%j.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# ============================================================
# FastQC quality control for cactus moth RADseq FASTQ files
#
# Input:
#   Raw_data/
#
# Output:
#   QC/fastqc_reports/
#
# The project directory is detected automatically from the
# location of this script.
# ============================================================

set -euo pipefail

# Define project directory
ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

# Define input and output directories
raw_data_dir="${ORIGINAL_DIR}/Raw_data"
output_dir="${OUTPUT_DIR}/QC/fastqc_reports"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Load FastQC
module purge
module load bioinfo/FastQC/0.12.1

# Move to raw data directory
cd "${raw_data_dir}" || exit 1

# Loop through all R1 FASTQ files
for r1_file in *_1.fastq.gz
do
    sample=$(basename "${r1_file}" _1.fastq.gz)
    r2_file="${sample}_2.fastq.gz"

    if [[ -f "${r2_file}" ]]
    then
        echo "Processing sample: ${sample}"
        echo "R1: ${r1_file}"
        echo "R2: ${r2_file}"

        fastqc \
            --threads 4 \
            "${r1_file}" "${r2_file}" \
            --outdir "${output_dir}"
    else
        echo "WARNING: R2 file not found for ${sample}, skipping ${r1_file}"
    fi
done

echo "FastQC analysis complete."
echo "Reports saved in: ${output_dir}"
