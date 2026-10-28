#!/bin/bash
#SBATCH --job-name=fastqc_filtered
#SBATCH --output=fastqc_filtered_%j.out
#SBATCH --error=fastqc_filtered_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# ============================================================
# FastQC quality control for filtered cactus moth RADseq reads
#
# Input:
#   Clean_data/process_radtags_adapters/
#
# Output:
#   QC/fastqc_filtered_reports/
#
# The project directory is detected automatically from the
# location of this script.
# ============================================================

set -euo pipefail

# Define project directory
ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

# Define input and output directories
raw_data_dir="${ORIGINAL_DIR}/Clean_data/process_radtags_adapters"
output_dir="${OUTPUT_DIR}/QC/fastqc_filtered_reports"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Load FastQC
module purge
module load bioinfo/FastQC/0.12.1

# Move to filtered data directory
cd "${raw_data_dir}" || exit 1

# Loop through all R1 filtered files
for r1_file in *_1.fq.gz
do
    # Derive sample name by removing _1.fq.gz
    sample=$(basename "${r1_file}" _1.fq.gz)

    # Define corresponding R2 file
    r2_file="${sample}_2.fq.gz"

    # Check that R2 exists
    if [[ -f "${r2_file}" ]]
    then
        echo "Processing sample: ${sample}"
        echo "R1: ${r1_file}"
        echo "R2: ${r2_file}"

        # Run FastQC on paired-end filtered files
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
