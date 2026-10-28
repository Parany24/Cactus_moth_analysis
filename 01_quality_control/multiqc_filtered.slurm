#!/bin/bash
#SBATCH --job-name=cactus_multiqc_filtered
#SBATCH --output=multiqc_filtered_%j.out
#SBATCH --error=multiqc_filtered_%j.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

reports_dir="${ORIGINAL_DIR}/QC/fastqc_filtered_reports"
output_dir="${OUTPUT_DIR}/QC/multiqc_filtered_reports"

if [[ ! -d "${reports_dir}" ]]; then
    echo "ERREUR : dossier FastQC introuvable : ${reports_dir}" >&2
    exit 1
fi

mkdir -p "${output_dir}"

module purge
module load bioinfo/MultiQC/1.35

multiqc "${reports_dir}" \
    --outdir "${output_dir}" \
    --filename cactus_moth_multiqc_filtered_report.html

echo "Rapport MultiQC enregistré dans : ${output_dir}/cactus_moth_multiqc_filtered_report.html"
