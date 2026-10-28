#!/bin/bash
#SBATCH --job-name=cactus_cov
#SBATCH --output=coverage_%A_%a.out
#SBATCH --error=coverage_%A_%a.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL

set -euo pipefail

ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

BAM_DIR="${ORIGINAL_DIR}/Mapping/bam"
SAMPLE_LIST="${ORIGINAL_DIR}/metadata/samples_bam.txt"
REF_FAI="${ORIGINAL_DIR}/genome_reference_cactorum/Ccactorum.fcs_gx_2026_FINAL.fa.fai"
OUT_DIR="${OUTPUT_DIR}/Coverage/per_sample"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERREUR : SLURM_ARRAY_TASK_ID n'est pas défini." >&2
    exit 1
fi

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERREUR : liste d'échantillons introuvable : ${SAMPLE_LIST}" >&2
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
if [[ -z "${SAMPLE}" ]]; then
    echo "ERREUR : aucun échantillon à la ligne ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
OUT="${OUT_DIR}/${SAMPLE}.coverage.tsv"

for f in "${BAM}" "${BAM}.bai" "${REF_FAI}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERREUR : fichier introuvable : ${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUT_DIR}"

module purge
module load bioinfo/samtools/1.20

GENOME_SIZE=$(awk '{sum += $2} END {print sum}' "${REF_FAI}")

echo "Sample: ${SAMPLE}"
echo "BAM: ${BAM}"
echo "Genome size: ${GENOME_SIZE}"

printf 'SRA\tgenome_size_bp\ttotal_depth_sum\tgenome_wide_mean_depth\tcovered_bases\tbreadth_covered_percent\tmean_depth_covered_sites\n' > "${OUT}"

samtools depth -aa "${BAM}" |
awk -v sample="${SAMPLE}" -v genome_size="${GENOME_SIZE}" '
BEGIN { OFS="\t" }
{
    depth_sum += $3
    if ($3 > 0) {
        covered_bases += 1
        covered_depth_sum += $3
    }
}
END {
    genome_wide_mean_depth = depth_sum / genome_size
    breadth_covered_percent = (covered_bases / genome_size) * 100
    if (covered_bases > 0) {
        mean_depth_covered_sites = covered_depth_sum / covered_bases
    } else {
        mean_depth_covered_sites = 0
    }
    print sample, genome_size, depth_sum, genome_wide_mean_depth, covered_bases, breadth_covered_percent, mean_depth_covered_sites
}' >> "${OUT}"

echo "Coverage finished for ${SAMPLE}"
