#!/bin/bash
#SBATCH --job-name=cactus_radtags
#SBATCH --output=process_radtags_adapters_%j.out
#SBATCH --error=process_radtags_adapters_%j.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=48:00:00

set -euo pipefail

ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

RAW_DIR="${ORIGINAL_DIR}/Raw_data"
CLEAN_DIR="${OUTPUT_DIR}/Clean_data/process_radtags_adapters"

if [[ ! -d "${RAW_DIR}" ]]; then
    echo "ERREUR : données brutes introuvables : ${RAW_DIR}" >&2
    exit 1
fi

mkdir -p "${CLEAN_DIR}"

module purge
module load bioinfo/Stacks/2.68

process_radtags \
    -P \
    -p "${RAW_DIR}" \
    -o "${CLEAN_DIR}" \
    --renz-1 nspI \
    --renz-2 dpnII \
    -c \
    -q \
    -r \
    --adapter-1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter-2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --adapter-mm 2 \
    --threads 4

echo "process_radtags terminé."
echo "Fichiers nettoyés : ${CLEAN_DIR}"
