#!/bin/bash
#SBATCH --job-name=admixture_K
#SBATCH --output=admixture_K_%A_%a.out
#SBATCH --error=admixture_K_%A_%a.err
#SBATCH --array=2-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

module load bioinfo/ADMIXTURE/1.3.0

cd /work/user/gdobigny/Cactus_moth_analysis/05_population_structure

K=${SLURM_ARRAY_TASK_ID}

echo "Running ADMIXTURE K=${K}"
date

admixture --cv /work/user/gdobigny/Cactus_moth_RADseq/Stacks/population_structure/admixture_LD_numeric.bed "${K}" \
    | tee "admixture_K${K}.log"

echo "Finished K=${K}"
date
