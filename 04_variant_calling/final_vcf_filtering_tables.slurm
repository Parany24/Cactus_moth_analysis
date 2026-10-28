#!/bin/bash
#SBATCH -J cactus_final_vcf
#SBATCH -o final_filtering_%j.out
#SBATCH -e final_filtering_%j.err
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --mail-user=paranymamionah@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================
# Final VCF filtering workflow for cactus moth ddRADseq data
#
# Input:
#   STACKS VCFs from populations
#
# Datasets:
#   1. all_snps
#   2. random_snp
#
# Filters:
#   1. SNPs only
#   2. biallelic SNPs only
#   3. max missing 0.80
#   4. MAF >= 0.05
#
# No heterozygosity filter is applied here.
#
# Outputs:
#   final filtered VCFs
#   Table_S2_filtering_summary.tsv
#   Table_S3_missing_data_summary.tsv
#   final .imiss and .lmiss files
# ============================================================

ORIGINAL_DIR="/work/user/gdobigny/Cactus_moth_RADseq"
OUTPUT_PROJECT_DIR="/work/user/gdobigny/Cactus_moth_analysis"

OUT_DIR="${OUTPUT_PROJECT_DIR}/Stacks/final_vcf_filtering"
mkdir -p "${OUT_DIR}"

module purge
module load bioinfo/VCFtools/0.1.16
module load bioinfo/Bcftools/1.9

# Input raw VCFs
VCF_ALL="${ORIGINAL_DIR}/Stacks/populations_all_snps/populations.snps.sorted.vcf.gz"
VCF_RANDOM="${ORIGINAL_DIR}/Stacks/populations_random_snp/populations.snps.sorted.vcf.gz"

# Filtering thresholds
MAX_MISSING="0.80"
MIN_MAF="0.05"

# Tables for manuscript
TABLE_S2="${OUT_DIR}/Table_S2_filtering_summary.tsv"
TABLE_S3="${OUT_DIR}/Table_S3_missing_data_summary.tsv"

# Initialize tables
echo -e "Dataset\tStep\tFiltering_step\tN_samples\tN_SNPs\tSNPs_removed_from_previous_step\tPercent_retained_from_previous_step\tPercent_retained_from_initial\tVCF_file" > "${TABLE_S2}"

echo -e "Dataset\tN_samples\tN_SNPs\tMean_individual_missing_percent\tSD_individual_missing_percent\tMin_individual_missing_percent\tMax_individual_missing_percent\tN_individuals_gt20_missing\tN_individuals_gt30_missing\tIMISS_file\tLMISS_file\tFinal_VCF" > "${TABLE_S3}"

# ------------------------------------------------------------
# Function to count samples and SNPs
# ------------------------------------------------------------

count_samples() {
    local vcf="$1"
    bcftools query -l "${vcf}" | wc -l
}

count_snps() {
    local vcf="$1"
    bcftools view -H "${vcf}" | wc -l
}

# ------------------------------------------------------------
# Function to add one filtering step to Table S2
# ------------------------------------------------------------

add_table_s2_row() {
    local dataset="$1"
    local step="$2"
    local filtering_step="$3"
    local vcf="$4"
    local previous_snps="$5"
    local initial_snps="$6"

    local n_samples
    local n_snps
    local removed
    local retained_previous
    local retained_initial

    n_samples=$(count_samples "${vcf}")
    n_snps=$(count_snps "${vcf}")

    if [[ "${step}" == "0" ]]; then
        removed=0
        retained_previous=100
        retained_initial=100
    else
        removed=$((previous_snps - n_snps))
        retained_previous=$(awk -v n="${n_snps}" -v p="${previous_snps}" 'BEGIN{if(p>0){print (n/p)*100}else{print "NA"}}')
        retained_initial=$(awk -v n="${n_snps}" -v i="${initial_snps}" 'BEGIN{if(i>0){print (n/i)*100}else{print "NA"}}')
    fi

    echo -e "${dataset}\t${step}\t${filtering_step}\t${n_samples}\t${n_snps}\t${removed}\t${retained_previous}\t${retained_initial}\t${vcf}" >> "${TABLE_S2}"
}

# ------------------------------------------------------------
# Function to summarize final missingness and add to Table S3
# ------------------------------------------------------------

add_table_s3_row() {
    local dataset="$1"
    local final_vcf="$2"
    local imiss_file="$3"
    local lmiss_file="$4"

    local n_samples
    local n_snps

    n_samples=$(count_samples "${final_vcf}")
    n_snps=$(count_snps "${final_vcf}")

    awk -v dataset="${dataset}" \
        -v nsamples="${n_samples}" \
        -v nsnps="${n_snps}" \
        -v imiss="${imiss_file}" \
        -v lmiss="${lmiss_file}" \
        -v finalvcf="${final_vcf}" '
    NR > 1 {
        miss = $5 * 100
        n += 1
        sum += miss
        sumsq += miss * miss

        if (n == 1 || miss < min) min = miss
        if (n == 1 || miss > max) max = miss

        if (miss > 20) gt20 += 1
        if (miss > 30) gt30 += 1
    }
    END {
        if (n > 1) {
            mean = sum / n
            sd = sqrt((sumsq - (sum * sum / n)) / (n - 1))
        } else if (n == 1) {
            mean = sum
            sd = 0
        } else {
            mean = "NA"
            sd = "NA"
            min = "NA"
            max = "NA"
        }

        print dataset "\t" nsamples "\t" nsnps "\t" mean "\t" sd "\t" min "\t" max "\t" gt20+0 "\t" gt30+0 "\t" imiss "\t" lmiss "\t" finalvcf
    }' "${imiss_file}" >> "${TABLE_S3}"
}

# ------------------------------------------------------------
# Main loop over datasets
# ------------------------------------------------------------

for DATASET in all_snps random_snp
do
    echo "============================================================"
    echo "Processing dataset: ${DATASET}"

    if [[ "${DATASET}" == "all_snps" ]]; then
        INPUT_VCF="${VCF_ALL}"
    else
        INPUT_VCF="${VCF_RANDOM}"
    fi

    if [[ ! -f "${INPUT_VCF}" ]]; then
        echo "ERROR: input VCF not found:"
        echo "${INPUT_VCF}"
        exit 1
    fi

    DATASET_DIR="${OUT_DIR}/${DATASET}"
    mkdir -p "${DATASET_DIR}"

    echo "Input VCF: ${INPUT_VCF}"
    echo "Output directory: ${DATASET_DIR}"

    # Make sure input VCF is indexed
    if [[ ! -f "${INPUT_VCF}.tbi" && ! -f "${INPUT_VCF}.csi" ]]; then
        echo "Indexing input VCF"
        tabix -p vcf "${INPUT_VCF}"
    fi

    # Step 0: initial raw STACKS VCF
    INITIAL_SNPS=$(count_snps "${INPUT_VCF}")
    PREVIOUS_SNPS="${INITIAL_SNPS}"

    add_table_s2_row \
        "${DATASET}" \
        "0" \
        "Initial STACKS VCF" \
        "${INPUT_VCF}" \
        "${PREVIOUS_SNPS}" \
        "${INITIAL_SNPS}"

    # Step 1: retain SNPs only
    STEP1_VCF="${DATASET_DIR}/${DATASET}.step1.snps_only.vcf.gz"

    bcftools view \
        -v snps \
        -Oz \
        -o "${STEP1_VCF}" \
        "${INPUT_VCF}"

    tabix -p vcf "${STEP1_VCF}"

    add_table_s2_row \
        "${DATASET}" \
        "1" \
        "Retained SNP variants only" \
        "${STEP1_VCF}" \
        "${PREVIOUS_SNPS}" \
        "${INITIAL_SNPS}"

    PREVIOUS_SNPS=$(count_snps "${STEP1_VCF}")

    # Step 2: retain biallelic SNPs only
    STEP2_VCF="${DATASET_DIR}/${DATASET}.step2.biallelic_snps.vcf.gz"

    bcftools view \
        -m2 \
        -M2 \
        -v snps \
        -Oz \
        -o "${STEP2_VCF}" \
        "${STEP1_VCF}"

    tabix -p vcf "${STEP2_VCF}"

    add_table_s2_row \
        "${DATASET}" \
        "2" \
        "Retained biallelic SNPs only" \
        "${STEP2_VCF}" \
        "${PREVIOUS_SNPS}" \
        "${INITIAL_SNPS}"

    PREVIOUS_SNPS=$(count_snps "${STEP2_VCF}")

    # Step 3: retain SNPs genotyped in at least 80% of individuals
    STEP3_PREFIX="${DATASET_DIR}/${DATASET}.step3.maxmissing_${MAX_MISSING}"
    STEP3_VCF="${STEP3_PREFIX}.recode.vcf"

    vcftools \
        --gzvcf "${STEP2_VCF}" \
        --max-missing "${MAX_MISSING}" \
        --recode \
        --recode-INFO-all \
        --out "${STEP3_PREFIX}"

    bgzip -f "${STEP3_VCF}"
    STEP3_VCF_GZ="${STEP3_VCF}.gz"
    tabix -p vcf "${STEP3_VCF_GZ}"

    add_table_s2_row \
        "${DATASET}" \
        "3" \
        "Retained SNPs genotyped in at least 80% of individuals" \
        "${STEP3_VCF_GZ}" \
        "${PREVIOUS_SNPS}" \
        "${INITIAL_SNPS}"

    PREVIOUS_SNPS=$(count_snps "${STEP3_VCF_GZ}")

    # Step 4: retain SNPs with MAF >= 0.05
    STEP4_PREFIX="${DATASET_DIR}/${DATASET}.step4.maf_${MIN_MAF}"
    STEP4_VCF="${STEP4_PREFIX}.recode.vcf"

    vcftools \
        --gzvcf "${STEP3_VCF_GZ}" \
        --maf "${MIN_MAF}" \
        --recode \
        --recode-INFO-all \
        --out "${STEP4_PREFIX}"

    bgzip -f "${STEP4_VCF}"
    STEP4_VCF_GZ="${STEP4_VCF}.gz"
    tabix -p vcf "${STEP4_VCF_GZ}"

    add_table_s2_row \
        "${DATASET}" \
        "4" \
        "Retained SNPs with MAF >= 0.05" \
        "${STEP4_VCF_GZ}" \
        "${PREVIOUS_SNPS}" \
        "${INITIAL_SNPS}"

    # Final VCF
    FINAL_VCF="${DATASET_DIR}/${DATASET}.final.filtered.vcf.gz"
    FINAL_TBI="${FINAL_VCF}.tbi"

    cp "${STEP4_VCF_GZ}" "${FINAL_VCF}"
    cp "${STEP4_VCF_GZ}.tbi" "${FINAL_TBI}"

    # Final missingness reports
    FINAL_PREFIX="${DATASET_DIR}/${DATASET}.final.filtered"

    vcftools \
        --gzvcf "${FINAL_VCF}" \
        --missing-indv \
        --out "${FINAL_PREFIX}"

    vcftools \
        --gzvcf "${FINAL_VCF}" \
        --missing-site \
        --out "${FINAL_PREFIX}"

    IMISS_FILE="${FINAL_PREFIX}.imiss"
    LMISS_FILE="${FINAL_PREFIX}.lmiss"

    add_table_s3_row \
        "${DATASET}" \
        "${FINAL_VCF}" \
        "${IMISS_FILE}" \
        "${LMISS_FILE}"

    echo "Finished dataset: ${DATASET}"
    echo "Final VCF: ${FINAL_VCF}"
    echo "Final SNP number: $(count_snps "${FINAL_VCF}")"
    echo "Final sample number: $(count_samples "${FINAL_VCF}")"
done

echo "============================================================"
echo "Final filtering workflow finished."
echo
echo "Table S2:"
echo "${TABLE_S2}"
cat "${TABLE_S2}"
echo
echo "Table S3:"
echo "${TABLE_S3}"
cat "${TABLE_S3}"
echo
echo "Final VCFs:"
ls -lh "${OUT_DIR}/all_snps/all_snps.final.filtered.vcf.gz"*
ls -lh "${OUT_DIR}/random_snp/random_snp.final.filtered.vcf.gz"*
