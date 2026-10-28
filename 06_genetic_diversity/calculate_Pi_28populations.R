#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Nucleotide diversity (Pi) for 28 populations
# No bcftools required
#
# Same input files and general workflow as the Ho script.
#
# For each common SNP and population:
#
# Pi = 2 * n_ref * n_alt / [n_chr * (n_chr - 1)]
#
# where:
#   n_ref = number of reference alleles
#   n_alt = number of alternate alleles
#   n_chr = number of called chromosomes
#
# Population summary:
#   Mean_Pi = mean per-SNP Pi over usable SNPs
#   SD_Pi   = standard deviation of per-SNP Pi
#
# IMPORTANT:
# The calculation averages Pi over the selected common SNPs.
# It is not a whole-genome, per-base estimate if the SNP list
# contains only selected polymorphic sites.
# ============================================================


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

project_dir <- "/work/user/gdobigny/Cactus_moth_RADseq"

vcf_file <- file.path(
  project_dir,
  "Stacks/final_vcf_filtering/all_snps/all_snps.final.filtered.vcf.gz"
)

snp_file <- file.path(
  project_dir,
  "06_results/heterozygosity_He/common_SNPs/common_SNP_list.tsv"
)

metadata_file <- file.path(
  project_dir,
  "06_results/heterozygosity_He/common_SNPs/amova_metadata.tsv"
)

output_file <- file.path(
  project_dir,
  "06_results/nucleotide_diversity_Pi/common_SNPs/Pi_28populations.tsv"
)


# ------------------------------------------------------------
# 2. Check input files
# ------------------------------------------------------------

for (f in c(vcf_file, snp_file, metadata_file)) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}


# ------------------------------------------------------------
# 3. Read the common SNP list
# ------------------------------------------------------------

cat("Reading common SNP list...\n")

common_snps <- fread(snp_file)

if (!"ID" %in% names(common_snps)) {
  stop("The common SNP list must contain a column named ID.")
}

common_ids <- unique(
  trimws(as.character(common_snps$ID))
)

common_ids <- common_ids[
  !is.na(common_ids) &
    common_ids != "" &
    common_ids != "."
]

if (length(common_ids) == 0) {
  stop("The common SNP list is empty.")
}

cat("SNPs in common list:", length(common_ids), "\n")


# ------------------------------------------------------------
# 4. Read metadata
# ------------------------------------------------------------

cat("Reading metadata...\n")

metadata <- fread(metadata_file)

if (!all(c("Sample", "Pop") %in% names(metadata))) {
  stop("Metadata must contain columns named Sample and Pop.")
}

metadata[, sample := trimws(as.character(Sample))]
metadata[, population := trimws(as.character(Pop))]

metadata <- metadata[
  !is.na(sample) & sample != "" &
    !is.na(population) & population != ""
]

if (anyDuplicated(metadata$sample)) {

  duplicated_samples <- unique(
    metadata$sample[duplicated(metadata$sample)]
  )

  stop(
    "Duplicate sample IDs in metadata: ",
    paste(head(duplicated_samples, 10), collapse = ", ")
  )
}

population_ids <- sort(unique(metadata$population))

cat("Individuals in metadata:", nrow(metadata), "\n")
cat("Populations in metadata:", length(population_ids), "\n")


# ------------------------------------------------------------
# 5. Read the VCF header and identify samples
# ------------------------------------------------------------

cat("Reading VCF header...\n")

vcf_connection <- gzfile(vcf_file, open = "rt")

header_fields <- NULL

repeat {

  current_line <- readLines(
    vcf_connection,
    n = 1,
    warn = FALSE
  )

  if (length(current_line) == 0) {
    close(vcf_connection)
    stop("Could not find the #CHROM header in the VCF.")
  }

  if (startsWith(current_line, "#CHROM")) {

    header_fields <- strsplit(
      current_line,
      "\t",
      fixed = TRUE
    )[[1]]

    break
  }
}

close(vcf_connection)

if (length(header_fields) < 10) {
  stop("The VCF header contains fewer than 10 columns.")
}

vcf_samples <- trimws(
  header_fields[10:length(header_fields)]
)

if (anyDuplicated(vcf_samples)) {
  stop("Duplicate sample names found in VCF header.")
}

missing_samples <- setdiff(
  metadata$sample,
  vcf_samples
)

if (length(missing_samples) > 0) {
  stop(
    "Metadata samples missing from VCF: ",
    paste(head(missing_samples, 20), collapse = ", ")
  )
}

cat("Samples found in VCF:", length(vcf_samples), "\n")


# ------------------------------------------------------------
# 6. Precompute the VCF sample indices for each population
# ------------------------------------------------------------

population_samples <- lapply(
  population_ids,
  function(p) {
    metadata$sample[metadata$population == p]
  }
)

names(population_samples) <- population_ids

population_indices <- lapply(
  population_samples,
  function(samples) {
    match(samples, vcf_samples)
  }
)

names(population_indices) <- population_ids

for (p in population_ids) {

  if (anyNA(population_indices[[p]])) {
    stop(
      "Could not map every sample to the VCF for population: ",
      p
    )
  }
}


# ------------------------------------------------------------
# 7. Initialize per-population counters
# ------------------------------------------------------------

# For each population and matched SNP:
#   pi     = per-SNP nucleotide diversity
#   called = number of individuals with valid diploid genotypes
#
# NA Pi values are retained when no valid genotype is available.

pi_counts <- lapply(
  population_ids,
  function(p) {
    list(
      pi = numeric(0),
      called = integer(0)
    )
  }
)

names(pi_counts) <- population_ids

total_matched <- 0L
seen_ids <- character(0)


# ------------------------------------------------------------
# 8. Process VCF records in chunks
# ------------------------------------------------------------

cat("Reading compressed VCF records...\n")

chunk_size <- 2000L


# ------------------------------------------------------------
# 8a. Process one chunk and retain common SNPs
# ------------------------------------------------------------

process_chunk <- function(lines) {

  if (length(lines) == 0) {
    return(NULL)
  }

  fields <- strsplit(
    lines,
    "\t",
    fixed = TRUE
  )

  field_lengths <- lengths(fields)

  if (any(field_lengths != length(header_fields))) {

    bad_record <- which(
      field_lengths != length(header_fields)
    )[1]

    stop(
      "Unexpected number of columns in a VCF record. ",
      "Expected ", length(header_fields),
      " but found ", field_lengths[bad_record], "."
    )
  }

  # The SNP identifier is in VCF column 3.
  ids <- vapply(
    fields,
    function(x) x[3],
    character(1)
  )

  keep <- ids %in% common_ids

  if (!any(keep)) {
    return(NULL)
  }

  fields <- fields[keep]
  ids <- ids[keep]

  # Check duplicate IDs both within the current chunk and
  # across chunks.
  duplicated_ids <- duplicated(ids) | ids %in% seen_ids

  if (any(duplicated_ids)) {
    stop(
      "Duplicate common SNP ID found in VCF: ",
      ids[which(duplicated_ids)[1]]
    )
  }

  seen_ids <<- c(seen_ids, ids)

  genotypes <- lapply(
    fields,
    function(x) x[10:length(x)]
  )

  list(
    ids = ids,
    genotypes = genotypes
  )
}


# ------------------------------------------------------------
# 8b. Calculate per-SNP Pi for each population
# ------------------------------------------------------------

process_records <- function(records) {

  if (length(records) == 0) {
    return(invisible(NULL))
  }

  processed <- process_chunk(records)

  if (is.null(processed)) {
    return(invisible(NULL))
  }

  total_matched <<- total_matched + length(processed$ids)

  # Process each common SNP.
  for (i in seq_along(processed$ids)) {

    gt_fields <- processed$genotypes[[i]]

    # Calculate Pi separately for each population.
    for (p in population_ids) {

      sample_indices <- population_indices[[p]]

      gt <- gt_fields[sample_indices]

      # Extract GT before any FORMAT separator.
      # Example: 0/1:25:... becomes 0/1.
      gt <- sub(":.*$", "", gt)

      # Convert phased genotypes such as 0|1 to 0/1.
      gt <- gsub("|", "/", gt, fixed = TRUE)

      # Accept valid diploid, biallelic genotypes only.
      valid <- gt %in% c(
        "0/0",
        "0/1",
        "1/0",
        "1/1"
      )

      gt_valid <- gt[valid]
      n_called <- length(gt_valid)

      pi_value <- NA_real_

      if (n_called > 0) {

        # Count alleles across called diploid individuals.
        n_ref <- (
          2 * sum(gt_valid == "0/0") +
            sum(gt_valid %in% c("0/1", "1/0"))
        )

        n_alt <- (
          2 * sum(gt_valid == "1/1") +
            sum(gt_valid %in% c("0/1", "1/0"))
        )

        n_chr <- n_ref + n_alt

        # Mean fraction of chromosome pairs differing at this SNP.
        if (n_chr > 1) {

          pi_value <- (
            2 * n_ref * n_alt
          ) / (
            n_chr * (n_chr - 1)
          )
        }
      }

      pop_counts <- pi_counts[[p]]

      pop_counts$pi <- c(
        pop_counts$pi,
        pi_value
      )

      pop_counts$called <- c(
        pop_counts$called,
        n_called
      )

      pi_counts[[p]] <<- pop_counts
    }
  }

  invisible(NULL)
}


# ------------------------------------------------------------
# 9. Stream the VCF without loading it entirely into memory
# ------------------------------------------------------------

vcf_connection <- gzfile(vcf_file, open = "rt")

# Skip VCF header lines.
repeat {

  current_line <- readLines(
    vcf_connection,
    n = 1,
    warn = FALSE
  )

  if (length(current_line) == 0) {
    close(vcf_connection)
    stop("No variant records found in VCF.")
  }

  if (!startsWith(current_line, "#")) {
    break
  }
}

# Process the first variant record.
process_records(current_line)

# Process all remaining records in chunks.
records_read <- 1L

repeat {

  records <- readLines(
    vcf_connection,
    n = chunk_size,
    warn = FALSE
  )

  if (length(records) == 0) {
    break
  }

  process_records(records)

  records_read <- records_read + length(records)

  if (records_read %% 100000L < chunk_size) {
    cat(
      "VCF records read:", records_read,
      "| common SNPs matched:", total_matched,
      "\n"
    )
  }
}

close(vcf_connection)

cat("\nVCF processing finished.\n")
cat("VCF records read:", records_read, "\n")
cat("Common SNPs matched in VCF:", total_matched, "\n")

if (total_matched == 0) {
  stop(
    "No common SNP IDs matched VCF column 3. ",
    "Check SNP identifier formats."
  )
}


# ------------------------------------------------------------
# 10. Calculate mean and SD of Pi for each population
# ------------------------------------------------------------

cat("Calculating population-level Pi summaries...\n")

results <- lapply(
  population_ids,
  function(p) {

    pi_values <- pi_counts[[p]]$pi
    called <- pi_counts[[p]]$called

    if (length(pi_values) != length(called)) {
      stop("Internal counter mismatch for population: ", p)
    }

    # Usable sites have a finite per-SNP Pi value.
    usable <- is.finite(pi_values)

    n_used <- sum(usable)

    mean_pi <- if (n_used > 0) {
      mean(pi_values[usable])
    } else {
      NA_real_
    }

    sd_pi <- if (n_used > 1) {
      sd(pi_values[usable])
    } else {
      NA_real_
    }

    data.table(
      Population = p,
      N_individuals = sum(metadata$population == p),
      N_SNPs_matched_VCF = length(pi_values),
      N_SNPs_used = n_used,
      Mean_Pi = mean_pi,
      SD_Pi = sd_pi
    )
  }
)

results <- rbindlist(results)

# ------------------------------------------------------------
# 11. Save results
# ------------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  results,
  file = output_file,
  sep = "\t",
  na = "NA"
)

cat("\nResults:\n")
print(results)

cat("\nSaved to:\n", output_file, "\n")

cat("\nDone.\n")
