#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Observed heterozygosity (Ho) for 28 populations
# No bcftools required
# ============================================================

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
  "06_results/homozygosity/common_SNPs/Ho_28populations.tsv"
)

# ------------------------------------------------------------
# 1. Check input files
# ------------------------------------------------------------

for (f in c(vcf_file, snp_file, metadata_file)) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
}

# ------------------------------------------------------------
# 2. Read the common SNP list
# ------------------------------------------------------------

cat("Reading common SNP list...\n")

common_snps <- fread(snp_file)

if (!"ID" %in% names(common_snps)) {
  stop("The common SNP list must contain a column named ID.")
}

common_ids <- unique(trimws(as.character(common_snps$ID)))
common_ids <- common_ids[!is.na(common_ids) & common_ids != ""]

if (length(common_ids) == 0) {
  stop("The common SNP list is empty.")
}

cat("SNPs in common list:", length(common_ids), "\n")

# ------------------------------------------------------------
# 3. Read metadata
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
# 4. Read the VCF header and identify samples
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

vcf_samples <- trimws(header_fields[10:length(header_fields)])

if (anyDuplicated(vcf_samples)) {
  stop("Duplicate sample names found in VCF header.")
}

missing_samples <- setdiff(metadata$sample, vcf_samples)

if (length(missing_samples) > 0) {
  stop(
    "Metadata samples missing from VCF: ",
    paste(head(missing_samples, 20), collapse = ", ")
  )
}

cat("Samples found in VCF:", length(vcf_samples), "\n")

# Precompute the VCF genotype-column indices for each population.
# This avoids repeatedly matching sample names for every SNP.

population_samples <- lapply(
  population_ids,
  function(p) metadata$sample[metadata$population == p]
)

names(population_samples) <- population_ids

population_indices <- lapply(
  population_samples,
  function(samples) match(samples, vcf_samples)
)

names(population_indices) <- population_ids

# Check that the sample indices are valid.
for (p in population_ids) {
  if (anyNA(population_indices[[p]])) {
    stop("Could not map every sample to the VCF for population: ", p)
  }
}

# ------------------------------------------------------------
# 5. Initialize counters
# ------------------------------------------------------------

# For each population and SNP:
# het    = number of heterozygous individuals
# called = number of individuals with a valid diploid genotype

ho_counts <- lapply(
  population_ids,
  function(p) {
    list(
      het = numeric(0),
      called = numeric(0)
    )
  }
)

names(ho_counts) <- population_ids

total_matched <- 0L
seen_ids <- character(0)

# ------------------------------------------------------------
# 6. Process VCF records
# ------------------------------------------------------------

cat("Reading compressed VCF records...\n")

chunk_size <- 2000L

process_chunk <- function(lines) {

  if (length(lines) == 0) {
    return(NULL)
  }

  # Split each VCF record into tab-separated fields.
  fields <- strsplit(lines, "\t", fixed = TRUE)

  # All VCF records should have the same number of columns as
  # the header.
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

  # SNP identifier is column 3.
  ids <- vapply(
    fields,
    function(x) x[3],
    character(1)
  )

  # Keep only SNPs in the common SNP list.
  keep <- ids %in% common_ids

  if (!any(keep)) {
    return(NULL)
  }

  fields <- fields[keep]
  ids <- ids[keep]

  # Ensure that each common SNP ID occurs only once in the VCF.
  duplicated_ids <- ids %in% seen_ids

  if (any(duplicated_ids)) {
    stop(
      "Duplicate common SNP ID found in VCF: ",
      ids[which(duplicated_ids)[1]]
    )
  }

  # Save matched IDs in the parent environment.
  seen_ids <<- c(seen_ids, ids)

  # Genotype/sample fields begin at column 10.
  genotypes <- lapply(
    fields,
    function(x) x[10:length(x)]
  )

  list(
    ids = ids,
    genotypes = genotypes
  )
}

process_records <- function(records) {

  if (length(records) == 0) {
    return(invisible(NULL))
  }

  processed <- process_chunk(records)

  if (is.null(processed)) {
    return(invisible(NULL))
  }

  total_matched <<- total_matched + length(processed$ids)

  # Loop over each matched SNP.
  for (i in seq_along(processed$ids)) {

    gt_fields <- processed$genotypes[[i]]

    # Calculate Ho separately for each population.
    for (p in population_ids) {

      sample_indices <- population_indices[[p]]

      # Select genotype fields corresponding to the population.
      gt <- gt_fields[sample_indices]

      # Extract GT, which appears before the first colon.
      # Examples: 0/1:25:... becomes 0/1
      gt <- sub(":.*$", "", gt)

      # Convert phased genotypes, e.g. 0|1, to 0/1.
      gt <- gsub("|", "/", gt, fixed = TRUE)

      # Keep valid diploid genotypes only.
      valid <- gt %in% c(
        "0/0",
        "0/1",
        "1/0",
        "1/1"
      )

      # Count heterozygotes.
      het <- gt %in% c("0/1", "1/0")

      n_het <- sum(het & valid, na.rm = TRUE)
      n_called <- sum(valid, na.rm = TRUE)

      # IMPORTANT:
      # Explicitly update the parent environment so the counters
      # are retained after this function finishes.
      pop_counts <- ho_counts[[p]]

      pop_counts$het <- c(
        pop_counts$het,
        n_het
      )

      pop_counts$called <- c(
        pop_counts$called,
        n_called
      )

      ho_counts[[p]] <<- pop_counts
    }
  }

  invisible(NULL)
}

# ------------------------------------------------------------
# 7. Stream the VCF without loading it entirely into memory
# ------------------------------------------------------------

vcf_connection <- gzfile(vcf_file, open = "rt")

# Skip header lines and locate the first variant record.
first_record <- NULL

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
    first_record <- current_line
    break
  }
}

# Process the first record.
process_records(first_record)

# Process the remaining records in chunks.
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
}

close(vcf_connection)

cat("Common SNPs matched in VCF:", total_matched, "\n")

if (total_matched == 0) {
  stop(
    "No common SNP IDs matched VCF column 3. ",
    "Check SNP identifier formats."
  )
}

# ------------------------------------------------------------
# 8. Calculate mean per-SNP observed heterozygosity
# ------------------------------------------------------------

cat("Calculating Ho for each population...\n")

results <- lapply(
  population_ids,
  function(p) {

    het <- ho_counts[[p]]$het
    called <- ho_counts[[p]]$called

    if (length(het) != length(called)) {
      stop("Internal counter mismatch for population: ", p)
    }

    # Per-SNP Ho = heterozygous individuals / called individuals.
    ho_per_snp <- rep(NA_real_, length(called))

    has_called <- called > 0

    ho_per_snp[has_called] <-
      het[has_called] / called[has_called]

    n_used <- sum(has_called)

    mean_ho <- if (n_used > 0) {
      mean(ho_per_snp[has_called])
    } else {
      NA_real_
    }

    sd_ho <- if (n_used > 1) {
      sd(ho_per_snp[has_called])
    } else {
      NA_real_
    }

    data.table(
      Population = p,
      N_individuals = sum(metadata$population == p),
      N_SNPs_matched_VCF = length(called),
      N_SNPs_used = n_used,
      Mean_Ho = mean_ho,
      SD_Ho = sd_ho
    )
  }
)

results <- rbindlist(results)

# ------------------------------------------------------------
# 9. Save results
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
