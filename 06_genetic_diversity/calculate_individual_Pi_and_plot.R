
suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# INDIVIDUAL NUCLEOTIDE DIVERSITY BY POPULATION
#
# For each individual:
#   Calculate its mean nucleotide divergence from the other
#   individuals in its population, using the common SNPs.
#
# Figures:
#   Boxplot by population, with the median shown by the
#   central line of each box.
#
# Ordering:
#   Native populations, sorted by increasing median,
#   followed by Invasive populations, sorted by increasing median.
#
# Note:
#   This individual-level measure represents mean divergence
#   from other members of the same population.
#   It is not the conventional population-level pi.
# ============================================================


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

base_dir <- "06_results/nucleotide_diversity_Pi/common_SNPs"

vcf_file <- "Stacks/final_vcf_filtering/all_snps/all_snps.final.filtered.vcf.gz"

common_snp_file <- file.path(
  "06_results/heterozygosity_He/common_SNPs",
  "common_SNP_list.tsv"
)

metadata_file <- file.path(
  "06_results/heterozygosity_He/common_SNPs",
  "amova_metadata.tsv"
)

population_pi_file <- file.path(
  base_dir,
  "Pi_bootstrap_mean_28populations.tsv"
)


# ------------------------------------------------------------
# 2. Output files
# ------------------------------------------------------------

individual_file <- file.path(
  base_dir,
  "Pi_individual_by_sample.tsv"
)

median_file <- file.path(
  base_dir,
  "Pi_individual_median_by_population.tsv"
)

plot_svg <- file.path(
  base_dir,
  "Pi_individual_28populations_boxplot_median.svg"
)

plot_pdf <- file.path(
  base_dir,
  "Pi_individual_28populations_boxplot_median.pdf"
)


# ------------------------------------------------------------
# 3. Check input files
# ------------------------------------------------------------

input_files <- c(
  vcf_file,
  common_snp_file,
  metadata_file,
  population_pi_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Input files not found:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 4. Read the list of common SNPs
# ------------------------------------------------------------

common_table <- read_tsv(
  common_snp_file,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

if (ncol(common_table) == 0) {
  stop("The common SNP file is empty.")
}

# Look for a recognized SNP identifier column.
candidate_cols <- which(
  tolower(names(common_table)) %in%
    c("snp_id", "id", "snp", "locus")
)

if (length(candidate_cols) > 0) {
  id_col <- names(common_table)[candidate_cols[1]]
} else {
  id_col <- names(common_table)[1]
  warning(
    "No recognized SNP_ID/ID column found; using ",
    id_col
  )
}

common_ids <- unique(
  trimws(as.character(common_table[[id_col]]))
)

common_ids <- common_ids[
  !is.na(common_ids) &
    common_ids != "" &
    common_ids != "."
]

if (length(common_ids) == 0) {
  stop("No usable SNP identifiers were found.")
}

cat("\n========================================\n")
cat("COMMON SNPs\n")
cat("========================================\n")
cat("Identifier column:", id_col, "\n")
cat("Number of common SNPs:", length(common_ids), "\n")


# ------------------------------------------------------------
# 5. Read individual metadata
# ------------------------------------------------------------

metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

if (!all(c("Sample", "Pop") %in% names(metadata))) {
  stop(
    "The metadata file must contain Sample and Pop. ",
    "Columns found: ",
    paste(names(metadata), collapse = ", ")
  )
}

metadata <- metadata %>%
  transmute(
    Sample = as.character(Sample),
    Population = as.character(Pop)
  ) %>%
  filter(
    !is.na(Sample),
    !is.na(Population),
    Sample != "",
    Population != ""
  )

if (anyDuplicated(metadata$Sample)) {
  stop("Some individuals appear more than once in the metadata.")
}


# ------------------------------------------------------------
# 6. Read population origins
# ------------------------------------------------------------

population_summary <- read_tsv(
  population_pi_file,
  show_col_types = FALSE
)

if (!all(c("Population", "Origin") %in% names(population_summary))) {
  stop(
    "The population-level Pi table must contain ",
    "the Population and Origin columns."
  )
}

population_summary <- population_summary %>%
  transmute(
    Population = as.character(Population),
    Origin = case_when(
      str_to_lower(str_trim(as.character(Origin))) == "native" ~
        "Native",
      str_to_lower(str_trim(as.character(Origin))) == "invasive" ~
        "Invasive",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Origin)) %>%
  distinct(Population, .keep_all = TRUE)

origin_lookup <- setNames(
  population_summary$Origin,
  population_summary$Population
)


# ------------------------------------------------------------
# 7. Read the VCF header
# ------------------------------------------------------------

con <- gzfile(vcf_file, open = "rt")

vcf_header_line <- NULL

repeat {
  line <- readLines(con, n = 1, warn = FALSE)

  if (length(line) == 0) {
    break
  }

  if (startsWith(line, "#CHROM")) {
    vcf_header_line <- line
    break
  }
}

close(con)

if (is.null(vcf_header_line)) {
  stop("Could not find the VCF #CHROM header line.")
}

vcf_header <- strsplit(
  vcf_header_line,
  "\t",
  fixed = TRUE
)[[1]]

if (length(vcf_header) < 10) {
  stop("The VCF does not appear to contain individual sample columns.")
}

vcf_samples <- vcf_header[10:length(vcf_header)]
n_samples <- length(vcf_samples)

cat("\n========================================\n")
cat("VCF\n")
cat("========================================\n")
cat("Number of individuals in the VCF:", n_samples, "\n")


# ------------------------------------------------------------
# 8. Match samples to populations
# ------------------------------------------------------------

missing_metadata <- setdiff(
  vcf_samples,
  metadata$Sample
)

if (length(missing_metadata) > 0) {
  stop(
    "VCF individuals missing from the metadata:\n",
    paste(missing_metadata, collapse = ", ")
  )
}

metadata_ordered <- metadata[
  match(vcf_samples, metadata$Sample),
]

population_for_sample <- metadata_ordered$Population

origin_for_sample <- unname(
  origin_lookup[population_for_sample]
)

if (anyNA(origin_for_sample)) {
  missing_origins <- unique(
    population_for_sample[is.na(origin_for_sample)]
  )

  stop(
    "Origin not found for the following populations: ",
    paste(missing_origins, collapse = ", ")
  )
}

cat(
  "Number of populations represented in the samples: ",
  length(unique(population_for_sample)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Build the genotype dosage matrix
#
# 0/0 = ALT dosage 0
# 0/1 or 1/0 = ALT dosage 1
# 1/1 = ALT dosage 2
# Other genotypes = missing
# ------------------------------------------------------------

dosage_matrix <- matrix(
  NA_real_,
  nrow = n_samples,
  ncol = length(common_ids),
  dimnames = list(vcf_samples, common_ids)
)

matched_ids <- logical(length(common_ids))
chunk_size <- 20000L

cat("\nReading the VCF in chunks...\n")

con <- gzfile(vcf_file, open = "rt")

# Skip VCF metadata lines until the #CHROM header.
repeat {
  line <- readLines(con, n = 1, warn = FALSE)

  if (length(line) == 0) {
    close(con)
    stop("End of VCF reached before the #CHROM header.")
  }

  if (startsWith(line, "#CHROM")) {
    break
  }
}

n_variants_read <- 0L
n_variants_matched <- 0L

repeat {

  lines <- readLines(
    con,
    n = chunk_size,
    warn = FALSE
  )

  if (length(lines) == 0) {
    break
  }

  # Lines after #CHROM should normally contain variants.
  data_lines <- lines[
    !startsWith(lines, "#")
  ]

  if (length(data_lines) == 0) {
    next
  }

  n_variants_read <- n_variants_read + length(data_lines)

  # Quickly extract the third column (variant ID).
  ids <- sub(
    "^(?:[^\\t]*\\t){2}([^\\t]*).*$",
    "\\1",
    data_lines,
    perl = TRUE
  )

  hit_positions <- which(ids %in% common_ids)

  if (length(hit_positions) == 0) {
    next
  }

  for (k in hit_positions) {

    fields <- strsplit(
      data_lines[k],
      "\t",
      fixed = TRUE
    )[[1]]

    if (length(fields) < 9 + n_samples) {
      next
    }

    snp_id <- fields[3]
    snp_index <- match(snp_id, common_ids)

    if (is.na(snp_index)) {
      next
    }

    # Individual sample columns in the VCF.
    sample_fields <- fields[10:(9 + n_samples)]

    # Extract GT before any additional fields such as DP or GQ.
    gt <- sub(
      ":.*$",
      "",
      sample_fields
    )

    # Normalize phased and unphased diploid genotype separators.
    gt <- gsub(
      "\\|",
      "/",
      gt
    )

    d <- rep(NA_real_, n_samples)

    d[gt == "0/0"] <- 0
    d[gt %in% c("0/1", "1/0")] <- 1
    d[gt == "1/1"] <- 2

    dosage_matrix[, snp_index] <- d
    matched_ids[snp_index] <- TRUE
    n_variants_matched <- n_variants_matched + 1L
  }

  if (n_variants_read %% 200000L < chunk_size) {
    cat(
      "Variants read:", n_variants_read,
      "| Common SNPs found:", sum(matched_ids),
      "\n"
    )
  }
}

close(con)

cat("\nVCF reading completed.\n")
cat("Variants read:", n_variants_read, "\n")
cat("Common SNPs found:", sum(matched_ids), "\n")
cat("Common SNPs requested:", length(common_ids), "\n")

if (sum(matched_ids) == 0) {
  stop(
    "No common SNPs were found in the VCF. ",
    "Please check the SNP identifiers."
  )
}

if (sum(!matched_ids) > 0) {
  warning(
    sum(!matched_ids),
    " SNPs from the list were not found in the VCF."
  )
}

# Keep only SNPs found in the VCF.
dosage_matrix <- dosage_matrix[
  ,
  matched_ids,
  drop = FALSE
]

cat(
  "SNPs actually used: ",
  ncol(dosage_matrix),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Calculate individual nucleotide divergence
#
# At each SNP:
#   p_i = ALT dosage of individual i / 2
#
# For individual i and its peers j, the expected divergence
# between two randomly sampled alleles is:
#
#   d_ij = p_i * (1 - p_j) + (1 - p_i) * p_j
#
# At each locus, calculate the mean divergence from the other
# non-missing individuals in the same population, excluding i.
#
# The final individual index is the mean of these divergences
# across loci. Loci with fewer than two valid genotypes in the
# population are excluded from the calculation.
# ------------------------------------------------------------

individual_total <- numeric(n_samples)
individual_n_loci <- integer(n_samples)

populations <- unique(population_for_sample)
n_loci <- ncol(dosage_matrix)

cat("\n========================================\n")
cat("INDIVIDUAL DIVERGENCE CALCULATION\n")
cat("========================================\n")

for (pop in populations) {

  sample_indices <- which(
    population_for_sample == pop
  )

  if (length(sample_indices) < 2) {
    warning(
      "Population ", pop,
      " has fewer than two individuals and will be skipped."
    )
    next
  }

  cat(
    "Population:", pop,
    "| Individuals:", length(sample_indices),
    "\n"
  )

  pop_matrix <- dosage_matrix[
    sample_indices,
    ,
    drop = FALSE
  ]

  for (locus in seq_len(n_loci)) {

    g <- pop_matrix[, locus]

    valid_local <- which(!is.na(g))
    n_valid <- length(valid_local)

    # Each individual must be comparable to at least
    # one other individual at this locus.
    if (n_valid < 2) {
      next
    }

    # ALT allele frequency carried by each individual.
    p <- g[valid_local] / 2

    sum_p <- sum(p)

    # Sum divergences against all valid individuals,
    # including self, then remove the self-comparison
    # and divide by the number of other valid individuals.
    divergence_sum <- (
      p * (n_valid - sum_p) +
        (1 - p) * sum_p -
        2 * p * (1 - p)
    )

    divergence_mean <- divergence_sum / (n_valid - 1)

    global_indices <- sample_indices[valid_local]

    individual_total[global_indices] <-
      individual_total[global_indices] + divergence_mean

    individual_n_loci[global_indices] <-
      individual_n_loci[global_indices] + 1L
  }
}


# ------------------------------------------------------------
# 11. Create the individual-level table
# ------------------------------------------------------------

individual_pi <- ifelse(
  individual_n_loci > 0,
  individual_total / individual_n_loci,
  NA_real_
)

individual_results <- tibble(
  Sample = vcf_samples,
  Population = population_for_sample,
  Origin = origin_for_sample,
  Pi_individual = individual_pi,
  N_loci_used = individual_n_loci
) %>%
  mutate(
    Origin = factor(
      Origin,
      levels = c("Native", "Invasive")
    )
  )

if (all(is.na(individual_results$Pi_individual))) {
  stop(
    "The calculation did not produce any valid individual values."
  )
}

write_tsv(
  individual_results %>%
    mutate(Origin = as.character(Origin)),
  individual_file
)

cat("\nIndividual-level table created:\n")
cat(" - ", individual_file, "\n", sep = "")

cat(
  "Individuals with a valid value: ",
  sum(!is.na(individual_results$Pi_individual)),
  " / ",
  nrow(individual_results),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Calculate the median for each population
#
# The median is calculated across individual indices
# for all individuals within each population.
# ------------------------------------------------------------

population_medians <- individual_results %>%
  group_by(
    Population,
    Origin
  ) %>%
  summarise(
    Median_Pi_individual = median(
      Pi_individual,
      na.rm = TRUE
    ),
    N_individuals = sum(
      !is.na(Pi_individual)
    ),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(Median_Pi_individual)
  ) %>%
  mutate(
    Origin = factor(
      Origin,
      levels = c("Native", "Invasive")
    )
  ) %>%
  arrange(
    Origin,
    Median_Pi_individual
  )

write_tsv(
  population_medians %>%
    mutate(Origin = as.character(Origin)),
  median_file
)

cat("\nPopulation medians:\n")
print(population_medians)

cat("\nPopulation median table created:\n")
cat(" - ", median_file, "\n", sep = "")


# ------------------------------------------------------------
# 13. Order populations by individual-level median
#
# Native populations first, followed by Invasive populations.
# Within each group, populations are sorted in ascending order.
# ------------------------------------------------------------

population_order <- population_medians %>%
  arrange(
    Origin,
    Median_Pi_individual
  ) %>%
  pull(Population)

plot_data <- individual_results %>%
  filter(
    !is.na(Pi_individual)
  ) %>%
  mutate(
    Population = factor(
      Population,
      levels = population_order
    ),
    Origin = factor(
      Origin,
      levels = c("Native", "Invasive")
    )
  )

n_native <- sum(
  population_medians$Origin == "Native"
)

cat("\nPopulation order in the figure:\n")
print(
  population_medians %>%
    arrange(
      Origin,
      Median_Pi_individual
    )
)


# ------------------------------------------------------------
# 14. Boxplot of individual values by population
#
# One box per population:
#   - central line: median of individual indices
#   - box: interquartile range
#   - whiskers: standard boxplot rule
#
# Individual points are not overlaid.
# No mean symbol is displayed.
# ------------------------------------------------------------

host_colors <- c(
  "Native" = "#BDBDBD",
  "Invasive" = "#E69F00"
)

p_individual <- ggplot(
  plot_data,
  aes(
    x = Population,
    y = Pi_individual,
    fill = Origin
  )
) +
  geom_boxplot(
    width = 0.72,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.4,
    fatten = 1.3
  ) +
  geom_vline(
    xintercept = n_native + 0.5,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = host_colors,
    breaks = c("Native", "Invasive"),
    name = "Origin"
  ) +
  scale_x_discrete(
    limits = population_order,
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    name = "Individual mean pairwise nucleotide divergence",
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    title = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),

    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),

    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      color = "black",
      face = "plain"
    ),

    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),

    axis.title.y = element_text(
      size = 11,
      color = "black"
    ),

    axis.ticks = element_line(
      color = "black"
    ),

    legend.position = c(0.03, 0.97),

    legend.justification = c(0, 1),

    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.3
    ),

    legend.title = element_text(
      size = 9
    ),

    legend.text = element_text(
      size = 9
    )
  )


# ------------------------------------------------------------
# 15. Export SVG and PDF figures
# ------------------------------------------------------------

ggsave(
  filename = plot_svg,
  plot = p_individual,
  width = 13,
  height = 6.5,
  units = "in"
)

ggsave(
  filename = plot_pdf,
  plot = p_individual,
  width = 13,
  height = 6.5,
  units = "in"
)


# ------------------------------------------------------------
# 16. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("CALCULATION COMPLETED\n")
cat("========================================\n")

cat("\nIndividual-level table:\n")
cat(" - ", individual_file, "\n", sep = "")

cat("\nPopulation median table:\n")
cat(" - ", median_file, "\n", sep = "")

cat("\nBoxplots:\n")
cat(" - ", plot_svg, "\n", sep = "")
cat(" - ", plot_pdf, "\n", sep = "")

cat("\nFinished without errors.\n")

