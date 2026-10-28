#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# DIVERSITE NUCLEOTIDIQUE INDIVIDUELLE PAR POPULATION
#
# Pour chaque individu :
#   Calcul de sa divergence nucleotidique moyenne par rapport
#   aux autres individus de sa population, sur les SNP communs.
#
# Figures :
#   Boxplot par population, mediane indiquee par la ligne
#   centrale de chaque boite.
#
# Ordre :
#   Native, triees par mediane croissante,
#   puis Invasive, triees par mediane croissante.
#
# Attention :
#   Cette mesure individuelle est une divergence moyenne
#   par rapport aux autres membres de la population.
#   Ce n'est pas le pi classique de population.
# ============================================================


# ------------------------------------------------------------
# 1. Fichiers d'entree
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
# 2. Fichiers de sortie
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
# 3. Verifier les fichiers
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
    "Fichiers introuvables :\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 4. Lire la liste des SNP communs
# ------------------------------------------------------------

common_table <- read_tsv(
  common_snp_file,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

if (ncol(common_table) == 0) {
  stop("Le fichier de SNP communs est vide.")
}

# Chercher une colonne d'identifiants connue.
candidate_cols <- which(
  tolower(names(common_table)) %in%
    c("snp_id", "id", "snp", "locus")
)

if (length(candidate_cols) > 0) {
  id_col <- names(common_table)[candidate_cols[1]]
} else {
  id_col <- names(common_table)[1]
  warning(
    "Aucune colonne SNP_ID/ID reconnue ; utilisation de ",
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
  stop("Aucun identifiant SNP utilisable.")
}

cat("\n========================================\n")
cat("SNP COMMUNS\n")
cat("========================================\n")
cat("Colonne des identifiants :", id_col, "\n")
cat("Nombre de SNP communs :", length(common_ids), "\n")


# ------------------------------------------------------------
# 5. Lire les metadonnees des individus
# ------------------------------------------------------------

metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

if (!all(c("Sample", "Pop") %in% names(metadata))) {
  stop(
    "Le fichier de metadonnees doit contenir Sample et Pop. ",
    "Colonnes trouvees : ",
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
  stop("Des individus apparaissent plusieurs fois dans les metadonnees.")
}


# ------------------------------------------------------------
# 6. Lire l'origine des populations
# ------------------------------------------------------------

population_summary <- read_tsv(
  population_pi_file,
  show_col_types = FALSE
)

if (!all(c("Population", "Origin") %in% names(population_summary))) {
  stop(
    "Le tableau de Pi par population doit contenir ",
    "les colonnes Population et Origin."
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
# 7. Lire l'en-tete du VCF
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
  stop("Impossible de trouver la ligne d'en-tete #CHROM du VCF.")
}

vcf_header <- strsplit(
  vcf_header_line,
  "\t",
  fixed = TRUE
)[[1]]

if (length(vcf_header) < 10) {
  stop("Le VCF ne semble pas contenir de colonnes individuelles.")
}

vcf_samples <- vcf_header[10:length(vcf_header)]
n_samples <- length(vcf_samples)

cat("\n========================================\n")
cat("VCF\n")
cat("========================================\n")
cat("Nombre d'individus dans le VCF :", n_samples, "\n")


# ------------------------------------------------------------
# 8. Associer les echantillons aux populations
# ------------------------------------------------------------

missing_metadata <- setdiff(
  vcf_samples,
  metadata$Sample
)

if (length(missing_metadata) > 0) {
  stop(
    "Individus du VCF absents des metadonnees :\n",
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
    "Origine non trouvee pour les populations : ",
    paste(missing_origins, collapse = ", ")
  )
}

cat(
  "Nombre de populations dans les echantillons : ",
  length(unique(population_for_sample)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Creer la matrice de genotypes
#
# 0/0 = dosage ALT 0
# 0/1 ou 1/0 = dosage ALT 1
# 1/1 = dosage ALT 2
# Autres genotypes = manquants
# ------------------------------------------------------------

dosage_matrix <- matrix(
  NA_real_,
  nrow = n_samples,
  ncol = length(common_ids),
  dimnames = list(vcf_samples, common_ids)
)

matched_ids <- logical(length(common_ids))
chunk_size <- 20000L

cat("\nLecture du VCF par blocs...\n")

con <- gzfile(vcf_file, open = "rt")

# Passer les lignes de metadonnees du VCF jusqu'a #CHROM.
repeat {
  line <- readLines(con, n = 1, warn = FALSE)

  if (length(line) == 0) {
    close(con)
    stop("Fin du VCF avant la ligne #CHROM.")
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

  # En principe, les blocs apres #CHROM contiennent les variants.
  data_lines <- lines[
    !startsWith(lines, "#")
  ]

  if (length(data_lines) == 0) {
    next
  }

  n_variants_read <- n_variants_read + length(data_lines)

  # Extraire rapidement la troisieme colonne (ID).
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

    # Colonnes individuelles du VCF
    sample_fields <- fields[10:(9 + n_samples)]

    # Extraire GT avant les eventuels champs DP, GQ, etc.
    gt <- sub(
      ":.*$",
      "",
      sample_fields
    )

    # Normaliser les separateurs diploides phasing/non-phasing
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
      "Variants lus :", n_variants_read,
      "| SNP communs trouves :", sum(matched_ids),
      "\n"
    )
  }
}

close(con)

cat("\nLecture du VCF terminee.\n")
cat("Variants lus :", n_variants_read, "\n")
cat("SNP communs retrouves :", sum(matched_ids), "\n")
cat("SNP communs demandes :", length(common_ids), "\n")

if (sum(matched_ids) == 0) {
  stop(
    "Aucun SNP commun n'a ete retrouve dans le VCF. ",
    "Verifier les identifiants de SNP."
  )
}

if (sum(!matched_ids) > 0) {
  warning(
    sum(!matched_ids),
    " SNP de la liste n'ont pas ete retrouves dans le VCF."
  )
}

# Garder uniquement les SNP retrouves dans le VCF.
dosage_matrix <- dosage_matrix[
  ,
  matched_ids,
  drop = FALSE
]

cat(
  "SNP effectivement utilises : ",
  ncol(dosage_matrix),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Calculer la divergence individuelle
#
# A chaque SNP :
#   p_i = dosage ALT de l'individu / 2
#
# Pour un individu i et ses pairs j, la divergence attendue
# entre deux alleles tires aleatoirement est :
#
#   d_ij = p_i * (1 - p_j) + (1 - p_i) * p_j
#
# A chaque locus, on calcule la moyenne envers les autres
# individus non manquants de la population, en excluant i.
#
# L'indice individuel final est la moyenne de ces divergences
# par locus. Les loci sans au moins deux genotypes valides
# dans la population ne contribuent pas au calcul.
# ------------------------------------------------------------

individual_total <- numeric(n_samples)
individual_n_loci <- integer(n_samples)

populations <- unique(population_for_sample)
n_loci <- ncol(dosage_matrix)

cat("\n========================================\n")
cat("CALCUL DE LA DIVERGENCE INDIVIDUELLE\n")
cat("========================================\n")

for (pop in populations) {

  sample_indices <- which(
    population_for_sample == pop
  )

  if (length(sample_indices) < 2) {
    warning(
      "Population ", pop,
      " avec moins de deux individus : ignoree."
    )
    next
  }

  cat(
    "Population :", pop,
    "| individus :", length(sample_indices),
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

    # Un individu doit pouvoir etre compare a au moins
    # un autre individu au locus.
    if (n_valid < 2) {
      next
    }

    # Frequence allelique ALT portee par chaque individu
    p <- g[valid_local] / 2

    sum_p <- sum(p)

    # Somme des divergences envers tous les individus valides,
    # y compris soi-meme, puis retrait de la comparaison avec
    # soi-meme et division par le nombre de pairs valides.
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
# 11. Tableau par individu
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
    "Le calcul n'a produit aucune valeur individuelle valide."
  )
}

write_tsv(
  individual_results %>%
    mutate(Origin = as.character(Origin)),
  individual_file
)

cat("\nTableau individuel cree :\n")
cat(" - ", individual_file, "\n", sep = "")

cat(
  "Individus avec une valeur valide : ",
  sum(!is.na(individual_results$Pi_individual)),
  " / ",
  nrow(individual_results),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Calculer la mediane par population
#
# La mediane est calculee entre les indices individuels
# des individus de chaque population.
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

cat("\nMediane par population :\n")
print(population_medians)

cat("\nTableau des medianes cree :\n")
cat(" - ", median_file, "\n", sep = "")


# ------------------------------------------------------------
# 13. Ordonner les populations selon leur mediane individuelle
#
# Native d'abord, puis Invasive.
# L'ordre au sein de chaque groupe est croissant.
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

cat("\nOrdre des populations sur la figure :\n")
print(
  population_medians %>%
    arrange(
      Origin,
      Median_Pi_individual
    )
)


# ------------------------------------------------------------
# 14. Boxplot des valeurs individuelles par population
#
# Une boite par population :
#   - ligne centrale : mediane des indices individuels
#   - boite : intervalle interquartile
#   - moustaches : regle standard du boxplot
#
# Pas de points individuels superposes.
# Pas de symbole de moyenne.
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
# 15. Exporter SVG et PDF
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
# 16. Bilan final
# ------------------------------------------------------------

cat("\n========================================\n")
cat("CALCUL TERMINE\n")
cat("========================================\n")

cat("\nTableau par individu :\n")
cat(" - ", individual_file, "\n", sep = "")

cat("\nMediane par population :\n")
cat(" - ", median_file, "\n", sep = "")

cat("\nBoxplot :\n")
cat(" - ", plot_svg, "\n", sep = "")
cat(" - ", plot_pdf, "\n", sep = "")

cat("\nTermine sans erreur.\n")
