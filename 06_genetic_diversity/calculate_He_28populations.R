
library(data.table)

outdir <- "06_results/heterozygosity_He/common_SNPs"

vcf <- "Stacks/final_vcf_filtering/all_snps/all_snps.final.filtered.vcf.gz"

# Official classification of the 28 populations
pops <- fread(file.path(outdir, "populations_manuscript.tsv"))
meta <- fread(file.path(outdir, "amova_metadata.tsv"))
common <- fread(file.path(outdir, "common_SNP_list.tsv"))

stopifnot(
  all(c("PopulationID", "host_category") %in% names(pops)),
  all(c("Sample", "Pop") %in% names(meta)),
  "ID" %in% names(common)
)

stopifnot(nrow(pops) == 28)
stopifnot(sum(pops$host_category == "Invasive") == 14)
stopifnot(sum(pops$host_category == "native") == 14)

# Assign each individual to its manuscript population
meta <- merge(
  meta[, .(Sample, Pop)],
  pops[, .(PopulationID, host_category)],
  by.x = "Pop", by.y = "PopulationID"
)

if (anyDuplicated(meta$Sample))
  stop("A sample appears more than once in the metadata.")

common_ids <- unique(common$ID)
n_pop <- nrow(pops)

# Accumulators for population-level diversity across SNPs
sum_He <- numeric(n_pop)
n_loci <- integer(n_pop)

con <- gzfile(vcf, open = "rt")
header <- NULL

repeat {
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  if (startsWith(line, "#CHROM")) {
    header <- strsplit(line, "\t", fixed = TRUE)[[1]]
    break
  }
}

if (is.null(header)) stop("VCF header not found.")

samples <- header[10:length(header)]
idx <- match(samples, meta$Sample)

if (anyNA(idx)) {
  close(con)
  stop("Some VCF samples are missing from the metadata.")
}

pop_index <- match(meta$Pop[idx], pops$PopulationID)

if (anyNA(pop_index)) {
  close(con)
  stop("A population is missing from the manuscript population table.")
}

gt_index <- integer(0)
n_matched <- 0L

repeat {
  lines <- readLines(con, n = 1000)
  if (length(lines) == 0) break

  for (ln in lines) {
    if (startsWith(ln, "#")) next
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(f) < 10 || !(f[3] %in% common_ids)) next

    fmt <- strsplit(f[9], ":", fixed = TRUE)[[1]]
    gi <- match("GT", fmt)
    if (is.na(gi)) next

    n_matched <- n_matched + 1L

    gt <- vapply(f[10:length(f)], function(x) {
      z <- strsplit(x, ":", fixed = TRUE)[[1]]
      if (length(z) < gi) "." else z[gi]
    }, character(1))

    for (p in seq_len(n_pop)) {
      inds <- which(pop_index == p)
      g <- gt[inds]
      valid <- grepl("^(0[/|]0|0[/|]1|1[/|]0|1[/|]1)$", g)

      if (!any(valid)) next

      g <- g[valid]
      alleles <- unlist(strsplit(g, "[/|]"))
      alleles <- as.integer(alleles)
      q <- mean(alleles)

      sum_He[p] <- sum_He[p] + 2 * q * (1 - q)
      n_loci[p] <- n_loci[p] + 1L
    }
  }
}

close(con)

if (n_matched == 0)
  stop("No common SNPs matched the VCF identifiers.")

result <- copy(pops)

result[, `:=`(
  He = sum_He / n_loci,
  N_SNP_used = n_loci
)]

if (any(result$N_SNP_used == 0))
  stop("At least one population has no usable SNPs.")

setorder(result, host_category, He)

fwrite(
  result,
  file.path(outdir, "He_28populations_manuscript.tsv"),
  sep = "\t"
)

# Plot one He value per population, grouped by host category
pdf(file.path(outdir, "He_28populations_by_host.pdf"),
    width = 12, height = 7)

par(mar = c(10, 5, 3, 1))

cols <- ifelse(result$host_category == "Invasive",
               "tomato", "steelblue")

plot(seq_len(nrow(result)), result$He,
     pch = 19, col = cols,
     xaxt = "n", xlab = "",
     ylab = "Expected heterozygosity (He)",
     main = "Expected heterozygosity by population and host category")

axis(1, at = seq_len(nrow(result)),
     labels = result$PopulationID, las = 2, cex.axis = 0.8)

abline(v = sum(result$host_category == "Invasive") + 0.5,
       lty = 2, col = "grey50")

legend("topleft",
       legend = c("Invasive", "native"),
       col = c("tomato", "steelblue"),
       pch = 19, bty = "n")

dev.off()

# Independent comparison between the two host categories
cat("Number of common SNPs used:", n_matched, "\n")
cat("Number of populations:", nrow(result), "\n")
print(result)

cat("\nUnpaired Wilcoxon test:\n")

print(wilcox.test(
  He ~ host_category,
  data = result,
  exact = FALSE
))

