#!/usr/bin/env Rscript
# Bin500 correlation: read 500 bp binned fragment BED files, merge by chrom/bin, compute Pearson
# correlation on log2(counts), and save a correlation heatmap + dendrogram (matches rmarkdowns/3-CUTandTag-Make_the_bed.Rmd).
# Usage: Rscript bin500_correlation.R <output.pdf> <bed1.bed> [bed2.bed ...]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript bin500_correlation.R <output.pdf> <bed1.bed> [bed2.bed ...]")
}
output_pdf <- args[1]
bed_paths <- args[-1]
bed_paths <- bed_paths[file.exists(bed_paths)]
if (length(bed_paths) == 0) stop("No existing BED files provided.")
if (length(bed_paths) < 2) {
  stop("Need at least 2 bin500 BED files for correlation analysis.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(corrplot)
})

# Extract sample name from path: .../SAMPLE.bowtie2.fragmentsCount.bin500.bed -> SAMPLE
sample_from_path <- function(p) {
  base <- basename(p)
  sub("\\.bowtie2\\.fragmentsCount\\.bin500\\.bed$", "", base)
}

dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)

frag_count <- NULL
for (i in seq_along(bed_paths)) {
  path <- bed_paths[i]
  sample_name <- sample_from_path(path)
  d <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  colnames(d) <- c("chrom", "bin", sample_name)
  if (is.null(frag_count)) {
    frag_count <- d
  } else {
    frag_count <- dplyr::full_join(frag_count, d, by = c("chrom", "bin"))
  }
}

# Count columns (all but chrom, bin)
count_cols <- seq(3, ncol(frag_count))
M <- frag_count[, count_cols] %>%
  as.matrix() %>%
  {. + 1} %>%
  log2() %>%
  cor(use = "pairwise.complete.obs")

# Hierarchical clustering (same as Rmd: 1 - cor, complete linkage); used for dendrogram panel.
hc <- hclust(as.dist(1 - M), method = "complete")

# Viridis-like palette (from Rmd)
pal <- colorRampPalette(rev(c(
  "#440154", "#482878", "#3E4989", "#31688E", "#26828E",
  "#1F9E89", "#35B779", "#6DCD59", "#B4DD2C", "#FDE725"
)))(100)

# Two-panel figure: (1) correlation heatmap reordered by clustering; (2) dendrogram of the same clustering.
# Heatmap shows pairwise correlations; dendrogram shows sample similarity as a tree. Order matches between panels.
pdf(output_pdf, width = 10, height = 12)
par(mfrow = c(2, 1), mar = c(5, 4, 4, 2))
# Color scale must cover the matrix range (corrplot errors otherwise); use data range so plot works 
col_lim <- range(M, na.rm = TRUE)
corrplot(M, method = "color",
         order = "hclust",
         hclust.method = "complete",
         type = "full",
         tl.col = "black",
         tl.cex = 0.6,
         tl.pos = "l",
         number.digits = 2,
         number.cex = 1,
         is.corr = FALSE,
         title = "CUT&Tag 500 bp bin pairwise correlations",
         col = pal,
         col.lim = col_lim,
         mar = c(0, 0, 2, 0))
plot(hc, main = "Hierarchical clustering dendrogram",
     xlab = "", ylab = "Distance", sub = "")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
dev.off()
