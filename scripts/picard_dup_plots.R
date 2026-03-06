#!/usr/bin/env Rscript
# Picard duplicate QC plots: read MarkDuplicates metrics, left-join with alignment summary, save merged RDS and plots.
# Usage: Rscript picard_dup_plots.R <sample_list.txt> <align_summary.rds> <align_dup_summary.rds> <duplication_rate.pdf> <estimated_library_size.pdf> <metrics_paths.txt>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript picard_dup_plots.R <sample_list> <align_summary_rds> <rds_out> <fig_dup_rate> <fig_lib_size> <metrics_paths_file>")
}
sample_list_path <- args[1]
align_summary_rds <- args[2]
rds_path <- args[3]
fig_dup_rate_path <- args[4]
fig_lib_size_path <- args[5]
metrics_manifest_path <- args[6]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

metrics_paths <- readLines(metrics_manifest_path)
metrics_paths <- trimws(metrics_paths)
metrics_paths <- metrics_paths[nzchar(metrics_paths)]
if (length(metrics_paths) == 0) stop("No paths in metrics manifest: ", metrics_manifest_path)

# Optional: order by sample list so output matches alignment_qc sample order
if (file.exists(sample_list_path)) {
  sampleList <- scan(file = sample_list_path, what = "character", quiet = TRUE)
  if (length(sampleList) > 0) {
    ordered <- character(0)
    for (s in sampleList) {
      suffix <- paste0(s, ".picard.rmDup.txt")
      match <- metrics_paths[endsWith(metrics_paths, suffix)]
      if (length(match) > 0) ordered <- c(ordered, match[1])
    }
    if (length(ordered) > 0) metrics_paths <- ordered
  }
}

dir.create(dirname(rds_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(fig_dup_rate_path), showWarnings = FALSE, recursive = TRUE)

# Parse one Picard MarkDuplicates metrics file; return one-row data frame or NULL.
parse_picard_metrics <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- read.table(path, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  if (nrow(x) == 0) return(NULL)
  data.frame(
    READ_PAIRS_EXAMINED = as.numeric(as.character(x$READ_PAIRS_EXAMINED[1])),
    PERCENT_DUPLICATION = as.numeric(as.character(x$PERCENT_DUPLICATION[1])),
    ESTIMATED_LIBRARY_SIZE = as.numeric(as.character(x$ESTIMATED_LIBRARY_SIZE[1])),
    stringsAsFactors = FALSE
  )
}

dupResult <- c()
for (path in metrics_paths) {
  row <- parse_picard_metrics(path)
  if (is.null(row)) next
  # Sample name from filename: <sample>.picard.rmDup.txt
  sample <- sub("\\.picard\\.rmDup\\.txt$", "", basename(path))
  parts <- strsplit(sample, "_")[[1]]
  Genotype <- if (length(parts) >= 1) parts[1] else sample
  Target   <- if (length(parts) >= 2) parts[2] else NA_character_
  Replicate <- if (length(parts) >= 3) parts[3] else NA_character_
  dupResult <- data.frame(
    Sample = sample,
    Genotype = Genotype,
    Target = if (is.na(Target)) sample else Target,
    Replicate = if (is.na(Replicate)) "" else Replicate,
    DuplicationRatePct = row$PERCENT_DUPLICATION * 100,
    READ_PAIRS_EXAMINED = row$READ_PAIRS_EXAMINED,
    ESTIMATED_LIBRARY_SIZE = row$ESTIMATED_LIBRARY_SIZE,
    stringsAsFactors = FALSE
  ) %>% rbind(dupResult, .)
}

if (is.null(dupResult) || nrow(dupResult) == 0) {
  stop("No valid Picard metrics found. Tried: ", paste(metrics_paths, collapse = ", "))
}
dupResult$Genotype <- factor(dupResult$Genotype)

# Left-join with alignment summary and save merged RDS to new file
alignDupSummary <- readRDS(align_summary_rds) %>%
  left_join(dupResult, by = c("Genotype", "Replicate", "Target", "Sample"))
saveRDS(alignDupSummary, rds_path)

theme_qc <- function() {
  theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
      plot.title = element_text(size = 12, hjust = 0.5, margin = margin(0, 0, 10, 0))
    )
}

# Duplication rate (%)
dupFig1a <- dupResult %>%
  ggplot(aes(x = Genotype, y = DuplicationRatePct)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Replicate), position = position_jitter(width = 0.15, height = 0)) +
  ylab("Duplication Rate (%)") +
  xlab("") +
  ggtitle("Picard duplicate rate") +
  theme_qc()
if (length(unique(dupResult$Target)) > 1L) dupFig1a <- dupFig1a + facet_wrap(~ Target)
ggsave(fig_dup_rate_path, dupFig1a, width = 8, height = 6, dpi = 300, bg = "white")

# Estimated library size
libsizeFig1a <- dupResult %>%
  ggplot(aes(x = Genotype, y = ESTIMATED_LIBRARY_SIZE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Replicate), position = position_jitter(width = 0.15, height = 0)) +
  ylab("Estimated Library Size") +
  xlab("") +
  ggtitle("Estimated library size") +
  theme_qc()
if (length(unique(dupResult$Target)) > 1L) libsizeFig1a <- libsizeFig1a + facet_wrap(~ Target, scales = "free_y")
ggsave(fig_lib_size_path, libsizeFig1a, width = 8, height = 6, dpi = 300, bg = "white")

message("Picard duplicate QC done: ", rds_path, " (alignment + duplication summary)")
