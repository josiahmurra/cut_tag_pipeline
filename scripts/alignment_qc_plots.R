#!/usr/bin/env Rscript
# Alignment QC plots: read bowtie2 summary files and plot alignment rate, mapped fragments, sequencing depth.
# Usage: Rscript alignment_qc_plots.R <sample_list.txt> <alignSummary.rds> <alignment_rate.pdf> <mapped_fragments.pdf> <sequencing_depth.pdf> <summary_paths.txt>
#   summary_paths.txt: one bowtie2 summary file path per line (avoids ARG_MAX and spaces-in-path issues).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript alignment_qc_plots.R <sample_list> <rds_out> <fig_align_rate> <fig_mapped> <fig_depth> <summary_paths_file>")
}
sample_list_path <- args[1]
rds_path <- args[2]
fig_align_rate_path <- args[3]
fig_mapped_path <- args[4]
fig_depth_path <- args[5]
summary_manifest_path <- args[6]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

summary_paths <- readLines(summary_manifest_path)
summary_paths <- trimws(summary_paths)
summary_paths <- summary_paths[nzchar(summary_paths)]
if (length(summary_paths) == 0) stop("No paths in summary manifest: ", summary_manifest_path)

dir.create(dirname(rds_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(fig_align_rate_path), showWarnings = FALSE, recursive = TRUE)

sampleList <- scan(file = sample_list_path, what = "character", quiet = TRUE)
if (length(sampleList) == 0) stop("Sample list is empty.")

# Parse one bowtie2 summary file; return a one-row data frame or NULL on error
parse_bowtie2_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- read.table(path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  # Total reads: first line usually "N reads; of these:"
  total_reads <- NA_real_
  if (nrow(x) >= 1 && length(x$V1) >= 1) {
    total_reads <- as.numeric(as.character(x$V1[1]))
  }
  if (is.na(total_reads)) return(NULL)
  # Overall alignment rate: line containing "overall alignment rate"
  align_rate <- NA_real_
  idx <- which(grepl("overall alignment rate", apply(x, 1L, paste, collapse = " ")))
  if (length(idx) > 0L) {
    pct <- sub("%.*", "", x$V1[idx[1L]]) # Strip "%" so we can convert to numeric 
    align_rate <- as.numeric(as.character(pct))
  }
  # Mapped fragment count: bowtie2 reports concordant + discordant; use V1[4]+V1[5] if present (common layout)
  mapped <- NA_real_
  if (nrow(x) >= 5 && all(!is.na(suppressWarnings(as.numeric(x$V1[4:5]))))) {
    mapped <- as.numeric(as.character(x$V1[4])) + as.numeric(as.character(x$V1[5]))
  }
  if (is.na(mapped) && !is.na(align_rate)) {
    mapped <- round(total_reads * (align_rate / 100))
  }
  data.frame(
    SequencingDepth = total_reads,
    MappedFragNum = mapped,
    AlignmentRate = align_rate,
    stringsAsFactors = FALSE
  )
}

alignResult <- c()
for (path in summary_paths) {
  row <- parse_bowtie2_summary(path)
  if (is.null(row)) next
  # Sample name from filename: <sample>.bowtie2.txt
  sample <- sub("\\.bowtie2\\.txt$", "", basename(path))
  parts <- strsplit(sample, "_")[[1]]
  Genotype <- if (length(parts) >= 1) parts[1] else sample
  Target   <- if (length(parts) >= 2) parts[2] else NA_character_
  Replicate <- if (length(parts) >= 3) parts[3] else NA_character_
  alignResult <- data.frame(
    Sample = sample,
    Genotype = Genotype,
    Target = if (is.na(Target)) sample else Target,
    Replicate = if (is.na(Replicate)) "" else Replicate,
    SequencingDepth = row$SequencingDepth,
    MappedFragNum = row$MappedFragNum,
    AlignmentRate = row$AlignmentRate,
    stringsAsFactors = FALSE
  ) %>% rbind(alignResult, .)
}

if (is.null(alignResult) || nrow(alignResult) == 0) {
  stop("No valid bowtie2 summary files found. Tried: ", paste(summary_paths, collapse = ", "),
       ". If paths are relative, R may be running in a different working directory; use absolute paths in the manifest.")
}
alignResult$Genotype <- factor(alignResult$Genotype)

saveRDS(alignResult, rds_path)

# Plotting helper
theme_qc <- function() {
  theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
      plot.title = element_text(size = 12, hjust = 0.5, margin = margin(0, 0, 10, 0))
    )
}

# Alignment rate
alignFig1a <- alignResult %>%
  ggplot(aes(x = Genotype, y = AlignmentRate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Replicate), position = position_jitter(width = 0.15, height = 0)) +
  ylab("Percentage of Reads Aligned") +
  xlab("") +
  ggtitle("Alignment Rate") +
  theme_qc()
if (length(unique(alignResult$Target)) > 1L) alignFig1a <- alignFig1a + facet_wrap(~ Target)
ggsave(fig_align_rate_path, alignFig1a, width = 8, height = 6, dpi = 300, bg = "white")

# Mapped fragments
alignFig2a <- alignResult %>%
  ggplot(aes(x = Genotype, y = MappedFragNum)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Replicate), position = position_jitter(width = 0.15, height = 0)) +
  ylab("Count") +
  xlab("") +
  ggtitle("Mapped Fragments") +
  theme_qc()
if (length(unique(alignResult$Target)) > 1L) alignFig2a <- alignFig2a + facet_wrap(~ Target)
ggsave(fig_mapped_path, alignFig2a, width = 8, height = 6, dpi = 300, bg = "white")

# Sequencing depth
alignFig3a <- alignResult %>%
  ggplot(aes(x = Genotype, y = SequencingDepth)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Replicate), position = position_jitter(width = 0.15, height = 0)) +
  ylab("Count") +
  xlab("") +
  ggtitle("Sequencing Depth") +
  theme_qc()
if (length(unique(alignResult$Target)) > 1L) alignFig3a <- alignFig3a + facet_wrap(~ Target)
ggsave(fig_depth_path, alignFig3a, width = 8, height = 6, dpi = 300, bg = "white")

message("Alignment QC done: ", rds_path)
