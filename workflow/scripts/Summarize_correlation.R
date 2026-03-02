#!/usr/bin/env Rscript

# Summarize Correlation: Compute pairwise Spearman/Pearson correlations and
# generate heatmap plots from per-base error correlation files.
#
# Usage:
#   Rscript Summarize_correlation.R <output_table.csv> <output_plot.pdf> <input1.txt> [input2.txt ...]
#
# Inputs:  Per-pair correlation txt files (output of Correlation.sh).
#          Both strands (for + rev) are pooled per sample pair before computing
#          correlation coefficients.
# Outputs: CSV table (Sample_a, Sample_b, N_points, Spearman, Pearson)
#          PDF with two side-by-side heatmaps (Spearman and Pearson).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript Summarize_correlation.R <out_table.csv> <out_plot.pdf> <input1.txt> ...\n")
  quit(status = 1)
}

out_table   <- args[1]
out_pdf     <- args[2]
input_files <- args[3:length(args)]

# ── Filename parsing ──────────────────────────────────────────────────────────
# Filename format: {sample_a}_vs_{sample_b}_{strand}.txt
parse_filename <- function(path) {
  fname <- sub("\\.txt$", "", basename(path))
  vs_idx   <- regexpr("_vs_", fname, fixed = TRUE)
  sample_a <- substr(fname, 1, vs_idx - 1)
  rest     <- substr(fname, vs_idx + 4, nchar(fname))
  strand   <- sub(".*_(for|rev)$", "\\1", rest)
  sample_b <- sub("_(for|rev)$",   "",    rest)
  list(sample_a = sample_a, sample_b = sample_b, strand = strand)
}

# ── Read all files, pool strands per pair ─────────────────────────────────────
pair_data <- list()   # key = "sample_a|||sample_b"

for (f in input_files) {
  if (!file.exists(f)) {
    warning("File not found, skipping: ", f)
    next
  }

  info <- parse_filename(f)
  dat  <- tryCatch(
    read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) == 0) next

  # Keep only the two error columns; coerce to numeric; drop incomplete rows
  e1 <- suppressWarnings(as.numeric(dat$Perbase_error_1))
  e2 <- suppressWarnings(as.numeric(dat$Perbase_error_2))
  ok <- !is.na(e1) & !is.na(e2)
  if (sum(ok) == 0) next

  chunk <- data.frame(e1 = e1[ok], e2 = e2[ok])

  key <- paste(info$sample_a, info$sample_b, sep = "|||")
  if (is.null(pair_data[[key]])) {
    pair_data[[key]] <- list(
      sample_a = info$sample_a,
      sample_b = info$sample_b,
      data     = chunk
    )
  } else {
    pair_data[[key]]$data <- rbind(pair_data[[key]]$data, chunk)
  }
}

if (length(pair_data) == 0) {
  stop("No valid correlation data found in the input files.")
}

# ── Compute Spearman and Pearson per pair ─────────────────────────────────────
results <- lapply(pair_data, function(x) {
  d <- x$data
  n <- nrow(d)
  if (n < 3) {
    spearman <- NA_real_
    pearson  <- NA_real_
  } else {
    spearman <- tryCatch(
      cor(d$e1, d$e2, method = "spearman"), error = function(e) NA_real_)
    pearson  <- tryCatch(
      cor(d$e1, d$e2, method = "pearson"),  error = function(e) NA_real_)
  }
  data.frame(
    Sample_a = x$sample_a,
    Sample_b = x$sample_b,
    N_points = n,
    Spearman = round(spearman, 6),
    Pearson  = round(pearson,  6),
    stringsAsFactors = FALSE
  )
})

cor_df <- do.call(rbind, results)
cor_df <- cor_df[order(cor_df$Sample_a, cor_df$Sample_b), ]
rownames(cor_df) <- NULL

# ── Write CSV table ───────────────────────────────────────────────────────────
dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)
write.csv(cor_df, out_table, row.names = FALSE)
cat("Table written to:", out_table, "\n")

# ── Shared theme: theme_classic, 6 pt black font ──────────────────────────────
base_theme <- theme_classic(base_size = 6) +
  theme(
    text            = element_text(size = 6, color = "black"),
    axis.text       = element_text(size = 6, color = "black"),
    axis.text.x     = element_text(size = 6, color = "black",
                                   angle = 45, hjust = 1, vjust = 1),
    axis.title      = element_text(size = 6, color = "black"),
    legend.text     = element_text(size = 6, color = "black"),
    legend.title    = element_text(size = 6, color = "black"),
    plot.title      = element_text(size = 6, color = "black", hjust = 0.5),
    strip.text      = element_text(size = 6, color = "black"),
    legend.position = "bottom",
    legend.key.width = unit(1, "cm")
  )

# ── Heatmap factory ───────────────────────────────────────────────────────────
make_heatmap <- function(df, fill_col, title) {
  df$fill_val <- df[[fill_col]]
  df$label    <- sprintf("%.2f", df$fill_val)

  ggplot(df, aes(x = Sample_a, y = Sample_b, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label), size = 6 / .pt, color = "black") +
    scale_fill_viridis_c(
      option   = "turbo",
      limits   = c(-1, 1),
      name     = title,
      na.value = "grey80"
    ) +
    labs(title = title, x = "Sample A", y = "Sample B") +
    base_theme
}

p_spearman <- make_heatmap(cor_df, "Spearman", "Spearman Correlation")
p_pearson  <- make_heatmap(cor_df, "Pearson",  "Pearson Correlation")

# ── Render PDF: 6 × 5 inches, two panels side by side ────────────────────────
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(out_pdf, width = 6, height = 5)
grid.arrange(p_spearman, p_pearson, ncol = 2)
dev.off()

cat("PDF written to:", out_pdf, "\n")
