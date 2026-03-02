#!/usr/bin/env Rscript

# Plot Summary Histograms: Generate a 2-column x 3-row PDF of alignment histograms
# Input:  tab-delimited histogram txt file (Var, Value, Count) from Make_histograms.sh
# Output: PDF with Read Length, Mapping Quality, Insertion Size,
#         Deletion Size, and Coverage histograms

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript Plot_summary_histograms.R <input_histograms.txt> <output.pdf>\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]

# ── Data ──────────────────────────────────────────────────────────────────────
dat <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ── Parse COV bracket notation "[1-1]" or "[>=1000]" → numeric ───────────────
parse_cov_value <- function(x) {
  x <- sub("\\[",  "", x)   # remove leading [
  x <- sub("\\]",  "", x)   # remove trailing ]
  x <- sub(">=",   "", x)   # remove >= prefix
  x <- sub("-.*",  "", x)   # keep only the first number before "-"
  suppressWarnings(as.numeric(x))
}

# ── Weighted 99 % winsorisation + 50-bin re-binning ──────────────────────────
winsorize_and_bin <- function(values, counts, n_bins = 50) {
  # Drop NAs and zero-count rows
  keep   <- !is.na(values) & counts > 0
  values <- values[keep]
  counts <- counts[keep]
  if (length(values) == 0) return(NULL)

  # Sort ascending
  ord    <- order(values)
  values <- values[ord]
  counts <- counts[ord]

  total      <- sum(counts)
  cum_frac   <- cumsum(counts) / total

  # Weighted quantiles at 1 % and 99 %
  q01_idx <- which(cum_frac >= 0.01)[1]
  q99_idx <- which(cum_frac >= 0.99)[1]
  if (is.na(q01_idx)) q01_idx <- 1
  if (is.na(q99_idx)) q99_idx <- length(values)

  q01 <- values[q01_idx]
  q99 <- values[q99_idx]

  # Degenerate case: all values identical
  if (q01 >= q99) {
    q01 <- min(values)
    q99 <- max(values)
  }
  if (q01 >= q99) {
    return(data.frame(x = values, count = counts))
  }

  # Winsorise: cap values (accumulate counts at the caps)
  values_win <- pmax(pmin(values, q99), q01)

  df <- data.frame(v = values_win, c = counts) %>%
    group_by(v) %>%
    summarise(c = sum(c), .groups = "drop")

  # 50 equal-width bins across [q01, q99]
  breaks   <- seq(q01, q99, length.out = n_bins + 1)
  bin_mids <- (breaks[-length(breaks)] + breaks[-1]) / 2

  bin_idx <- cut(df$v, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin_idx[is.na(bin_idx)] <- 1L   # safety: pin any edge case to first bin

  data.frame(x = bin_mids[bin_idx], count = df$c) %>%
    group_by(x) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    arrange(x)
}

# ── Shared theme (8 pt, black, theme_classic) ─────────────────────────────────
base_theme <- theme_classic(base_size = 8) +
  theme(
    text         = element_text(size = 8, color = "black"),
    axis.text    = element_text(size = 8, color = "black"),
    axis.title   = element_text(size = 8, color = "black"),
    legend.text  = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8, color = "black"),
    plot.title   = element_text(size = 8, color = "black", hjust = 0.5),
    strip.text   = element_text(size = 8, color = "black"),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black")
  )

# ── Per-panel plot factory ────────────────────────────────────────────────────
make_panel <- function(var_code, title, xlab) {
  sub <- dat[dat$Var == var_code, ]

  if (nrow(sub) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste("No data:", var_code),
                 size = 8 / .pt, color = "black") +
        theme_void() +
        ggtitle(title) +
        base_theme
    )
  }

  # Parse values
  values <- if (var_code == "COV") {
    sapply(sub$Value, parse_cov_value)
  } else {
    suppressWarnings(as.numeric(sub$Value))
  }

  bd <- winsorize_and_bin(values, sub$Count, n_bins = 50)

  if (is.null(bd) || nrow(bd) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Insufficient data",
                 size = 8 / .pt, color = "black") +
        theme_void() +
        ggtitle(title) +
        base_theme
    )
  }

  bar_width <- if (nrow(bd) > 1) diff(range(bd$x)) / 50 else 1

  ggplot(bd, aes(x = x, y = count)) +
    geom_col(width = bar_width, fill = "black", color = NA) +
    labs(title = title, x = xlab, y = "Count") +
    base_theme
}

# ── Build the five panels ─────────────────────────────────────────────────────
panels <- list(
  make_panel("RL",   "Read Length Distribution",   "Read Length (bp)"),
  make_panel("MAPQ", "Mapping Quality Distribution","Mapping Quality (MAPQ)"),
  make_panel("INS",  "Insertion Size Distribution", "Insertion Size (bp)"),
  make_panel("DEL",  "Deletion Size Distribution",  "Deletion Size (bp)"),
  make_panel("COV",  "Coverage Distribution",       "Coverage (x)")
)

# ── Render to PDF: 6 × 7 inches, 2-column × 3-row grid ───────────────────────
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
pdf(output_file, width = 6, height = 7)
do.call(grid.arrange, c(panels, list(ncol = 2, nrow = 3)))
dev.off()

cat("PDF saved to:", output_file, "\n")
