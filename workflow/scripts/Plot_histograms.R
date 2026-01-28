#!/usr/bin/env Rscript

# Plot Histograms: Generate PDF histograms from histogram data
# Usage: Rscript Plot_histograms.R <input_histograms.txt> <output_histograms.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript Plot_histograms.R <input_histograms.txt> <output_histograms.pdf>\n")
  quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]

# Read histogram data
data <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Function to parse COV bracket values like "[1-1]" -> extract first number
parse_cov_value <- function(x) {
  if (grepl("^\\[", x)) {
    # Extract the first number from bracket format [n-m]
    as.numeric(gsub("\\[(\\d+)-.*", "\\1", x))
  } else {
    as.numeric(x)
  }
}

# Function to create histogram with 50 bins and outlier handling
create_histogram <- function(df, var_name, var_label) {
  if (nrow(df) == 0) {
    return(NULL)
  }

  # Parse values (handle COV bracket format)
  if (var_name == "COV") {
    df$Value_numeric <- sapply(df$Value, parse_cov_value)
  } else {
    df$Value_numeric <- as.numeric(df$Value)
  }

  # Remove NA values
  df <- df[!is.na(df$Value_numeric), ]

  if (nrow(df) == 0) {
    return(NULL)
  }

  # Expand data: replicate each value by its count for proper histogram
  # For large counts, we'll use weighted histogram approach instead

  # Calculate quantiles for outlier handling (weighted by count)
  total_count <- sum(df$Count)
  df <- df[order(df$Value_numeric), ]
  df$cumsum <- cumsum(df$Count)
  df$cum_frac <- df$cumsum / total_count

  # Find 0.01 and 0.99 quantile values
  q01_idx <- which(df$cum_frac >= 0.01)[1]
  q99_idx <- which(df$cum_frac >= 0.99)[1]

  if (is.na(q01_idx)) q01_idx <- 1
  if (is.na(q99_idx)) q99_idx <- nrow(df)

  q01_val <- df$Value_numeric[q01_idx]
  q99_val <- df$Value_numeric[q99_idx]

  # Handle edge case where q01 == q99
  if (q01_val >= q99_val) {
    q01_val <- min(df$Value_numeric)
    q99_val <- max(df$Value_numeric)
  }

  # Create 50 bins between q01 and q99
  bin_breaks <- seq(q01_val, q99_val, length.out = 51)

  # Assign values to bins (outliers go to first/last bin)
  df$bin <- cut(df$Value_numeric,
                breaks = c(-Inf, bin_breaks[-c(1, length(bin_breaks))], Inf),
                labels = FALSE,
                include.lowest = TRUE)

  # Aggregate counts by bin
  bin_data <- df %>%
    group_by(bin) %>%
    summarise(
      Count = sum(Count),
      Value_mid = mean(Value_numeric),
      .groups = "drop"
    )

  # Calculate bin midpoints for x-axis
  bin_mids <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2
  bin_data$x <- bin_mids[pmin(bin_data$bin, length(bin_mids))]

  # Create plot
  p <- ggplot(bin_data, aes(x = x, y = Count)) +
    geom_col(fill = "gray40", color = "black", linewidth = 0.2) +
    labs(x = var_label, y = "Count", title = var_label) +
    theme_classic(base_size = 6) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black", hjust = 0.5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )

  return(p)
}

# Define variable labels
var_labels <- c(
  "RL" = "Read Length",
  "MAPQ" = "Mapping Quality",
  "INS" = "Insertion Size",
  "DEL" = "Deletion Size",
  "COV" = "Coverage"
)

# Get unique variables in the data
vars_present <- unique(data$Var)
vars_to_plot <- intersect(names(var_labels), vars_present)

# Create plots for each variable
plots <- list()
for (var in vars_to_plot) {
  df_var <- data[data$Var == var, ]
  p <- create_histogram(df_var, var, var_labels[var])
  if (!is.null(p)) {
    plots[[var]] <- p
  }
}

# Calculate PDF dimensions
# 2 columns, 3 inches wide each = 6 inches total width
# Height based on number of plots (2 per row)
n_plots <- length(plots)
n_rows <- ceiling(n_plots / 2)
plot_height <- 2  # inches per plot
total_height <- n_rows * plot_height + 0.5  # extra space for margins

# Save PDF
pdf(output_file, width = 6, height = total_height)

if (length(plots) > 0) {
  # Arrange plots in 2 columns
  do.call(grid.arrange, c(plots, ncol = 2))
}

dev.off()

cat("PDF saved to:", output_file, "\n")
