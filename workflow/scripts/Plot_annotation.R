#!/usr/bin/env Rscript

# Plot Annotation: Plot mean coverage, per-base error, and reactivity
# around genomic features from averaged annotation files.
#
# Reads one forward-strand and one reverse-strand averaged annotation file,
# then produces a 3-panel PDF (Coverage | Per-base Error | Reactivity) arranged
# in a single row.  Lines are coloured by Sample (Treatment/Control) and use
# solid (forward) or dotted (reverse) line type for the genome strand.
#
# Usage:
#   Rscript Plot_annotation.R <forward.txt.gz> <reverse.txt.gz> <output.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript Plot_annotation.R <forward.txt.gz> <reverse.txt.gz> <output.pdf>\n")
  quit(status = 1)
}

for_file <- args[1]
rev_file <- args[2]
out_pdf  <- args[3]

# ── Read both files and tag with genome strand ────────────────────────────────
read_strand_file <- function(path, genome_strand) {
  tryCatch({
    con <- gzfile(path, "rt")
    on.exit(close(con))
    dat <- read.delim(con, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    dat$Genome_strand <- genome_strand
    dat
  }, error = function(e) {
    warning("Could not read: ", path, " — ", conditionMessage(e))
    NULL
  })
}

dat_for <- read_strand_file(for_file, "forward")
dat_rev <- read_strand_file(rev_file, "reverse")

dat <- do.call(rbind, Filter(Negate(is.null), list(dat_for, dat_rev)))
if (is.null(dat) || nrow(dat) == 0) stop("No data loaded from input files.")

# ── Coerce to numeric ─────────────────────────────────────────────────────────
dat$Distance      <- suppressWarnings(as.numeric(dat$Distance))
dat$Coverage      <- suppressWarnings(as.numeric(dat$Coverage))
dat$Perbase_error <- suppressWarnings(as.numeric(dat$Perbase_error))
dat$Reactivity    <- suppressWarnings(as.numeric(dat$Reactivity))

# Order genome strand so forward appears first in the legend
dat$Genome_strand <- factor(dat$Genome_strand, levels = c("forward", "reverse"))

# ── Aggregate across BED feature strand ("+"/"-") ─────────────────────────────
# The "Strand" column inside each file is the BED feature strand; we average
# it away so each (Distance, Sample, Genome_strand) has a single value.
dat_cov <- dat %>%
  filter(!is.na(Distance), !is.na(Coverage)) %>%
  group_by(Distance, Sample, Genome_strand) %>%
  summarise(Coverage = mean(Coverage, na.rm = TRUE), .groups = "drop")

dat_err <- dat %>%
  filter(!is.na(Distance), !is.na(Perbase_error)) %>%
  group_by(Distance, Sample, Genome_strand) %>%
  summarise(Perbase_error = mean(Perbase_error, na.rm = TRUE), .groups = "drop")

# Reactivity is only non-empty for Treatment rows
dat_react <- dat %>%
  filter(!is.na(Distance), !is.na(Reactivity)) %>%
  group_by(Distance, Sample, Genome_strand) %>%
  summarise(Reactivity = mean(Reactivity, na.rm = TRUE), .groups = "drop")

# ── Shared aesthetics ─────────────────────────────────────────────────────────
sample_colors  <- c("Treatment" = "black", "Control"   = "grey50")
strand_types   <- c("forward"   = "solid", "reverse"   = "dotted")

base_theme <- theme_classic(base_size = 6) +
  theme(
    text            = element_text(size = 6, color = "black"),
    axis.text       = element_text(size = 6, color = "black"),
    axis.title      = element_text(size = 6, color = "black"),
    legend.text     = element_text(size = 6, color = "black"),
    legend.title    = element_text(size = 6, color = "black"),
    plot.title      = element_text(size = 6, color = "black", hjust = 0.5),
    strip.text      = element_text(size = 6, color = "black"),
    axis.line       = element_line(color = "black"),
    axis.ticks      = element_line(color = "black"),
    legend.position = "none"
  )

# ── Panel 1: Mean coverage ────────────────────────────────────────────────────
p1 <- ggplot(dat_cov,
             aes(x = Distance, y = Coverage,
                 color = Sample, linetype = Genome_strand,
                 group = interaction(Sample, Genome_strand))) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_linetype_manual(values = strand_types, name = "Strand") +
  labs(title = "Mean Coverage", x = "Distance (bp)", y = "Coverage") +
  base_theme

# ── Panel 2: Per-base error ───────────────────────────────────────────────────
p2 <- ggplot(dat_err,
             aes(x = Distance, y = Perbase_error,
                 color = Sample, linetype = Genome_strand,
                 group = interaction(Sample, Genome_strand))) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_linetype_manual(values = strand_types, name = "Strand") +
  labs(title = "Per-Base Error", x = "Distance (bp)", y = "Per-base error") +
  base_theme

# ── Panel 3: Reactivity ───────────────────────────────────────────────────────
p3 <- ggplot(dat_react,
             aes(x = Distance, y = Reactivity,
                 color = Sample, linetype = Genome_strand,
                 group = interaction(Sample, Genome_strand))) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_linetype_manual(values = strand_types, name = "Strand") +
  labs(title = "Reactivity", x = "Distance (bp)", y = "Reactivity") +
  base_theme

# Move legend to right side of the Reactivity panel only (single vertical column)
p3 <- p3 + theme(
  legend.position  = "right",
  legend.direction = "vertical"
)

# ── Render PDF: 7 × 3 inches, three panels in one row ────────────────────────
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(out_pdf, width = 7, height = 3)
grid.arrange(p1, p2, p3, ncol = 3, nrow = 1, widths = c(1, 1, 1.4))
dev.off()

cat("PDF written to:", out_pdf, "\n")
