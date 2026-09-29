#!/usr/bin/env Rscript
# ===========================================================================
# 08_render_figures.R  — render the fig*.rds objects to PAPER-READY figures.
# Reads Plot/robustness_check/NLSY/fig*_{equalized,unequalized}.rds (styled
# ggplot objects from 07) and writes vector PDF (for LaTeX) + 300-dpi PNG.
# Headless-safe: PDF via cairo_pdf, PNG via cairo. Run locally on the fig-rds
# folder brought back from Midway (or right after --stage=results).
# ===========================================================================
HPC      <- TRUE
dir_root <- if (HPC) "/home/weiqiw/synthetic_dynasty_analysis" else "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)
suppressPackageStartupMessages({library(ggplot2)})

dir_fig <- file.path(dir_root, "Plot", "robustness_check", "NLSY")
dir_pdf <- file.path(dir_fig, "pdf"); dir_png <- file.path(dir_fig, "png")
dir.create(dir_pdf, showWarnings = FALSE); dir.create(dir_png, showWarnings = FALSE)

# per-figure page size (inches), matching 07's create-figX dims
dims <- function(nm) {
  if (grepl("UPDOWN|_ADJ", nm))                 c(7.6, 4.8)
  else if (grepl("RANKRANK", nm))               c(9.0, 4.6)
  else if (grepl("^fig[58]_", nm))              c(9.0, 4.6)   # class-by-survey panels
  else if (grepl("_MTE_CLS|_IM_CLS", nm))       c(10.0, 4.0)  # class-compare panels
  else                                          c(7.6, 4.6)
}
has_cairo <- capabilities("cairo")

figs <- list.files(dir_fig, pattern = "^fig.*_(equalized|unequalized)\\.rds$", full.names = TRUE)
if (!length(figs)) stop("[08] no fig*_{equalized,unequalized}.rds in ", dir_fig, " -- run --stage=results first")
n <- 0L
for (f in figs) {
  nm <- sub("\\.rds$", "", basename(f)); p <- readRDS(f); wh <- dims(nm)
  ggsave(file.path(dir_pdf, paste0(nm, ".pdf")), p, width = wh[1], height = wh[2],
         device = if (has_cairo) cairo_pdf else "pdf")
  ggsave(file.path(dir_png, paste0(nm, ".png")), p, width = wh[1], height = wh[2],
         dpi = 300, type = if (has_cairo) "cairo" else NULL)
  n <- n + 1L
}
cat(sprintf("[08] rendered %d figures -> %s (PDF) and %s (PNG, 300 dpi)\n", n, dir_pdf, dir_png))
