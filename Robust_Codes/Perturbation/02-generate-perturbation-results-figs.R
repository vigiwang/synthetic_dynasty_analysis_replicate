#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "Perturbation"
title <- "generate_perturbation_results_figs"

HPC <- TRUE # Set FALSE to run locally with downloaded Data/main_results inputs

if (HPC) {
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
} else {
  dir_root <- "~/Desktop/synthetic_dynasty_analysis"
}

dir_root <- path.expand(dir_root)
setwd(dir_root)

dir_log <- file.path(dir_root, "code", "_LOGS", section)
log_path <- file.path(
  dir_log,
  paste0(subsection, "_", title, "_log.txt")
)

dir_perturbation_data <- file.path(dir_root, "Data", "Perturbation_estimation")
dir_perturbation_plot <- file.path(dir_root, "Plot", "Perturbation")
dir_perturbation_fun <- file.path(dir_root, "perturbation_fun")

create_dir_if_missing <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message("Created directory: ", dir)
  } else {
    message("Directory already exists: ", dir)
  }
}

create_dir_if_missing(dir_root)
create_dir_if_missing(dir_log)
create_dir_if_missing(dir_perturbation_plot)

options(warn = -1)


#-------------------------------------------------------------------------------
# U.S. Occupational Mobility — Replication Files
#
# Project:     U.S. Occupational Mobility Analysis
# Repository:  https://github.com/synthetic_dynasty_analysis
#
# Script:      02-generate-perturbation-results-figs.R
#
# Inputs:
#   - Canonical unrounded estimation object:
#       Data/Perturbation_estimation/perturbation_results.rds
#   - Paper figure style:
#       Functions/plot.R
#   - Perturbation helpers:
#       perturbation_fun/perturbation_inputs.R
#       perturbation_fun/perturbation_figures.R
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/Perturbation_generate_perturbation_results_figs_log.txt
#   - Figures (ggplot RDS always; vector PDF local runs only). Each of the
#     three figure families is split into four measure-block sub-figures
#     (a = movement + SSM, b = MTE Class 1-5, c = AIM + IM Class 1-5,
#     d = Class-4 targets; AMTE is tables-only) so every sub-figure is
#     readable at text width:
#       Plot/Perturbation/figure_K1{a,b,c,d}_trend_changes_*.{rds,pdf}
#       Plot/Perturbation/figure_K2{a,b,c,d}_slope_calibration_*.{rds,pdf}
#       Plot/Perturbation/figure_K3{a,b,c,d}_cohort_contributions_*.{rds,pdf}
#
# Description:
#   Renders the appendix perturbation figures (manuscript lettering K)
#   from the canonical unrounded perturbation results. Family K.1:
#   deterministic Delta beta with paired-bootstrap 95% intervals (both
#   windows, both mechanisms). Family K.2: the K.1 marks over the signed
#   total-system and targeted benchmark bands. Family K.3: cohort-specific
#   OLS slope contributions in four window-mechanism columns. Grayscale
#   only; style follows Functions/plot.R. The canonical rds keeps the
#   original figure_J*_data key names.
#
#   PDF rendering is LOCAL-ONLY: Midway nodes ship no serif font, so an
#   HPC run saves the ggplot rds objects only. Download the rds files and
#   render the PDFs on the local machine (rerun this script locally, or
#   print the downloaded rds objects to cairo_pdf at the sizes below).
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES                       #
#-------------------------------------------------#
if (HPC) {
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("dplyr")
  library("ggplot2")
  library("scales")
} else {
  packages <- c("dplyr", "ggplot2", "scales")
  installed <- packages %in% rownames(installed.packages())
  if (any(!installed)) install.packages(packages[!installed])
  invisible(lapply(packages, library, character.only = TRUE))
}

t_start <- Sys.time()

#-------------------------------------------------#
#  SOURCE PAPER STYLE AND HELPERS                  #
#-------------------------------------------------#
source(file.path(dir_root, "Functions", "plot.R"))   # visual-style authority
source(file.path(dir_perturbation_fun, "perturbation_inputs.R"))
source(file.path(dir_perturbation_fun, "perturbation_figures.R"))

#-------------------------------------------------#
#  SECTION 1: LOAD CANONICAL RESULTS               #
#-------------------------------------------------#
cfg <- get_perturbation_config()
res <- readRDS(file.path(dir_perturbation_data, cfg$artifacts$results_rds))
dictionary <- res$measure_dictionary

#-------------------------------------------------#
#  SECTION 2: SUB-FIGURE SPECIFICATION             #
#-------------------------------------------------#
# measure blocks and per-sub-figure canvas sizes (inches):
#   movement = all movement measures + SSM (pp per decade)
#   mte      = MTE Class 1-5 (generations per decade; AMTE tables-only)
#   memory   = AIM + IM Class 1-5 (index units per decade)
#   class4   = steady-state share Class 4 + MFPT into Class 4
fig_groups <- perturbation_figure_groups()
group_letter <- c(movement = "a", mte = "b", memory = "c", class4 = "d")

# all K1/K2 sub-figures share ONE canvas width; combined with the
# equalized y-label column (see .figure_ytext_margin) this makes the
# panels identical in position and length across sub-figures, so they
# align when stacked in the manuscript
fig_dims <- list(
  K1 = list(movement = c(15, 8.5), mte = c(15, 8), memory = c(15, 8.5),
            class4 = c(15, 7)),
  K2 = list(movement = c(15, 8.5), mte = c(15, 8), memory = c(15, 8.5),
            class4 = c(15, 7)),
  K3 = list(movement = c(15, 9), mte = c(15, 8), memory = c(15, 9),
            class4 = c(15, 5))
)

# align_axes: K1/K2 PDFs render through perturbation_align_axis_grob(),
# which fixes every y-axis label column to one width so the panels sit
# at identical positions across sub-figures (the rds keeps the plain
# ggplot object)
save_subfigure <- function(fig, stem, dims, align_axes = FALSE) {
  saveRDS(fig, file.path(dir_perturbation_plot, paste0(stem, ".rds")))
  if (!HPC) {
    pdf_path <- file.path(dir_perturbation_plot, paste0(stem, ".pdf"))
    if (align_axes) {
      g <- perturbation_align_axis_grob(fig)
      cairo_pdf(pdf_path, width = dims[1], height = dims[2])
      grid::grid.draw(g)
      dev.off()
    } else {
      ggplot2::ggsave(pdf_path, fig, width = dims[1], height = dims[2],
                      device = cairo_pdf, limitsize = FALSE)
    }
  }
}

#-------------------------------------------------#
#  SECTION 3: FIGURE K.1 (TREND CHANGES)           #
#-------------------------------------------------#
for (g in names(fig_groups)) {
  key <- paste0("K1", group_letter[[g]])
  fig <- make_figure_K1_trend_changes(
    res$figure_J1_data, dictionary, measure_ids = fig_groups[[g]])
  save_subfigure(fig, cfg$artifacts$figures[[key]], fig_dims$K1[[g]],
                 align_axes = TRUE)
}

#-------------------------------------------------#
#  SECTION 4: FIGURE K.2 (SLOPE CALIBRATION)       #
#-------------------------------------------------#
for (g in names(fig_groups)) {
  key <- paste0("K2", group_letter[[g]])
  fig <- make_figure_K2_slope_calibration(
    res$figure_J2_data, dictionary, measure_ids = fig_groups[[g]])
  save_subfigure(fig, cfg$artifacts$figures[[key]], fig_dims$K2[[g]],
                 align_axes = TRUE)
}

#-------------------------------------------------#
#  SECTION 5: FIGURE K.3 (COHORT CONTRIBUTIONS)    #
#-------------------------------------------------#
for (g in names(fig_groups)) {
  key <- paste0("K3", group_letter[[g]])
  fig <- make_figure_K3_cohort_contributions(
    res$figure_J3_data, dictionary, measure_ids = fig_groups[[g]],
    divider_full = cfg$contribution_divider_full,
    divider_late = cfg$contribution_divider_late)
  save_subfigure(fig, cfg$artifacts$figures[[key]], fig_dims$K3[[g]])
}

#-------------------------------------------------#
#  LOG                                             #
#-------------------------------------------------#
sink(log_path, split = TRUE)
cat("Script: 02-generate-perturbation-results-figs.R\n")
cat("Mode: ", ifelse(HPC, "HPC (Midway)", "local"), "\n")
cat("Input: ", file.path(dir_perturbation_data, cfg$artifacts$results_rds), "\n")
cat("Outputs:\n", paste0("  ", cfg$artifacts$figures, collapse = "\n"),
    ifelse(HPC, "\n(.rds only; download and render PDFs locally)\n",
           "\n(.rds + .pdf)\n"))
cat("Measures: ", nrow(dictionary), "\n")
cat("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " minutes\n", sep = "")
cat("Perturbation figures completed successfully.\n")
sink()
