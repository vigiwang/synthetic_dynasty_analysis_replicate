#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "Perturbation_estimation"
title <- "calculate_perturbation_results_tables"

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
create_dir_if_missing(dir_perturbation_data)

options(warn = -1)


#-------------------------------------------------------------------------------
# U.S. Occupational Mobility — Replication Files
#
# Project:     U.S. Occupational Mobility Analysis
# Repository:  https://github.com/synthetic_dynasty_analysis
#
# Script:      03-calculate-perturbation-results-tables.R
#
# Inputs:
#   - Canonical unrounded estimation object:
#       Data/Perturbation_estimation/perturbation_results.rds
#   - Perturbation helpers:
#       perturbation_fun/perturbation_inputs.R
#       perturbation_fun/perturbation_tables.R
#
# Outputs:
#   - Log file:
#       code/_LOGS/Data/Perturbation_estimation_calculate_perturbation_results_tables_log.txt
#   - LaTeX tables:
#       Data/Perturbation_estimation/table_J1_sparse_cell_design.tex
#       Data/Perturbation_estimation/table_J2_cohort_trend_results.tex
#       Data/Perturbation_estimation/table_J3_slope_bootstrap_calibration.tex
#
# Description:
#   Renders the three Appendix J LaTeX tables from the canonical unrounded
#   perturbation results. Table J.1: sparse-cell design and stationary-TV /
#   bootstrap calibration by mechanism and era. Table J.2: cohort-trend
#   results for all 20 paper measures under both mechanisms (full and late
#   panels). Table J.3: perturbation-induced trend changes relative to
#   total-system and targeted bootstrap slope variation.
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES                       #
#-------------------------------------------------#
if (HPC) {
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("dplyr")
} else {
  packages <- c("dplyr")
  installed <- packages %in% rownames(installed.packages())
  if (any(!installed)) install.packages(packages[!installed])
  invisible(lapply(packages, library, character.only = TRUE))
}

t_start <- Sys.time()

#-------------------------------------------------#
#  SOURCE HELPERS                                  #
#-------------------------------------------------#
source(file.path(dir_perturbation_fun, "perturbation_inputs.R"))
source(file.path(dir_perturbation_fun, "perturbation_tables.R"))

#-------------------------------------------------#
#  SECTION 1: LOAD CANONICAL RESULTS               #
#-------------------------------------------------#
cfg <- get_perturbation_config()
res <- readRDS(file.path(dir_perturbation_data, cfg$artifacts$results_rds))
dictionary <- res$measure_dictionary

#-------------------------------------------------#
#  SECTION 2: TABLE J.1                            #
#-------------------------------------------------#
table_J1 <- prepare_table_J1_sparse_cell_design(res$table_J1_data)
write_table_J1_tex(table_J1,
  file.path(dir_perturbation_data, cfg$artifacts$table_J1))

#-------------------------------------------------#
#  SECTION 3: TABLE J.2                            #
#-------------------------------------------------#
table_J2 <- prepare_table_J2_cohort_trend_results(res$table_J2_data, dictionary)
write_table_J2_tex(table_J2,
  file.path(dir_perturbation_data, cfg$artifacts$table_J2))

#-------------------------------------------------#
#  SECTION 4: TABLE J.3                            #
#-------------------------------------------------#
table_J3 <- prepare_table_J3_slope_bootstrap_calibration(
  res$table_J3_data, dictionary, cfg$table_J3_mechanism_order)
write_table_J3_tex(table_J3,
  file.path(dir_perturbation_data, cfg$artifacts$table_J3))

#-------------------------------------------------#
#  LOG                                             #
#-------------------------------------------------#
sink(log_path, split = TRUE)
cat("Script: 03-calculate-perturbation-results-tables.R\n")
cat("Mode: ", ifelse(HPC, "HPC (Midway)", "local"), "\n")
cat("Input: ", file.path(dir_perturbation_data, cfg$artifacts$results_rds), "\n")
cat("Outputs: ", paste(c(cfg$artifacts$table_J1, cfg$artifacts$table_J2,
                         cfg$artifacts$table_J3), collapse = ", "), "\n")
cat("Measures: ", nrow(dictionary), "\n")
cat("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " minutes\n", sep = "")
cat("Perturbation tables completed successfully.\n")
sink()
