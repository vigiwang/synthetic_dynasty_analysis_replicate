#----------Preliminaries----------#
rm(list = ls())
section    <- "NLSY_Validation"
subsection <- "bootstrap"
title      <- "nlsy_mobility_bootstrap_estimation"

HPC <- TRUE # Set FALSE to run locally and test the baseline results

# Project root (the analyst's project). Under it:
#   Functions/            : shared original functions (utils.R, CopulaFunctions1.R)
#   NLSY/                 : NLSY matrix estimator + adapter helpers + this script
#   Data/                 : source inputs (nlsy*_input.rds)
#   Data/NLSY_estimation/  : results (a new folder parallel to main_results, robustness_check)
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
}

if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root) # Set working directory
dir_fun  <- file.path(dir_root, "Functions")             # shared original functions
dir_nlsy <- file.path(dir_root, "NLSY_fun")                  # NLSY matrix estimator + adapter
dir_data <- file.path(dir_root, "Data", "NLSY_estimation", "inputs") # source inputs
dir_out  <- file.path(dir_root, "Data", "NLSY_estimation")          # results
dir_log  <- file.path(dir_out, "_LOGS")
log_path <- paste0(dir_log, "/", subsection, "_", title, "_log.txt")

create_dir_if_missing <- function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message("Created directory: ", dir)
  } else {
    message("Directory already exists: ", dir)
  }
}

create_dir_if_missing(dir_out)
create_dir_if_missing(dir_log)
options(warn = -1)
#-------------------------------------------------------------------------------
# NLSY79/NLSY97 Occupational & Income Mobility — Bootstrap of the Analyst's Measures
#
# Project:     U.S. Mobility Analysis — NLSY validation
#
# Script:      05_hpc_bootstrap.R
#
# Inputs:
#   - Data/  : analytic dyads (self-contained, like gss_fc_occ10_5class.rds)
#       nlsy79_income_input.rds / nlsy79_occ_input.rds / nlsy97_income_input.rds
#   - Functions/ : shared original functions (utils.R, CopulaFunctions1.R)
#   - NLSY/  : the MATRIX single_simulation() estimator (synthetic_dynasty_estimation.R)
#       and the adapter (nlsy_data_helpers.R, nlsy_transition_helpers.R,
#       nlsy_measure_bootstrap.R : three resamplers + weighted-matrix construction)
#
# Outputs (Data/NLSY_estimation/):
#   - <SAMPLE>_baseline.rds : point-estimate measures per cohort (+ "Full")
#   - <SAMPLE>_boot.rds      : list of `times` replicate results; each a named
#                              list (cohort -> full single_simulation() output)
#   - <SAMPLE>_ci_summary.csv: percentile CIs for headline measures per cohort
#   - se_comparison_design_vs_individual.csv : design-based vs naive individual
#       bootstrap SEs (se_ratio > 1 = naive understates the SE)
#
# Description:
#   Loops over the three matched-cohort NLSY analytic samples, each run TWICE:
#   with its design-based resampler (NLSY79 household cluster; NLSY97 survey
#   Fay-BRR) and with the naive individual bootstrap for SE comparison. Each
#   replicate rebuilds the WEIGHTED transition-matrix series (recomputing pooled
#   income cut-points inside each replicate; EGP fixed; NO smoothing) and applies
#   the analyst's single_simulation() unchanged via mclapply. Mirrors
#   01_generate_main_results_estimation.R (times = 2000).
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES                      #
#-------------------------------------------------#

if(HPC){
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("dplyr")
  library("tidyverse")
  library("survey")
  library("parallel")
  library("markovchain")
  library("logmult")
  library("assertthat")
  library("pracma")
  library("ipfr")
  library("combinat")
  library("readr")
}else{
  packages <- c(
    "dplyr",
    "tidyverse",
    "survey",
    "parallel",
    "markovchain",
    "logmult",
    "assertthat",
    "pracma",
    "ipfr",
    "combinat",
    "readr"
  )

  install_and_load <- function(pkg_list) {
    for (pkg in pkg_list) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing missing package: ", pkg)
        install.packages(pkg, dependencies = TRUE)
      }
      library(pkg, character.only = TRUE)
    }
  }

  install_and_load(packages)
}


#-------------------------------------------------#
#  LOAD THE FUNCTIONS                             #
#-------------------------------------------------#

# Shared original functions (from Functions/)
source(file.path(dir_fun, "utils.R"))
source(file.path(dir_fun, "CopulaFunctions1.R"))

# NLSY estimator: the MATRIX-based single_simulation(tm, size, tol, weighted).
# This DIFFERS from the GSS dyad-based single_simulation in Functions/ (which
# smooths and takes individual dyads), so it is kept in NLSY/ to avoid clobbering
# the GSS pipeline. Adapter helper functions are here too (NEW; separate from the
# original Functions scripts).
source(file.path(dir_nlsy, "synthetic_dynasty_estimation.R"))
source(file.path(dir_nlsy, "nlsy_data_functions.R"))
source(file.path(dir_nlsy, "nlsy_estimation_functions.R"))


#-------------------------------------------------------#
#     LOAD DATA AND SET GLOBAL PARAMETERS:
#-------------------------------------------------------#

# Global parameters (shared by all three samples):
n_cores <- detectCores() - 2
times   <- if (HPC) 2000L else 20L   # full run on HPC; small smoke value locally
size    <- 5
tol     <- 0.01
set.seed(980625)

# The analytic samples. Each appears TWICE: once with its DESIGN-BASED scheme
# (household cluster for NLSY79; Fay-BRR survey replicates for NLSY97) and once
# with the NAIVE individual bootstrap (rows resampled i.i.d., ignoring sibling
# clustering / the survey design). Comparing the two shows how much the naive
# bootstrap UNDERSTATES the standard errors (see the SE comparison at the end).
# One HPC submission loops over all of them; subset SAMPLES below to run fewer.
sample_configs <- list(
  # ---- design-based (main) ----
  nlsy79_income     = list(input = "nlsy79_income_input.rds", scheme = "household",
                           type = "income", cohorts = 1961:1964,
                           parent = "parent", child = "child", weight = "weight", hhid = "hhid"),
  nlsy79_occ        = list(input = "nlsy79_occ_input.rds", scheme = "household",
                           type = "occ", cohorts = 1961:1964,
                           parent = "status_p", child = "status_c", weight = "weight", hhid = "hhid"),
  nlsy97_income     = list(input = "nlsy97_income_input.rds", scheme = "survey",
                           type = "income", cohorts = 1980:1984,
                           parent = "parent", child = "child", weight = "weight", hhid = NULL),
  # ---- family-size UNADJUSTED income (parallel analysis) ----
  # Same dyads / weights / cluster ids as the primary runs above; the only
  # difference is parent/child hold size-undivided real$ income (no sqrt-size
  # denominator). Lets us bootstrap the mobility measures without the
  # family-size equivalization, to isolate the adjustment's effect.
  nlsy79_income_unadj = list(input = "nlsy79_income_unadj_input.rds", scheme = "household",
                             type = "income", cohorts = 1961:1964,
                             parent = "parent", child = "child", weight = "weight", hhid = "hhid"),
  nlsy97_income_unadj = list(input = "nlsy97_income_unadj_input.rds", scheme = "survey",
                             type = "income", cohorts = 1980:1984,
                             parent = "parent", child = "child", weight = "weight", hhid = NULL),
  # ---- individual bootstrap (naive comparison) ----
  nlsy79_income_ind = list(input = "nlsy79_income_input.rds", scheme = "individual",
                           type = "income", cohorts = 1961:1964,
                           parent = "parent", child = "child", weight = "weight", hhid = NULL),
  nlsy79_occ_ind    = list(input = "nlsy79_occ_input.rds", scheme = "individual",
                           type = "occ", cohorts = 1961:1964,
                           parent = "status_p", child = "status_c", weight = "weight", hhid = NULL),
  nlsy97_income_ind = list(input = "nlsy97_income_input.rds", scheme = "individual",
                           type = "income", cohorts = 1980:1984,
                           parent = "parent", child = "child", weight = "weight", hhid = NULL),
  # ---- family-size UNADJUSTED income, individual bootstrap (naive comparison) ----
  nlsy79_income_unadj_ind = list(input = "nlsy79_income_unadj_input.rds", scheme = "individual",
                                 type = "income", cohorts = 1961:1964,
                                 parent = "parent", child = "child", weight = "weight", hhid = NULL),
  nlsy97_income_unadj_ind = list(input = "nlsy97_income_unadj_input.rds", scheme = "individual",
                                 type = "income", cohorts = 1980:1984,
                                 parent = "parent", child = "child", weight = "weight", hhid = NULL)
)
SAMPLES <- names(sample_configs)   # 10: 3 design main + 2 unadj design + 3 individual + 2 unadj individual

# Headline scalar measures for the CI summary; the FULL single_simulation output is
# retained in the saved *_boot.rds for any further post-processing (curves, matrices).
accessors <- list(
  AMTE      = function(r) r$AMTE,
  lambda2   = function(r) r$lambda2,
  Altham    = function(r) r$altham_index,
  GeenensD  = function(r) r$GeenensD,
  Hellinger = function(r) r$HellingerDep
)

# Open log (captures the per-sample progress for the whole submission):
sink(log_path, split = TRUE)
cat("NLSY mobility bootstrap | HPC =", HPC, "| cores =", n_cores,
    "| base times =", times, "\n\n")


#-------------------------------------------------------#
#     RUN EACH SAMPLE: baseline -> bootstrap -> CI
#-------------------------------------------------------#

for (SAMPLE in SAMPLES) {

  cfg <- sample_configs[[SAMPLE]]
  cat("==================  ", SAMPLE, " (", cfg$scheme, ", ", cfg$type,
      ")  ==================\n")

  input_data <-
    readr::read_rds(
      file.path(
        dir_data,
        cfg$input
      )
    )

  # For the NLSY97 survey design, build the Fay-BRR replicate weights (a FIXED
  # replicate set) and use its column count as the replicate count.
  if (cfg$scheme == "survey") {
    design97 <-
      svydesign(
        ids     = ~vpsu,
        strata  = ~vstrat,
        weights = ~weight,
        data    = input_data,
        nest    = TRUE
      )
    repdesign97 <-
      as.svrepdesign(
        design97,
        type    = "Fay",
        fay.rho = 0.5,
        mse     = TRUE
      )
    fay_repweights <- weights(repdesign97, "analysis")
    times_use      <- ncol(fay_repweights)
  } else {
    fay_repweights <- NULL
    times_use      <- times
  }

  #-------------------------------------------------#
  #     CALCULATE THE BASELINE RESULTS:             #
  #-------------------------------------------------#

  baseline_rst <-
    baseline_measures(
      data       = input_data,
      parent_col = cfg$parent,
      child_col  = cfg$child,
      weight_col = cfg$weight,
      cohort_col = "cohort",
      cohorts    = cfg$cohorts,
      type       = cfg$type,
      size       = size,
      tol        = tol
    )

  write_rds(
    baseline_rst,
    file.path(
      dir_out,
      paste0(SAMPLE, "_baseline.rds"))
  )
  cat("   baseline finished.\n")

  #-------------------------------------------------#
  #     CALCULATE AND PROCESS BOOTSTRAP RESULTS     #
  #-------------------------------------------------#

  boot_rst_lst <-
    mclapply(
      1:times_use,
      function(i){
        set.seed(i)
        single_bootstrap_replicate(
          data          = input_data,
          scheme        = cfg$scheme,
          parent_col    = cfg$parent,
          child_col     = cfg$child,
          weight_col    = cfg$weight,
          cohort_col    = "cohort",
          cohorts       = cfg$cohorts,
          type          = cfg$type,
          size          = size,
          tol           = tol,
          hhid_col      = cfg$hhid,
          repweight_vec = if (cfg$scheme == "survey") fay_repweights[, i] else NULL
        )
      },
      mc.cores = n_cores
    )

  write_rds(
    boot_rst_lst,
    file.path(
      dir_out,
      paste0(SAMPLE, "_boot.rds")),
    compress = "xz"
  )
  cat("   bootstrap finished (", times_use, "replicates).\n")

  #-------------------------------------------------#
  #     PROCESS: PERCENTILE CI SUMMARY              #
  #-------------------------------------------------#

  ci_summary <-
    do.call(
      rbind,
      lapply(
        c("Full", as.character(cfg$cohorts)),
        function(co)
          do.call(
            rbind,
            lapply(
              names(accessors),
              function(m)
                cbind(
                  measure = m,
                  boot_ci(
                    boot     = boot_rst_lst,
                    accessor = accessors[[m]],
                    cohort   = co,
                    point    = baseline_rst,
                    fay_rho  = if (cfg$scheme == "survey") 0.5 else NULL
                  )
                )
            )
          )
      )
    )

  write_csv(
    as.data.frame(ci_summary),
    file.path(
      dir_out,
      paste0(SAMPLE, "_ci_summary.csv"))
  )
  cat("   CI summary finished. Saved baseline/boot/ci to:", dir_out, "\n\n")
}

#-------------------------------------------------#
#     SE COMPARISON: DESIGN-BASED vs INDIVIDUAL   #
#-------------------------------------------------#

# For each sample, join the design-based CI summary to its individual-bootstrap
# counterpart and compute se_ratio = design SE / individual SE. Ratios > 1 show
# how much the naive individual bootstrap UNDERSTATES the design-based SE
# (sibling clustering for NLSY79; stratified/clustered Fay-BRR for NLSY97).
design_to_ind <- c(
  nlsy79_income = "nlsy79_income_ind",
  nlsy79_occ    = "nlsy79_occ_ind",
  nlsy97_income = "nlsy97_income_ind"
)

se_comparison <-
  do.call(
    rbind,
    lapply(
      names(design_to_ind),
      function(s) {
        des <- read_csv(
          file.path(dir_out, paste0(s, "_ci_summary.csv")),
          show_col_types = FALSE
        )
        ind <- read_csv(
          file.path(dir_out, paste0(design_to_ind[[s]], "_ci_summary.csv")),
          show_col_types = FALSE
        )
        des %>%
          select(measure, cohort, estimate, se_design = se) %>%
          left_join(
            ind %>% select(measure, cohort, se_individual = se),
            by = c("measure", "cohort")
          ) %>%
          mutate(
            sample   = s,
            se_ratio = round(se_design / se_individual, 3)
          ) %>%
          select(sample, measure, cohort, estimate, se_design, se_individual, se_ratio)
      }
    )
  )

write_csv(
  se_comparison,
  file.path(
    dir_out,
    "se_comparison_design_vs_individual.csv")
)

cat("SE comparison (design vs individual) finished.\n\n")
cat("All NLSY mobility bootstrap tasks completed successfully.\n")
sink()
