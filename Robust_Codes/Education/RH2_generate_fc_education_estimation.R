#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "robustness_check"
title <- "fc_education_estimation"

HPC <- TRUE # Set False if you wish to run locally and test the baseline results

# Specify the root directory:
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "~/Desktop/synthetic_dynasty_analysis"
}

setwd(dir_root) # Set working directory
# Define subdirectories for logs and figures:
dir_log <- paste0(dir_root, "/code/", "_LOGS","/", section)
log_path <- paste0(dir_log, "/", subsection, "_", title, "_log.txt")
dir_data <- paste0(dir_root, "/" ,section, "/", subsection)
# Ensure all necessary directories exist under your root folder
# if not, the function will create folders for you

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
create_dir_if_missing(dir_data)
options(warn = -1)

#-------------------------------------------------------------------------------
# U.S. Occupational Mobility — Replication Files
#
# Project:     U.S. Occupational Mobility Analysis
# Repository:  https://github.com/synthetic_dynasty_analysis_replicate
#
# Script:      RH2_generate_fc_education_estimation.R
#
# Inputs:
#   - GSS father–child educational pairs (5-category years of schooling):
#       Data/robustness_check/gss_fc_edu_5class.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/synthetic_dynasty_estimation.R
#           (Markov-based mobility estimators and related calculations)
#       Functions/process_results.R
#           (post-processing of baseline and bootstrap results, including
#            construction of bias-corrected estimators)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Data/robustness_check_fc_education_estimation_log.txt
#   - Estimation output:
#       Data/robustness_check/5class_fc_edu_bc.rds
#
# Description:
#   This script generates the robustness results using father–child
#   educational mobility pairs constructed from the General Social Survey (GSS).
#   The analysis restricts the sample to respondents born between 1945 and 1990
#   and observed at ages 25–55 at the time of survey.
#
#   The script first computes baseline mobility measures and then provides
#   code to implement the bootstrap procedure. Due to the computational
#   intensity of the bootstrap, users may wish to run this step on
#   high-performance computing (HPC) resources (e.g., Midway).
#
#   For full replication without access to HPC resources, the repository
#   also includes the authors’ precomputed bootstrap outputs, which are
#   used to construct the final bias-corrected estimators.
#
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
if(HPC){
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("here")
  library("mgcv")
  library("tidyverse")
  library("dplyr")
  library("haven")
  library("tidyr")
  library("naniar")
  library("stringr")
  library("gt")
  library("janitor")
  library("kableExtra")
  library("tm")
  library("openxlsx")
  library("scales")
  library("assertthat")
  library("readxl")
  library("labelled")
  library("knitr")
  library("parallel")
  library("readr")
  library("combinat")
  library("ipfr")
  library("ggplot2")
  library("pracma")
  library("markovchain")
  library("logmult")
}else{
  packages <- c(
    "here",
    "mgcv",
    "tidyverse",
    "dplyr",
    "haven",
    "tidyr",
    "naniar",
    "stringr",
    "gt",
    "janitor",
    "kableExtra",
    "tm",
    "openxlsx",
    "scales",
    "assertthat",
    "readxl",
    "labelled",
    "parallel",
    "readr",
    "combinat",
    "ipfr",
    "pracma",
    "markovchain",
    "logmult"
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

# Load helper functions
source(
  file.path(
    dir_root,
    "Functions",
    "utils.R"
  )
)

# Load main mobility estimators
source(
  file.path(
    dir_root,
    "Functions",
    "synthetic_dynasty_estimation.R"
  )
)

# Load copula-based mobility estimators (Prof. Coleman)
source(
  file.path(
    dir_root,
    "Functions",
    "CopulaFunctions1.R"
  )
)

# Load functions for processing estimation results and bias correction
source(
  file.path(
    dir_root,
    "Functions",
    "process_results.R"
  )
)


#--------------------------------------------------------#
#     LOAD DATA AND SET GLOBAL PARAMETERS:
#-------------------------------------------------------#

gss_data_5class_fc_edu <- readr::read_rds(
  file.path(
    dir_data,
    "gss_fc_edu_5class.rds"
  )
)

n_cores <- detectCores() - 2
times <- 2500
tol <- 0.01
size <- 5
gamma_par <- 1
set.seed(980625)


#-------------------------------------------------------#
#     CALCULATE THE BASELINE RESULTS:
#-------------------------------------------------------#

# The two steps of calculate_baseline() are kept separate here so that
# the baseline cohort coverage can be recorded: for the sparse
# father-child sample, early cohorts may yield non-primitive transition
# matrices and be dropped by process_results(), so the baseline is not
# guaranteed to span the full 1945-1990 range.

baseline_raw <-
  mclapply(
    1,
    function(i){
      rst <-
        single_simulation_lst_baseline(
          gss_data = gss_data_5class_fc_edu,
          size     = 5,     # Number of classes in the EGP typology
          age_l    = 25,    # Lower bound of respondents' age at survey
          age_u    = 55,    # Upper bound of respondents' age at survey
          year_l   = 1920,  # Lower bound of birth cohorts used for GAM training
          # (wider than the main sample to capture more information)
          year_u   = 1997,  # Upper bound of birth cohorts used for GAM training
          # (wider than the main sample to capture more information)
          est_l = "1945",
          est_u = "1990",
          tol = 0.01,
          gamma_par = 1)
    },
    mc.cores = n_cores
  )

# Record the cohorts with a valid baseline result (mirrors the drop rule
# used inside process_results()):

baseline_cohort_vec <-
  as.integer(names(baseline_raw[[1]]))[
    !sapply(baseline_raw[[1]], function(x) all(is.na(x[[1]])))
  ]

baseline_rst <-
  process_results(
    boot_rst = baseline_raw,
    times = 1,
    size = size,
    n_cores = n_cores,
    reshape = TRUE
  )

assert_that(
  length(baseline_rst) == length(baseline_cohort_vec)
)

print("Baseline result of father-child sample has finished!")


#-------------------------------------------------#
#     CALCULATE AND PROCESS BOOTSTRAP RESULTS     #
#-------------------------------------------------#

boot_rst_lst <-
  mclapply(
    1:times,
    function(i){
      set.seed(i)
      rst <-
        single_simulation_lst(
          gss_data = gss_data_5class_fc_edu,
          size = 5,
          age_l = 25,
          age_u = 55,
          year_l = 1920,
          year_u = 1997,
          est_l = "1945",
          est_u = "1990",
          tol = 0.01,
          gamma_par = 1)
    },
    mc.cores = n_cores
  )

# valid_boot() requires the cohorts passing the 5% non-primitive screen
# to form a contiguous range; in the sparse father-child sample an interior
# cohort can fail the screen while its neighbours pass. Restrict the
# bootstrap to the longest contiguous block of passing cohorts before
# the validation:

cohort_names <- names(boot_rst_lst[[1]])

non_primitive_rate <-
  sapply(
    cohort_names,
    function(cohort){
      mean(
        sapply(
          boot_rst_lst,
          function(boot){
            length(boot[[cohort]]) == 1 && all(is.na(boot[[cohort]]))
          }
        )
      )
    }
  )

valid_flag <- non_primitive_rate < 0.05
run_lst <- rle(valid_flag)
run_end <- cumsum(run_lst$lengths)
run_start <- run_end - run_lst$lengths + 1
valid_runs <- which(run_lst$values)

assert_that(
  length(valid_runs) > 0
)

longest_run <- valid_runs[which.max(run_lst$lengths[valid_runs])]
cohort_keep <- cohort_names[run_start[longest_run]:run_end[longest_run]]

print(paste("Cohorts kept for the bootstrap validation:",
            min(cohort_keep), "to", max(cohort_keep)))

boot_rst_lst <-
  lapply(
    boot_rst_lst,
    function(boot){
      boot[cohort_keep]
    }
  )

valid_boot_rst <-
  valid_boot(
    boot_rst = boot_rst_lst,
    valid_rate = 0.05,
    boot_times = times,
    valid_times = 2000
  )

# Keep only the bootstrap replicates that survived the validation:
# valid_boot() subsets to valid_times replicates and pads with NULL when
# fewer remain, which the sparse father-child sample can trigger.

valid_boot_rst$boot_rst <-
  Filter(Negate(is.null), valid_boot_rst$boot_rst)

print(paste("The number of valid bootstrap replicates is",
            length(valid_boot_rst$boot_rst)))

gc()


boot_rst <-
  process_results(
    boot_rst = valid_boot_rst$boot_rst,
    times = times,
    size = size,
    n_cores = n_cores,
    reshape = TRUE
  )

print("Bootstrap result of father-child sample has finished!")



#-------------------------------------------------#
#           CALCULATE BIAS-CORRECTED RESULTS      #
#-------------------------------------------------#

# Align the baseline cohorts with the cohorts surviving the bootstrap
# validation. Both lists are contiguous but may cover different ranges
# (valid_boot() drops sparse bootstrap cohorts; the baseline itself may
# drop early cohorts), so the bias correction is restricted to their
# common cohort range, computed from the data rather than hardcoded:

valid_cohort_vec <- as.integer(names(valid_boot_rst$boot_rst[[1]]))

assert_that(
  length(boot_rst) == length(valid_cohort_vec)
)

common_cohort_vec <- intersect(valid_cohort_vec, baseline_cohort_vec)

assert_that(
  all(diff(common_cohort_vec) == 1)
)

print(paste("The bias-corrected result covers birth cohorts",
            min(common_cohort_vec), "to", max(common_cohort_vec)))

bias_corrected_rst <-
  generate_validation_data(
    process_boot_rst_lst = boot_rst[match(common_cohort_vec, valid_cohort_vec)],
    process_baseline_lst = baseline_rst[match(common_cohort_vec, baseline_cohort_vec)]
  )

write_rds(
  bias_corrected_rst,
  file.path(
    dir_data,
    "5class_fc_edu_bc.rds")
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved processed data to: ", dir_data, "\n")
cat("Father-child baseline and bias-corrected estimation have been completed successfully.\n")
cat("Bias-corrected results cover birth cohorts ",
    min(common_cohort_vec), "to", max(common_cohort_vec), "\n")

sink()
