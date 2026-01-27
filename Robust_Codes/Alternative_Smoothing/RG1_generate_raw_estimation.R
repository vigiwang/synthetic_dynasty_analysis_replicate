#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "robustness_check"
title <- "raw_estimation"

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
# Script:      RG1_generate_raw_estimation.R
#
# Inputs:
#   - GSS father–child occupational pairs (5-class EGP):
#       Data/main_results/gss_fc_occ10_5class.rds
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
#       code/_LOGS/Data/robustness_check_raw_estimation_log.txt
#   - Estimation output:
#       Data/robustness_check/raw_bc.rds
#
# Description:
#   This script generates the robustness results using no-smoothing father–child
#   occupational mobility pairs constructed from the General Social Survey (GSS).
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

gss_data_main <- readr::read_rds(
  file.path(
    dir_root,
    "Data",
    "main_results",
    "gss_fc_occ10_5class.rds"
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


#-------------------------------------------------#
#           CALCULATE BASELINE RESULTS            #
#-------------------------------------------------#

baseline_rst <- 
  calculate_baseline_newform(
    gss_data = gss_data_main,
    size     = 5,     # Number of classes in the EGP typology
    age_l    = 25,    # Lower bound of respondents' age at survey
    age_u    = 55,    # Upper bound of respondents' age at survey
    year_l   = 1945,  # Lower bound of birth cohorts used for GAM training
    # (wider than the main sample to capture more information)
    year_u   = 1990,  # Upper bound of birth cohorts used for GAM training
    # (wider than the main sample to capture more information)
    bin_width = NULL,
    method = "raw"
  )

print("Baseline result of raw sample has finished!")

#-------------------------------------------------#
#     CALCULATE AND PROCESS BOOTSTRAP RESULTS     #
#-------------------------------------------------#

boot_rst_lst <- 
  mclapply(
    1:times,
    function(i){
      set.seed(i)
      rst <-
        single_simulation_lst_newform(
          gss_data = gss_data_main, 
          age_l = 25, 
          age_u = 55, 
          year_l = 1945, 
          year_u = 1990, 
          boot = TRUE,
          bin_width = NULL,
          method = "raw"
        )
    },
    mc.cores = n_cores
  )


valid_boot_rst <-
  valid_boot(
    boot_rst = boot_rst_lst,
    valid_rate = 0.05,
    boot_times = times,
    valid_times = 2000
  )

# Process the bootstrap results:

boot_rst <- 
  process_results(
    boot_rst = valid_boot_rst$boot_rst,
    times = times,
    size = size,
    n_cores = n_cores,
    reshape = TRUE
  )

print("Bootstrap result of raw sample has finished!")

#-------------------------------------------------#
#           CALCULATE BIAS-CORRECTED RESULTS      #
#-------------------------------------------------#


bias_corrected_rst <-
  generate_validation_data(
    process_boot_rst_lst = boot_rst, 
    process_baseline_lst = baseline_rst
  )

write_rds(
  bias_corrected_rst,
  file.path(
    dir_data,
    "raw_bc.rds")
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved processed data to: ", dir_data, "\n")
cat("Without smoothing baseline and bias-corrected estimation  have been completed successfully.\n")

sink()


