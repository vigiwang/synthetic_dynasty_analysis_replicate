#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "main_results"
title <- "occ_pairs_generation"

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
# Repository:  https://github.com/synthetic_dynasty_analysis
#
# Script:      00_generate_main_results_data.R
#
# Inputs:
#   - GSS microdata (via the {gssr} package): gssr::data(gss_all)
#   - Occupational crosswalk (Morgan, 2017):
#       data/raw/occ10-to-egp-class-crosswalk.csv
#
# Outputs:
#   - Log file:
#       code/main_results/_LOGS/data_generation.txt
#   - Analytical dataset:
#       data/main_results/raw/gss_fc_occ10_5class.rds
#
# Description:
#   This script constructs the main analytical sample from the General Social
#   Survey (GSS). It codes both respondents' and their fathers' occupations into 
#   the 2010 Census Occupational Classification (COC). These COC codes are then 
#   mapped to the 5-class EGP typology using the Morgan (2017) occupational 
#   crosswalk.
#
# Notes:
#   - The GSS input is accessed programmatically through {gssr}. For exact
#     replication, ensure the installed {gssr} version (and its bundled GSS
#     release) matches the version used in the analysis. 
#
#-------------------------------------------------------------------------------

#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#
if(HPC){
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("tidyverse")
  library("readr")
  library("haven")
  library("dplyr")
}else{
  packages <-
    c(
      "tidyverse",
      "readr",
      "dplyr",
      "haven"
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
#     LOAD GSS DATA AND OCCUPATIONAL CROSSWALK
#-------------------------------------------------#

# Load GSS microdata (via the {gssr} package)
# Install from r-universe if not already available
if(HPC){
  gss_all <- read_dta(
    file.path(
      dir_root,
      "Data",
      "gss7224_r2.dta"
    )
  )
  
}else{
  if (!requireNamespace("gssr", quietly = TRUE)) {
    install.packages(
      "gssr",
      repos = c(
        "https://kjhealy.r-universe.dev",
        "https://cloud.r-project.org"
      )
    )
  }
  library(gssr)
  data(gss_all)
}


# Load the 2010 Census Occupation → 5-Class EGP crosswalk (Morgan, 2017)

crosswalk <- readr::read_csv(
  file.path(
    dir_root,
    "Data",
    "occ10-to-egp-class-crosswalk.csv"
  )
)

#-------------------------------------------------#
#     CALCULATE THE 5-CLASS TYPOLOGY
#-------------------------------------------------#

# Link the occ10 code with EGP class:

gss_fc_occ10_egp <-
  gss_all %>%
  mutate(
    uniqid = row_number()
  ) %>%
  dplyr::select(
    year,
    paocc10, # fathers' occ10 codes
    occ10, # respondents' occ10 codes
    age,
    cohort,
    uniqid,
    sample,
    wtssps,
    race
  ) %>%
  filter(
    !is.na(occ10) & !is.na(paocc10) # only keep pairs with occupation information
  ) %>%
  filter(
    !is.na(cohort) & cohort != 9999 # only keep children with birth cohort information
  ) %>%
  merge(
    crosswalk %>%
      dplyr::select(
        -title
      ),
    by.x = "occ10",
    by.y = "occ10"
  ) %>%
  rename(
    status_c = egp10_10,
    label_c = egp_label
  ) %>%
  merge(
    crosswalk %>%
      dplyr::select(
        -title
      ),
    by.x = "paocc10",
    by.y = "occ10"
  ) %>%
  rename(
    status_p = egp10_10,
    label_p = egp_label
  )

# Code the EGP class into 5 Class Typology:

gss_fc_occ10_5class <-
  gss_fc_occ10_egp %>%
  filter(
    status_c %in% c(1, 2, 3, 4, 7, 8, 9, 10, 11) & 
    status_p %in% c(1, 2, 3, 4, 7, 8, 9, 10, 11)
  ) %>%
  mutate(
    status_c = case_when(
      status_c == 1 ~ 1,
      status_c == 2 ~ 2,
      status_c %in% c(3,8) ~ 3,
      status_c %in% c(9,7) ~ 4,
      status_c %in% c(4,10,11) ~ 5
    ),
    status_p = case_when(
      status_p == 1 ~ 1,
      status_p == 2 ~ 2,
      status_p %in% c(3,8) ~ 3,
      status_p %in% c(9,7) ~ 4,
      status_p %in% c(4,10,11) ~ 5
    )
  )

write_rds(
  gss_fc_occ10_5class,
  file.path(
    dir_data,
    "gss_fc_occ10_5class.rds"),
  compress = "xz"
  )

print("5-Class EGP Typology Father-Child Pairs have been created!")

# Open log：
sink(log_path, split = TRUE)

width_curr <- getOption("width")
options(width = 300)
print(
  head(gss_fc_occ10_5class) %>%
    dplyr::select(
      occ10,
      paocc10,
      status_c,
      status_p,
      age,
      cohort
    )
  )

cat("Saved processed data to: ", dir_data, "\n")
cat("Main results raw data generation completed successfully.\n")
sink()









