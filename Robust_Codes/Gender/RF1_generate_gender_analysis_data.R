#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "robustness_check"
title <- "gender"

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
# Script:      RF1_generate_gender_analysis_data.R
#
# Inputs:
#   - GSS microdata (via the {gssr} package): gssr::data(gss_all)
#   - Occupational crosswalk (Morgan, 2017):
#       Data/occ10-to-egp-class-crosswalk.csv
#
# Outputs:
#   - Log file:
#       code/_LOGS/Data/robustness_check_gender_log.txt
#   - Analytical dataset:
#       - Data/robustness_check/5class_mc_bc.rds
#       - Data/robustness_check/5class_fd_bc.rds
#       - Data/robustness_check/5class_fs_bc.rds
#       - Data/robustness_check/5class_pc_bc.rds
#
# Description:
#   This script constructs the main analytical sample for class typology
#   robustness check from the General Social Survey (GSS). It codes both respondents' 
#   and their parents' occupations into the 2010 Census Occupational Classification (COC). 
#   These COC codes are then mapped to the 5-class EGP typology using the Morgan (2017) 
#   occupational crosswalk.
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
#     GENERATE MOTHER-CHILD PAIRS:
#-------------------------------------------------#

gss_mc_occ10_egp <-
  gss_all %>%
  mutate(
    uniqid = row_number()
  ) %>%
  dplyr::select(
    year,
    maocc10, # mothers' occ10 codes
    occ10, # respondents' occ10 codes
    age,
    cohort,
    uniqid
  ) %>%
  filter(
    !is.na(occ10) & !is.na(maocc10) # only keep pairs with occupation information
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
    by.x = "maocc10",
    by.y = "occ10"
  ) %>%
  rename(
    status_p = egp10_10,
    label_p = egp_label
  ) 

gss_mc_occ10_5class <-
  gss_mc_occ10_egp %>%
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

#-------------------------------------------------#
#    PARENT OF HIGHER STATUS - CHILD PAIRS:
#-------------------------------------------------#

# Code both parents' occupation information into EGP class:

gss_pc_occ10_egp <-
  gss_all %>%
  mutate(
    uniqid = row_number()
  ) %>%
  dplyr::select(
    year,
    maocc10, # mothers' occ10 codes
    paocc10, # father's occ10 codes
    occ10, # respondents' occ10 codes
    age,
    cohort,
    uniqid
  ) %>%
  filter(
    !is.na(occ10) & !is.na(maocc10) & !is.na(paocc10) # only keep pairs with occupation information
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
    status_f = egp10_10,
    label_f = egp_label
  ) %>%
  merge(
    crosswalk %>%
      dplyr::select(
        -title
      ),
    by.x = "maocc10",
    by.y = "occ10"
  ) %>%
  rename(
    status_m = egp10_10,
    label_m = egp_label
  ) 

# Use 5-class typology:

gss_pc_occ10_5class <-
  gss_pc_occ10_egp %>%
  filter(
    status_c %in% c(1, 2, 3, 4, 7, 8, 9, 10, 11) & 
      status_f %in% c(1, 2, 3, 4, 7, 8, 9, 10, 11) &
      status_m %in% c(1, 2, 3, 4, 7, 8, 9, 10, 11) 
  ) %>%
  mutate(
    status_c = case_when(
      status_c == 1 ~ 1,
      status_c == 2 ~ 2,
      status_c %in% c(3,8) ~ 3,
      status_c %in% c(9,7) ~ 4,
      status_c %in% c(4,10,11) ~ 5
    ),
    status_f = case_when(
      status_f == 1 ~ 1,
      status_f == 2 ~ 2,
      status_f %in% c(3,8) ~ 3,
      status_f %in% c(9,7) ~ 4,
      status_f %in% c(4,10,11) ~ 5
    ),
    status_m = case_when(
      status_m == 1 ~ 1,
      status_m == 2 ~ 2,
      status_m %in% c(3,8) ~ 3,
      status_m %in% c(9,7) ~ 4,
      status_m %in% c(4,10,11) ~ 5
    )
  ) %>%
  mutate(
    status_p = ifelse(
      status_f < status_m,
      status_f,
      status_m
    )
  )

#-------------------------------------------------#
#             ONLY KEEP FATHER-SON PAIRS:
#-------------------------------------------------#

gss_fc_occ10_5class <- readr::read_rds(
  file.path(
    dir_root,
    "Data",
    "main_results",
    "gss_fc_occ10_5class.rds"
  )
)

gss_fs_occ10_5class <-
  gss_fc_occ10_5class %>%
  left_join(
    gss_all %>%
      mutate(
        uniqid = row_number()
      ) %>%
      dplyr::select(
        uniqid,
        sex
      ),
    by = "uniqid"
  ) %>%
  filter(
    sex == 1
  )


#-------------------------------------------------#
#           ONLY KEEP FARTHER-DAUGHTER PAIRS:
#-------------------------------------------------#

gss_fd_occ10_5class <-
  gss_fc_occ10_5class %>%
  left_join(
    gss_all %>%
      mutate(
        uniqid = row_number()
      ) %>%
      dplyr::select(
        uniqid,
        sex
      ),
    by = "uniqid"
  ) %>%
  filter(
    sex == 2
  )




#-----------------------------#
#  Save the Data:
#-----------------------------#

write_rds(
  gss_fs_occ10_5class,
  file.path(
    dir_data,
    "gss_fs_occ10_5class.rds")
)

write_rds(
  gss_fd_occ10_5class,
  file.path(
    dir_data,
    "gss_fd_occ10_5class.rds")
)

write_rds(
  gss_mc_occ10_5class,
  file.path(
    dir_data,
    "gss_mc_occ10_5class.rds")
)


write_rds(
  gss_pc_occ10_5class,
  file.path(
    dir_data,
    "gss_pc_occ10_5class.rds")
)

# Open log：
sink(log_path, split = TRUE)

cat("Saved processed data to: ", dir_data, "\n")
cat("Occupational pairs under alternative gender related restrictions
    for robustness check has been generated.\n")

sink()









