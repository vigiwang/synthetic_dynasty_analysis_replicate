#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "robustness_check"
title <- "education"

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
# Script:      RF15_generate_education_analysis_data.R
#
# Inputs:
#   - GSS microdata (via the {gssr} package): gssr::data(gss_all)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Data/robustness_check_education_log.txt
#   - Analytical dataset:
#       - Data/robustness_check/gss_fc_edu_5class.rds
#       - Data/robustness_check/gss_mc_edu_5class.rds
#       - Data/robustness_check/gss_ms_edu_5class.rds
#       - Data/robustness_check/gss_md_edu_5class.rds
#       - Data/robustness_check/gss_pc_edu_5class.rds
#
# Description:
#   This script constructs the analytical samples for the educational
#   mobility robustness check from the General Social Survey (GSS). It
#   codes respondents' and their parents' highest completed years of
#   schooling (educ, paeduc, maeduc) into a 5-category typology:
#
#       Category 1: more than 16 years  (>16)
#       Category 2: exactly 16 years    (16)
#       Category 3: 13 to 15 years      (13-15)
#       Category 4: exactly 12 years    (12)
#       Category 5: fewer than 12 years (<12)
#
#   Category 1 is the highest level, mirroring the orientation of the
#   5-class EGP typology (Class 1 = highest) so that all mobility
#   measures keep the same interpretation across the occupational and
#   educational analyses.
#
#   Five dyad samples are generated: father-child, mother-child,
#   mother-son, mother-daughter, and highest-educated-parent-child.
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
#     LOAD GSS DATA
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


#-------------------------------------------------#
#     CODE YEARS OF SCHOOLING INTO 5 CATEGORIES
#-------------------------------------------------#

# Category 1 (>16) is the highest level, mirroring the 5-class EGP
# typology where Class 1 is the highest class:

gss_edu <-
  gss_all %>%
  mutate(
    uniqid = row_number()
  ) %>%
  dplyr::select(
    year,
    educ,   # respondents' highest year of school completed
    paeduc, # fathers' highest year of school completed
    maeduc, # mothers' highest year of school completed
    sex,
    age,
    cohort,
    uniqid
  ) %>%
  mutate(
    educ = as.numeric(zap_labels(educ)),
    paeduc = as.numeric(zap_labels(paeduc)),
    maeduc = as.numeric(zap_labels(maeduc)),
    year = as.numeric(zap_labels(year)),
    age = as.numeric(zap_labels(age)),
    cohort = as.numeric(zap_labels(cohort)),
    sex = as.numeric(zap_labels(sex))
  ) %>%
  filter(
    !is.na(cohort) & cohort != 9999 # only keep children with birth cohort information
  ) %>%
  mutate(
    status_c = case_when(
      educ > 16 ~ 1,
      educ == 16 ~ 2,
      educ >= 13 & educ <= 15 ~ 3,
      educ == 12 ~ 4,
      educ < 12 ~ 5
    ),
    label_c = case_when(
      status_c == 1 ~ ">16",
      status_c == 2 ~ "16",
      status_c == 3 ~ "13-15",
      status_c == 4 ~ "12",
      status_c == 5 ~ "<12"
    ),
    status_f = case_when(
      paeduc > 16 ~ 1,
      paeduc == 16 ~ 2,
      paeduc >= 13 & paeduc <= 15 ~ 3,
      paeduc == 12 ~ 4,
      paeduc < 12 ~ 5
    ),
    label_f = case_when(
      status_f == 1 ~ ">16",
      status_f == 2 ~ "16",
      status_f == 3 ~ "13-15",
      status_f == 4 ~ "12",
      status_f == 5 ~ "<12"
    ),
    status_m = case_when(
      maeduc > 16 ~ 1,
      maeduc == 16 ~ 2,
      maeduc >= 13 & maeduc <= 15 ~ 3,
      maeduc == 12 ~ 4,
      maeduc < 12 ~ 5
    ),
    label_m = case_when(
      status_m == 1 ~ ">16",
      status_m == 2 ~ "16",
      status_m == 3 ~ "13-15",
      status_m == 4 ~ "12",
      status_m == 5 ~ "<12"
    )
  )


#-------------------------------------------------#
#     GENERATE FATHER-CHILD PAIRS:
#-------------------------------------------------#

gss_fc_edu_5class <-
  gss_edu %>%
  filter(
    !is.na(status_c) & !is.na(status_f) # only keep pairs with education information
  ) %>%
  mutate(
    status_p = status_f,
    label_p = label_f
  ) %>%
  dplyr::select(
    year,
    educ,
    paeduc,
    age,
    cohort,
    uniqid,
    label_c,
    status_c,
    label_p,
    status_p
  )


#-------------------------------------------------#
#     GENERATE MOTHER-CHILD PAIRS:
#-------------------------------------------------#

gss_mc_edu_5class <-
  gss_edu %>%
  filter(
    !is.na(status_c) & !is.na(status_m) # only keep pairs with education information
  ) %>%
  mutate(
    status_p = status_m,
    label_p = label_m
  ) %>%
  dplyr::select(
    year,
    educ,
    maeduc,
    sex,
    age,
    cohort,
    uniqid,
    label_c,
    status_c,
    label_p,
    status_p
  )


#-------------------------------------------------#
#             ONLY KEEP MOTHER-SON PAIRS:
#-------------------------------------------------#

gss_ms_edu_5class <-
  gss_mc_edu_5class %>%
  filter(
    sex == 1
  )


#-------------------------------------------------#
#           ONLY KEEP MOTHER-DAUGHTER PAIRS:
#-------------------------------------------------#

gss_md_edu_5class <-
  gss_mc_edu_5class %>%
  filter(
    sex == 2
  )


#-------------------------------------------------#
#    PARENT OF HIGHER EDUCATION - CHILD PAIRS:
#-------------------------------------------------#

gss_pc_edu_5class <-
  gss_edu %>%
  filter(
    !is.na(status_c) & !is.na(status_f) & !is.na(status_m) # only keep pairs with education information
  ) %>%
  mutate(
    status_p = ifelse(
      status_f < status_m,
      status_f,
      status_m
    ),
    label_p = ifelse(
      status_f < status_m,
      label_f,
      label_m
    )
  ) %>%
  dplyr::select(
    year,
    educ,
    paeduc,
    maeduc,
    age,
    cohort,
    uniqid,
    label_c,
    status_c,
    label_p,
    status_p
  )


#-----------------------------#
#  Save the Data:
#-----------------------------#

write_rds(
  gss_fc_edu_5class,
  file.path(
    dir_data,
    "gss_fc_edu_5class.rds")
)

write_rds(
  gss_mc_edu_5class,
  file.path(
    dir_data,
    "gss_mc_edu_5class.rds")
)

write_rds(
  gss_ms_edu_5class,
  file.path(
    dir_data,
    "gss_ms_edu_5class.rds")
)

write_rds(
  gss_md_edu_5class,
  file.path(
    dir_data,
    "gss_md_edu_5class.rds")
)

write_rds(
  gss_pc_edu_5class,
  file.path(
    dir_data,
    "gss_pc_edu_5class.rds")
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved processed data to: ", dir_data, "\n")
cat("Educational pairs under the 5-category years-of-schooling typology
    for robustness check has been generated.\n")

sink()
