#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "robustness_check"
title <- "alternative_class_typology_generation"

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
# Script:      RB1_generate_class_typology_data.R
#
# Inputs:
#   - GSS microdata (via the {gssr} package): gssr::data(gss_all)
#   - Occupational crosswalk (Morgan, 2017):
#       Data/occ10-to-egp-class-crosswalk.csv
#
# Outputs:
#   - Log file:
#       code/_LOGS/Data/robustness_check_alternative_class_typology_generation_log.txt
#   - Analytical dataset:
#       - Data/robustness_check/gss_fc_occ10_6class.rds
#       - Data/robustness_check/gss_fc_occ10_7class.rds
#       - Data/robustness_check/gss_fc_occ10_5class_noself.rds
#       - Data/robustness_check/gss_fc_occ10_5class_nosfarmer.rds
#
# Description:
#   This script constructs the analytical sample for alternative class typology
#   robustness check from the General Social Survey (GSS). It codes both respondents' 
#   and their fathers' occupations into the 2010 Census Occupational Classification (COC). These COC codes are then 
#   mapped to the 6-class EGP typology(Pfeffer and Hertel (2015)), 7-class EGP typology(Morgan (2017)), benchmark 
#   5-class typology without self-employed individuals(Erikson and Goldthorpe (2010)),benchmark 5-class typology 
#   without farmers, using the Morgan (2017) occupational crosswalk.
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
#     IDENTIFY THE SELF_EMPLOYED INDIVIDUALS
#-------------------------------------------------#

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
    wrkstat,
    wrkslf,
    pawrkslf
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

gss_fc_occ10_egp_self <-
  gss_fc_occ10_egp %>%
  mutate(
    status_c_11 = case_when(
      !status_c %in% c(1,2,7,12,NA) & wrkslf == 1 ~ 5,
      .default = status_c
    ),
    status_p_11 = case_when(
      !status_p %in% c(1,2,7,12,NA) & pawrkslf == 1 ~ 5,
      .default = status_p
    )
  )

#-------------------------------------------------#
#             CODE 7-CLASS TYPOLOGY:
#-------------------------------------------------#

gss_fc_occ10_7class <-
  gss_fc_occ10_egp_self %>%
  filter(
    status_c %in% c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11) & 
      status_p %in% c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)
  ) %>%
  mutate(
    status_c = case_when(
      status_c_11 == 1 ~ 1,
      status_c_11 == 2 ~ 2,
      status_c_11 %in% c(3,4) ~ 3,
      status_c_11 %in% c(5,7) ~ 4,
      status_c_11 == 8 ~ 5,
      status_c_11 == 9 ~ 6,
      status_c_11 %in% c(10,11) ~ 7
    ),
    status_p = case_when(
      status_p_11 == 1 ~ 1,
      status_p_11 == 2 ~ 2,
      status_p_11 %in% c(3,4) ~ 3,
      status_p_11 %in% c(5,7) ~ 4,
      status_p_11 == 8 ~ 5,
      status_p_11 == 9 ~ 6,
      status_p_11 %in% c(10,11) ~ 7
    )
  )


#-------------------------------------------------#
#             CODE 6-CLASS TYPOLOGY:
#-------------------------------------------------#

gss_fc_occ10_6class <-
  gss_fc_occ10_egp_self %>%
  filter(
    status_c_11 %in% c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11) &   
      status_p_11 %in% c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11)
  ) %>%
  mutate(
    status_c_11 = case_when(
      status_c_11 == 1 ~ 1,
      status_c_11 == 2 ~ 2,
      status_c_11 %in% c(3,4) ~ 3,
      status_c_11 %in% c(5,7) ~ 4,
      status_c_11 %in% c(8,9) ~ 5,
      status_c_11 %in% c(10,11) ~ 6
    ),
    status_p_11 = case_when(
      status_p_11 == 1 ~ 1,
      status_p_11 == 2 ~ 2,
      status_p_11 %in% c(3,4) ~ 3,
      status_p_11 %in% c(5,7) ~ 4,
      status_p_11 %in% c(8,9) ~ 5,
      status_p_11 %in% c(10,11) ~ 6
    )
  ) %>%
  dplyr::select(
    -status_c,
    -status_p
  ) %>%
  rename(
    status_c = status_c_11,
    status_p = status_p_11
  )

#-------------------------------------------------#
#             REMOVE ALL SELF-EMPLOYED:
#-------------------------------------------------#

gss_fc_occ10_5class_noself <-
  gss_fc_occ10_egp_self %>%
  filter(
    status_c_11 %in% c(1, 2, 3, 4, 8, 9, 10, 11) & 
      status_p_11 %in% c(1, 2, 3, 4, 8, 9, 10, 11)
  ) %>%
  mutate(
    status_c = case_when(
      status_c_11 == 1 ~ 1,
      status_c_11== 2 ~ 2,
      status_c_11 %in% c(3,8) ~ 3,
      status_c_11 %in% c(9) ~ 4,
      status_c_11 %in% c(4,10,11) ~ 5
    ),
    status_p = case_when(
      status_p_11 == 1 ~ 1,
      status_p_11 == 2 ~ 2,
      status_p_11 %in% c(3,8) ~ 3,
      status_p_11 %in% c(9) ~ 4,
      status_p_11 %in% c(4,10,11) ~ 5
    )
  )

assert_that(
  all(
    gss_fc_occ10_5class_noself %>%
      filter(
        wrkslf == 1
      ) %>%
      pull(
        status_c
      )) %in%
    c(1,2),
  
  all(
    gss_fc_occ10_5class_noself %>%
      filter(
        pawrkslf == 1
      ) %>%
      pull(
        status_p
      )) %in%
    c(1,2)
  
)


#-------------------------------------------------#
#             REMOVE FARMERS ONLY:
#-------------------------------------------------#

gss_fc_occ10_5class_nofarmer <-
  gss_fc_occ10_egp %>%
  filter(
    status_c %in% c(1, 2, 3, 4, 8, 9, 10, 11) & 
      status_p %in% c(1, 2, 3, 4, 8, 9, 10, 11)
  ) %>%
  mutate(
    status_c = case_when(
      status_c == 1 ~ 1,
      status_c == 2 ~ 2,
      status_c %in% c(3,8) ~ 3,
      status_c %in% c(9) ~ 4,
      status_c %in% c(4,10,11) ~ 5
    ),
    status_p = case_when(
      status_p == 1 ~ 1,
      status_p == 2 ~ 2,
      status_p %in% c(3,8) ~ 3,
      status_p %in% c(9) ~ 4,
      status_p %in% c(4,10,11) ~ 5
    )
  )
assert_that(
  sum(gss_fc_occ10_5class_nofarmer$label_p == "IVc") == 0 &
    sum(gss_fc_occ10_5class_nofarmer$label_c == "IVc") == 0
)


#-----------------------------#
#  Save the Data:
#-----------------------------#

write_rds(
  gss_fc_occ10_6class,
  file.path(
    dir_data,
    "gss_fc_occ10_6class.rds")
)

write_rds(
  gss_fc_occ10_7class,
  file.path(
    dir_data,
    "gss_fc_occ10_7class.rds")
)

write_rds(
  gss_fc_occ10_5class_nofarmer,
  file.path(
    dir_data,
    "gss_fc_occ10_5class_nofarmer.rds")
)

write_rds(
  gss_fc_occ10_5class_noself,
  file.path(
    dir_data,
    "gss_fc_occ10_5class_noself.rds")
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved processed data to: ", dir_data, "\n")
cat("Father-Child occupational pair under alternative class typologies
    for robustness check has been generated.\n")

sink()









