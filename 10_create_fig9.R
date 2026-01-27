#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "main_results"
title <- "plot_fig9"

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
dir_data <- paste0(dir_root, "/" ,"Data", "/", subsection)
dir_plot <- paste0(dir_root, "/" ,section, "/", subsection)
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
create_dir_if_missing(dir_plot)
options(warn = -1)


#-------------------------------------------------------------------------------
# U.S. Occupational Mobility — Replication Files
#
# Project:     U.S. Occupational Mobility Analysis
# Repository:  https://github.com/synthetic_dynasty_analysis_replicate
#
# Script:      10_create_fig9.R
#
# Inputs:
#   - GSS father–child occupational pairs (5-class EGP):
#       Data/main_results/gss_fc_occ10_5class.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/main_results_plot_fig9_log.txt
#   - Plot:
#       Plot/main_results/fig9_bin5_LML.rds
#       Plot/main_results/fig9_bin10_LML.rds
#
# Description:
#   This script recreate Figure 9: Log-multiplicative layer effects
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#

if(HPC){
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("tidyverse")
  library("dplyr")
  library("plotly")
  library("parallel")
  library("logmult")
}else{
  packages <-
    c(
      "tidyverse",
      "dplyr",
      "plotly",
      "parallel",
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
#           LOAD BIAS_CORRECTED RESULTS           #
#-------------------------------------------------#

gss_data_main <- readr::read_rds(
  file.path(
    dir_data,
    "gss_fc_occ10_5class.rds"
  )
)
#-------------------------------------------------#
#          SOURCE FUNCTIONS FOR ANALYSIS          #
#-------------------------------------------------#

# Load helper functions:
source(
  file.path(
    dir_root,
    "Functions",
    "utils.R"
  )
)

# Load plot functions:
source(
  file.path(
    dir_root,
    "Functions",
    "plot.R"
  )
)

bins_vec_5 <- c(
  "[1945,1950)",
  "[1950,1955)",
  "[1955,1960)",
  "[1960,1965)",
  "[1965,1970)",
  "[1970,1975)",
  "[1975,1980)",
  "[1980,1985)",
  "[1985,1990)"
)

bins_vec_10 <- c(
  "[1945,1954]",
  "[1955,1964]",
  "[1965,1974]",
  "[1975,1990]"
)

#---------------------------------------------------------------#
#     PLOT LOG-MULTIPLICATIVE LAYER EFFECT WITH BINNED DATA     #
#--------------------------------------------------------------#

#-----------------------------#
#  Bin = 5:
#-----------------------------#

bin5_LML <-
  plot_log_multiplicative(
    data = gss_data_main,
    bin_width = 5,
    year_l = 1945,
    year_u = 1990,
    age_l = 25,
    age_u = 55,
    cohort_bins = bins_vec_5,
    tick_vals = seq(0.55,1.25,0.1),
    range = c(0.55, 1.25)
    )

#-----------------------------#
#  Bin = 10:
#-----------------------------#

bin10_LML <-
  plot_log_multiplicative(
    data = gss_data_main,
    bin_width = 10,
    year_l = 1945,
    year_u = 1990,
    age_l = 25,
    age_u = 55,
    cohort_bins = bins_vec_10,
    tick_vals = seq(0.55,1.25,0.1),
    range = c(0.55, 1.25)
  )

#-----------------------------#
#  Save the Plot:
#-----------------------------#
# 1. Save Bin 5 LML
saveRDS(
  bin5_LML,
  file = path.expand(file.path(dir_plot, "fig9_bin5_LML.rds"))
)

# 2. Save Bin 10 LML
saveRDS(
  bin10_LML,
  file = path.expand(file.path(dir_plot, "fig9_bin10_LML.rds"))
)

# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure 9 of the paper’s main results have been successfully replicated.\n")

sink()







