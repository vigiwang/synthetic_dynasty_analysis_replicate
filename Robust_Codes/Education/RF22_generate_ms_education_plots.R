#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "robustness_check"
subsubsection <- "education"
title <- "plot_figF28_F30"

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
dir_plot <- paste0(dir_root, "/" ,section, "/", subsection, "/", subsubsection)
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
# Script:      RF22_generate_ms_education_plots.R
#
# Inputs:
#   - Mother-son education bias-corrected results:
#       Data/robustness_check/5class_ms_edu_bc.rds
#   - Mother-son educational pairs (5-category years of schooling):
#       Data/robustness_check/gss_ms_edu_5class.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/robustness_check_plot_figF28_F30_log.txt
#   - Plot:
#       Plot/robustness_check/education/figF28_ms_edu_AIM.rds
#       Plot/robustness_check/education/figF29_ms_edu_SSM.rds
#       Plot/robustness_check/education/figF30_ms_edu_bin5_LML.rds

# Description:
#   This script creates Figure F28: Aggregate intergenerational memory at
#   generation t = 1; Figure F29: Steady-state mobility (overall
#   mobility at the terminal synthetic generation, with bias-corrected
#   confidence intervals); Figure F30:
#   Log-multiplicative layer effect, for the mother-son educational
#   mobility sample.
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
#                LOAD GSS DATA
#-------------------------------------------------#

gss_data_bc <- readr::read_rds(
  file.path(
    dir_data,
    "5class_ms_edu_bc.rds"
  )
)$validation_lst


gss_data_main <- readr::read_rds(
  file.path(
    dir_data,
    "gss_ms_edu_5class.rds"
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


cohort_bins <-
  tibble(
    cohort_index = seq(1:length(gss_data_bc)),
    cohort_bins_vec = seq(1945,1990)
  )
cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)
year_vec <- as.character(seq(1945,1990,2))
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

#=======================#
#      AIM
#=======================#


AIM <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "memory_df",
    attribute_col = "class",
    metric_name = "AIM",
    y_title = "Memory",
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.11, 0.43, by = 0.02),
    x_range = c(1944,1991),
    y_range = c(0.11,0.43),
    y_tick_text = to_text(seq(0.11, 0.43, by = 0.02))
  ) %>%
  layout(
    width = 1000,
    height = 600
  )


#=======================#
#      SSM
#=======================#

# Steady-state mobility = overall mobility at the terminal synthetic
# generation (the distribution has converged under the production stop
# rule), taken from movement_df with bias-corrected CIs:

SSM <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = seq(1, 10),
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Historical Mobility",
    y_title = "Steady state mobility",
    steady_state = TRUE,
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.49, 0.75, by = 0.02),
    x_range = c(1944,1991),
    y_range = c(0.49,0.75),
    y_tick_text = to_text(seq(0.49, 0.75, by = 0.02))
  ) %>%
  layout(
    width = 1000,
    height = 600
  )


#===========================================#
#   Bin = 5 log-multiplicative layer effect
#==========================================#

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
#  Save the Plot:
#-----------------------------#

saveRDS(
  AIM,
  file = path.expand(
    file.path(
      dir_plot,
      "figF28_ms_edu_AIM.rds"
    )
  )
)

saveRDS(
  SSM,
  file = path.expand(
    file.path(
      dir_plot,
      "figF29_ms_edu_SSM.rds"
    )
  )
)

saveRDS(
  bin5_LML,
  file = path.expand(
    file.path(
      dir_plot,
      "figF30_ms_edu_bin5_LML.rds"
    )
  )
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure H4-H3 have been successfully replicated.\n")

sink()
