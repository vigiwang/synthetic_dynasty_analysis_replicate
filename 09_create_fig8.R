#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "main_results"
title <- "plot_fig8"

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
# Script:      09_create_fig8.R
#
# Inputs:
#   - GSS father–child occupational pairs (5-class EGP):
#       Data/main_results/main_rst_bc.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/main_results_plot_fig8_log.txt
#   - Plot:
#       Plot/main_results/fig8_IM_1945_t.rds
#       Plot/main_results/fig8_IM_1990_t.rds
#       Plot/main_results/fig8_AIM_1945_t.rds
#       Plot/main_results/fig8_AIM_1990_t.rds
#
# Description:
#   This script recreate Figure 8: Memory across synthetic generations in the 1945 and 1990 
#   birth cohorts
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
}else{
  packages <-
    c(
      "tidyverse",
      "dplyr",
      "plotly",
      "parallel"
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

gss_data_bc <- readr::read_rds(
  file.path(
    dir_data,
    "main_rst_bc.rds"
  )
)$validation_lst

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

year_vec <- 
  as.character(seq(1945,1990,3))


cohort_bins <-
  tibble(
    cohort_index = seq(1:length(gss_data_bc)),
    cohort_bins_vec = seq(1945,1990)
  )

cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)

#-----------------------------------------------------------------#
#                   PLOT THE IM BY CLASS                         #
#----------------------------------------------------------------#

#-----------------------------#
#  1945 Cohort: IM over Time
#-----------------------------#

IM_1945 <- gss_data_bc[[1]]$memory_df
IM_1945_t <-
  plot_memory_contrast(
    data = IM_1945,
    initial = TRUE,
    t_num = c(1,2,3,4),
    y_tick_vals = seq(0,0.22,0.02),
    y_range = c(-0.01,0.22),
    y_tick_text = to_text(seq(0,0.22,0.02)),
    legend_x = 0.85,
    legend_y = 0.9,
    individual = TRUE,
    x_tick_vals = seq(1,4)
    )

#-----------------------------#
#  1990 Cohort: IM over Time
#-----------------------------#

IM_1990 <- gss_data_bc[[46]]$memory_df
IM_1990_t <-
  plot_memory_contrast(
    data = IM_1990,
    initial = TRUE,
    t_num = c(1,2,3,4),
    y_tick_vals = seq(0,0.22,0.02),
    y_range = c(-0.01,0.22),
    y_tick_text = to_text(seq(0,0.22,0.02)),
    legend_x = 0.85,
    legend_y = 0.9,
    individual = TRUE,
    x_tick_vals = seq(1,4)
  ) 



#-----------------------------#
#  1945 Cohort: AIM over Time
#-----------------------------#

AIM_1945_t <-
  plot_memory_contrast(
    data = IM_1945,
    initial = TRUE,
    t_num = c(1,2,3,4),
    y_tick_vals = seq(0,0.22,0.02),
    y_range = c(-0.01,0.22),
    y_tick_text = to_text(seq(0,0.22,0.02)),
    legend_x = 0.9,
    legend_y = 0.9,
    individual = FALSE,
    x_tick_vals = seq(1,4)
  ) 


#-----------------------------#
#  1990 Cohort: AIM over Time
#-----------------------------#

AIM_1990_t <-
  plot_memory_contrast(
    data = IM_1990,
    initial = TRUE,
    t_num = c(1,2,3,4),
    y_tick_vals = seq(0,0.22,0.02),
    y_range = c(-0.01,0.22),
    y_tick_text = to_text(seq(0,0.22,0.02)),
    legend_x = 0.9,
    legend_y = 0.9,
    individual = FALSE,
    x_tick_vals = seq(1,4)
  ) 


#-----------------------------#
#  Save the Plot:
#-----------------------------#

# 1. Save IM 1945
saveRDS(
  IM_1945_t,
  file = path.expand(file.path(dir_plot, "fig8_IM_1945_t.rds"))
)

# 2. Save AIM 1945
saveRDS(
  AIM_1945_t,
  file = path.expand(file.path(dir_plot, "fig8_AIM_1945_t.rds"))
)

# 3. Save IM 1990
saveRDS(
  IM_1990_t,
  file = path.expand(file.path(dir_plot, "fig8_IM_1990_t.rds"))
)

# 4. Save AIM 1990
saveRDS(
  AIM_1990_t,
  file = path.expand(file.path(dir_plot, "fig8_AIM_1990_t.rds"))
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure 8 of the paper’s main results have been successfully replicated.\n")

sink()




