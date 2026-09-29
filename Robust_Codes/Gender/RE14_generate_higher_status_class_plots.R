#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "robustness_check"
subsubsection <- "gender"
title <- "plot_figE25_E27"

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
# Script:      RE14_generate_higher_status_class_plots.R
#
# Inputs:
#   - Parent of higher status-child bias-corrected results:
#       Data/robustness_check/5class_pc_bc.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/robustness_check_plot_figE25_E27_log.txt
#   - Plot:
#       Plot/robustness_check/gender/figE25_pc_MTE.rds
#       Plot/robustness_check/gender/figE26_pc_IM.rds
#       Plot/robustness_check/gender/figE27_pc_MFPT.rds

# Description:
#   This script creates Figure E25: Mean time to exit from each class;
#   Figure E26: Intergenerational memory at generation t = 1 by class;
#   Figure E27: Mean first passage time between classes, for the parent
#   of higher status-child sample.
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
  library("ggplot2")
}else{
  packages <-
    c(
      "tidyverse",
      "dplyr",
      "plotly",
      "parallel",
      "ggplot2"
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
    "5class_pc_bc.rds"
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
#              PLOT THE MEAN TIME TO EXIT BY CLASS                #
#----------------------------------------------------------------#

MTE_class <-
  plot_measure_by_class(
    data = gss_data_bc,
    measure = "MTE_df",
    attribute_col = "class",
    metric_name = paste0("Class", seq(1,5)),
    y_title = "Number of generations",
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.7,2.7,0.2),
    x_range = c(1944,1991),
    y_range = c(0.7,2.7),
    y_tick_text = sprintf("%.1f", seq(0.7,2.7,0.2)),
    legend_x = 0.985,
    legend_y = 0.95
  ) %>%
  layout(
    width = 1000,
    height = 600
  )


#-----------------------------------------------------------------#
#                   PLOT THE IM BY CLASS                         #
#----------------------------------------------------------------#

IM_class <-
  plot_measure_by_class(
    data = gss_data_bc,
    measure = "memory_df",
    attribute_col = "class",
    metric_name = paste0("Class", seq(1,5)),
    y_title = "Memory",
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(-0.02,0.38,0.04),
    x_range = c(1944,1991),
    y_range = c(-0.02,0.38),
    y_tick_text = to_text(seq(-0.02,0.38,0.04)),
    legend_x = 0.985,
    legend_y = 0.95
  ) %>%
  layout(
    width = 1000,
    height = 600
  )


#===========================================================#
#           PLOT THE MEAN FIRST PASSAGE TIME                #
#===========================================================#

MFP_updown_new <-
  MFP_class_comparison_new(gss_data_bc, "mean", 5, 10, 15, FALSE, NULL)$updown +
  labs(title = "")


#-----------------------------#
#  Save the Plot:
#-----------------------------#

saveRDS(
  MTE_class,
  file = path.expand(
    file.path(
      dir_plot,
      "figE25_pc_MTE.rds"
    )
  )
)

saveRDS(
  IM_class,
  file = path.expand(
    file.path(
      dir_plot,
      "figE26_pc_IM.rds"
    )
  )
)

saveRDS(
  MFP_updown_new,
  file = path.expand(
    file.path(
      dir_plot,
      "figE27_pc_MFPT.rds"
    )
  )
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure E25-E27 have been successfully replicated.\n")

sink()
