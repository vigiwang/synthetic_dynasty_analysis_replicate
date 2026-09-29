#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "robustness_check"
subsubsection <- "alternative_smoothing"
title <- "plot_figH4_H6"

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
# Script:      RG4_generate_bin5_plots.R
#
# Inputs:
#   - Bin = 5 Bias-corrected results:
#       Data/robustness_check/bin5_bc.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/robustness_check_plot_figH4_H6_log.txt
#   - Plot:
#       Plot/robustness_check/alternative_smoothing/figH4_bin5_OM.rds
#       Plot/robustness_check/alternative_smoothing/figH4_bin5_SM.rds
#       Plot/robustness_check/alternative_smoothing/figH4_bin5_EM.rds
#       Plot/robustness_check/alternative_smoothing/figH4_bin5_UP_DOWN.rds
#       Plot/robustness_check/alternative_smoothing/figH5_bin5_OM_by_t.rds
#       Plot/robustness_check/alternative_smoothing/figH6_bin5_AIM.rds

# Description:
#   This script recreate Figure H4: Overall Mobility(t = 1), 
#   Structural Mobility(t = 1), and Exchange Mobility(t = 1),
#   Upward-Downward Mobility(t = 1); Figure H5: Overall Mobility 
#   across synthetic generations.in paper; Figure H6: Aggregate 
#   intergenerational memory at generation t = 1; 
#
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
#                LOAD GSS DATA
#-------------------------------------------------#

gss_data_bc <- readr::read_rds(
  file.path(
    dir_data,
    "bin5_bc.rds"
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


cohort_bins <-
  tibble(
    cohort_index = seq(1:length(gss_data_bc)),
    cohort_bins_vec = c(
      "[1945,1950)",
      "[1950,1955)",
      "[1955,1960)",
      "[1960,1965)",
      "[1965,1970)",
      "[1970,1975)",
      "[1975,1980)",
      "[1980,1985)",
      "[1985,1990]"
    )
  )
cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)
year_vec <- c(
  "[1945,1950)",
  "[1950,1955)",
  "[1955,1960)",
  "[1960,1965)",
  "[1965,1970)",
  "[1970,1975)",
  "[1975,1980)",
  "[1980,1985)",
  "[1985,1990]"
)
#-----------------------------------------------------------------#
#   PLOT THE FOUR PLOTS: OM1, SM1, EM1, UP-DOWN1                 #
#----------------------------------------------------------------#

#=================#
# Overall Mobility
#================#


overall_mobility <- 
  plot_overall_mobility(
    data = gss_data_bc,
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Historical Mobility",
    label_name = "Overall mobility",
    t_num = c(1),
    initial = TRUE,
    x_tick_vals = seq(0,8,1),
    y_tick_vals = seq(0.55, 0.75, by = 0.05),
    y_tick_text = to_text(seq(0.55, 0.75, by = 0.05)),
    x_range = c(-0.5, 8.5),
    y_range = c(0.55, 0.75),
    legend_x = 0.7,
    legend_y = 0.8,
    legend_size = 20
  )  %>% 
  plotly::layout(
    width = 1000, 
    height = 600
  )

#====================#
# Structural Mobility
#===================#

structural_mobility <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "movement_df",
    y_title = "Probability to move",
    attribute_col = "category",
    metric_name = "Structural Mobility",
    x_tick_vals = seq(0,8,1),
    y_tick_vals = seq(0.08, 0.28, by = 0.02),
    x_range = c(-0.5, 8.5),
    y_range = c(0.08,0.28),
    y_tick_text = to_text(seq(0.08, 0.28, by = 0.02))
  ) %>% 
  plotly::layout(
    width = 1000, 
    height = 600
  )

#====================#
# Exchange Mobility
#===================#

exchange_mobility <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Exchange Mobility",
    y_title = "Probability to move",
    x_tick_vals = seq(0,8,1),
    y_tick_vals = seq(0.4, 0.6, by = 0.02),
    x_range = c(-0.5, 8.5),
    y_range = c(0.4,0.6),
    y_tick_text = to_text(seq(0.4, 0.6, by = 0.02))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )

#=========================#
# Upward-Downward Mobility
#========================#

updown_mobility <-
  plot_updown_mobility(
    data = gss_data_bc,
    t_num = c(1),
    x_tick_vals = seq(0,8,1),
    y_tick_vals = seq(0.24, 0.46, 0.02),
    x_range = c(-0.5, 8.5),
    y_range = c(0.24,0.46),
    y_tick_text = to_text(seq(0.24, 0.46, 0.02)),
    legend_x = 0.7,
    legend_y= 0.15
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )

#=======================#
# Overall Mobility by t
#=======================#


overall_mobility_synthetic <- 
  plot_overall_mobility(
    data = gss_data_bc,
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Historical Mobility",
    label_name = "Overall mobility",
    t_num = c(1,2,3,4,5),
    initial = FALSE,
    x_tick_vals = seq(0,8,1),
    y_tick_vals = seq(0.55, 0.75, by = 0.02),
    y_tick_text = to_text(seq(0.55, 0.75, by = 0.02)),
    x_range = c(-0.5, 8.5),
    y_range = c(0.55, 0.75),
    legend_x = 0.78,
    legend_y = 0.36,
    legend_size = 12
  ) %>% 
  layout(
    width = 1000, 
    height = 600
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
    x_tick_vals = seq(0,8,1),
    y_tick_vals = seq(0, 0.22, by = 0.02),
    x_range = c(-0.5, 8.5),
    y_range = c(0,0.22),
    y_tick_text = to_text(seq(0, 0.22, by = 0.02))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )


#-----------------------------#
#  Save the Plot:
#-----------------------------#

saveRDS(
  updown_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figH4_bin5_UP_DOWN.rds"
    )
  )
)

saveRDS(
  overall_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figH4_bin5_OM.rds"
    )
  )
)

saveRDS(
  structural_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figH4_bin5_SM.rds"
    )
  )
)

saveRDS(
  exchange_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figH4_bin5_EM.rds"
    )
  )
)

saveRDS(
  overall_mobility_synthetic,
  file = path.expand(
    file.path(
      dir_plot,
      "figH5_bin5_OM_by_t.rds"
    )
  )
)

saveRDS(
  AIM,
  file = path.expand(
    file.path(
      dir_plot,
      "figH6_bin5_AIM.rds"
    )
  )
)
# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure H4-G6 have been successfully replicated.\n")

sink()



