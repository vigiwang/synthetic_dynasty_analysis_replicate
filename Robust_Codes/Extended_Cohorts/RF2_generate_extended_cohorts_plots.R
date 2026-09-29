#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "robustness_check"
subsubsection <- "extended_cohorts"
title <- "plot_figJ1_J4"

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
# Script:      RF2_generate_extended_cohorts.R
#
# Inputs:
#   - Extended-cohorts Bias-Corrected estimation results:
#       Data/robustness_check/more_cohorts_bc.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/robustness_check_plot_figJ1_J4_log.txt
#   - Plot:
#       Plot/robustness_check/extended_cohorts/figJ1_extended_cohorts_OM.rds
#       Plot/robustness_check/extended_cohorts/figJ1_extended_cohorts_SM.rds
#       Plot/robustness_check/extended_cohorts/figJ1_extended_cohorts_EM.rds
#       Plot/robustness_check/extended_cohorts/figJ1_extended_cohorts_UP_DOWN.rds
#       Plot/robustness_check/extended_cohorts/figJ2_extended_cohorts_OM_by_t.rds
#       Plot/robustness_check/extended_cohorts/figJ3_extended_cohorts_AIM.rds
#       Plot/robustness_check/extended_cohorts/figJ4_extended_cohorts_bin5_LML.rds

# Description:
#   This script recreate Figure J1: Overall Mobility(t = 1), 
#   Structural Mobility(t = 1), and Exchange Mobility(t = 1),
#   Upward-Downward Mobility(t = 1); Figure J2: Overall Mobility 
#   across synthetic generations.in paper; Figure J3: Aggregate 
#   intergenerational memory at generation t = 1; Figure J4: Log-multiplicative
#   Layer effect.
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
    "more_cohorts_bc.rds"
  )
)$validation_lst


gss_data_main <- readr::read_rds(
  file.path(
    dir_root,
    "Data",
    "main_results",
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


cohort_bins <-
  tibble(
    cohort_index = seq(1:length(gss_data_bc)),
    cohort_bins_vec = seq(1925,1990)
  )
cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)
year_vec <- as.character(seq(1925,1990,4))
bins_vec_5 <- c(
  "[1925,1930)",
  "[1930,1935)",
  "[1935,1940)",
  "[1940,1945)",
  "[1945,1950)",
  "[1950,1955)",
  "[1955,1960)",
  "[1960,1965)",
  "[1965,1970)",
  "[1970,1975)",
  "[1975,1980)",
  "[1980,1990)",
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
    x_tick_vals = seq(1925, 1990, by = 4),
    y_tick_vals = seq(0.55, 0.75, by = 0.02),
    y_tick_text = to_text(seq(0.55, 0.75, by = 0.02)),
    x_range = c(1924, 1991),
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
    x_tick_vals = seq(1925,1990,4),
    y_tick_vals = seq(0.1, 0.38, by = 0.02),
    x_range = c(1924,1991),
    y_range = c(0.1,0.38),
    y_tick_text = to_text(seq(0.1, 0.38, by = 0.02))
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
    x_tick_vals = seq(1925,1990,4),
    y_tick_vals = seq(0.3, 0.6, by = 0.02),
    x_range = c(1924,1991),
    y_range = c(0.3,0.6),
    y_tick_text = to_text(seq(0.3, 0.6, by = 0.02))
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
    x_tick_vals = seq(1925,1990,4),
    y_tick_vals = seq(0.22, 0.46, 0.02),
    x_range = c(1924,1991),
    y_range = c(0.22,0.46),
    y_tick_text = to_text(seq(0.22, 0.46, 0.02)),
    legend_x = 0.7,
    legend_y= 0.2
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
    x_tick_vals = seq(1925, 1990, by = 4),
    y_tick_vals = seq(0.55, 0.75, by = 0.02),
    y_tick_text = to_text(seq(0.55, 0.75, by = 0.02)),
    x_range = c(1924, 1991),
    y_range = c(0.55, 0.75),
    legend_x = 0.85,
    legend_y = 0.98,
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
    x_tick_vals = seq(1925,1990,4),
    y_tick_vals = seq(0.06, 0.26, by = 0.02),
    x_range = c(1924,1991),
    y_range = c(0.06,0.26),
    y_tick_text = to_text(seq(0.06, 0.26, by = 0.02))
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
    year_l = 1925,
    year_u = 1990,
    age_l = 25,
    age_u = 55,
    cohort_bins = bins_vec_5,
    tick_vals = seq(0.45,1.35,0.1),
    range = c(0.45, 1.35)
  )


#-----------------------------#
#  Save the Plot:
#-----------------------------#

saveRDS(
  updown_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ1_extended_cohorts_UP_DOWN.rds"
    )
  )
)

saveRDS(
  overall_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ1_extended_cohorts_OM.rds"
    )
  )
)

saveRDS(
  structural_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ1_extended_cohorts_SM.rds"
    )
  )
)

saveRDS(
  exchange_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ1_extended_cohorts_EM.rds"
    )
  )
)

saveRDS(
  overall_mobility_synthetic,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ2_extended_cohorts_OM_by_t.rds"
    )
  )
)

saveRDS(
  AIM,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ3_extended_cohorts_AIM.rds"
    )
  )
)

saveRDS(
  bin5_LML,
  file = path.expand(
    file.path(
      dir_plot,
      "figJ4_extended_cohorts_bin5_LML.rds"
    )
  )
)

# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure J1-H4 have been successfully replicated.\n")

sink()



