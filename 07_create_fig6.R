#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "main_results"
title <- "plot_fig6"

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
# Repository:  https://github.com/synthetic_dynasty_analysis
#
# Script:      07_create_fig6.R
#
# Inputs:
#   - Baseline Result Data:
#       data/main_results/estimation/main_rst_baseline.rds
#   - Bootstrap Result Data:
#       data/main_results/estimation/main_rst_boot.rds
#
#   - Helper functions:
#       _Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       _Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/main_results/_LOGS/plot_fig6.txt
#   - Plot:
#       data/main_results/plots/fig6_AIM_delta.jpg
#
# Description:
#   This script recreate Figure 6: Rate of decay in aggregate intergenerational  
#   memory between generation t = 1 and t = 2.
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


#--------------------------------------------------------#
#           LOAD BASELINE AND BOOTSTRAP RESULTS         
#-------------------------------------------------------#


# Bootstrap Results:

gss_data_boot <- readr::read_rds(
  file.path(
    dir_data,
    "main_rst_boot.rds"
  )
)

# Baseline Results:

gss_data_baseline <- readr::read_rds(
  file.path(
    dir_data,
    "main_rst_baseline.rds"
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

# Load data process functions:
source(
  file.path(
    dir_root,
    "Functions",
    "process_results.R"
  )
)

year_vec <- 
  as.character(seq(1945,1990,3))


cohort_bins <-
  tibble(
    cohort_index = seq(1:length(gss_data_baseline)),
    cohort_bins_vec = seq(1945,1990)
  )

cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)

#------------------------------------------------------------------#
#                     PLOT THE DELTA AIM                          #
#----------------------------------------------------------------#

#------------------------------------------------------#
#     Calculate the measure in bootstrap samples:
#------------------------------------------------------#

# Calculate the delta AIM change from the bootstrap results:

delta_AIM_boot_rst <-
lapply(
AIM_info <-
  lapply(
    modify_cohort_lsts(gss_data_boot, qt = 0.95),
    function(cohort_boot_lst){
      AIM_lst <-
        lapply(
          cohort_boot_lst,
          function(boot){
            AIM_df <-
              boot$`Memory Curve` %>%
              dplyr::select(
                AIM,
                t
              ) %>%
              mutate(
                AIM_change = AIM - lead(AIM,default = 0),
                AIM_change_r = AIM_change / AIM
              ) %>%
              drop_na()
            return(AIM_df)
          }
        )
      return(AIM_lst)
    }
  ),
  function(cohort){
    AIM <-
      lapply(
       as.list(seq(1,3)), 
        function(t){map_dfr(AIM_info[[1]], ~ .x[t, ])})
       return(AIM)
      }
)


# Calculate the CI from the bootstrap results:

delta_AIM_boot_stats <-
  bind_rows(
    lapply(
      delta_AIM_boot_rst ,
      function(cohort_lsts){
        stat_df <-
          bind_rows(
            lapply(
              cohort_lsts,
              function(AIM_df){
                AIM_df <-
                  AIM_df %>%
                  mutate(
                    mean_boot = mean(AIM_change),
                    CI_np_l = 
                      quantile(AIM_change, 0.025, na.rm = TRUE),
                    CI_np_u = 
                      quantile(AIM_change, 0.975, na.rm = TRUE),
                    mean_r_boot = mean(AIM_change_r),
                    CI_np_l_r = 
                      quantile(AIM_change_r, 0.025, na.rm = TRUE),
                    CI_np_u_r = 
                      quantile(AIM_change_r, 0.975, na.rm = TRUE)
                  ) %>%
                  dplyr::select(
                    t,
                    mean_boot,
                    CI_np_l,
                    CI_np_u,
                    mean_r_boot,
                    CI_np_l_r,
                    CI_np_u_r
                  ) %>%
                  unique()
              }
            )
          )
      }
    ),
    .id = "cohort_index"
  )

#---------------------------------------------------#
#      Calculate the measure in baseline sample:
#--------------------------------------------------#

# Calculate the delta AIM measure in baseline results:

delta_AIM_baseline <-
  bind_rows(
    lapply(
      gss_data_baseline,
      function(cohort){
        AIM_info <-
          cohort$`Memory Measure` %>%
          filter(
            class == "AIM"
          ) %>%
          dplyr::select(
            mean,
            t
          ) %>%
          mutate(
            AIM_change_baseline = 
              mean - lead(mean, default = 0),
            AIM_change_r_baseline = AIM_change_baseline / mean
          )
      }
    ),
    .id = "index"
  ) %>%
  filter(
    mean != 0
  ) %>%
  mutate(
    cohort_index = as.character(rep(seq(1945,1990), each = 3))
  )

#----------------------------------------------#
#  Calculate the bias-corrected result:
#---------------------------------------------#

delta_AIM_bc_df <-
  delta_AIM_baseline %>%
  dplyr::select(
    -mean
  ) %>%
  left_join(
    delta_AIM_boot_stats,
    by = c("t","cohort_index")
  ) %>%
  mutate(
    mean_bc = 2 * AIM_change_baseline - mean_boot,
    CI_bc_l = 2 * AIM_change_baseline - CI_np_u,
    CI_bc_u = 2 * AIM_change_baseline - CI_np_l,
    
    mean_bc_r = 2 * AIM_change_r_baseline - mean_r_boot,
    CI_bc_l_r = 2 * AIM_change_r_baseline - CI_np_u_r,
    CI_bc_u_r = 2 * AIM_change_r_baseline - CI_np_l_r
  )

delta_AIM_bc_lst<-
  split(
    delta_AIM_bc_df,
    delta_AIM_bc_df$t
  )

#-----------------------------#
#  Generate the Plot:
#-----------------------------#

config(
 AIM_change_r <-
  plot_AIM_change(
    plot_data = delta_AIM_bc_lst,
    t = 1,
    y_title = "$\\Delta_{AIM}\\;\\text{(proportion of memory eliminated by }t=2\\text{)}$", 
    x_tick_vals = seq(1945,1990,3), 
    mean_col = "mean_bc_r", 
    CI_u_col = "CI_bc_u_r",
    CI_l_col = "CI_bc_l_r",
    y_tick_vals = seq(0.65,0.9,0.025), 
    x_range = c(1944,1991), 
    y_range = c(0.65, 0.9), 
    y_tick_text = to_text(seq(0.65,0.9,0.025))
  ),mathjax = "cdn") %>% 
  layout(
    width = 1000, 
    height = 600
    )

#-----------------------------#
#  Save the Plot:
#-----------------------------#

saveRDS(
  AIM_change_r,
  file = path.expand(
    file.path(
      dir_plot,
      "fig6_AIM_delta.rds"  
    )
  )
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure 6 of the paper’s main results have been successfully replicated.\n")

sink()




