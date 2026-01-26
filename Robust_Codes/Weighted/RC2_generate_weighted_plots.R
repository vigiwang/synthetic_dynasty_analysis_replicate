#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "robustness_check"
subsubsection <- "weighted"
title <- "plot_figC1_C7"

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
# Repository:  https://github.com/synthetic_dynasty_analysis
#
# Script:      R6b_create_figC1_C7.R
#
# Inputs:
#   - GSS father–child occupational pairs (5-class EGP):
#       data/robust/class_typology/weighted_bc.rds
#   - Helper functions:
#       _Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       _Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/robust/_LOGS/plot_figC1_4.txt
#   - Plot:
#       data/robust/plots/age/figC1_weighted_OM.png
#       data/robust/plots/age/figC1_weighted_SM.png
#       data/robust/plots/age/figC1_weighted_EM.png
#       data/robust/plots/age/figC1_weighted_UP_DOWN.png
#       data/robust/plots/age/figC2_weighted_OM_by_t.png
#       data/robust/plots/age/figC3_weighted_AIM.png
#       data/robust/plots/age/figC4_weighted_bin5_LML.png

# Description:
#   This script recreate Figure C1: Overall Mobility(t = 1), 
#   Structural Mobility(t = 1), and Exchange Mobility(t = 1),
#   Upward-Downward Mobility(t = 1); Figure C2: Overall Mobility 
#   across synthetic generations.in paper; Figure C3: Aggregate 
#   intergenerational memory at generation t = 1; Figure C4:
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
    "weighted_bc.rds"
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

gss_data_weighted <- 
  bind_rows(
    lapply(
      split(gss_data_main,gss_data_main$year),
      function(df){
        df <-
          df %>%
          mutate(
            wn = wtssps / mean(wtssps)
          ) %>%
          mutate(
            flag = mean(wn)
          )
        assert_that(
          abs(unique(df$flag) -  1) < 1e-8
        )
        return(df)
      }
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
year_vec <- as.character(seq(1945,1990,3))
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
#=========================================================#
#   Fig C1: PLOT THE FOUR PLOTS: OM1, SM1, EM1, UP-DOWN1  #
#=========================================================#

#---------------------#
# Overall Mobility
#--------------------#

overall_mobility <- 
  plot_overall_mobility(
    data = gss_data_bc,
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Historical Mobility",
    label_name = "Overall mobility",
    t_num = c(1),
    initial = TRUE,
    x_tick_vals = seq(1945, 1990, by = 3),
    y_tick_vals = seq(0.55, 0.75, by = 0.02),
    y_tick_text = to_text(seq(0.55, 0.75, by = 0.02)),
    x_range = c(1944, 1991),
    y_range = c(0.55, 0.75),
    legend_x = 0.7,
    legend_y = 0.8,
    legend_size = 20
  )  %>% 
  plotly::layout(
    width = 1000, 
    height = 600
  )


#---------------------#
# Structural Mobility
#--------------------#


structural_mobility <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "movement_df",
    y_title = "Probability to move",
    attribute_col = "category",
    metric_name = "Structural Mobility",
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.1, 0.3, by = 0.02),
    x_range = c(1944,1991),
    y_range = c(0.1,0.3),
    y_tick_text = to_text(seq(0.1, 0.3, by = 0.02))
  ) %>% 
  plotly::layout(
    width = 1000, 
    height = 600
  )


#---------------------#
# Exchange Mobility
#--------------------#


exchange_mobility <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Exchange Mobility",
    y_title = "Probability to move",
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.4, 0.6, by = 0.02),
    x_range = c(1944,1991),
    y_range = c(0.4,0.6),
    y_tick_text = to_text(seq(0.4, 0.6, by = 0.02))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )

#--------------------------#
# Upward-Downward Mobility
#--------------------------#


updown_mobility <-
  plot_updown_mobility(
    data = gss_data_bc,
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.25, 0.45, 0.02),
    x_range = c(1944,1991),
    y_range = c(0.25,0.45),
    y_tick_text = to_text(seq(0.25, 0.45, 0.02)),
    legend_x = 0.7,
    legend_y= 0.2
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )


#=========================================================#
#   Fig C2: PLOT THE OVERALL MOBILITY BY T
#=========================================================#


overall_mobility_synthetic <- 
  plot_overall_mobility(
    data = gss_data_bc,
    measure = "movement_df",
    attribute_col = "category",
    metric_name = "Historical Mobility",
    label_name = "Overall mobility",
    t_num = c(1,2,3,4),
    initial = FALSE,
    x_tick_vals = seq(1945, 1990, by = 3),
    y_tick_vals = seq(0.55, 0.75, by = 0.02),
    y_tick_text = to_text(seq(0.55, 0.75, by = 0.02)),
    x_range = c(1944, 1991),
    y_range = c(0.55, 0.75),
    legend_x = 0.7,
    legend_y = 0.98,
    legend_size = 18
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )



#=========================================================#
#  Fig C3: MEAN TIME TO EXIT FROM EACH ORIGIN CLASS(MTEi)
#=========================================================#


MTE_class <- 
  plot_measure_by_class(
    data = gss_data_bc,
    measure = "MTE_df",
    attribute_col = "class",
    metric_name = paste0("Class", seq(1,5)),
    y_title = "Number of generations",
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(1,2.1,0.1),
    x_range = c(1944,1991),
    y_range = c(1,2.1),
    y_tick_text = c("1.0",seq(1.1,1.9,0.1),"2.0","2.1"),
    legend_x = 0.985,
    legend_y = 0.95
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )


#================================================================#
#   Fig C4: AGGREGATE INTERGENERATIONAL MEMORY AT GENERATION T = 1
#================================================================#


AIM <- 
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "memory_df",
    attribute_col = "class",
    metric_name = "AIM",
    y_title = "Memory",
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0, 0.22, by = 0.02),
    x_range = c(1944,1991),
    y_range = c(0,0.22),
    y_tick_text = to_text(seq(0, 0.22, by = 0.02))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )

#========================================================================#
#   Fig C5: INTERGENERATIONAL MEMORY AT GENERATION T = 1 by ORIGIN CLASS
#========================================================================#

IM_class <-
  plot_measure_by_class(
    data = gss_data_bc,
    measure = "memory_df",
    attribute_col = "class",
    metric_name = paste0("Class", seq(1,5)),
    y_title = "Memory",
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.02,0.24,0.02),
    x_range = c(1944,1991),
    y_range = c(0.02,0.24),
    y_tick_text = to_text(seq(0.02,0.24,0.02)),
    legend_x = 0.985,
    legend_y = 0.95
  ) %>% 
  layout(
    width = 1000, 
    height = 600
  )


#========================================================================#
#   Fig C6: MEMORY ACROSS SYNTHETIC GENERATIONS in 1945 & 1990 COHORTS
#========================================================================#

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


#========================================================================#
#   Fig C7: LOG-MULTIPLICATIVE LAYER EFFECT WITH BINNED DATA  
#========================================================================#


#-----------------------------#
#  Bin = 5:
#-----------------------------#

bin5_LML <-
  plot_log_multiplicative_weight(
    data = gss_data_weighted,
    bin_width = 5,
    year_l = 1945,
    year_u = 1990,
    age_l = 25,
    age_u = 55,
    cohort_bins = bins_vec_5
  )

#-----------------------------#
#  Bin = 10:
#-----------------------------#

bin10_LML <-
  plot_log_multiplicative_weight(
    data = gss_data_weighted,
    bin_width = 10,
    year_l = 1945,
    year_u = 1990,
    age_l = 25,
    age_u = 55,
    cohort_bins = bins_vec_10
  )


#-----------------------------#
#  Save the Plot:
#-----------------------------#

# 1. Up/Down Mobility (Weighted)
saveRDS(
  updown_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figC1_weighted_UP_DOWN.rds"
    )
  )
)

# 2. Overall Mobility (Weighted)
saveRDS(
  overall_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figC1_weighted_OM.rds"
    )
  )
)

# 3. Structural Mobility (Weighted)
saveRDS(
  structural_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figC1_weighted_SM.rds"
    )
  )
)

# 4. Exchange Mobility (Weighted)
saveRDS(
  exchange_mobility,
  file = path.expand(
    file.path(
      dir_plot,
      "figC1_weighted_EM.rds"
    )
  )
)

# 5. Overall Mobility Synthetic (Weighted)
saveRDS(
  overall_mobility_synthetic,
  file = path.expand(
    file.path(
      dir_plot,
      "figC2_weighted_OM_by_t.rds"
    )
  )
)

# 6. MTE Class (Weighted)
saveRDS(
  MTE_class,
  file = path.expand(
    file.path(
      dir_plot,
      "figC3_weighted_MTE.rds"
    )
  )
)

# 7. AIM (Weighted)
saveRDS(
  AIM,
  file = path.expand(
    file.path(
      dir_plot,
      "figC4_weighted_AIM.rds"
    )
  )
)

# 8. IM Class (Weighted)
saveRDS(
  IM_class,
  file = path.expand(
    file.path(
      dir_plot,
      "figC5_weighted_IM.rds"
    )
  )
)

# 9. IM 1945 (Weighted)
saveRDS(
  IM_1945_t,
  file = path.expand(
    file.path(
      dir_plot,
      "figC6_weighted_IM_1945_t.rds"
    )
  )
)

# 10. IM 1990 (Weighted)
saveRDS(
  IM_1990_t,
  file = path.expand(
    file.path(
      dir_plot,
      "figC6_weighted_IM_1990_t.rds"
    )
  )
)

# 11. AIM 1945 (Weighted)
saveRDS(
  AIM_1945_t,
  file = path.expand(
    file.path(
      dir_plot,
      "figC6_weighted_AIM_1945_t.rds"
    )
  )
)

# 12. AIM 1990 (Weighted)
saveRDS(
  AIM_1990_t,
  file = path.expand(
    file.path(
      dir_plot,
      "figC6_weighted_AIM_1990_t.rds"
    )
  )
)

# 13. Bin5 LML (Weighted)
saveRDS(
  bin5_LML,
  file = path.expand(
    file.path(
      dir_plot,
      "figC7_weighted_bin5_LML.rds"
    )
  )
)

# 14. Bin10 LML (Weighted)
saveRDS(
  bin10_LML,
  file = path.expand(
    file.path(
      dir_plot,
      "figC7_weighted_bin10_LML.rds"
    )
  )
)

# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure C1-C7 have been successfully replicated.\n")

sink()



