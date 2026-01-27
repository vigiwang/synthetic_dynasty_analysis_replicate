#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "main_results"
title <- "plot_fig1"

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
# Script:      02_create_fig1.R
#
# Inputs:
#   - GSS father–child occupational pairs (5-class EGP):
#       data/main_results/gss_fc_occ10_5class.rds
#   - Helper functions:
#       _Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       _Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/plot/_LOGS/plot_fig1.txt
#   - Plot:
#       Plot/main_results/fig1_father.pdf
#       Plot/main_results/fig1_child.pdf
# Description:
#   This script recreate Figure 1: Sample sizes by class origin, destination, 
#   and birth cohort in paper.
#
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES AND CMED R PACKAGE   #
#-------------------------------------------------#

if(HPC){
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("tidyverse")
  library("dplyr")
  library("haven")
  library("ggplot2")
  library("patchwork")
}else{
  packages <-
    c(
      "tidyverse",
      "dplyr",
      "ggplot2",
      "haven",
      "patchwork"
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

year_vec <- 
  as.character(seq(1945,1990,3))

#-----------------------------------------------------------------#
#   PLOT THE MARGINAL DISTRIBUTION FOR CHILDREN ACROSS COHORTS   #
#----------------------------------------------------------------#


child_5_org <- 
  plot_marginal(
    gss_data_main %>%
      mutate(
        across(where(is.labelled), as.numeric)
      ) %>%
      filter(
        age >= 25 & age <= 55
      ) %>%
      filter(
        cohort >= 1945 & cohort <= 1990
      ) %>%
      mutate(
        status_c = case_when(
          status_c == 1 ~ "Class1",
          status_c == 2 ~ "Class2",
          status_c == 3 ~ "Class3",
          status_c == 4 ~ "Class4",
          status_c == 5 ~ "Class5"
        )
      ), 
    status_c, 
    c("Class1","Class2","Class3","Class4","Class5"),
    5)


#-----------------------------------------------------------------#
#   PLOT THE MARGINAL DISTRIBUTION FOR FATHERS ACROSS COHORTS    #
#----------------------------------------------------------------#


father_5_org <- 
  plot_marginal(
    gss_data_main %>%
      mutate(
        across(where(is.labelled), as.numeric)
      ) %>%
      filter(
        age >= 25 & age <= 55
      ) %>%
      filter(
        cohort >= 1945 & cohort <= 1990
      ) %>%
      mutate(
        status_p = case_when(
          status_p == 1 ~ "Class1",
          status_p == 2 ~ "Class2",
          status_p == 3 ~ "Class3",
          status_p == 4 ~ "Class4",
          status_p == 5 ~ "Class5"
        )
      ), 
    status_p, 
    c("Class1","Class2","Class3","Class4","Class5"),
    5)


#-----------------------------#
#  Save the Plot:
#-----------------------------#

ggsave(
  filename = file.path(
    dir_plot,
    "fig1_child.pdf"
  ),
  plot     = child_5_org$p,
  width    = 17,
  height   = 12,
  units    = "in",
  device   = cairo_pdf,    
  bg       = "white"
)


ggsave(
  filename = file.path(
    dir_plot,
    "fig1_father.pdf"
  ),
  plot     = father_5_org$p,
  width    = 17,
  height   = 12,
  units    = "in",
  device   = cairo_pdf,    
  bg       = "white"
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure 1 of the paper’s main results have been successfully replicated.\n")

sink()






