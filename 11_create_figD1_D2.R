#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "main_results"
title <- "plot_figD1_D2"

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
# Script:      11_create_figD1_D2.R
#
# Inputs:
#   - Bias-corrected main results:
#       Data/main_results/main_rst_bc.rds
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/main_results_plot_figD1_D2_log.txt
#   - Plot:
#       Plot/main_results/figD1_MTE_class_CI.rds
#       Plot/main_results/figD2_IM_class_CI.rds
#
# Description:
#   This script recreates Figure D.1: Mean time to exit from each origin
#   class with 95% confidence intervals, and Figure D.2: Intergenerational
#   memory at generation t = 1 by origin class with 95% confidence
#   intervals. Each class is drawn in its own stacked panel (shared birth
#   cohort axis) so the bias-corrected confidence bands remain readable.
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


cohort_bins <-
  tibble(
    cohort_index = seq(1:length(gss_data_bc)),
    cohort_bins_vec = seq(1945,1990)
  )
cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)


#-------------------------------------------------#
#   STACKED BY-CLASS PANELS WITH CONFIDENCE BANDS #
#-------------------------------------------------#

# One panel per class (shared cohort axis) so the 95% bias-corrected
# bands stay readable; class styles mirror plot_measure_by_class().
# Defined here so Functions/plot.R stays untouched.

plot_measure_by_class_CI <- function(
    data, measure, attribute_col, metric_name, y_title,
    t_num, x_tick_vals, x_range, y_tick_lst){
  plot_data <-
    cross_cohort_plot_data(
      data,
      measure,
      attribute_col,
      metric_name,
      t_num
    )
  class_colors  <- rev(c("#000000", "#1F1F1F", "#3D3D3D", "#5C5C5C", "#7A7A7A", "#999999", "#B7B7B7"))[seq_along(metric_name)]
  class_dashes  <- c("dash", "dot", "dashdot", "solid", "dash", "longdash", "longdashdot")[seq_along(metric_name)]
  class_symbols <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "circle", "square")[seq_along(metric_name)]

  panel_lst <- list()
  for(i in seq_along(metric_name)){
    m_i <- metric_name[i]
    dat_i <- plot_data %>% dplyr::filter(measure == m_i)
    tick_vals_i <- y_tick_lst[[i]]
    tick_text_i <- gsub("^(-?)0\\.", "\\1.", sprintf("%.2f", tick_vals_i))
    panel_lst[[i]] <-
      plot_ly() %>%
      add_ribbons(
        data = dat_i,
        x = as.numeric(cohort_bins_vec),
        ymin = ~CI_bc_l,
        ymax = ~CI_bc_u,
        fillcolor = "rgba(0,0,0,0.1)",
        line = list(color = "transparent"),
        name = "Confidence Interval",
        showlegend = FALSE,
        hoverinfo = "skip"
      ) %>%
      add_trace(
        data = dat_i,
        x = as.numeric(cohort_bins_vec),
        y = ~mean_bc,
        type = "scatter",
        mode = "lines+markers",
        name = m_i,
        showlegend = TRUE,
        line   = list(color = class_colors[i], dash = class_dashes[i], width = 2),
        marker = list(color = class_colors[i], symbol = class_symbols[i], size = 7,
                      line = list(color = "#000000", width = 0.5))
      ) %>%
      layout(
        yaxis = list(
          tickvals = tick_vals_i,
          ticktext = tick_text_i,
          range = range(tick_vals_i),
          tickes = "outside", ticklen = 5,
          tickwidth = 1, showline = TRUE,
          tickfont = list(size = 15, family = "Times New Roman", color = "black"),
          showgrid = FALSE, zeroline = FALSE
        )
      )
  }

  fig <-
    subplot(
      panel_lst,
      nrows = length(panel_lst),
      shareX = TRUE,
      titleY = FALSE,
      margin = 0.025
    ) %>%
    layout(
      font = list(family = "Times New Roman", color = "black"),
      xaxis = list(
        title = list(
          text = "Birth cohorts",
          font = list(
            size = 26,
            family = "Times New Roman",
            color = "black")),
        tickes = "outside", ticklen = 5,
        tickwidth = 1, showline = TRUE,
        tickfont = list(size = 15, family = "Times New Roman", color = "black"),
        showgrid = FALSE,
        tickangle = 0,
        tickmode = "array",
        tickvals = x_tick_vals,
        range = x_range,
        tickangle = 0),
      annotations = list(
        list(
          text = y_title,
          textangle = -90,
          x = -0.072,
          y = 0.5,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "middle",
          showarrow = FALSE,
          font = list(
            size = 26,
            family = "Times New Roman",
            color = "black"))),
      legend = list(
        orientation = "v",
        x = 1.02,
        y = 0.98,
        font = list(
          size = 20, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      ),
      margin = list(l = 110, b = 30, r = 120)
    )
  return(fig)
}


#-----------------------------------------------------------------#
#     FIG D1: MEAN TIME TO EXIT BY CLASS WITH 95% CI              #
#----------------------------------------------------------------#

MTE_class_CI <-
  plot_measure_by_class_CI(
    data = gss_data_bc,
    measure = "MTE_df",
    attribute_col = "class",
    metric_name = paste0("Class", seq(1,5)),
    y_title = "Number of generations",
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    x_range = c(1944,1991),
    y_tick_lst = list(
      seq(1.08, 1.56, 0.12), # Class 1
      seq(1.20, 1.68, 0.12), # Class 2
      seq(1.32, 1.80, 0.12), # Class 3
      seq(0.84, 1.32, 0.12), # Class 4
      seq(1.56, 2.04, 0.12)  # Class 5
    )
  ) %>%
  layout(
    width = 1000,
    height = 600
  )


#-----------------------------------------------------------------#
#     FIG D2: IM AT t = 1 BY CLASS WITH 95% CI                    #
#----------------------------------------------------------------#

IM_class_CI <-
  plot_measure_by_class_CI(
    data = gss_data_bc,
    measure = "memory_df",
    attribute_col = "class",
    metric_name = paste0("Class", seq(1,5)),
    y_title = "Memory",
    t_num = c(1),
    x_tick_vals = seq(1945,1990,3),
    x_range = c(1944,1991),
    y_tick_lst = list(
      seq(0.06, 0.30, 0.06),  # Class 1
      seq(0.00, 0.24, 0.06),  # Class 2
      seq(-0.06, 0.18, 0.06), # Class 3
      seq(0.00, 0.24, 0.06),  # Class 4
      seq(0.00, 0.24, 0.06)   # Class 5
    )
  ) %>%
  layout(
    width = 1000,
    height = 600
  )


#-----------------------------#
#  Save the Plot:
#-----------------------------#

saveRDS(
  MTE_class_CI,
  file = path.expand(
    file.path(
      dir_plot,
      "figD1_MTE_class_CI.rds"
    )
  )
)

saveRDS(
  IM_class_CI,
  file = path.expand(
    file.path(
      dir_plot,
      "figD2_IM_class_CI.rds"
    )
  )
)


# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure D.1-D.2 have been successfully replicated.\n")

sink()
