#----------Preliminaries----------#
rm(list = ls())
section <- "Plot"
subsection <- "robustness_check"
subsubsection <- "alternative_measures"
title <- "plot_figE1_E5"

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
# Script:      RE_generate_alternative_measure_plots.R
#
# Inputs:
#   - Baseline Result Data:
#       Data/main_results/main_rst_baseline.rds
#   - Bootstrap Result Data:
#       Data/main_results//main_rst_boot.rds
#
#   - Helper functions:
#       Functions/utils.R
#           (general-purpose helper functions used throughout the analysis)
#       Functions/plot.R
#           (codes used to process data and draw plots)
#
# Outputs:
#   - Log file:
#       code/_LOGS/Plot/robustness_check_plot_figE1_E5_log.txt
#   - Plot:
#       Plot/robustness_check/alternative_measure/figE1_outflow_t1.rds
#       Plot/robustness_check/alternative_measure/figE1_outflow_ss.rds
#       Plot/robustness_check/alternative_measure/figE2_MFPT.rds
#       Plot/robustness_check/alternative_measure/figE3_dprime_t.rds
#       Plot/robustness_check/alternative_measure/figE3_lambda2.rds
#       Plot/robustness_check/alternative_measure/figE4_altham_index.rds
#       Plot/robustness_check/alternative_measure/figE5_HD.rds
#
# Description:
#   This script recreate Figure E1: Class outflows at generation t = 1 
#   and at the steady state; Figure E2: Mean first passage time; Figure E3:
#   Alternative measures of memory; Figure E4:Alternative measures of dependence.
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



#--------------------------------------------------------#
#           LOAD BASELINE AND BOOTSTRAP RESULTS         
#-------------------------------------------------------#

# Bootstrap Results:

gss_data_boot <- 
  readr::read_rds(
  file.path(
    dir_root,
    "Data",
    "main_results",
    "main_rst_boot.rds"
  )
)

# Baseline Results:

gss_data_baseline <- 
  readr::read_rds(
  file.path(
    dir_root,
    "Data",
    "main_results",
    "main_rst_baseline.rds"
  )
)

# Bias-corrected Results:

gss_data_bc <- 
  readr::read_rds(
  file.path(
    dir_root,
    "Data",
    "main_results",
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
    cohort_index = seq(1:length(gss_data_bc)),
    cohort_bins_vec = seq(1945,1990)
  )

cohort_bins_vec <- as.character(cohort_bins$cohort_bins_vec)

#===========================================================#
#   Fig D1: PLOT THE CLASS OUTFLOW AT t = ! & STEADY STATE  #
#===========================================================#

#--------------------------------------------------#
# Calculate the outflow measures in bootstrap data:
#--------------------------------------------------#

boot_outflow <-
  lapply(
    modify_cohort_lsts(gss_data_boot, qt = 0.95),
    function(cohort_lst){
      lst <-
        lapply(
          cohort_lst,
          function(boot_lst){
            M <- matrix(1, nrow = 5, ncol = 5)
            diag(M) <- 0
            # Process the Original tm_num(t = 1):
            # -- Convert it to tm_perc:
            tm_num <- 
              boot_lst$tm_num %>% 
              dplyr::select(
                -rowsum,
                -validation
              )
            tm_perc <- tm_num / rowSums(tm_num)
            outflow_initial <- diag(M %*% t(as.matrix(tm_perc)))
            # Process the Steady State TM:
            steady_state = boot_lst$`Steady State TM`
            outflow_ss <- diag(M %*% t(as.matrix(steady_state)))
            return(
              list(
                outflow_initial = outflow_initial,
                outflow_ss = outflow_ss
              )
            )
          }
        )
    }
  )

#--------------------------------------------------#
# Calculate the outflow statistics in bootstrap data:
#--------------------------------------------------#

boot_outflow_stats <-
  lapply(
    boot_outflow,
    function(cohort){
      # Process Initial:
      initial_df <-
        as.data.frame(
          do.call(
            rbind,
            lapply(
              cohort,
              function(lst){
                lst$outflow_initial
              }
            )
          )) 
      
      colnames(initial_df) <- c("Class1","Class2","Class3","Class4","Class5")
      
      initial_df <-
        initial_df %>%
        summarise(across(
          everything(),
          list(
            mean_boot   = ~ mean(.x, na.rm = TRUE),
            CI_np_u_boot = ~ quantile(.x, 0.975, na.rm = TRUE),
            CI_np_l_boot = ~ quantile(.x, 0.025, na.rm = TRUE)
          ),
          .names = "{.col}_{.fn}"
        )) %>%
        pivot_longer(
          cols = matches("^Class\\d+_"),  
          names_to = c("class", "stat"),
          names_pattern = "^(Class\\d+)_(.+)$",   
          values_to = "value"
        )  %>%
        pivot_wider(
          names_from  = stat,                            
          values_from = value
        ) 
      
      # Process steady state:
      ss_df <-
        as.data.frame(
          do.call(
            rbind,
            lapply(
              cohort,
              function(lst){
                lst$outflow_ss
              }
            )
          )) 
      
      colnames(ss_df) <- c("Class1","Class2","Class3","Class4","Class5")
      
      ss_df <-
        ss_df %>%
        summarise(across(
          everything(),
          list(
            mean_boot   = ~ mean(.x, na.rm = TRUE),
            CI_np_u_boot = ~ quantile(.x, 0.975, na.rm = TRUE),
            CI_np_l_boot = ~ quantile(.x, 0.025, na.rm = TRUE)
          ),
          .names = "{.col}_{.fn}"
        )) %>%
        pivot_longer(
          cols = matches("^Class\\d+_"),  
          names_to = c("class", "stat"),
          names_pattern = "^(Class\\d+)_(.+)$",   
          values_to = "value"
        )  %>%
        pivot_wider(
          names_from  = stat,                            
          values_from = value
        ) 
      return(
        list(
          initial = initial_df,
          ss = ss_df
        )
      )
    }
  )

#--------------------------------------------------#
# Calculate the outflow measures in baseline data:
#--------------------------------------------------#

baseline_outflow_initial <-
  lapply(
    gss_data_baseline,
    function(cohort){
      tm_num <-
        cohort$tm_num %>%
        mutate(Pair = str_replace_all(Pair, "[–—]", "-")) %>%    
        separate(Pair, into = c("status_p", "status_c"), 
                 sep = "-", convert = TRUE) %>%
        mutate(
          status_p = factor(status_p, levels = 1:5),
          status_c = factor(status_c, levels = 1:5)
        ) %>%
        group_by(status_p, status_c) %>%
        summarise(value = first(mean), .groups = "drop") %>%      
        complete(status_p, status_c, fill = list(value = 0)) %>%
        pivot_wider(names_from = status_c, values_from = value) %>%
        arrange(status_p) %>%
        column_to_rownames("status_p") %>%
        as.matrix()
      
      # Calculate outflow:
      tm_perc <- tm_num / rowSums(tm_num)
      M <- matrix(1, nrow = 5, ncol = 5)
      diag(M) <- 0
      outflow_initial <- 
        tibble(
          class = c("Class1","Class2","Class3","Class4","Class5"),
          mean_baseline = diag(M %*% t(as.matrix(tm_perc)))
        )
      return(outflow_initial)
    }
  )

baseline_outflow_ss <-
  lapply(
    gss_data_baseline,
    function(cohort){
      tm_num <-
        cohort$steady_state_TM %>%
        mutate(Pair = str_replace_all(Pair, "[–—]", "-")) %>%    
        separate(Pair, into = c("status_p", "status_c"), 
                 sep = "-", convert = TRUE) %>%
        mutate(
          status_p = factor(status_p, levels = 1:5),
          status_c = factor(status_c, levels = 1:5)
        ) %>%
        group_by(status_p, status_c) %>%
        summarise(value = first(mean), .groups = "drop") %>%      
        complete(status_p, status_c, fill = list(value = 0)) %>%
        pivot_wider(names_from = status_c, values_from = value) %>%
        arrange(status_p) %>%
        column_to_rownames("status_p") %>%
        as.matrix()
      
      # Calculate outflow:
      tm_perc <- tm_num / rowSums(tm_num)
      M <- matrix(1, nrow = 5, ncol = 5)
      diag(M) <- 0
      outflow_ss <- 
        tibble(
          class = c("Class1","Class2","Class3","Class4","Class5"),
          mean_baseline = diag(M %*% t(as.matrix(tm_perc)))
        )
      return(outflow_ss)
    }
  )


#--------------------------------------------------#
# Calculate the outflow  bias-corrected measures:
#--------------------------------------------------#

outflow_initial_bc <-
  lapply(
    seq_len(length(boot_outflow_stats)),
    function(i){
      boot_df <-
        boot_outflow_stats[[i]]$initial %>%
        left_join(
          baseline_outflow_initial[[i]],
          by = "class"
        ) %>%
        mutate(
          mean_bc = 2 * mean_baseline - mean_boot,
          CI_bc_l = 2 * mean_baseline - CI_np_u_boot,
          CI_bc_u = 2 * mean_baseline - CI_np_l_boot
        ) 
    }
  )

outflow_ss_bc <-
  lapply(
    seq_len(length(boot_outflow_stats)),
    function(i){
      boot_df <-
        boot_outflow_stats[[i]]$ss %>%
        left_join(
          baseline_outflow_ss[[i]],
          by = "class"
        ) %>%
        mutate(
          mean_bc = 2 * mean_baseline - mean_boot,
          CI_bc_l = 2 * mean_baseline - CI_np_u_boot,
          CI_bc_u = 2 * mean_baseline - CI_np_l_boot
        ) 
    }
  )

outflow_bc <-
  list(
    initial = outflow_initial_bc,
    ss = outflow_ss_bc
  )

#---------------------------------#
#  Generate the outflow at t = 1:
#---------------------------------#

outflow_initial <-
  plot_outflow_contrast(
    data = bind_rows(
      outflow_bc$initial,
      .id = "cohort_index"
    ),
    y_tick_vals = seq(0.45,1,0.05),
    y_range = c(0.45,1) ,
    y_tick_text = c(to_text(seq(0.45,0.95,0.05)),"1"), 
    legend_x = 0.95, 
    legend_y = 0.9, 
    x_tick_vals = seq(1945,1990,3) , 
    x_range = c(1944,1991)) %>% 
  layout(
    width = 1000, 
    height = 600
    )

#---------------------------------#
#  Generate the outflow at t = 1:
#---------------------------------#

outflow_ss <-
  plot_outflow_contrast(
    data = bind_rows(
      outflow_bc$ss,
      .id = "cohort_index"
    ),
    y_tick_vals = seq(0.45,1,0.05),
    y_range = c(0.45,1) ,
    y_tick_text = c(to_text(seq(0.45,0.95,0.05)),"1"), 
    legend_x = 0.95, 
    legend_y = 0.9, 
    x_tick_vals = seq(1945,1990,3) , 
    x_range = c(1944,1991)) %>% 
  layout(
    width = 1000, 
    height = 600
    )

#===========================================================#
#      Fig D2: PLOT THE MEAN FIRST PASSAGE TIME             #
#===========================================================#

MFP_updown_new <-
  MFP_class_comparison_new(gss_data_bc, "mean", 5, 10, 15, FALSE, NULL)$updown +
  labs(title = "")

#===========================================================#
#   Fig D3: PLOT THE ALTERNATIVE MEASURES OF MEMORY         #
#===========================================================#

#---------------------------------#
# Maximum distance to stationarity:
#---------------------------------#

dprime_t <-
  plot_d(
    data = gss_data_bc,
    t_num = 1,
    size = 5,
    y_title = "d",
    y_tick_vals = seq(0,0.4,0.04),
    x_tick_vals = seq(1945,1990,3),
    x_range = c(1944,1991),
    y_range = c(0,0.4),
    y_tick_text = to_text(seq(0,0.4,0.04))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
    )


#---------------------------------#
# Second largest eigenvalue modulus:
#---------------------------------#

config(
lambda2 <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "lambda2_df",
    attribute_col = "Attribute",
    metric_name = "lambda2",
    y_title = TeX("\\lvert \\lambda_{2} \\rvert"),
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.1, 0.5, by = 0.04),
    x_range = c(1944, 1991),
    y_range = c(0.1,0.5),
    y_tick_text = to_text(seq(0.1, 0.5, by = 0.04))
  ), mathjax = "cdn"
) %>% 
  layout(
    width = 1000, 
    height = 600
    )

#===========================================================#
#   Fig D4: ALTERNATIVE MEASURES OF DEPENDENCE
#===========================================================#

#---------------------------------#
# Altham Index  :
#---------------------------------#

altham_index <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "altham",
    attribute_col = "Attribute",
    metric_name = "altham",
    y_title = "Altham index",
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(10, 30, by = 3),
    x_range = c(1944, 1991),
    y_range = c(10,28),
    y_tick_text = to_text(seq(10, 30, by = 3))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
    )

#---------------------------------#
# Hellinger Dependence :
#---------------------------------#


Hellinger_Dependece <-
  plot_single_measure(
    data = gss_data_bc,
    t_num = c(1),
    measure = "HellingersDep",
    attribute_col = "Attribute",
    metric_name = "HellingerDep",
    y_title = "Hellinger dependence index",
    x_tick_vals = seq(1945,1990,3),
    y_tick_vals = seq(0.10, 0.50, by = 0.04),
    x_range = c(1944, 1991),
    y_range = c(0.10,0.50),
    y_tick_text = to_text(seq(0.10, 0.50, by = 0.04))
  ) %>% 
  layout(
    width = 1000, 
    height = 600
    )

#-----------------------------#
#  Save the Plot:
#-----------------------------#

# 1. Outflow Initial
saveRDS(
  outflow_initial,
  file = path.expand(
    file.path(
      dir_plot,
      "figE1_outflow_initial.rds"
    )
  )
)

# 2. Outflow SS
saveRDS(
  outflow_ss,
  file = path.expand(
    file.path(
      dir_plot,
      "figE1_outflow_ss.rds"
    )
  )
)

# 3. MFPT (Kept as ggsave/PDF)
ggsave(
  filename = file.path(
    dir_plot,
    "figE2_MFPT.pdf"
  ),
  plot     = MFP_updown_new,
  width    = 17,
  height   = 12,
  units    = "in",
  device   = cairo_pdf,
  bg       = "white"
)

# 4. dprime_t
saveRDS(
  dprime_t,
  file = path.expand(
    file.path(
      dir_plot,
      "figE3_dprime_t.rds"
    )
  )
)

# 5. Lambda2
saveRDS(
  lambda2,
  file = path.expand(
    file.path(
      dir_plot,
      "figE3_lambda2.rds"
    )
  )
)

# 6. Altham Index
saveRDS(
  altham_index,
  file = path.expand(
    file.path(
      dir_plot,
      "figE4_altham_index.rds"
    )
  )
)

# 7. Hellinger Dependence
saveRDS(
  Hellinger_Dependece,
  file = path.expand(
    file.path(
      dir_plot,
      "figE5_HD.rds"
    )
  )
)
# Open log：
sink(log_path, split = TRUE)

cat("Saved plots to: ", dir_plot ,"\n")
cat("Figure E1-D4 of the paper’s main results have been successfully replicated.\n")

sink()
