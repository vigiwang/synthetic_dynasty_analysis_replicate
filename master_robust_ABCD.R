# Set Working Directory:
rm(list = ls())
HPC <- TRUE
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "~/Desktop/synthetic_dynasty_analysis" 
}
setwd(dir_root)
source("Robust_Codes/Age/RA1_generate_age2564_estimation.R")
print("Age 25-64 Robustness Check Estimation Finished!")

source("Robust_Codes/Age/RA2_generate_age2564_plots.R")
print("Age 25-64 Robustness Check Plots Generation Finished!")

source("Robust_Codes/Age/RA3_generate_age3064_estmation.R")
print("Age 30-64 Robustness Check Estimation Finished!")

source("Robust_Codes/Age/RA4_generate_age3064_plots.R")
print("Age 30-64 Robustness Check Plots Generation Finished!")

source("Robust_Codes/Age/RA5_generate_age3040_estimation.R")
print("Age 30-40 Robustness Check Estimation Finished!")

source("Robust_Codes/Age/RA6_generate_age3040_plot.R")
print("Age 30-40 Robustness Check Plots Generation Finished!")

source("Robust_Codes/Class_Typology/RB1_generate_class_typology_data.R")
print("Alternative Class Typology Data Generation Finished!")

source("Robust_Codes/Class_Typology/RB2_generate_6class_estimation.R")
print("6 Class Estimation Finished!")

source("Robust_Codes/Class_Typology/RB3_generate_6class_plots.R")
print("6 Class Plots Generation Finished!")

source("Robust_Codes/Class_Typology/RB4_generate_7class_estimation.R")
print("7 Class Estimation Finished!")

source("Robust_Codes/Class_Typology/RB5_generate_7class_plots.R")
print("7 Class Plots Generation Finished!")

source("Robust_Codes/Class_Typology/RB6_generate_5class_noself_estimation.R")
print("5 Class(without self-employed) Estimation Finished!")

source("Robust_Codes/Class_Typology/RB7_generate_noself_plots.R")
print("5 Class No-self Plots Generation Finished!")

source("Robust_Codes/Class_Typology/RB8_generate_5class_nofarmer_estimation.R")
print("5 Class(without farmer) Finished!")

source("Robust_Codes/Class_Typology/RB9_generate_nofarmer_plots.R")
print("5 Class No-farmer Plots Generation Finished!")

source("Robust_Codes/Weighted/RC1_generate_weighted_estimation.R")
print("Weighted-sample estimation finished!")

source("Robust_Codes/Weighted/RC2_generate_weighted_plots.R")
print("Weighted-sample Plots Generation Finished!")

source("Robust_Codes/Alternative_Measure/RD_generate_alternative_measure_plots.R")
print("Alternative-Measure Plots Generation Finished!")

print("Main Results Estimation has finished!")