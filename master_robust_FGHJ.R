# Set Working Directory:
rm(list = ls())
HPC <- TRUE
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "~/Desktop/synthetic_dynasty_analysis" 
}
setwd(dir_root)
source("Robust_Codes/Gender/RF1_generate_gender_analysis_data.R")
print("Gender Robustness Check Data Generation Finished!")

source("Robust_Codes/Gender/RF2_generate_father_daughter_estimation.R")
print("Father-Daughter Estimation Finished!")

source("Robust_Codes/Gender/RF3_generate_father_daughter_plots.R")
print("Father-Daughter Plots Generated!")

source("Robust_Codes/Gender/RF4_generate_father_son_estimation.R")
print("Father-Son Estimation Finished!")

source("Robust_Codes/Gender/RF5_generate_father_son_plots.R")
print("Father-Son Plots Generated!")

source("Robust_Codes/Gender/RF6_generate_mc_estimation.R")
print("Mother-Child Estimation Finished!")

source("Robust_Codes/Gender/RF7_generate_mother_child_plots.R")
print("Mother-Child Plots Generated!")

source("Robust_Codes/Gender/RF8_generate_higher_status_child_estimation.R")
print("Higher-Status Parent-Child Estimation Finished!")

source("Robust_Codes/Gender/RF9_generate_higher_status_child_plots.R")
print("Higher Status Parent-Child Plots Generated!")

source("Robust_Codes/Alternative_Smoothing/RH1_generate_raw_estimation.R")
print("No-Smoothing Estimation Finished!")

source("Robust_Codes/Alternative_Smoothing/RH2_generate_raw_plots.R")
print("No-smoothing Plots Generated!")

source("Robust_Codes/Alternative_Smoothing/RH3_generate_bin5_estimation.R")
print("Bin = 5 Estimation Finished!")

source("Robust_Codes/Alternative_Smoothing/RH4_generate_bin5_plots.R")
print("Bin = 5 Plots Generated!")

source("Robust_Codes/Alternative_Smoothing/RH5_generate_bin10_estimation.R")
print("Bin = 10 Estimation Finished!")

source("Robust_Codes/Alternative_Smoothing/RH6_generate_bin10_plots.R")
print("Bin = 10 Plots Generated!")

source("Robust_Codes/Race/RG1_generate_white_estimation.R")
print("White-only Estimation Finished!")

source("Robust_Codes/Race/RG2_generate_white_plots.R")
print("White-only sample Plots Generated!")
 
source("Robust_Codes/Extended_Cohorts/RJ1_generate_extended_cohorts_estimation.R")
print("Extended-Cohorts Estimation Finished!")

source("Robust_Codes/Extended_Cohorts/RJ2_generate_extended_cohorts_plots.R")
print("Extended-Cohorts sample Plots Generated!")

