# Set Working Directory:
rm(list = ls())
HPC <- TRUE
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "~/Desktop/synthetic_dynasty_analysis"
}
setwd(dir_root)
source("Robust_Codes/Education/RF15_generate_education_analysis_data.R")
print("Education Robustness Check Data Generation Finished!")

source("Robust_Codes/Education/RF16_generate_fc_education_estimation.R")
print("Father-Child Education Estimation Finished!")

source("Robust_Codes/Education/RF17_generate_mc_education_estimation.R")
print("Mother-Child Education Estimation Finished!")

source("Robust_Codes/Education/RF18_generate_ms_education_estimation.R")
print("Mother-Son Education Estimation Finished!")

source("Robust_Codes/Education/RF19_generate_md_education_estimation.R")
print("Mother-Daughter Education Estimation Finished!")

source("Robust_Codes/Education/RF20_generate_pc_education_estimation.R")
print("Higher-Education Parent-Child Estimation Finished!")

# The mother-son and mother-daughter occupational estimations below can be
# commented out if their outputs (5class_ms_bc.rds, 5class_md_bc.rds) have
# already been generated in a previous run:

source("Robust_Codes/Gender/RF1b_generate_mother_analysis_data.R")
print("Gender Robustness Check (Mother-Son/Daughter) Data Generation Finished!")

source("Robust_Codes/Gender/RF10_generate_mother_son_estimation.R")
print("Mother-Son Estimation Finished!")

source("Robust_Codes/Gender/RF11_generate_mother_daughter_estimation.R")
print("Mother-Daughter Estimation Finished!")

source("Robust_Codes/Gender/RF12_generate_mother_son_plots.R")
print("Mother-Son Plots Generated!")

source("Robust_Codes/Gender/RF13_generate_mother_daughter_plots.R")
print("Mother-Daughter Plots Generated!")

source("Robust_Codes/Education/RF21_generate_mc_education_plots.R")
print("Mother-Child Education Plots Generated!")

source("Robust_Codes/Education/RF22_generate_ms_education_plots.R")
print("Mother-Son Education Plots Generated!")

source("Robust_Codes/Education/RF23_generate_md_education_plots.R")
print("Mother-Daughter Education Plots Generated!")

source("Robust_Codes/Education/RF24_generate_pc_education_plots.R")
print("Higher-Education Parent-Child Education Plots Generated!")
