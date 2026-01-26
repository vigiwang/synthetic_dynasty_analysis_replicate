# Set Working Directory:
rm(list = ls())
HPC <- TRUE
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "~/Desktop/synthetic_dynasty_analysis" 
}
setwd(dir_root)
source("00_generate_main_results_data.R")
print("EGP 5-Class Typology Father-Child Pairs Generated!")
source("01_generate_main_results_estimation.R")
print("Benchmark Estimation Finished!")
source("02_create_fig1.R")
print("Fig 1 Generated!")
source("03_create_fig2.R")
print("Fig 2 Generated!")
source("04_create_fig3.R")
print("Fig 3 Generated!")
source("05_create_fig4.R")
print("Fig 4 Generated!")
source("06_create_fig5.R")
print("Fig 5 Generated!")
source("07_create_fig6.R")
print("Fig 6 Generated!")
source("08_create_fig7.R")
print("Fig 7 Generated!")
source("09_create_fig8.R")
print("Fig 8 Generated!")
source("10_create_fig9.R")
print("Fig 9 Generated!")
print("Main Results Estimation has finished!")
