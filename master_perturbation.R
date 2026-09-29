# Set Working Directory:
rm(list = ls())
HPC <- TRUE
if(HPC){
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
}else{
  dir_root <- "~/Desktop/synthetic_dynasty_analysis"
}
setwd(path.expand(dir_root))
# Note: each numbered script carries its own HPC flag in its
# preliminaries and resets the workspace; keep those flags consistent
# with the flag set here. Relative paths below resolve from dir_root,
# which every numbered script restores via setwd().
source("Robust_Codes/Perturbation/01-calculate-perturbation-exercise.R")
print("Perturbation Exercise Calculated (Appendix J estimation)!")
source("Robust_Codes/Perturbation/02-generate-perturbation-results-figs.R")
print("Fig K.1-K.3 Sub-figures Generated!")
source("Robust_Codes/Perturbation/03-calculate-perturbation-results-tables.R")
print("Table J.1-J.3 Generated!")
print("Perturbation Replication has finished!")
