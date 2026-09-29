# ============================================================================
# master_nlsy_pipeline.R — NLSY validation pipeline (synthetic_dynasty layout)
# ----------------------------------------------------------------------------
# Deploy: place Robust_Codes/NLSY/ (these scripts) at dir_root, the raw-data
# folder at dir_root/Data/NLSY/, and run:   Rscript master_nlsy_pipeline.R
# GSS benchmark is read from Data/main_results/main_rst_bc.rds (shipped).
#
#   00  validate canonical inputs (fail loudly on a missing/incomplete input)
#   01  Data/NLSY raw -> complete occupation dyads          [Data/NLSY_estimation/dyads]
#   02  Data/NLSY raw -> complete income dyads (baseline)    [Data/NLSY_estimation/dyads]
#   02b NLSY97 non-relative clean-household screen (cleans dyads BEFORE weighting)
#   03  samples + ID lists + attach weights -> weighted finals
#   04  finals -> 5 bootstrap inputs                        [Data/NLSY_estimation/inputs]
#   05  bootstrap (HPC; single multi-core mclapply job — submit via your SLURM)
#   06  process bootstrap -> bias-corrected CIs + precomputes
#   07  Section-2/3 figures (plot.R cosmetics) -> Plot/robustness_check/NLSY/fig*.rds
#   08  render fig*.rds -> paper-ready PDF + 300-dpi PNG
# ============================================================================
HPC <- TRUE                 # TRUE on the cluster
if (HPC) {
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
} else {
  dir_root <- "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
}
# NLSY_PROJECT_ROOT overrides both (portable across machines / the Midway path).
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT")
setwd(dir_root)
Sys.setenv(NLSY_PROJECT_ROOT = dir_root)   # child scripts inherit the same root

RUN_00 <- TRUE;  RUN_01  <- TRUE;  RUN_02 <- TRUE;  RUN_02B <- TRUE
RUN_03 <- TRUE;  RUN_04  <- TRUE
RUN_05 <- TRUE  # bootstrap re-run (hours; usually via your SLURM job)
RUN_06 <- TRUE;  RUN_07  <- TRUE;  RUN_08 <- FALSE

run_r <- function(f) { message("\n########## ", f, " ##########")
  if (system2("Rscript", c("--vanilla", shQuote(file.path(dir_root, "Robust_Codes", "NLSY", f)))) != 0)
    stop("stage failed: ", f) }

if (RUN_00)  run_r("00_validate_inputs.R")
if (RUN_01)  run_r("01_generate_complete_occ_dyads.R")
if (RUN_02)  run_r("02_generate_complete_income_dyads.R")
if (RUN_02B) run_r("02b_apply_nlsy97_nonrelative_clean_household.R")
if (RUN_03)  run_r("03_attach_weights.R")
if (RUN_04)  run_r("04_generate_bootstrap_inputs.R")
if (RUN_05)  run_r("05_hpc_bootstrap.R")
if (RUN_06)  run_r("06_process_bootstrap.R")
if (RUN_07)  run_r("07_generate_results.R")
if (RUN_08)  run_r("08_render_figures.R")
message("\nNLSY pipeline finished. Figures: Plot/robustness_check/NLSY/ (fig*.rds, pdf/, png/)")
