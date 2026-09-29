#----------Preliminaries----------#
rm(list = ls())
section <- "Data"
subsection <- "Perturbation_estimation"
title <- "calculate_perturbation_exercise"

HPC <- TRUE # Set FALSE to run locally with downloaded Data/main_results inputs

if (HPC) {
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
} else {
  dir_root <- "~/Desktop/synthetic_dynasty_analysis"
}

dir_root <- path.expand(dir_root)
setwd(dir_root)

dir_log <- file.path(dir_root, "code", "_LOGS", section)
log_path <- file.path(
  dir_log,
  paste0(subsection, "_", title, "_log.txt")
)

dir_main_results <- file.path(dir_root, "Data", "main_results")
dir_perturbation_data <- file.path(dir_root, "Data", "Perturbation_estimation")
dir_perturbation_plot <- file.path(dir_root, "Plot", "Perturbation")
dir_perturbation_fun <- file.path(dir_root, "perturbation_fun")

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
create_dir_if_missing(dir_perturbation_data)
create_dir_if_missing(dir_perturbation_plot)

options(warn = -1)


#-------------------------------------------------------------------------------
# U.S. Occupational Mobility — Replication Files
#
# Project:     U.S. Occupational Mobility Analysis
# Repository:  https://github.com/synthetic_dynasty_analysis
#
# Script:      01-calculate-perturbation-exercise.R
#
# Inputs:
#   - Stored Midway main-analysis objects (read-only):
#       Data/main_results/gss_fc_occ10_5class.rds
#       Data/main_results/main_rst_baseline.rds
#       Data/main_results/main_rst_boot.rds
#       Data/main_results/main_rst_bc.rds
#   - Paper estimation functions:
#       Functions/synthetic_dynasty_estimation.R
#   - Perturbation production helpers:
#       perturbation_fun/perturbation_inputs.R
#       perturbation_fun/perturbation_mechanisms.R
#       perturbation_fun/perturbation_measures.R
#       perturbation_fun/perturbation_trends_bootstrap.R
#
# Outputs:
#   - Log file:
#       code/_LOGS/Data/Perturbation_estimation_calculate_perturbation_exercise_log.txt
#   - Canonical unrounded estimation object:
#       Data/Perturbation_estimation/perturbation_results.rds
#
# Description:
#   Appendix J (Sparsity and Estimation Uncertainty): every substantive
#   deterministic and bootstrap calculation behind Tables J.1-J.3 and
#   Figures J.1-J.3. A fixed sparse-cell mask (raw dyad count <= 5) is
#   built once from the deterministic production sample; the two paper
#   perturbation mechanisms (within-row redistribution and augmented
#   pseudo-dyad table) are applied at kappa = 1; all 20 paper measures
#   are recomputed with the production estimators; full- and late-window
#   per-decade OLS trends are estimated for the baseline and both
#   mechanisms; the 2,000 stored GAM-smoothed bootstrap replicates are
#   processed with replicate identity preserved and paired within
#   replicate; stationary-TV design calibration (Table J.1), slope
#   calibration (Table J.3 / Figure J.2), and cohort-specific OLS slope
#   contributions (Figure J.3) are computed; concise paper-benchmark
#   checks are run; and one canonical unrounded result object is saved.
#-------------------------------------------------------------------------------


#-------------------------------------------------#
#  INSTALL/LOAD DEPENDENCIES                       #
#-------------------------------------------------#
if (HPC) {
  .libPaths(c("/home/weiqiw/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  library("dplyr")
  library("tidyr")
  library("tibble")
  library("haven")
  library("assertthat")
} else {
  packages <- c("dplyr", "tidyr", "tibble", "haven", "assertthat")
  installed <- packages %in% rownames(installed.packages())
  if (any(!installed)) install.packages(packages[!installed])
  invisible(lapply(packages, library, character.only = TRUE))
}

t_start <- Sys.time()

#-------------------------------------------------#
#  SOURCE PRODUCTION FUNCTIONS AND HELPERS         #
#-------------------------------------------------#
source(file.path(dir_root, "Functions", "synthetic_dynasty_estimation.R"))

source(file.path(dir_perturbation_fun, "perturbation_inputs.R"))
source(file.path(dir_perturbation_fun, "perturbation_mechanisms.R"))
source(file.path(dir_perturbation_fun, "perturbation_measures.R"))
source(file.path(dir_perturbation_fun, "perturbation_trends_bootstrap.R"))

#-------------------------------------------------#
#  SECTION 1: LOAD AND VALIDATE INPUTS             #
#-------------------------------------------------#
cfg <- get_perturbation_config()
dictionary <- get_perturbation_measure_dictionary()
inputs <- load_perturbation_inputs(cfg, dir_main_results)
std <- standardize_perturbation_inputs(inputs, cfg)
message("Inputs validated: n = ", std$n_total, "; cohort n range = [",
        min(std$n_by_cohort), ", ", max(std$n_by_cohort), "]")

#-------------------------------------------------#
#  SECTION 2: FIXED SPARSE-CELL MASK               #
#-------------------------------------------------#
mask_obj <- build_fixed_sparse_mask(std$raw_mats, cfg)
masks <- mask_obj$masks

#-------------------------------------------------#
#  SECTION 3: DETERMINISTIC PERTURBATIONS +        #
#             ALL 20 PAPER MEASURES                #
#-------------------------------------------------#
message("Deterministic perturbations and paper measures (kappa = ", cfg$kappa, ")")
cohort_objects <- lapply(names(std$smoothed), function(nm) {
  s <- std$smoothed[[nm]]; r <- std$raw_mats[[nm]]; A <- masks[[nm]]
  res <- list(
    Baseline = list(
      transition_matrix = s$P,
      father_origin_marginal = s$mu0,
      child_destination_marginal = as.numeric(s$mu0 %*% s$P),
      stationary_distribution = stationary_distribution_exact(s$P),
      measures = calculate_reported_paper_measures(s$P, s$mu0, cfg$steady_state_tol))
  )
  for (mech in cfg$mechanisms) {
    p <- apply_perturbation_to_cohort(mech, s$P, s$mu0, r$n_ic, r$N_c, A, cfg$kappa)
    stopifnot(!p$capped)  # zero capped cohorts at kappa = 1 (paper contract)
    res[[mech]] <- list(
      transition_matrix = p$transition_matrix,
      father_origin_marginal = p$father_origin_marginal,
      child_destination_marginal = p$child_destination_marginal,
      stationary_distribution = stationary_distribution_exact(p$transition_matrix),
      measures = calculate_reported_paper_measures(
        p$transition_matrix, p$father_origin_marginal, cfg$steady_state_tol))
  }
  list(cohort = s$cohort, results = res)
})
names(cohort_objects) <- names(std$smoothed)
measures_long <- calculate_cohort_measure_results(cohort_objects, dictionary, cfg$kappa)

#-------------------------------------------------#
#  SECTION 4: DETERMINISTIC TRENDS                 #
#-------------------------------------------------#
deterministic_trends <- calculate_deterministic_trends(measures_long, dictionary, cfg)

#-------------------------------------------------#
#  SECTION 5: PAIRED BOOTSTRAP (2,000 REPLICATES)  #
#-------------------------------------------------#
message("Paired bootstrap trends (", cfg$n_boot, " stored replicates)")
bootstrap_trends <- calculate_paired_bootstrap_trends(
  inputs$boot, std, masks, dictionary, cfg, deterministic_trends)
message("Targeted slope bootstrap")
targeted <- calculate_targeted_slope_bootstrap(
  inputs$boot, std, masks, dictionary, cfg, deterministic_trends)
bootstrap_trends <- dplyr::left_join(
  bootstrap_trends, targeted,
  by = c("replicate_id", "window", "measure_id", "mechanism"))
paired_summary <- summarize_paired_delta_beta(bootstrap_trends)
stopifnot(all(paired_summary$n_replicates == cfg$n_boot))

#-------------------------------------------------#
#  SECTION 6: TABLE J.1 DESIGN CALIBRATION         #
#-------------------------------------------------#
message("Stationary-TV design calibration (Table J.1)")
table_J1_data <- calculate_stationary_design_calibration(
  std, masks, mask_obj$inventory, inputs$boot, cfg)

#-------------------------------------------------#
#  SECTION 7: TABLE J.2 / J.3 / FIGURE DATA        #
#-------------------------------------------------#
table_J2_data <- dplyr::left_join(deterministic_trends, paired_summary,
  by = c("window", "measure_id", "mechanism"))

table_J3_data <- calculate_slope_calibration(bootstrap_trends, deterministic_trends)
table_J3_data <- dplyr::left_join(table_J3_data,
  paired_summary[, c("window", "measure_id", "mechanism",
                     "delta_beta_lo", "delta_beta_hi")],
  by = c("window", "measure_id", "mechanism"))

fig_base <- dplyr::left_join(table_J2_data,
  dictionary[, c("measure_id", "decade_factor")], by = "measure_id")
figure_J1_data <- dplyr::mutate(fig_base,
  delta_beta_display = delta_beta * decade_factor,
  delta_beta_lo_display = delta_beta_lo * decade_factor,
  delta_beta_hi_display = delta_beta_hi * decade_factor)

figure_J2_data <- dplyr::left_join(figure_J1_data,
  table_J3_data[, c("window", "measure_id", "mechanism", "q95_total", "q95_target",
                    "u_total_lo", "u_total_hi", "u_target_lo", "u_target_hi")],
  by = c("window", "measure_id", "mechanism"))
figure_J2_data <- dplyr::mutate(figure_J2_data,
  q95_total_display = q95_total * decade_factor,
  q95_target_display = q95_target * decade_factor,
  u_total_lo_display = u_total_lo * decade_factor,
  u_total_hi_display = u_total_hi * decade_factor,
  u_target_lo_display = u_target_lo * decade_factor,
  u_target_hi_display = u_target_hi * decade_factor)

contributions <- calculate_cohort_contributions(
  measures_long, dictionary, cfg, deterministic_trends)
figure_J3_data <- dplyr::mutate(
  dplyr::left_join(contributions,
    dictionary[, c("measure_id", "decade_factor")], by = "measure_id"),
  q_delta_display = q_delta * decade_factor)

#-------------------------------------------------#
#  SECTION 8: PAPER-BENCHMARK CHECKS               #
#-------------------------------------------------#
# The values hard-coded below are BENCHMARKS, not computational inputs:
# they are the rounded numbers printed in Appendix J of the paper, coded
# here so that every rerun of this script is checked against the
# published results. Nothing upstream reads these values — Sections 1-7
# compute everything at full precision from the raw Data/main_results/
# inputs alone. Only AFTER all computation is finished are the fresh
# results rounded to the paper's display precision and compared with
# these anchors. A mismatch stops the script BEFORE the canonical result
# object is saved, signalling that something upstream (input vintage,
# sample filter, mechanism, measure function, trend fit, or bootstrap
# pairing) no longer reproduces the paper — it is a tripwire that
# detects replication failure; it cannot create replication.
benchmark_failures <- character(0)
chk <- function(label, value, expected, digits) {
  value <- unname(value)   # drop names (e.g., from fac[m]) before comparison
  if (!isTRUE(all.equal(round(value, digits), expected, tolerance = 1e-9))) {
    benchmark_failures <<- c(benchmark_failures,
      sprintf("%s: got %s, expected %s", label,
              format(round(value, digits), nsmall = digits), format(expected)))
  }
}
g_tr <- function(w, m, mech, col) {
  d <- deterministic_trends
  d[[col]][d$window == w & d$measure_id == m & d$mechanism == mech]
}
fac <- stats::setNames(dictionary$decade_factor, dictionary$measure_id)

## Table J.2 full-period anchors: the baseline / redistribution /
## augmentation slopes exactly as printed in the paper's Table J.2
## (Panel A), used only to benchmark the freshly computed results
j2_anchors <- list(
  historical = c(-0.530, -0.477, -0.411), structural = c(-2.111, -2.432, -2.426),
  exchange = c(1.581, 1.955, 2.015), upward = c(-1.206, -1.014, -1.137),
  downward = c(0.676, 0.537, 0.725), AMTE = c(0.009, 0.008, 0.007),
  AIM = c(0.002, 0.002, 0.001), SSM = c(-0.367, -0.122, -0.146),
  ss_C4_share = c(-0.699, -0.339, -0.388), mfp_to_C4 = c(2.127, 0.666, 0.786))
for (m in names(j2_anchors)) {
  a <- j2_anchors[[m]]
  chk(paste0(m, " base"), g_tr("full", m, "Redistribution", "beta_base") * fac[m], a[1], 3)
  chk(paste0(m, " redist"), g_tr("full", m, "Redistribution", "beta_pert") * fac[m], a[2], 3)
  chk(paste0(m, " augment"), g_tr("full", m, "Augmentation", "beta_pert") * fac[m], a[3], 3)
}
## IM Class 1 full-period anchors
chk("IM_C1 base", g_tr("full", "IM_C1", "Redistribution", "beta_base") * fac["IM_C1"], 0.002, 3)
chk("IM_C1 redist", g_tr("full", "IM_C1", "Redistribution", "beta_pert") * fac["IM_C1"], -0.002, 3)
chk("IM_C1 augment", g_tr("full", "IM_C1", "Augmentation", "beta_pert") * fac["IM_C1"], -0.001, 3)
pdir <- function(m, mech) paired_summary$p_direction_retained[
  paired_summary$window == "full" & paired_summary$measure_id == m &
  paired_summary$mechanism == mech]
chk("IM_C1 P(dir) redist", pdir("IM_C1", "Redistribution"), 0.41, 2)
chk("IM_C1 P(dir) augment", pdir("IM_C1", "Augmentation"), 0.46, 2)
## slope-calibration R_total anchors
rtt <- function(m, mech) table_J3_data$R_total[
  table_J3_data$window == "full" & table_J3_data$measure_id == m &
  table_J3_data$mechanism == mech]
chk("ss_C4 R_total redist", rtt("ss_C4_share", "Redistribution"), 1.35, 2)
chk("ss_C4 R_total augment", rtt("ss_C4_share", "Augmentation"), 1.16, 2)
chk("mfp_C4 R_total redist", rtt("mfp_to_C4", "Redistribution"), 1.09, 2)
chk("mfp_C4 R_total augment", rtt("mfp_to_C4", "Augmentation"), 1.00, 2)
## late-period anchors: the sensitive Class-4 long-run results exactly
## as printed in the paper's Table J.2 (Panel B)
chk("ss_C4 late base", g_tr("late", "ss_C4_share", "Redistribution", "beta_base") * 100, -0.410, 3)
chk("ss_C4 late redist", g_tr("late", "ss_C4_share", "Redistribution", "beta_pert") * 100, 0.704, 3)
chk("ss_C4 late augment", g_tr("late", "ss_C4_share", "Augmentation", "beta_pert") * 100, 0.452, 3)
chk("ss_C4 late ddelta redist", g_tr("late", "ss_C4_share", "Redistribution", "delta_beta") * 100, 1.113, 3)
chk("ss_C4 late ddelta augment", g_tr("late", "ss_C4_share", "Augmentation", "delta_beta") * 100, 0.861, 3)
chk("mfp late base", g_tr("late", "mfp_to_C4", "Redistribution", "beta_base"), 1.284, 3)
chk("mfp late redist", g_tr("late", "mfp_to_C4", "Redistribution", "beta_pert"), -1.980, 3)
chk("mfp late augment", g_tr("late", "mfp_to_C4", "Augmentation", "beta_pert"), -1.530, 3)
chk("mfp late ddelta redist", g_tr("late", "mfp_to_C4", "Redistribution", "delta_beta"), -3.264, 3)
chk("mfp late ddelta augment", g_tr("late", "mfp_to_C4", "Augmentation", "delta_beta"), -2.815, 3)
## Table J.1 era anchors: the four design/calibration rows exactly as
## printed in the paper's Table J.1, used only to benchmark the freshly
## computed medians
j1 <- dplyr::summarise(dplyr::group_by(table_J1_data, mechanism, era),
  cells = stats::median(n_sparse_cells), rows = stats::median(n_affected_rows),
  tv = stats::median(stationary_tv_k1), q50 = stats::median(total_q50),
  q95 = stats::median(total_q95), R = stats::median(R_total_stationary),
  kb = stats::median(kappa_boot95_total, na.rm = TRUE),
  tvmu0 = stats::median(tv_mu0), tvmu1 = stats::median(tv_mu1),
  capped = sum(capped), .groups = "drop")
g_j1 <- function(mech, e, col) j1[[col]][j1$mechanism == mech & j1$era == e]
j1_anchors <- list(
  list("Redistribution", "early", 2, 2, 0.0053, 0.0105, 0.0191, 0.269, 3.80, 0.0000, 0.0038),
  list("Redistribution", "late", 6, 4, 0.0195, 0.0159, 0.0292, 0.702, 1.42, 0.0000, 0.0165),
  list("Augmentation", "early", 2, 2, 0.0044, 0.0105, 0.0191, 0.238, 4.56, 0.0026, 0.0032),
  list("Augmentation", "late", 6, 4, 0.0154, 0.0159, 0.0292, 0.538, 1.92, 0.0063, 0.0149))
for (a in j1_anchors) {
  pre <- paste("J1", a[[1]], a[[2]])
  chk(paste(pre, "cells"), g_j1(a[[1]], a[[2]], "cells"), a[[3]], 0)
  chk(paste(pre, "rows"), g_j1(a[[1]], a[[2]], "rows"), a[[4]], 0)
  chk(paste(pre, "tv"), g_j1(a[[1]], a[[2]], "tv"), a[[5]], 4)
  chk(paste(pre, "q50"), g_j1(a[[1]], a[[2]], "q50"), a[[6]], 4)
  chk(paste(pre, "q95"), g_j1(a[[1]], a[[2]], "q95"), a[[7]], 4)
  chk(paste(pre, "R"), g_j1(a[[1]], a[[2]], "R"), a[[8]], 3)
  chk(paste(pre, "kb"), g_j1(a[[1]], a[[2]], "kb"), a[[9]], 2)
  chk(paste(pre, "tvmu0"), g_j1(a[[1]], a[[2]], "tvmu0"), a[[10]], 4)
  chk(paste(pre, "tvmu1"), g_j1(a[[1]], a[[2]], "tvmu1"), a[[11]], 4)
  chk(paste(pre, "capped"), g_j1(a[[1]], a[[2]], "capped"), 0, 0)
}
## headline-direction check: none of the ten headline full-window trend
## directions may change under either mechanism
headline <- c("historical", "structural", "exchange", "upward", "downward",
              "AMTE", "AIM", "SSM", "ss_C4_share", "mfp_to_C4")
for (m in headline) for (mech in cfg$mechanisms) {
  if (sign(g_tr("full", m, mech, "beta_pert")) !=
      sign(g_tr("full", m, mech, "beta_base"))) {
    benchmark_failures <- c(benchmark_failures,
      sprintf("headline direction changed: %s under %s", m, mech))
  }
}

benchmark_pass <- length(benchmark_failures) == 0
if (!benchmark_pass) {
  message("PAPER-BENCHMARK FAILURES:")
  for (f in benchmark_failures) message("  ", f)
  stop("Paper-benchmark checks failed (", length(benchmark_failures),
       " mismatches). The canonical result object was NOT saved. Identify the ",
       "first upstream mismatch (input sample, sparse mask, mechanism, ",
       "production measure, trend fit, bootstrap pairing, or formatting) ",
       "before replacing paper-facing artifacts.")
}
message("Paper-benchmark checks: PASS")

#-------------------------------------------------#
#  SECTION 9: SAVE CANONICAL RESULT OBJECT         #
#-------------------------------------------------#
cohort_results <- list(
  measures = measures_long,
  objects = cohort_objects
)
perturbation_results <- list(
  metadata = list(
    created_at = Sys.time(),
    sample_age_min = cfg$age_min, sample_age_max = cfg$age_max,
    sample_n = std$n_total,
    class_order = cfg$class_labels,
    cohorts = cfg$cohorts,
    early_era = cfg$early_era, late_era = cfg$late_era,
    full_window = cfg$full_window, late_window = cfg$late_window,
    kappa = cfg$kappa, sparse_threshold = cfg$sparse_threshold,
    bootstrap_replicates = cfg$n_boot,
    input_files = inputs$input_paths,
    benchmark_pass = benchmark_pass
  ),
  measure_dictionary = dictionary,
  sparse_mask = list(masks = masks, inventory = mask_obj$inventory),
  cohort_results = cohort_results,
  deterministic_trends = deterministic_trends,
  bootstrap_trends = bootstrap_trends,
  table_J1_data = table_J1_data,
  table_J2_data = table_J2_data,
  table_J3_data = table_J3_data,
  figure_J1_data = figure_J1_data,
  figure_J2_data = figure_J2_data,
  figure_J3_data = figure_J3_data
)
saveRDS(perturbation_results,
        file.path(dir_perturbation_data, "perturbation_results.rds"))

#-------------------------------------------------#
#  LOG                                             #
#-------------------------------------------------#
sink(log_path, split = TRUE)
cat("Script: 01-calculate-perturbation-exercise.R\n")
cat("Mode: ", ifelse(HPC, "HPC (Midway)", "local"), "\n")
cat("Inputs: ", paste(basename(inputs$input_paths), collapse = ", "), "\n")
cat("Output: ", file.path(dir_perturbation_data, "perturbation_results.rds"), "\n")
cat("Cohorts: ", length(cfg$cohorts), " (", min(cfg$cohorts), "-", max(cfg$cohorts), ")\n", sep = "")
cat("Classes: ", cfg$n_classes, "\n")
cat("Bootstrap replicates processed: ", cfg$n_boot, "\n")
cat("Measures: ", nrow(dictionary), "\n")
cat("Elapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " minutes\n", sep = "")
cat("Paper-benchmark checks: ", ifelse(benchmark_pass, "PASS", "FAIL"), "\n")
cat("Perturbation exercise completed successfully.\n")
sink()
