# ============================================================================
# 06_process_bootstrap.R
# ----------------------------------------------------------------------------
# Processes ALL bootstrap results: raw <sample>_{baseline,boot}.rds (10 samples)
# -> percentile stats for EVERY measure per cohort + "Full" -> the analyst's
# bias-corrected CIs (mean_bc = 2*base - boot_mean; reflected CI; Fay-BRR
# survey runs rescaled by k = 1/(1-rho) = 2), via the UNCHANGED
# Functions/process_results.R -> Data/NLSY_estimation/processed/<sample>_processed.rds
# (cached; delete a cache to force reprocessing) + scalar/SE summary CSVs.
# Then runs the nine report precomputes -> Data/NLSY_estimation/intermediate/.
# This is a script port of analysis/NLSY_bootstrap_process.Rmd (same code path).
# ============================================================================
#----------Preliminaries----------#
HPC <- TRUE   # set FALSE for a local test run
if (HPC) {
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
} else {
  dir_root <- "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
}
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)
source(file.path(dir_root, "NLSY_fun", "nlsy_config.R"))
suppressPackageStartupMessages({library(parallel); library(janitor)})  # tabyl() used by Functions/process_results.R
dir_boot <- file.path(root, "Data/NLSY_estimation")
dir_proc <- file.path(dir_boot, "processed")
if (!dir.exists(dir_proc)) dir.create(dir_proc)
n_cores  <- max(1L, detectCores() - 2L)
size     <- 5L
source(file.path(root, "Functions", "process_results.R"))   # analyst original, UNCHANGED

# ---- stale-cache invalidation: 05 writes results in place (Data/NLSY_estimation);
# a processed cache older than its raw boot results is dropped and rebuilt ----
for (f in list.files(dir_boot, pattern = "_boot\\.rds$")) {
  pc <- file.path(dir_proc, sub("_boot\\.rds$", "_processed.rds", f))
  if (file.exists(pc) && file.mtime(file.path(dir_boot, f)) > file.mtime(pc)) {
    file.remove(pc); message("[06] stale cache dropped (fresh bootstrap results): ", basename(pc))
  }
}

# ---- adapters (identical to the Rmd; analyst functions NOT modified) --------
patch_result <- function(r) {
  if (!is.list(r)) return(r)
  if (is.null(r$rho)) r$rho <- NA_real_
  tmn <- as.data.frame(r$tm_num)
  if (!"rowsum" %in% names(tmn)) tmn$rowsum <- tmn$validation
  r$tm_num <- tmn
  r
}
FAY_K <- 1
valid_bootstrap_results <- function(boot_rst_df, baseline_rst_df, merge_vec) {
  left_join(boot_rst_df, baseline_rst_df, by = merge_vec,
            suffix = c("_boot", "_baseline")) %>%
    mutate(
      mean_bc = mean_baseline - FAY_K * (mean_boot - mean_baseline),
      CI_bc_l = mean_baseline - FAY_K * (CI_np_u_boot - mean_baseline),
      CI_bc_u = mean_baseline - FAY_K * (CI_np_l_boot - mean_baseline),
      se_adj  = FAY_K * sd_boot
    ) %>%
    mutate(valid = ((mean_bc >= CI_bc_l - 1e-8) & (mean_bc <= CI_bc_u + 1e-8) |
                      mean_baseline == 1)) %>%
    dplyr::select(mean_baseline, mean_bc, CI_bc_u, CI_bc_l, everything())
}
load_run <- function(sample) {
  base <- readRDS(file.path(dir_boot, paste0(sample, "_baseline.rds")))
  boot <- readRDS(file.path(dir_boot, paste0(sample, "_boot.rds")))
  ok <- vapply(boot, function(rep) all(vapply(rep, is.list, logical(1))), logical(1))
  if (any(!ok)) message(sample, ": dropped ", sum(!ok), " replicate(s) with non-ergodic cohorts")
  list(base = lapply(base, patch_result),
       boot = lapply(boot[ok], function(rep) lapply(rep, patch_result)))
}
process_domain <- function(base, boot) {
  cohort_names <- setdiff(names(base), "Full")
  boot_stats_c <- process_results(lapply(boot, function(rep) rep[cohort_names]),
                                  length(boot), size, n_cores, TRUE)
  base_stats_c <- process_results(list(base[cohort_names]), 1, size, n_cores, TRUE)
  boot_stats_f <- process_results(lapply(boot, function(rep) setNames(rep["Full"], "9999")),
                                  length(boot), size, n_cores, TRUE)
  base_stats_f <- process_results(list(setNames(base["Full"], "9999")),
                                  1, size, n_cores, TRUE)
  boot_stats <- c(boot_stats_c, boot_stats_f)
  base_stats <- c(base_stats_c, base_stats_f)
  names(boot_stats) <- names(base_stats) <- c(cohort_names, "Full")
  list(boot_stats = boot_stats, baseline_stats = base_stats)
}

# ---- the ten runs -----------------------------------------------------------
runs <- tribble(
  ~sample,             ~domain,           ~design,                     ~fay_k,
  "nlsy79_income",     "NLSY79 income",   "household cluster",          1,
  "nlsy79_occ",        "NLSY79 occ",      "household cluster",          1,
  "nlsy97_income",     "NLSY97 income",   "survey Fay-BRR (rho = 0.5)", 2,
  "nlsy79_income_ind", "NLSY79 income",   "individual (naive)",         1,
  "nlsy79_occ_ind",    "NLSY79 occ",      "individual (naive)",         1,
  "nlsy97_income_ind", "NLSY97 income",   "individual (naive)",         1,
  "nlsy79_income_unadj",     "NLSY79 income (unadj)", "household cluster",          1,
  "nlsy97_income_unadj",     "NLSY97 income (unadj)", "survey Fay-BRR (rho = 0.5)", 2,
  "nlsy79_income_unadj_ind", "NLSY79 income (unadj)", "individual (naive)",         1,
  "nlsy97_income_unadj_ind", "NLSY97 income (unadj)", "individual (naive)",         1
)

processed <- list()
for (i in seq_len(nrow(runs))) {
  s     <- runs$sample[i]
  cache <- file.path(dir_proc, paste0(s, "_processed.rds"))
  if (file.exists(cache)) { processed[[s]] <- readRDS(cache); next }
  message("processing ", s, " ...")
  run   <- load_run(s)
  stats <- process_domain(run$base, run$boot)
  FAY_K <<- runs$fay_k[i]
  val   <- generate_validation_data(stats$boot_stats, stats$baseline_stats)
  FAY_K <<- 1
  names(val$validation_lst) <- names(stats$boot_stats)
  processed[[s]] <- list(boot_stats = stats$boot_stats, baseline_stats = stats$baseline_stats,
                         validation = val, fay_k = runs$fay_k[i], n_reps = length(run$boot))
  saveRDS(processed[[s]], cache)
  rm(run, stats, val); invisible(gc())
}

# ---- scalar summary + design-vs-individual SE comparison (CSV outputs) ------
scalar_from_validation <- function(v) {
  bind_rows(
    v$MTE_df        %>% transmute(measure = class, estimate = mean_baseline, mean_bc, CI_bc_l, CI_bc_u, se = se_adj),
    v$lambda2_df    %>% transmute(measure = "lambda2", estimate = mean_baseline, mean_bc, CI_bc_l, CI_bc_u, se = se_adj),
    v$GeenensD      %>% transmute(measure = "GeenensD", estimate = mean_baseline, mean_bc, CI_bc_l, CI_bc_u, se = se_adj),
    v$HellingersDep %>% transmute(measure = "HellingerDep", estimate = mean_baseline, mean_bc, CI_bc_l, CI_bc_u, se = se_adj),
    v$altham        %>% transmute(measure = "Altham", estimate = mean_baseline, mean_bc, CI_bc_l, CI_bc_u, se = se_adj))
}
scalar_bc <- bind_rows(lapply(names(processed), function(s) {
  vl <- processed[[s]]$validation$validation_lst
  bind_rows(lapply(names(vl), function(co)
    scalar_from_validation(vl[[co]]) %>% mutate(cohort = co))) %>% mutate(sample = s)
})) %>%
  left_join(runs, by = "sample") %>%
  select(domain, design, sample, cohort, measure, estimate, mean_bc, CI_bc_l, CI_bc_u, se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write_csv(scalar_bc, file.path(dir_proc, "nlsy_bc_scalar_summary.csv"))

pair_map <- c(nlsy79_income = "nlsy79_income_ind", nlsy79_occ = "nlsy79_occ_ind",
              nlsy97_income = "nlsy97_income_ind",
              nlsy79_income_unadj = "nlsy79_income_unadj_ind",
              nlsy97_income_unadj = "nlsy97_income_unadj_ind")
se_compare <- bind_rows(lapply(names(pair_map), function(s) {
  des <- scalar_bc %>% filter(sample == s) %>%
    select(domain, cohort, measure, estimate, se_design = se,
           CI_bc_l_design = CI_bc_l, CI_bc_u_design = CI_bc_u)
  ind <- scalar_bc %>% filter(sample == pair_map[[s]]) %>%
    select(cohort, measure, se_ind = se, CI_bc_l_ind = CI_bc_l, CI_bc_u_ind = CI_bc_u)
  des %>% left_join(ind, by = c("cohort", "measure")) %>%
    mutate(se_ratio = round(se_design / se_ind, 3),
           ci_width_ratio = round((CI_bc_u_design - CI_bc_l_design) /
                                    (CI_bc_u_ind - CI_bc_l_ind), 3))
}))
write_csv(se_compare, file.path(dir_proc, "nlsy_bc_se_comparison.csv"))

# ---- validity record: bc estimates outside their own bc CI (0 = clean) ------
rec <- bind_rows(lapply(names(processed), function(s) {
  processed[[s]]$validation$record_df %>%
    mutate(cohort = names(processed[[s]]$validation$validation_lst), sample = s) }))
n_viol <- sum(rowSums(as.matrix(rec[sapply(rec, is.numeric)]), na.rm = TRUE))
cat("[06] validity violations across all runs x cohorts:", n_viol, "(0 = clean)\n")

# ---- the nine report-input computations (consolidated) ---------------------
message("[06] computing report inputs (NLSY_fun/nlsy_report_inputs.R) ...")
source(file.path(root, "NLSY_fun", "nlsy_report_inputs.R"))
cat("[06] bootstrap processing + report precomputes complete\n")
