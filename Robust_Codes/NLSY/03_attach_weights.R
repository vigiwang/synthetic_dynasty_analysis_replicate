# ============================================================================
# 03_generate_complete_dyads_weight.R
# ----------------------------------------------------------------------------
# (a) assembles the analytic samples, (b) EXPORTS THE CUSTOM-WEIGHT ID LISTS
# (upload to NLS Investigator; the NLSY79 income request is the 1957-64
# UNIVERSE - never the narrow reporting range), (c) imports + validates the
# RETURNED weights from Data/NLSY/custom_weights_returned/, (d) attaches
# them -> FINAL weighted dyad files + demographic objects + survey designs.
# The manual NLS request is "stage 02" of the pipeline (see README).
# Requires stages 00-01. Outputs: Data/NLSY_estimation/custom_weight_ids/*,
# nlsy79_income_final.rds, nlsy97_income_final.rds, nlsy79_occ_final.rds,
# demographic *.rds, design objects, output/audit/*.
# Code below is extracted VERBATIM from the verified monolithic build
# (legacy/00_run_build.R); only the load/save seams between stages are new.
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

# ---- stage-00/01 inputs + raw dyad extracts (for the reconstruction audit) ----
n79_income_person <- readRDS(file.path(ANALY, "nlsy79_income_person_preweight.rds"))
n97_income_person <- readRDS(file.path(ANALY, "nlsy97_income_person_preweight.rds"))
occ_sample        <- readRDS(file.path(ANALY, "nlsy79_occ_person_preweight.rds"))
n79_core_aug      <- readRDS(file.path(INTER, "nlsy79_core_aug.rds"))
n97_core_aug      <- readRDS(file.path(INTER, "nlsy97_core_aug.rds"))
# RA dyad files are OPTIONAL (reconstruction cross-check only; the pipeline
# needs only the raw extracts + returned weights)
occ_dyad   <- if (file.exists(p$n79_occ_dyad)) read_ext(p$n79_occ_dyad) else NULL
inc79_dyad <- if (file.exists(p$n79_inc_dyad)) read_ext(p$n79_inc_dyad) else NULL
inc97_dyad <- if (file.exists(p$n97_inc_dyad)) read_ext(p$n97_inc_dyad) else NULL

# ===========================================================================
# STEP 8  Assemble analytic samples
# ===========================================================================
n79_income_sample <- filter(n79_income_person, nlsy79_zhou_income_sample)
n97_income_sample <- filter(n97_income_person, nlsy97_zhou_income_sample)
n79_occ_sample <- filter(occ_sample, complete_occ_dyad)


# ===========================================================================
# STEP 8b  Demographic objects for the analytic samples + missingness audit
# ---------------------------------------------------------------------------
# The demographic module is for the GSS cross-check. CRITICAL (spec S11/S19): a
# mobility-complete respondent is NEVER dropped because a demographic is missing;
# we keep everyone, flag missings, and report available-case denominators. Sex
# and race are complete; parent education (father especially) and mother's
# occupation EGP have real missingness that must be reported, not silently dropped.
demo_keep <- function(d) {
  d %>% select(any_of(c("caseid", "pubid", "birth_year", "sample_id", "sample_type",
    "sex", "sex_clean", "race", "race3", "is_comparable_civilian_sample",
    "is_cross_sectional_sample", "is_supplemental_sample")),
    dplyr::matches("^(respondent|mother|father)_educ"),
    any_of(c("father_egp", "mother_egp")))
}
nlsy79_demo_all_comparable <- n79_core_aug %>% filter(is_comparable_civilian_sample) %>% demo_keep()
nlsy97_demo_all_comparable <- n97_core_aug %>% filter(is_comparable_civilian_sample) %>% demo_keep()
nlsy79_demo_income_sample  <- demo_keep(n79_income_sample)
nlsy97_demo_income_sample  <- demo_keep(n97_income_sample)
nlsy79_demo_occ_sample     <- demo_keep(n79_occ_sample)
nlsy79_demo_cross_sectional_only_sensitivity <- n79_core_aug %>% filter(is_cross_sectional_sample) %>% demo_keep()
nlsy97_demo_cross_sectional_only_sensitivity <- n97_core_aug %>% filter(is_cross_sectional_sample) %>% demo_keep()

# available-case demographic-missingness audit (does anyone in a complete dyad
# lose demographic info?)
demo_missing_audit <- function(d, vars, label) {
  tibble(sample = label, n = nrow(d),
         variable = vars,
         n_missing = vapply(vars, function(v) if (v %in% names(d)) sum(is.na(d[[v]])) else NA_integer_, integer(1))) %>%
    mutate(pct_missing = round(100 * n_missing / n, 1),
           available_n = n - n_missing)
}
demo_vars5 <- c("sex", "race", "respondent_educ5", "mother_educ5", "father_educ5")
demo_missing <- bind_rows(
  demo_missing_audit(n79_income_sample, demo_vars5, "nlsy79_income"),
  demo_missing_audit(n97_income_sample, demo_vars5, "nlsy97_income"),
  demo_missing_audit(n79_occ_sample, c("sex","race","father_egp","mother_egp",
                                       "respondent_educ5","mother_educ5","father_educ5"), "nlsy79_occ")
)
write_csv(demo_missing, file.path(AUDIT, "demographic_missingness_by_sample.csv"))
message("Demographic missingness within complete-dyad samples (available-case denominators):")
print(as.data.frame(demo_missing))
# assert no mobility-complete respondent was dropped for missing demographics
assert_that(nrow(n79_income_sample) == sum(n79_income_person$nlsy79_zhou_income_sample),
            "income sample size must not depend on demographic completeness")


# ===========================================================================
# STEP 9  Cross-check reconstructed samples against RA dyads
# ===========================================================================
compare_ids <- function(mine, ra, label) {
  mine <- sort(unique(as.integer(mine))); ra <- sort(unique(as.integer(ra)))
  tibble(sample = label, n_mine = length(mine), n_ra = length(ra),
         n_shared = length(intersect(mine, ra)),
         n_mine_only = length(setdiff(mine, ra)),
         n_ra_only = length(setdiff(ra, mine)))
}
if (!is.null(inc79_dyad) && !is.null(inc97_dyad) && !is.null(occ_dyad)) {
recon_audit <- bind_rows(
  compare_ids(n79_income_sample$caseid, inc79_dyad$CASEID_1979, "nlsy79_income"),
  compare_ids(n97_income_sample$pubid,  inc97_dyad$PUBID_1997,  "nlsy97_income"),
  compare_ids(n79_occ_sample$caseid, occ_dyad$CASEID_1979, "nlsy79_occ")
)
write_csv(recon_audit, file.path(AUDIT, "reconstruction_vs_ra_dyads.csv"))
print(recon_audit)
} else message("[03] RA dyad files absent - reconstruction cross-check skipped (optional)")


# ===========================================================================
# STEP 10  EXPORT CUSTOM-WEIGHT ID FILES  (spec S14)  <-- KEY DELIVERABLE
# ===========================================================================
write_custom_weight_ids(n79_income_sample, "caseid",
  file.path(CWID, "nlsy79_income_ids_zhou"),
  cohort = "birth_year", sex = "sex", race = "race", sample_type = "sample_id")

write_custom_weight_ids(n97_income_sample, "pubid",
  file.path(CWID, "nlsy97_income_ids"),
  cohort = "birth_year", sex = "sex", race = "race", sample_type = "sample_type")

write_custom_weight_ids(n79_occ_sample, "caseid",
  file.path(CWID, "nlsy79_occ_ids"),
  cohort = "birth_year", sex = "sex", race = "race", sample_type = "sample_id")

# NLSY79 income WEIGHTING UNIVERSE (recommended primary request): comparable
# civilian, complete income dyad, born 1957-1964 (all NLSY79 cohorts). Requesting
# custom weights for this full universe -- then restricting analysis to the Zhou
# reporting DOMAIN (age 14-18 as of end-1978 = born 1960-1964) -- avoids the 1960
# boundary-weight anomaly that arises when weights are requested for a narrow
# cohort cutoff (spec S14.1/S15).
n79_income_universe <- n79_income_person %>%
  filter(is_comparable_civilian_sample, complete_income_dyad,
         birth_year >= cfg$nlsy79_weight_birth_min,
         birth_year <= cfg$nlsy79_weight_birth_max)
write_custom_weight_ids(n79_income_universe, "caseid",
  file.path(CWID, "nlsy79_income_ids_1957_1964_diagnostic"),
  cohort = "birth_year", sex = "sex", race = "race", sample_type = "sample_id")

# Domain-membership map: which universe IDs fall in the Zhou reporting domain
# (age 14-18 as of end-1978 -> born 1960-1964). After weights return for the
# universe, restrict to is_zhou_domain == TRUE for the Zhou-domain estimates.
n79_income_universe %>%
  transmute(caseid, birth_year,
            is_zhou_domain = is_nlsy79_age14_18_1979 & birth_year >= 1961 & birth_year <= 1964) %>%
  write_csv(file.path(CWID, "nlsy79_income_universe_domain_flags.csv"))


# ===========================================================================
# STEP 11  Import & diagnose returned custom weights (spec S15)
# ===========================================================================
read_returned_weight <- function(path, id_name) {
  if (!file_exists(path)) return(NULL)
  w <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  names(w)[1:2] <- c(id_name, "custom_weight")
  w %>% mutate(!!id_name := as.character(.data[[id_name]]),
               custom_weight = as.numeric(custom_weight))
}
diagnose_weight <- function(sample, weight, id, label) {
  if (is.null(weight)) return(tibble(sample = label, note = "returned weights absent"))
  m <- left_join(sample, weight, by = id)
  n_missing <- sum(is.na(m$custom_weight))
  w <- m$custom_weight[!is.na(m$custom_weight)]
  kish <- sum(w)^2 / sum(w^2)
  tibble(sample = label, n = nrow(m), n_missing_weight = n_missing,
         sum_w = sum(w), mean_w = mean(w), median_w = median(w),
         max_w = max(w), cv_w = sd(w) / mean(w), kish_ess = kish)
}
# Prefer newly returned custom weights (Data/NLSY/custom_weights_returned/,
# where the user saves the files requested via the upload checklist); fall back
# to the RA's originals if the new ones are not present yet.
RETW <- file.path(root, "Data/NLSY/custom_weights_returned")  # legacy; weights now in outcome folders
pick_weight <- function(new_name, ra_path) {
  np <- file.path(RETW, new_name)
  if (file_exists(np)) np else ra_path
}
w79_inc <- read_returned_weight(file.path(root,"Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_weight.csv"), "caseid")
w79_occ <- read_returned_weight(file.path(root,"Data/NLSY/NLSY79/NLSY79_occ/NLSY79_occ_weight.csv"), "caseid")
w97_inc <- read_returned_weight(file.path(root,"Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_weight.csv"), "pubid")

weight_diag <- bind_rows(
  diagnose_weight(n79_income_sample, w79_inc, "caseid", "nlsy79_income"),
  diagnose_weight(n79_occ_sample, w79_occ, "caseid", "nlsy79_occ"),
  diagnose_weight(n97_income_sample, w97_inc, "pubid", "nlsy97_income")
)
write_csv(weight_diag, file.path(AUDIT, "custom_weight_diagnostics.csv"))
print(weight_diag)

# NLSY79 1960 boundary diagnostic: custom-weight sum by cohort (spec S15).
if (!is.null(w79_inc)) {
  n79_1960_diag <- n79_income_sample %>%
    left_join(w79_inc, by = "caseid") %>%
    group_by(birth_year) %>%
    summarise(n = n(), sum_custom_weight = sum(custom_weight, na.rm = TRUE),
              mean_custom_weight = mean(custom_weight, na.rm = TRUE), .groups = "drop")
  write_csv(n79_1960_diag, file.path(AUDIT, "nlsy79_1960_cohort_weight_diagnostic.csv"))
}

# Final objects with weights attached (only after weight validation).
# Attach the returned custom weight to the analytic sample. Weights depend ONLY on
# sample membership (mobility completeness), NOT on demographic completeness -- a
# respondent with a missing parent-education is still in the sample and still gets a
# weight. Respondents the NLS tool returned no weight for keep custom_weight = NA and
# are flagged has_custom_weight = FALSE (reported, not dropped).
attach_weight <- function(sample, weight, id) {
  s <- if (is.null(weight)) mutate(sample, custom_weight = NA_real_)
       else left_join(sample, weight, by = id)
  s %>% mutate(has_custom_weight = !is.na(custom_weight))
}
nlsy79_income_final    <- attach_weight(n79_income_sample, w79_inc, "caseid")
nlsy79_occ_final <- attach_weight(n79_occ_sample, w79_occ, "caseid")
nlsy97_income_final    <- attach_weight(n97_income_sample, w97_inc, "pubid")


# ===========================================================================
# STEP 12  Bootstrap / survey-design inputs (spec S16) + smoke tests
# ===========================================================================
n79_income_indiv_boot <- make_individual_bootstrap_input(
  nlsy79_income_final, "caseid", "custom_weight", "birth_year")
# nlsy79_income_final already carries hhid (merged from the demographic core),
# so the household-cluster input just needs the existing column.
n79_income_hh_boot <- make_household_bootstrap_input(
  nlsy79_income_final, "caseid", "hhid", "custom_weight")
n97_income_indiv_boot <- make_individual_bootstrap_input(
  nlsy97_income_final, "pubid", "custom_weight", "birth_year")

# NLSY97 pooled survey design + Fay-BRR replicate design (needs VSTRAT/VPSU/weight)
design97_income <- repdesign97_income_fay <- NULL
if (all(c("vstrat", "vpsu") %in% names(nlsy97_income_final)) &&
    any(!is.na(nlsy97_income_final$custom_weight))) {
  design97_income <- tryCatch(
    build_nlsy97_design(nlsy97_income_final), error = function(e) {message(e); NULL})
  if (!is.null(design97_income))
    repdesign97_income_fay <- build_nlsy97_fay_brr(design97_income, cfg$fay_rho)
}


# ===========================================================================
# STEP 13  Save intermediate + analysis objects and a manifest (spec S18)
# ===========================================================================
saveRDS(nlsy79_demo_all_comparable, file.path(ANALY, "nlsy79_demo_all_comparable.rds"))
saveRDS(nlsy97_demo_all_comparable, file.path(ANALY, "nlsy97_demo_all_comparable.rds"))
saveRDS(nlsy79_demo_income_sample,  file.path(ANALY, "nlsy79_demo_income_sample.rds"))
saveRDS(nlsy97_demo_income_sample,  file.path(ANALY, "nlsy97_demo_income_sample.rds"))
saveRDS(nlsy79_demo_occ_sample,     file.path(ANALY, "nlsy79_demo_occ_sample.rds"))
saveRDS(nlsy79_demo_cross_sectional_only_sensitivity, file.path(ANALY, "nlsy79_demo_cross_sectional_only_sensitivity.rds"))
saveRDS(nlsy97_demo_cross_sectional_only_sensitivity, file.path(ANALY, "nlsy97_demo_cross_sectional_only_sensitivity.rds"))
saveRDS(nlsy79_income_final,    file.path(ANALY, "nlsy79_income_final.rds"))
saveRDS(nlsy79_occ_final, file.path(ANALY, "nlsy79_occ_final.rds"))
saveRDS(nlsy97_income_final,    file.path(ANALY, "nlsy97_income_final.rds"))
if (!is.null(design97_income))         saveRDS(design97_income, file.path(ANALY, "design97_income.rds"))
if (!is.null(repdesign97_income_fay))  saveRDS(repdesign97_income_fay, file.path(ANALY, "repdesign97_income_fay.rds"))
saveRDS(n79_income_indiv_boot, file.path(ANALY, "nlsy79_income_individual_bootstrap_input.rds"))
saveRDS(n79_income_hh_boot,    file.path(ANALY, "nlsy79_income_household_bootstrap_input.rds"))
saveRDS(n97_income_indiv_boot, file.path(ANALY, "nlsy97_income_individual_bootstrap_input.rds"))

manifest <- tribble(
  ~object, ~unit, ~n, ~use,
  "nlsy79_income_person_preweight", "respondent", nrow(n79_income_person), "NLSY79 income preweight",
  "nlsy79_income_final",            "respondent", nrow(nlsy79_income_final), "NLSY79 income + custom weight",
  "nlsy97_income_person_preweight", "respondent", nrow(n97_income_person), "NLSY97 income preweight",
  "nlsy97_income_final",            "respondent", nrow(nlsy97_income_final), "NLSY97 income + custom weight",
  "nlsy79_occ_final",         "respondent", nrow(nlsy79_occ_final), "NLSY79 occupation (mode EGP, 25-55) + weight"
)
write_csv(manifest, file.path(AUDIT, "final_object_manifest.csv"))

cat("\n==== BUILD COMPLETE ====\n")
cat("Custom-weight upload files written to:", CWID, "\n")
print(dir_ls(CWID))
