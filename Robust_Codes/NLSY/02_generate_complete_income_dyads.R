# ============================================================================
# 01_generate_complete_income_dyads.R
# ----------------------------------------------------------------------------
# Raw NLSY79/NLSY97 income extracts -> long income observations -> generation
# tagging (parent window w/ lived-at-home screen; adult ages 27-38 by income
# reference year) -> per-year PCE deflation + sqrt(household size)
# equivalization -> within-person averages -> COMPLETE INCOME DYADS (pre-weight).
# NLSY97 parent income: 1997 parent-reported primary, mean of youth-reported
# 1998-2001 fallback. Requires stage 00 (nlsy79_core_aug.rds).
# Outputs: Data/NLSY_estimation/intermediate/nlsy79_income_long.rds, nlsy97_income_long.rds,
#          nlsy97_core_aug.rds; Data/NLSY_estimation/dyads/nlsy79_income_person_preweight.rds,
#          nlsy97_income_person_preweight.rds
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

# ---- raw reads + stage-00 core ----
n79_inc_raw  <- read_ext(p$n79_inc)
n79_demo_raw <- read_ext(p$n79_demo)
n97_inc_raw  <- read_ext(p$n97_inc)
n97_demo_raw <- read_ext(p$n97_demo)
n79_core_aug <- readRDS(file.path(INTER, "nlsy79_core_aug.rds"))
n79_core     <- n79_core_aug            # superset of the core columns

# ===========================================================================
# STEP 3  NLSY97 respondent core + sample flags (spec S2, S9.2)
# ===========================================================================
n97_core <- n97_inc_raw %>%
  transmute(
    pubid_raw   = .data[["R0000100"]],
    sex         = as.integer(.data[["R0536300"]]),
    birth_year  = as.integer(.data[["R0536402"]]),
    sample_type = as.integer(.data[["R1235800"]]),
    race        = as.integer(.data[["R1482600"]])
  ) %>%
  standardize_id("pubid_raw", "pubid") %>%
  mutate(
    is_cross_sectional_sample = sample_type == cfg$nlsy97_cross_sectional_code,
    is_supplemental_sample    = sample_type == cfg$nlsy97_supplemental_code,
    is_comparable_civilian_sample = is_cross_sectional_sample | is_supplemental_sample
  )
assert_unique_id(n97_core, "pubid", "nlsy97_core")
# Harmonized race (NLSY97 KEY!RACE 1 Black/2 Hisp/3 Mixed/4 NBNH -> common scheme;
# "Mixed" folds into NBNH) and clean sex (NLSY97 0 = No Information -> NA).
n97_core <- n97_core %>% mutate(
  race3 = dplyr::case_when(race == 1L ~ "Black", race == 2L ~ "Hispanic",
                           race %in% c(3L, 4L) ~ "NBNH", TRUE ~ NA_character_),
  sex_clean = ifelse(sex %in% c(1L, 2L), sex, NA_integer_))

# NLSY97 demographic augmentation: SIDCODE, VSTRAT, VPSU, parent edu, HH size
n97_demo_sel <- n97_demo_raw %>%
  transmute(
    pubid_raw      = .data[["R0000100"]],
    sidcode        = as.character(.data[["R1193000"]]),
    vstrat         = as.integer(.data[["R1489700"]]),
    vpsu           = as.integer(.data[["R1489800"]]),
    hh_size_1997   = nls_clean_numeric(.data[["R1205400"]]),
    hgc_mother_raw = .data[["R1302500"]],
    hgc_father_raw = .data[["R1302400"]],
    degree_respondent_raw = .data[["Z9083900"]],   # CVC_HIGHEST_DEGREE_EVER (respondent)
    weight_cc_1997 = nls_clean_numeric(.data[["R1236101"]])
  ) %>% standardize_id("pubid_raw", "pubid") %>% select(-pubid_raw)

# NLSY97 respondent YEARS of schooling (CVC_HGC_EVER), consolidated from the
# NLSY97_edu_year extract into the demographic folder. This makes the respondent
# years-based like NLSY79, so educ5 is symmetric across cohorts and matches GSS educ.
n97_resp_hgc <- read_ext(p$n97_inc) %>%
  transmute(pubid = as.character(R0000100), hgc_respondent_raw = CVC_HGC_EVER)
assert_unique_id(n97_resp_hgc, "pubid", "nlsy97_resp_hgc")

n97_core_aug <- left_join_with_audit(
  base = n97_core, augment = n97_demo_sel, by = "pubid",
  label = "nlsy97_core_demo") %>%
  left_join(n97_resp_hgc, by = "pubid") %>%
  prepare_parent_education("hgc_mother_raw", "hgc_father_raw") %>%
  # PRIMARY respondent education is now YEARS (CVC_HGC_EVER); parents are years too
  add_harmonized_education("hgc_respondent_raw", "respondent", "years") %>%
  add_harmonized_education("hgc_mother_raw",     "mother",     "years") %>%
  add_harmonized_education("hgc_father_raw",     "father",     "years") %>%
  # keep the degree-based classification as a secondary/robustness measure
  mutate(respondent_educ5_degree = nlsy97_degree_to_educ5(degree_respondent_raw),
         respondent_educ4_degree = nlsy97_degree_to_educ4(degree_respondent_raw))

# ===========================================================================
# STEP 5  NLSY79 income: long -> generation tagging -> person summary (S13.1/3/4/6)
# ===========================================================================
by79 <- n79_core %>% select(caseid, birth_year)

# family size per survey year (need-adjustment denominator), from demographic extract
n79_size_long <- n79_demo_raw %>% standardize_id("R0000100", "caseid") %>%
  reshape_size_long("caseid", nlsy79_famsize_map)

# BDZ footnote-20 lived-at-home screen (parent-report "Version A"): TNFI counts
# as PARENTAL income only in survey rounds where the household record was
# parent-reported, i.e. HHI-24 == 1 (respondent living in the parental
# household). HHI-24 was NOT fielded in the 1981 survey, so income year 1980
# cannot be verified and is excluded from the parental average (same rule as
# the prior construction). Zhou/BDZ 2018: "We average all years in which
# parents reported their income" (fn20: "parents provided income reports on a
# special survey version").
n79_home_long <- read_ext(p$n79_inc_lab) %>%
  transmute(caseid = as.character(.data[["CASEID_1979"]]),
            `1979` = .data[["HHI-24_1979"]], `1980` = .data[["HHI-24_1980"]],
            `1982` = .data[["HHI-24_1982"]], `1983` = .data[["HHI-24_1983"]]) %>%
  pivot_longer(-caseid, names_to = "survey_year", values_to = "hhi24") %>%
  transmute(caseid, survey_year = as.integer(survey_year),
            lived_at_home = !is.na(hhi24) & hhi24 == 1)

# long income -> generation windows -> lived-at-home screen -> merge family size
# -> PCE-deflate & equivalize.
# Real income = nominal * pce_base / pce_incomeyear (BEA PCE, 2017=100, base 2023).
# Two income measures are stored: income_real (size-UNadjusted, robustness) and
# income_equiv = income_real / sqrt(family_size) (primary, need-adjusted).
n79_income_long <- n79_inc_raw %>%
  standardize_id("R0000100", "caseid") %>%
  reshape_income_long("caseid", nlsy79_tnfi_map) %>%
  tag_income_generation(by79, "caseid",
                        parent_survey_years = cfg$parent_rounds_79,
                        adult_age_min = cfg$adult_income_age_min,
                        adult_age_max = cfg$adult_income_age_max) %>%
  left_join(n79_home_long %>% mutate(caseid = caseid),
            by = c("caseid", "survey_year")) %>%
  mutate(parent_source_ok = coalesce(lived_at_home, FALSE),
         is_parent_obs    = is_valid_income & is_parent_window & parent_source_ok) %>%
  left_join(n79_size_long, by = c("caseid", "survey_year")) %>%
  add_real_equiv("income_nominal", "income_year", "family_size", pce)

message(sprintf(
  "NLSY79 lived-at-home screen: parent obs %s (screened) across income years {%s}",
  sum(n79_income_long$is_parent_obs, na.rm = TRUE),
  paste(sort(unique(n79_income_long$income_year[n79_income_long$is_parent_obs])), collapse = ",")))

n79_income_person <- summarise_income_person_real(n79_income_long, "caseid",
                       min_parent = cfg$min_parent_income_obs,
                       min_adult  = cfg$min_adult_income_obs) %>%
  left_join(n79_core_aug, by = "caseid")

# reporting-sample flag = comparable civilian & complete dyad, restricted to the
# reporting birth domain (born >= nlsy79_report_birth_min). The age-14-18-in-1978
# screen marks the comparable NLSY79 sample; the reporting floor (1961, = BDZ's
# age-17-in-1978 ceiling) drops the 1960 cohort from the analysis domain while it
# stays in the universe as weighting scaffolding (like 1957-59).
n79_income_person <- n79_income_person %>%
  mutate(nlsy79_zhou_income_sample =
           is_comparable_civilian_sample & is_nlsy79_age14_18_1979 & complete_income_dyad &
           birth_year >= cfg$nlsy79_report_birth_min & birth_year <= cfg$nlsy79_report_birth_max)

# ===========================================================================
# STEP 6  NLSY97 income  (SPECIAL parent-income handling)
# ---------------------------------------------------------------------------
# The NLSY97 parent questionnaire was fielded ONLY in Round 1 (1997), so 1997 is
# the sole genuinely parent-reported household income; 1998-2001 are
# youth-reported and noisy (8% don't-know vs 0% for parents; ~0.44 within-person
# correlation with the 1997 figure; mostly valid-skipped). See
# analysis/investigate_nlsy97_parent_income.R.
#
# We therefore do NOT average 1996-2000 (the RA's original approach). Instead:
#   * PARENT income = 1997 parent-reported value (primary), coalesced to
#     1998->1999->2000->2001 ONLY when 1997 is missing (build_nlsy97_parent_income);
#     the parent generation needs >=1 such value.
#   * ADULT income  = CV_INCOME_FAMILY, age 27-32, >=2 obs (unchanged).
# ===========================================================================
by97 <- n97_core %>% select(pubid, birth_year)
n97_inc_std <- n97_inc_raw %>% standardize_id("R0000100", "pubid")

# household size per survey year (need-adjustment denominator); NLSY97 uses
# CV_HH_SIZE as the documented household-size proxy for family size.
n97_size_long <- n97_demo_raw %>% standardize_id("R0000100", "pubid") %>%
  reshape_size_long("pubid", nlsy97_hhsize_map)
n97_size_by_year <- rename(n97_size_long, household_size = family_size)

# --- adult side: long -> age-27-32 window -> deflate & equivalize -> counts/means ---
n97_income_long <- n97_inc_std %>%
  reshape_income_long("pubid", nlsy97_income_map) %>%
  tag_income_generation(by97, "pubid",
                        parent_survey_years = cfg$parent_rounds_97,
                        adult_age_min = cfg$adult_income_age_min,
                        adult_age_max = cfg$adult_income_age_max) %>%
  left_join(rename(n97_size_by_year, family_size = household_size),
            by = c("pubid", "survey_year")) %>%
  add_real_equiv("income_nominal", "income_year", "family_size", pce)
n97_adult <- n97_income_long %>%
  group_by(pubid) %>%
  summarise(n_adult_income_obs      = sum(is_adult_obs, na.rm = TRUE),
            mean_adult_income_real  = mean(income_real[is_adult_obs],  na.rm = TRUE),
            mean_adult_income_equiv = mean(income_equiv[is_adult_obs], na.rm = TRUE),
            adult_income_years      = paste(sort(income_year[is_adult_obs]), collapse = ";"),
            .groups = "drop") %>%
  mutate(across(starts_with("mean_"), ~ ifelse(is.nan(.x), NA_real_, .x)))

# --- parent side, measure A (default): 1997 parent-reported PRIMARY; when 1997
#     is missing, impute with the AVERAGE of the available 1998-2001 youth-reported
#     years (each deflated to real$ and size-equivalized AT ITS OWN survey year
#     BEFORE averaging), rather than the first available single year. Provenance
#     flags (is_parent_reported / is_imputed / has_parent_income) come from
#     build_nlsy97_parent_income(); the equivalized/real VALUES come from
#     build_nlsy97_parent_income_avg(). For the 1997 parent-reported dyads the two
#     agree exactly; only the imputed dyads with >=2 fallback years change.
n97_parent <- build_nlsy97_parent_income(n97_inc_std, id_col = "pubid") %>%
  select(-any_of(c("parent_income_real", "parent_income_equiv"))) %>%
  left_join(build_nlsy97_parent_income_avg(n97_inc_std, n97_size_by_year, pce,
                                           id_col = "pubid"),
            by = "pubid")

# --- parent side, measure B (Zhou strict): average ALL valid 1996-2000 years ---
# Stored alongside measure A for comparison/robustness (Zhou BDZ 2018 p.1226).
n97_parent_zhou <- n97_income_long %>%
  group_by(pubid) %>%
  summarise(
    n_parent_obs_zhou           = sum(is_parent_obs, na.rm = TRUE),
    parent_income_real_zhou     = mean(income_real[is_parent_obs],  na.rm = TRUE),
    parent_income_equiv_zhou    = mean(income_equiv[is_parent_obs], na.rm = TRUE),
    .groups = "drop") %>%
  mutate(across(starts_with("parent_income_"), ~ ifelse(is.nan(.x), NA_real_, .x)))

# --- cross-validation of the 1998-2001 imputation sources vs 1997 ---
n97_parent_xval <- crossvalidate_nlsy97_parent_income(n97_inc_std)
write_csv(n97_parent_xval, file.path(AUDIT, "nlsy97_parent_income_crossvalidation.csv"))
message("NLSY97 parent-income cross-validation (fallback vs 1997 parent-reported):")
print(n97_parent_xval)

# provenance audit: how many parent values are primary (1997) vs imputed
n97_parent_prov <- n97_parent %>%
  summarise(n_any_parent = sum(has_parent_income),
            n_parent_reported_1997 = sum(parent_income_is_parent_reported),
            n_imputed_1998_2001    = sum(parent_income_is_imputed))
write_csv(n97_parent_prov, file.path(AUDIT, "nlsy97_parent_income_provenance.csv"))

# --- assemble NLSY97 person-level income record ---
# complete_income_dyad depends on cfg$nlsy97_parent_income_rule:
#   "coalesce_1997_primary" (default): >=1 coalesced parent value (measure A)
#   "zhou_avg_1996_2000_ge2"          : >=2 parent-reported years (measure B, Zhou)
n97_income_person <- n97_core_aug %>%
  left_join(n97_adult,       by = "pubid") %>%
  left_join(n97_parent,      by = "pubid") %>%
  left_join(n97_parent_zhou, by = "pubid") %>%
  mutate(
    n_adult_income_obs = coalesce(n_adult_income_obs, 0L),
    n_parent_obs_zhou  = coalesce(n_parent_obs_zhou, 0L),
    parent_dyad_ok = if (cfg$nlsy97_parent_income_rule == "zhou_avg_1996_2000_ge2")
      n_parent_obs_zhou >= cfg$min_parent_income_obs else has_parent_income,
    complete_income_dyad = coalesce(parent_dyad_ok & n_adult_income_obs >= cfg$min_adult_income_obs, FALSE),
    nlsy97_zhou_income_sample = is_comparable_civilian_sample & complete_income_dyad
  )


# ---- stage outputs ----
saveRDS(n79_income_long,   file.path(INTER, "nlsy79_income_long.rds"))
saveRDS(n97_income_long,   file.path(INTER, "nlsy97_income_long.rds"))
saveRDS(n97_core_aug,      file.path(INTER, "nlsy97_core_aug.rds"))
saveRDS(n79_income_person, file.path(ANALY, "nlsy79_income_person_preweight.rds"))
saveRDS(n97_income_person, file.path(ANALY, "nlsy97_income_person_preweight.rds"))
cat("\n[01] complete income dyads: NLSY79",
    sum(n79_income_person$complete_income_dyad, na.rm=TRUE), "| NLSY97",
    sum(n97_income_person$complete_income_dyad, na.rm=TRUE), "\n")
