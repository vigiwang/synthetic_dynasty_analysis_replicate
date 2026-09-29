# ============================================================================
# 00_generate_complete_occ_dyads.R
# ----------------------------------------------------------------------------
# Raw NLSY79 extracts -> respondent core + sample/cohort flags -> father/mother
# EGP (Census-1970 -> Dorn -> IPUMS -> Occ10EGP -> 5-class) -> child EGP
# life-trajectory + person measures (closest-to-age-40 primary, ages 25-55)
# -> COMPLETE OCCUPATION DYADS (pre-weight).
# Outputs: Data/NLSY_estimation/intermediate/nlsy79_core_aug.rds, nlsy79_child_egp_trajectory.rds
#          Data/NLSY_estimation/dyads/nlsy79_occ_person_preweight.rds
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

# ---- raw reads needed for this stage ----
n79_inc_raw  <- read_ext(p$n79_occ_raw)
n79_demo_raw <- read_ext(p$n79_occ_raw)

# ===========================================================================
# STEP 2  NLSY79 respondent core + sample/cohort flags (spec S2, S9.1)
# ===========================================================================
n79_core <- n79_inc_raw %>%
  transmute(
    caseid_raw = .data[["R0000100"]],
    birth_year = as.integer(.data[["R0000500"]]) + 1900L,   # 2-digit -> 4-digit
    sample_id  = as.integer(.data[["R0173600"]]),
    race       = as.integer(.data[["R0214700"]]),
    sex        = as.integer(.data[["R0214800"]])
  ) %>%
  standardize_id("caseid_raw", "caseid") %>%
  mutate(
    # NLSY79 age convention: age as of Dec 31, 1978 (the screener date), matching
    # the RA code and Zhou. age 14-18 as of end-1978 -> birth years 1960-1964
    # (1978-1964=14 ... 1978-1960=18). Using 1979-birthyear would wrongly drop 1960.
    age_end_1978 = 1978L - birth_year,
    # named sample-membership flags BEFORE any filtering (spec S2)
    is_cross_sectional_sample = sample_id %in% cfg$nlsy79_cross_sectional_codes,
    is_supplemental_sample    = sample_id %in% cfg$nlsy79_bh_supplement_codes,
    is_excluded_disadvantaged_supplement = sample_id %in% cfg$nlsy79_excluded_disadv_codes,
    is_excluded_military_sample           = sample_id %in% cfg$nlsy79_excluded_military,
    is_comparable_civilian_sample =
      is_cross_sectional_sample | is_supplemental_sample,
    is_nlsy79_age14_18_1979 = age_end_1978 >= cfg$nlsy79_baseline_age_min &
                              age_end_1978 <= cfg$nlsy79_baseline_age_max,
    is_nlsy79_birth1960_1964 = birth_year >= cfg$nlsy79_report_birth_min &
                               birth_year <= cfg$nlsy79_report_birth_max
  )
assert_unique_id(n79_core, "caseid", "nlsy79_core")
# Harmonized race/ethnicity (NLSY79 SAMPLE_RACE 1 Hisp/2 Black/3 NBNH) and clean
# sex, to the common Black / Hispanic / NBNH scheme used across both cohorts & GSS.
n79_core <- n79_core %>% mutate(
  race3 = dplyr::case_when(race == 1L ~ "Hispanic", race == 2L ~ "Black",
                           race == 3L ~ "NBNH", TRUE ~ NA_character_),
  sex_clean = ifelse(sex %in% c(1L, 2L), sex, NA_integer_))

# age-14-18 vs birth-1960-1964 discrepancy table (spec S9.1)
n79_cohort_discrepancy <- n79_core %>%
  filter(is_comparable_civilian_sample) %>%
  count(is_nlsy79_age14_18_1979, is_nlsy79_birth1960_1964, name = "n")
write_csv(n79_cohort_discrepancy, file.path(AUDIT, "nlsy79_age14_18_vs_birth1960_1964.csv"))

# ===========================================================================
# STEP 4  Audited augmentation merges (demographic -> income base) (spec S10)
# ===========================================================================
# NLSY79 demographic augmentation: respondent ID, HHID, parent education, resp edu
n79_demo_sel <- n79_demo_raw %>%
  transmute(
    caseid_raw     = .data[["R0000100"]],
    hhid           = as.character(.data[["R0000149"]]),
    hgc_mother_raw = .data[["R0006500"]],
    hgc_father_raw = .data[["R0007900"]],
    hgc_respondent_raw = .data[["T9900000"]]
  ) %>% standardize_id("caseid_raw", "caseid") %>% select(-caseid_raw)

n79_core_aug <- left_join_with_audit(
  base = n79_core, augment = n79_demo_sel, by = "caseid",
  label = "nlsy79_core_demo") %>%
  prepare_parent_education("hgc_mother_raw", "hgc_father_raw") %>%
  mutate(respondent_education = nls_clean_numeric(hgc_respondent_raw),
         respondent_education_missing_reason = nls_missing_reason(hgc_respondent_raw)) %>%
  # GSS-comparable 5-category harmonized education (all NLSY79 sources are years)
  add_harmonized_education("hgc_respondent_raw", "respondent", "years") %>%
  add_harmonized_education("hgc_mother_raw",     "mother",     "years") %>%
  add_harmonized_education("hgc_father_raw",     "father",     "years")

# ===========================================================================
# STEP 7  NLSY79 occupation: BUILD EGP OURSELVES via the crosswalk chain (spec S12)
# ---------------------------------------------------------------------------
# We construct EGP from raw occupation codes rather than consuming the RA dyad,
# so we can (a) use the COMPARABLE-CIVILIAN population like the income module
# (not cross-sectional only), (b) add MOTHER EGP (FAMOCC-19), and (c) create two
# measures the RA dyad lacks: child EGP CLOSEST TO AGE 40 and LONGEST-DURATION
# EGP (exposure-weighted). Chain: raw Census -> OCC1990dd (Dorn) -> OCC2010
# (IPUMS) -> EGP-11 (Occ10EGP) -> 5-class (map5). Validated: our father EGP
# reproduces the RA dyad 100%; child mode/best agree ~95-97% (window 25-55 vs
# the RA's 26-55). Adult career window = cfg$occ_age_min..occ_age_max.
max_duration_egp_available <- TRUE   # we now build a genuine longest-duration measure
cw <- load_occ_crosswalks(root)

occ_lab <- read_ext(p$n79_occ_raw)
occ_lab$caseid     <- as.character(occ_lab[["CASEID_1979"]])
occ_lab$birth_year <- as.integer(occ_lab[["Q1-3_A~Y_1979"]]) + 1900L
occ_lab$sample_id  <- as.integer(occ_lab[["SAMPLE_ID_1979"]])
occ_lab$race       <- as.integer(occ_lab[["SAMPLE_RACE_78SCRN"]])
occ_lab$sex        <- as.integer(occ_lab[["SAMPLE_SEX_1979"]])
# mother occupation (FAMOCC-19 = R0006900), Census 1970, from the Demographic extract
occ_lab <- left_join(occ_lab,
  transmute(n79_demo_raw, caseid = as.character(.data[["R0000100"]]),
            mother_occ_1970 = .data[["R0006900"]]), by = "caseid")

occ_cc <- occ_lab %>%
  mutate(is_cross_sectional_sample = sample_id %in% cfg$nlsy79_cross_sectional_codes,
         is_supplemental_sample    = sample_id %in% cfg$nlsy79_bh_supplement_codes,
         is_comparable_civilian_sample = is_cross_sectional_sample | is_supplemental_sample) %>%
  filter(is_comparable_civilian_sample) %>%
  mutate(father_egp = census_to_egp5(.data[["FAMOCC-26_1979"]], "1970", cw),
         mother_egp = census_to_egp5(mother_occ_1970,            "1970", cw))

# child EGP life-trajectory (one row per respondent-survey-year) + person measures
occ_traj <- build_child_egp_trajectory(
  occ_cc[, c("caseid", "birth_year", grep("^CPSOCC80|^OCCALL", names(occ_cc), value = TRUE))],
  "caseid", cw)
occ_child <- summarise_child_egp(occ_traj, "caseid",
               age_min = cfg$occ_age_min, age_max = cfg$occ_age_max,
               target_age = cfg$occ_target_age)
saveRDS(occ_traj, file.path(INTER, "nlsy79_child_egp_trajectory.rds"))

occ_sample <- occ_cc %>%
  select(caseid, birth_year, sample_id, race, sex,
         is_cross_sectional_sample, is_supplemental_sample, is_comparable_civilian_sample,
         father_egp, mother_egp) %>%
  left_join(occ_child, by = "caseid") %>%
  # attach hhid (sibling-cluster bootstrap) + harmonized race/sex + demographic education
  left_join(n79_core_aug %>%
              select(caseid, hhid, race3, sex_clean,
                     dplyr::matches("^(respondent|mother|father)_educ")),
            by = "caseid") %>%
  mutate(
    parent_egp = father_egp,                      # father-based dyad (RA convention)
    is_birth1960_1964 = birth_year >= cfg$nlsy79_report_birth_min &
                        birth_year <= cfg$nlsy79_report_birth_max,
    child_egp_primary = child_egp_age40,          # PRIMARY child measure = closest to age 40
    complete_occ_dyad = !is.na(parent_egp) & !is.na(child_egp_primary),
    complete_occ_dyad_mother = !is.na(mother_egp) & !is.na(child_egp_primary)
  )

# ---- stage outputs ----
saveRDS(n79_core_aug, file.path(INTER, "nlsy79_core_aug.rds"))
saveRDS(occ_sample,   file.path(ANALY, "nlsy79_occ_person_preweight.rds"))
cat("\n[00] complete occupation dyads:", sum(occ_sample$complete_occ_dyad),
    "of", nrow(occ_sample), "comparable-civilian respondents\n")
