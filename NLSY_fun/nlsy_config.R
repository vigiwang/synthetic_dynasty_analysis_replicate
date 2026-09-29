# ============================================================================
# nlsy_config.R — shared configuration for the numbered NLSY pipeline
# ----------------------------------------------------------------------------
# Sourced by every numbered script. Root = this folder (anchored by .here).
# ALL substantive choices live in the cfg list below (spec Section 6).
# ============================================================================
suppressPackageStartupMessages({
  library(fs); library(readr); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble); library(janitor)
  library(survey)
})

if (!exists("dir_root")) dir_root <- getwd()   # scripts setwd(dir_root) before sourcing
root <- dir_root

source(file.path(root, "NLSY_fun", "nlsy_data_functions.R"))
source(file.path(root, "NLSY_fun", "nlsy_estimation_functions.R"))

AUDIT <- file.path(root, "Data", "NLSY_estimation", "audit")
INTER <- file.path(root, "Data/NLSY_estimation/intermediate")
ANALY <- file.path(root, "Data/NLSY_estimation/dyads")
CWID  <- file.path(ANALY, "custom_weight_ids")
for (d in c(AUDIT, INTER, ANALY, CWID)) dir_create(d)

# ---------------------------------------------------------------------------
# Configuration (spec Section 6). Single source of substantive choices.
# ---------------------------------------------------------------------------
cfg <- list(
  nlsy79_baseline_age_min = 14L, nlsy79_baseline_age_max = 18L,
  parent_rounds_79 = 1979:1983, parent_rounds_97 = 1997:2001,
  # Adult income window 27-38: extends Zhou's 27-32 now that data run through 2023.
  # Ceiling 38 = the max age the youngest NLSY97 cohort (1984) can reach, so the
  # window is FULLY comparable across both surveys/all cohorts; it adds respondents
  # and centers nearer the age-40 lifetime-income ideal (less Haider-Solon bias).
  adult_income_age_min = 27L, adult_income_age_max = 38L,
  min_parent_income_obs = 2L, min_adult_income_obs = 2L,
  occ_target_age = 40L,
  occ_age_min = 25L, occ_age_max = 55L,   # adult career window for child EGP measures
  # NLSY97 parent-income rule. Zhou (BDZ 2018) averages ALL parent-reported years
  # 1996-2000 and requires >=2 obs in both generations. Our data-quality audit
  # showed only 1997 is genuinely parent-reported (1998-2001 are noisy youth
  # reports), so the DEFAULT here is the 1997-primary coalesce -- the one
  # documented deviation from Zhou (the "prominent issue" exception). Set to
  # "zhou_avg_1996_2000_ge2" to follow Zhou strictly. Both measures are stored.
  nlsy97_parent_income_rule = "coalesce_1997_primary",
  nlsy79_cross_sectional_codes  = 1:8,
  nlsy79_bh_supplement_codes    = c(10L, 11L, 13L, 14L),   # Black/Hispanic supplement
  nlsy79_excluded_disadv_codes  = c(9L, 12L),              # econ-disadvantaged non-B/H
  nlsy79_excluded_military      = 15:20,
  nlsy97_cross_sectional_code   = 1L,
  nlsy97_supplemental_code      = 0L,
  nlsy79_weight_birth_min = 1957L, nlsy79_weight_birth_max = 1964L,
  nlsy79_report_birth_min = 1961L, nlsy79_report_birth_max = 1964L,
  base_price_year = 2023L,
  # Income-class definition (project scheme, NOT Zhou within-cohort midranks):
  # pool father + child EQUIVALIZED real income WITHIN each survey (pooled across
  # cohorts), take 5 weighted-quantile cutoffs, and use those same cutoffs to
  # assign both the parent and adult income class. Cutoffs are recomputed inside
  # every bootstrap replicate (never frozen). This is a deliberate, user-chosen
  # departure from BDZ 2018's within-generation-and-cohort mid-ranks.
  income_class_scheme = "pooled_parent_child_per_survey",
  n_income_classes = 5L,
  income_class_value = "equiv",   # "equiv" = size-adjusted (primary); "real" = unadjusted robustness
  bootstrap_seed = 20260722L, fay_rho = 0.5,
  bootstrap_smoke_B = 20L
)
set.seed(cfg$bootstrap_seed)

# small helper: rename a set of R-number columns to friendly names via a named vector
rename_refnums <- function(df, mapping) {
  present <- mapping[mapping %in% names(df)]
  df %>% rename(!!!setNames(present, names(present)))
}
read_ext <- function(path) read_csv(path, show_col_types = FALSE,
                                    name_repair = "minimal", progress = FALSE)


# ---------------------------------------------------------------------------
# Raw-extract paths + PCE deflator
# ---------------------------------------------------------------------------
p <- list(
  n79_inc  = "Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_raw.csv",
  n79_inc_lab = "Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_raw.csv",  # labelled twin (carries HHI-24)
  n79_demo = "Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_raw.csv"
  , n79_occ_raw = "Data/NLSY/NLSY79/NLSY79_occ/NLSY79_occ_raw.csv",
  n79_occ_dyad = "Data/NLSY/NLSY79/NLSY79CPSJobs/NLSY79OccupationDyadsv2.csv",
  n79_inc_dyad = "Data/NLSY/NLSY79/NLSY79HHINC/NLSY79HHIncDyads.csv",
  # NLSY97: single consolidated raw input (HHINC + Demographic + family-relationship
  # roster merged on PUBID). Income and demographic reads point to the same file.
  n97_inc  = "Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_raw.csv",
  n97_demo = "Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_raw.csv",
  n97_famrel = "Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_raw.csv",
  n97_inc_dyad = "Data/NLSY/NLSY97/HHINC Data/NLSY97HHIncDyads.csv",
  # returned custom weights (already generated by the RA)
  w79_inc = "Data/NLSY/NLSY79/NLSY79HHINC/IncomeWeights79/NLSY79IncomeWeights.csv",
  w79_occ = "Data/NLSY/NLSY79/NLSY79CPSJobs/WeightsForOcc/OccWeights.csv",
  w97_inc = "Data/NLSY/NLSY97/HHINC Data/IncomeWeights/NLSY97IncomeWeights.csv"
)
p <- lapply(p, function(x) file.path(root, x))

# PCE chain-type price index (BEA DPCERG3A086NBEA, 2017=100), for real deflation.
# Source recorded in data/pce_index.csv header. Base year = cfg$base_price_year.
pce_path <- file.path(root, "Data", "NLSY", "NLSY97", "NLSY97_inc", "DPCERG3A086NBEA.csv")
if (!file_exists(pce_path)) stop("[build] PCE file missing: Data/NLSY/DPCERG3A086NBEA.csv (download from https://fred.stlouisfed.org/series/DPCERG3A086NBEA)")
pce <- load_pce_index(pce_path, cfg$base_price_year)
message(sprintf("PCE loaded: base year %d, base index %.3f (%s)",
                cfg$base_price_year, attr(pce, "base_index"), attr(pce, "source")))
