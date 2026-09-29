# ============================================================================
# 04_generate_bootstrap_inputs.R
# ----------------------------------------------------------------------------
# Weight-attached dyads -> the FIVE self-contained bootstrap input files,
# written to Data/NLSY_estimation/inputs/ AND staged (verified) to Data/ for
# 05_hpc_bootstrap.R:
#   nlsy79_income_input.rds        (universe 1957-64, zhou domain flag, equiv$)
#   nlsy79_occ_input.rds           nlsy97_income_input.rds
#   nlsy79_income_unadj_input.rds  nlsy97_income_unadj_input.rds (no sqrt-size)
# Includes the weight validation: coverage, ~4M persons per cohort, anomaly
# flag at >2.5x the median cohort total. Requires stage 03.
# ============================================================================

# ===== part A: the three primary inputs (verbatim from the verified script) ==
# ============================================================================
# analysis/06_make_boot_inputs.R
# ----------------------------------------------------------------------------
# Build the three self-contained bootstrap input files (uploaded to Midway):
#   Data/NLSY_estimation/inputs/nlsy79_income_input.rds  id,hhid,cohort,zhou,parent,child,weight
#   Data/NLSY_estimation/inputs/nlsy97_income_input.rds  id,vstrat,vpsu,cohort,parent,child,weight
#   Data/NLSY_estimation/inputs/nlsy79_occ_input.rds     id,hhid,cohort,status_p,status_c,
#                                                      status_c_mode,status_c_dur,weight
# NLSY79 income follows the universe-then-domain design: the file carries the
# full 1957-64 complete-dyad universe (weights were requested on this universe,
# avoiding the 1960 boundary-weight anomaly); `zhou` flags the BDZ 1961-64
# domain. Incomes are the PRIMARY equivalized real measures (PCE 2023$,
# / sqrt(family size)); the NLSY79 parental average now applies the BDZ fn20
# lived-at-home screen (HHI-24 == 1; income years 1978/79/81/82).
# ============================================================================
#----------Preliminaries----------#
HPC <- TRUE   # set FALSE for a local test run
if (HPC) {
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
} else {
  dir_root <- "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
}
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)
suppressPackageStartupMessages({library(dplyr); library(readr)})
out_dir <- "Data/NLSY_estimation/inputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- NLSY79 income: universe (1957-64) + returned custom weights ----------
pp <- readRDS("Data/NLSY_estimation/dyads/nlsy79_income_person_preweight.rds")
stopifnot("caseid" %in% names(pp))
w79 <- read_csv("Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_weight.csv",
                show_col_types = FALSE)
names(w79)[1:2] <- c("caseid", "custom_weight")
w79 <- w79 %>% transmute(caseid = as.character(caseid),
                         custom_weight = as.numeric(custom_weight))

inc79 <- pp %>%
  mutate(caseid = as.character(caseid)) %>%
  filter(is_comparable_civilian_sample,          # cross-section + B/H oversamples
         complete_income_dyad, birth_year >= 1957, birth_year <= 1964) %>%
  inner_join(w79, by = "caseid") %>%
  transmute(id = caseid, hhid = as.character(hhid), cohort = birth_year,
            zhou = nlsy79_zhou_income_sample,
            parent = mean_parent_income_equiv, child = mean_adult_income_equiv,
            weight = custom_weight) %>%
  filter(!is.na(parent), !is.na(child))          # drop rows whose EQUIV income is NA
stopifnot(!any(is.na(inc79$hhid)), !any(is.na(inc79$weight)))

# ---- weight diagnostics (rerun whenever a fresh weights CSV is dropped in) --
# 1. coverage: every universe dyad must receive a weight
n_universe <- pp %>% filter(is_comparable_civilian_sample, complete_income_dyad,
                            birth_year >= 1957, birth_year <= 1964) %>% nrow()
cat(sprintf("[weights] coverage: %d of %d universe dyads weighted (%s)\n",
            nrow(inc79), n_universe,
            ifelse(nrow(inc79) == n_universe, "OK", "CHECK: some IDs missing weights")))
# 2. cohort person-totals: weights carry two implied decimals (/100 = persons);
#    each single-year birth cohort should total ~4M persons -- flags any repeat
#    of the 1960 boundary-weight anomaly (delivered 1960 weights ~4-5x scale).
wt_tot <- inc79 %>% group_by(cohort) %>%
  summarise(persons_M = round(sum(weight) / 100 / 1e6, 2),
            mean_wt = round(mean(weight)), .groups = "drop")
cat("[weights] cohort person-totals (millions; ~4M each expected, watch 1960):\n")
print(as.data.frame(wt_tot))
if (max(wt_tot$persons_M) > 2.5 * stats::median(wt_tot$persons_M))
  cat("[weights] *** ANOMALY: a cohort total is >2.5x the median -- inspect before use ***\n")

# ---- NLSY97 income --------------------------------------------------------
f97 <- readRDS("Data/NLSY_estimation/dyads/nlsy97_income_final.rds")
inc97 <- f97 %>%
  filter(has_custom_weight, complete_income_dyad) %>%
  transmute(id = as.character(pubid), vstrat, vpsu, cohort = birth_year,
            parent = parent_income_equiv, child = mean_adult_income_equiv,
            weight = custom_weight) %>%
  filter(!is.na(parent), !is.na(child))          # drop rows whose EQUIV income is NA

# ---- NLSY79 occupation ----------------------------------------------------
focc <- readRDS("Data/NLSY_estimation/dyads/nlsy79_occ_final.rds")
occ79 <- focc %>%
  filter(has_custom_weight, !is.na(father_egp), !is.na(child_egp_age40)) %>%
  transmute(id = as.character(caseid), hhid = as.character(hhid),
            cohort = birth_year,
            status_p = father_egp, status_c = child_egp_age40,
            status_c_mode = child_egp_mode, status_c_dur = child_egp_longest_dur,
            weight = custom_weight)

# ---- optional validation vs archived / pre-existing inputs (skipped if absent) ----
arch79 <- "Data/NLSY_estimation/inputs/archive_pre_screen/nlsy79_income_input.rds"
old79 <- if (file.exists(arch79)) readRDS(arch79) else NULL
old97 <- if (file.exists(file.path(out_dir, "nlsy97_income_input.rds")))
           readRDS(file.path(out_dir, "nlsy97_income_input.rds")) else NULL
oldoc <- if (file.exists(file.path(out_dir, "nlsy79_occ_input.rds")))
           readRDS(file.path(out_dir, "nlsy79_occ_input.rds")) else NULL

if (!is.null(old79)) {
cat("=== NLSY79 income: old vs new (screened) ===\n")
cat(sprintf("rows: %d -> %d | zhou domain: %d -> %d\n",
            nrow(old79), nrow(inc79), sum(old79$zhou), sum(inc79$zhou)))
cat("new ids subset of old ids (weights valid):",
    all(inc79$id %in% as.character(old79$id)), "\n")
cat("per-cohort dyads (old | new):\n")
print(full_join(old79 %>% count(cohort, name = "old"),
                inc79 %>% count(cohort, name = "new"), by = "cohort") %>%
        arrange(cohort) %>% as.data.frame())
}

samecheck <- function(new, old, label) {
  new2 <- new %>% arrange(id) %>% as.data.frame()
  old2 <- old %>% mutate(id = as.character(id)) %>% arrange(id) %>% as.data.frame()
  ok <- isTRUE(all.equal(new2[names(old2)], old2, tolerance = 1e-10,
                         check.attributes = FALSE))
  cat(sprintf("%s unchanged vs current file: %s (rows %d vs %d)\n",
              label, ok, nrow(new2), nrow(old2)))
  ok
}
ok97 <- if (!is.null(old97)) samecheck(inc97, old97, "NLSY97 income") else FALSE
okoc <- if (!is.null(oldoc)) samecheck(occ79, oldoc, "NLSY79 occ")    else FALSE

# ---- write ----------------------------------------------------------------
saveRDS(inc79, file.path(out_dir, "nlsy79_income_input.rds"))
if (!ok97) saveRDS(inc97, file.path(out_dir, "nlsy97_income_input.rds"))
if (!okoc) saveRDS(occ79, file.path(out_dir, "nlsy79_occ_input.rds"))
cat("\nwrote", file.path(out_dir, "nlsy79_income_input.rds"), "\n")

# ---- quick substantive check on the new NLSY79 income input ----------------
source("NLSY_fun/nlsy_estimation_functions.R")
ud <- function(d) { cp <- pooled_weighted_cutpoints(d$parent, d$child, d$weight, 5)
  pt <- assign_tile(d$parent, cp); ct <- assign_tile(d$child, cp); W <- sum(d$weight)
  sprintf("UP=%.3f DOWN=%.3f imm=%.3f n=%d", sum(d$weight[ct > pt]) / W,
          sum(d$weight[ct < pt]) / W, sum(d$weight[ct == pt]) / W, nrow(d)) }
cat("\nnew NLSY79 income, full universe :", ud(inc79), "\n")
cat("new NLSY79 income, Zhou domain   :", ud(inc79 %>% filter(zhou)), "\n")
if (!is.null(old79)) cat("old NLSY79 income, Zhou domain   :", ud(old79 %>% filter(zhou)), "\n")


# ===== part B: the size-unadjusted twins (verbatim) ==========================
# ============================================================================
# analysis/_make_income_unadjusted_inputs.R
# ----------------------------------------------------------------------------
# Build the FAMILY-SIZE-UNADJUSTED (real, size-undivided) income bootstrap
# inputs for the parallel analysis. These are exact clones of the primary
# equivalized inputs -- SAME dyads, weights, cluster ids (hhid / vstrat+vpsu),
# cohort, and (NLSY79) zhou domain flag -- with ONLY the parent/child income
# columns swapped from income_EQUIV (real / sqrt(family_size)) to income_REAL
# (the same PCE-2023$ real income WITHOUT the sqrt-size denominator).
#
# "Real" here = mean over the qualifying years of the per-year real income
# (mean_*_income_real), i.e. the denominator removed at the per-year level then
# averaged -- NOT mean_equiv * sqrt(size) (size varies by year).
#
# Writes to BOTH the repo boot_inputs/ and the Midway staging Data/ folder:
#   nlsy79_income_unadj_input.rds   id,hhid,cohort,zhou,parent,child,weight
#   nlsy97_income_unadj_input.rds   id,vstrat,vpsu,cohort,parent,child,weight
# ============================================================================
suppressPackageStartupMessages(library(dplyr))
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)

repo_dir  <- "Data/NLSY_estimation/inputs"
dir.create(repo_dir, recursive = TRUE, showWarnings = FALSE)
midway    <- repo_dir   # single location: Data/NLSY_estimation/inputs (read directly by 05)

# --- real (size-unadjusted) lookups, keyed by id -----------------------------
r79 <- readRDS("Data/NLSY_estimation/dyads/nlsy79_income_person_preweight.rds") %>%
  transmute(id = as.character(caseid),
            parent_real = mean_parent_income_real,
            child_real  = mean_adult_income_real)
r97 <- readRDS("Data/NLSY_estimation/dyads/nlsy97_income_final.rds") %>%
  transmute(id = as.character(pubid),
            parent_real = parent_income_real,
            child_real  = mean_adult_income_real)

swap_to_real <- function(boot, lookup, label) {
  out <- boot %>% left_join(lookup, by = "id") %>%
    mutate(parent = parent_real, child = child_real) %>%
    select(-parent_real, -child_real)
  stopifnot(identical(names(out), names(boot)),          # schema unchanged
            nrow(out) == nrow(boot),                     # same dyads
            !any(is.na(out$parent)), !any(is.na(out$child)))
  cat(sprintf("%s: %d dyads | parent/child now size-UNADJUSTED real$ | weights & ids identical\n",
              label, nrow(out)))
  out
}

b79 <- readRDS(file.path(repo_dir, "nlsy79_income_input.rds"))
b97 <- readRDS(file.path(repo_dir, "nlsy97_income_input.rds"))
u79 <- swap_to_real(b79, r79, "NLSY79")
u97 <- swap_to_real(b97, r97, "NLSY97")

for (d in c(repo_dir, midway)) {
  saveRDS(u79, file.path(d, "nlsy79_income_unadj_input.rds"))
  saveRDS(u97, file.path(d, "nlsy97_income_unadj_input.rds"))
  cat("wrote:", file.path(d, "nlsy79_income_unadj_input.rds"),
      "and nlsy97_income_unadj_input.rds\n")
}

# --- quick sanity: unadjusted vs adjusted up/down on the reporting domain -----
source("NLSY_fun/nlsy_estimation_functions.R")
ud <- function(d, cohs) { d <- d %>% filter(cohort %in% cohs)
  cp <- pooled_weighted_cutpoints(d$parent, d$child, d$weight, 5)
  pt <- assign_tile(d$parent, cp); ct <- assign_tile(d$child, cp); W <- sum(d$weight)
  sprintf("UP=%.3f DOWN=%.3f imm=%.3f", sum(d$weight[ct>pt])/W, sum(d$weight[ct<pt])/W, sum(d$weight[ct==pt])/W) }
cat("\nNLSY79 1961-64  adjusted  :", ud(b79, 1961:1964), "\n")
cat("NLSY79 1961-64  UNadjusted:", ud(u79, 1961:1964), "\n")
cat("NLSY97 1980-84  adjusted  :", ud(b97, 1980:1984), "\n")
cat("NLSY97 1980-84  UNadjusted:", ud(u97, 1980:1984), "\n")


# ===== part C: confirm all five inputs are in place for 05_hpc_bootstrap =====
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)
for (f in c("nlsy79_income_input.rds", "nlsy79_occ_input.rds", "nlsy97_income_input.rds",
            "nlsy79_income_unadj_input.rds", "nlsy97_income_unadj_input.rds"))
  stopifnot(file.exists(file.path("Data/NLSY_estimation/inputs", f)))
cat("\n[04] five bootstrap inputs written to Data/NLSY_estimation/inputs (read directly by 05)\n")
