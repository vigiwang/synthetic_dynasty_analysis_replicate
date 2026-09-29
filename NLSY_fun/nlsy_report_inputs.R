# ============================================================================
# NLSY_fun/nlsy_report_inputs.R
# ----------------------------------------------------------------------------
# The nine post-bootstrap "report input" computations, consolidated into one
# file SOURCED by stage 06 (working directory = dir_root):
#   rank-rank (+ BDZ-exact replication), steady-state exchange (adjusted +
#   unadjusted), MTE delta-method SE three-way check, previous-construction
#   up/down decomposition (OPTIONAL - skips without RA dyads), adjusted-vs-
#   unadjusted up/down, family-size table, calendar income windows.
# All outputs -> Data/NLSY_estimation/intermediate/.
# ============================================================================
## ===================== [precompute_rankrank] =====================
# ============================================================================
# NLSY/_precompute_rankrank.R
# ----------------------------------------------------------------------------
# Rank-rank slope (BDZ 2018 positional persistence measure) for the current
# construction, with the survey design (household clusters for NLSY79; VSTRAT /
# VPSU for NLSY97), plus the previous RA construction for a sample comparison.
# Incomes ranked into weighted mid-percentiles WITHIN GENERATION AND COHORT;
# slope from a survey-weighted regression of child rank on parent rank.
# The rank-rank is tail-robust, so it is unaffected by the topcode / Pareto
# issue that makes the log-income elasticity non-comparable.
# Saves:
#   Data/NLSY_estimation/intermediate/rankrank_by_cohort.rds  (ours: survey, cohort, rr, rr_se)
#   Data/NLSY_estimation/intermediate/rankrank_compare.rds    (RA vs ours vs BDZ, per survey)
# ============================================================================
suppressPackageStartupMessages({library(dplyr); library(survey)})
options(survey.lonely.psu = "adjust")

wrank <- function(x, w) { o <- order(x); r <- numeric(length(x))
  r[o] <- (cumsum(w[o]) - 0.5 * w[o]) / sum(w); r }
addrank <- function(d) d %>% group_by(cohort) %>%
  mutate(rp = wrank(parent, weight), rc = wrank(child, weight)) %>% ungroup()

# ---- ours: per-cohort + Full rank-rank with design SE ----------------------
d79 <- readRDS("Data/NLSY_estimation/inputs/nlsy79_income_input.rds") %>%
  filter(zhou, parent > 0, child > 0) %>% transmute(cohort, parent, child, weight, hhid) %>% addrank()
d97 <- readRDS("Data/NLSY_estimation/inputs/nlsy97_income_input.rds") %>%
  filter(parent > 0, child > 0) %>% transmute(cohort, parent, child, weight, vstrat, vpsu) %>% addrank()
des79 <- svydesign(ids = ~hhid, weights = ~weight, data = d79)
des97 <- svydesign(ids = ~vpsu, strata = ~vstrat, weights = ~weight, data = d97, nest = TRUE)

rr_series <- function(des, cohorts, label) {
  fit <- function(sub) { m <- svyglm(rc ~ rp, design = sub)
    tibble::tibble(rr = unname(coef(m)[2]), rr_se = unname(SE(m)[2]), n = nrow(sub)) }
  bind_rows(
    lapply(cohorts, function(cc) tibble::tibble(survey = label, cohort = as.character(cc)) %>%
             bind_cols(fit(subset(des, cohort == cc)))),
    list(tibble::tibble(survey = label, cohort = "Full") %>% bind_cols(fit(des))))
}
rr_cohort <- bind_rows(rr_series(des79, 1961:1964, "NLSY79"),
                       rr_series(des97, 1980:1984, "NLSY97"))
saveRDS(rr_cohort, "Data/NLSY_estimation/intermediate/rankrank_by_cohort.rds")

# ---- RA previous version (household income), Full slope --------------------
rap <- "Data/NLSY/previous_construction"
ra_slope <- function(f, idcohort) {
  d <- read.csv(file.path(rap, f)) %>%
    transmute(cohort = .data[[idcohort]], parent = Lifetime_ParInc_Real,
              child = Lifetime_ChildInc_Real, weight = Weights / 100) %>%
    filter(parent > 0, child > 0, is.finite(weight), weight > 0) %>% addrank()
  unname(coef(lm(rc ~ rp, data = d, weights = weight))[2])
}
ra_ok <- file.exists(file.path(rap, "NLSY79HHIncDyads.csv"))
ra79 <- if (ra_ok) ra_slope("NLSY79HHIncDyads.csv", "BirthYear") else NA_real_
ra97 <- if (ra_ok) ra_slope("NLSY97HHIncDyads.csv", "BirthYear") else NA_real_
if (!ra_ok) message("[precompute] RA dyads absent - 'RA previous' column set NA (optional)")

full <- rr_cohort %>% filter(cohort == "Full")
rr_compare <- tibble::tibble(
  Survey        = c("NLSY79", "NLSY97"),
  `RA previous` = c(ra79, ra97),
  Ours          = full$rr,
  ours_lo       = full$rr - 1.96 * full$rr_se,
  ours_hi       = full$rr + 1.96 * full$rr_se,
  `BDZ Table 1` = c(0.426, 0.424),
  bdz_lo        = c(0.395, 0.400),
  bdz_hi        = c(0.453, 0.446))
rr_compare$inside_bdz <- rr_compare$Ours >= rr_compare$bdz_lo & rr_compare$Ours <= rr_compare$bdz_hi
saveRDS(rr_compare, "Data/NLSY_estimation/intermediate/rankrank_compare.rds")
print(as.data.frame(rr_compare), digits = 3)


## ===================== [precompute_rankrank_bdzexact] =====================
# ============================================================================
# NLSY/_precompute_rankrank_bdzexact.R
# ----------------------------------------------------------------------------
# The FULLY BDZ-aligned rank-rank replication: every Bloome-Dyer-Zhou 2018
# choice matched at once —
#   * cohorts NLSY79 born 1961-64, NLSY97 born 1980-84 (BDZ's exact ranges);
#   * adult income averaged over ages 27-32 (BDZ window), >=2 obs;
#   * family-size-adjusted (equivalized) real income (BDZ /sqrt(family size));
#   * BDZ parental-income rule = the multi-year AVERAGE of the parent-reported
#     rounds (NLSY79 income yrs 1978-82 lived-at-home; NLSY97 1996-2000), the
#     "average all years in which parents reported their income" of BDZ fn20.
# Incomes are ranked into weighted mid-percentiles WITHIN generation AND cohort
# (same machinery as _precompute_rankrank.R); slope from a weighted regression.
# NLSY79 uses the 1957-64 universe weights restricted to the 1961-64 domain
# (avoids the low-boundary custom-weight inflation); NLSY97 the on-file weights.
# We report the parent >=2 rule (BDZ literal) and >=1 (effective) because BDZ's
# reported NLSY97 value aligns with the latter (parents reported only in 1997).
# Saves Data/NLSY_estimation/intermediate/rankrank_bdzexact.rds.
# ============================================================================
suppressPackageStartupMessages(library(dplyr))

l79 <- readRDS("Data/NLSY_estimation/intermediate/nlsy79_income_long.rds") %>% mutate(id = as.character(caseid))
l97 <- readRDS("Data/NLSY_estimation/intermediate/nlsy97_income_long.rds") %>% mutate(id = as.character(pubid))
w79 <- readRDS("Data/NLSY_estimation/dyads/nlsy79_income_final.rds") %>% transmute(id = as.character(caseid), w = custom_weight, by = birth_year)
w97 <- readRDS("Data/NLSY_estimation/dyads/nlsy97_income_final.rds") %>% transmute(id = as.character(pubid), w = custom_weight, by = birth_year)

par_bdz <- function(l, minobs) l %>% filter(is_parent_obs) %>% group_by(id) %>%
  summarise(parent = mean(income_equiv, na.rm = TRUE), np = n(), .groups = "drop") %>% filter(np >= minobs)
adult2732 <- function(l) l %>% filter(is_valid_income, age_in_income_year >= 27, age_in_income_year <= 32) %>%
  group_by(id) %>% summarise(child = mean(income_equiv, na.rm = TRUE), nc = n(), .groups = "drop") %>% filter(nc >= 2)

wrank <- function(x, w) { o <- order(x); r <- numeric(length(x)); r[o] <- (cumsum(w[o]) - 0.5 * w[o]) / sum(w); r }
rr <- function(l, w, cohs, minpar) {
  d <- par_bdz(l, minpar) %>% inner_join(adult2732(l), by = "id") %>% inner_join(w, by = "id") %>%
    filter(by %in% cohs, parent > 0, child > 0, !is.na(parent), !is.na(child), w > 0) %>%
    group_by(by) %>% mutate(rp = wrank(parent, w), rc = wrank(child, w)) %>% ungroup()
  list(b = unname(coef(lm(rc ~ rp, d, weights = w))[2]), n = nrow(d)) }

a2 <- rr(l79, w79, 1961:1964, 2); a1 <- rr(l79, w79, 1961:1964, 1)
b2 <- rr(l97, w97, 1980:1984, 2); b1 <- rr(l97, w97, 1980:1984, 1)

rankrank_bdzexact <- tibble::tibble(
  Survey    = c("NLSY79 (1961–64)", "NLSY97 (1980–84)"),
  beta_ge2  = c(a2$b, b2$b), n_ge2 = c(a2$n, b2$n),
  beta_ge1  = c(a1$b, b1$b), n_ge1 = c(a1$n, b1$n),
  bdz       = c(0.426, 0.424), bdz_lo = c(0.395, 0.400), bdz_hi = c(0.453, 0.446))
saveRDS(rankrank_bdzexact, "Data/NLSY_estimation/intermediate/rankrank_bdzexact.rds")
print(as.data.frame(rankrank_bdzexact), digits = 3)


## ===================== [precompute_ss_exchange] =====================
# ============================================================================
# NLSY/_precompute_ss_exchange.R
# ----------------------------------------------------------------------------
# Exchange mobility AT THE STEADY STATE, with a bias-corrected 95% CI.
# The estimator returns each movement component as a vector over generations t
# (synthetic_dynasty_estimation.R). At the steady state (t = t_max) the marginal
# distribution has converged, so structural mobility -> 0 and exchange mobility
# equals total mobility: the long-run "circulation" mobility of the copula, net
# of any marginal shift. We take the last element of the Exchange Mobility vector
# from the baseline and from every bootstrap replicate, then apply the SAME
# bias-correction the report uses (mean_bc = 2*base - boot_mean; CI_bc from the
# reflected boot percentiles), with Fay-BRR rescaling (rho = 0.5) for NLSY97.
# Verified to reproduce the processed t=1 exchange CIs exactly.
# Saves Data/NLSY_estimation/intermediate/ss_exchange.rds (validation_lst-style, keyed by cohort).
# ============================================================================
suppressPackageStartupMessages(library(dplyr))

ex_ss <- function(r) { v <- r[["Movement Measure"]][["Exchange Mobility"]]
  as.numeric(v[length(v)]) }                       # steady-state (t_max) value

ss_build <- function(baseline_f, boot_f, cohs, rho = NULL) {
  base <- readRDS(baseline_f); boot <- readRDS(boot_f)
  bv <- sapply(cohs, function(c) ex_ss(base[[c]]))
  bm <- t(sapply(boot, function(rep) sapply(cohs, function(c) {
    r <- rep[[c]]; if (is.list(r)) ex_ss(r) else NA_real_ })))
  if (!is.null(rho))                               # Fay-BRR: rescale deviations 1/(1-rho)
    bm <- sweep(bm, 2, bv) / (1 - rho) + matrix(bv, nrow(bm), ncol(bm), byrow = TRUE)
  mb  <- colMeans(bm, na.rm = TRUE)
  qlo <- apply(bm, 2, quantile, 0.975, na.rm = TRUE)
  qhi <- apply(bm, 2, quantile, 0.025, na.rm = TRUE)
  setNames(lapply(seq_along(cohs), function(i) list(ss_df = tibble::tibble(
    mean_bc = 2 * bv[i] - mb[i],
    CI_bc_l = 2 * bv[i] - qlo[i],
    CI_bc_u = 2 * bv[i] - qhi[i]))), cohs)
}

ss79 <- ss_build("Data/NLSY_estimation/nlsy79_income_baseline.rds",
                 "Data/NLSY_estimation/nlsy79_income_boot.rds",
                 c("Full", as.character(1961:1964)), rho = NULL)
ss97 <- ss_build("Data/NLSY_estimation/nlsy97_income_baseline.rds",
                 "Data/NLSY_estimation/nlsy97_income_boot.rds",
                 c("Full", as.character(1980:1984)), rho = 0.5)
saveRDS(list(ss79 = ss79, ss97 = ss97), "Data/NLSY_estimation/intermediate/ss_exchange.rds")

cat("Steady-state exchange mobility (mean_bc [CI]):\n")
show <- function(ss, cohs) for (c in cohs) { d <- ss[[c]]$ss_df
  cat(sprintf("  %-5s %.3f [%.3f, %.3f]\n", c, d$mean_bc, d$CI_bc_l, d$CI_bc_u)) }
cat("NLSY79:\n"); show(ss79, c("Full", as.character(1961:1964)))
cat("NLSY97:\n"); show(ss97, c("Full", as.character(1980:1984)))


## ===================== [precompute_ss_exchange_unadj] =====================
# ============================================================================
# NLSY/_precompute_ss_exchange_unadj.R
# ----------------------------------------------------------------------------
# Steady-state exchange mobility (bias-corrected 95% CI) for the FAMILY-SIZE
# UNADJUSTED income runs, mirroring NLSY/_precompute_ss_exchange.R exactly
# but reading the *_unadj_* bootstrap objects. Saves
# Data/NLSY_estimation/intermediate/ss_exchange_unadj.rds (list ss79/ss97, keyed by cohort).
# ============================================================================
suppressPackageStartupMessages(library(dplyr))

ex_ss <- function(r) { v <- r[["Movement Measure"]][["Exchange Mobility"]]
  as.numeric(v[length(v)]) }                       # steady-state (t_max) value

ss_build <- function(baseline_f, boot_f, cohs, rho = NULL) {
  base <- readRDS(baseline_f); boot <- readRDS(boot_f)
  bv <- sapply(cohs, function(c) ex_ss(base[[c]]))
  bm <- t(sapply(boot, function(rep) sapply(cohs, function(c) {
    r <- rep[[c]]; if (is.list(r)) ex_ss(r) else NA_real_ })))
  if (!is.null(rho))                               # Fay-BRR: rescale deviations 1/(1-rho)
    bm <- sweep(bm, 2, bv) / (1 - rho) + matrix(bv, nrow(bm), ncol(bm), byrow = TRUE)
  mb  <- colMeans(bm, na.rm = TRUE)
  qlo <- apply(bm, 2, quantile, 0.975, na.rm = TRUE)
  qhi <- apply(bm, 2, quantile, 0.025, na.rm = TRUE)
  setNames(lapply(seq_along(cohs), function(i) list(ss_df = tibble::tibble(
    mean_bc = 2 * bv[i] - mb[i],
    CI_bc_l = 2 * bv[i] - qlo[i],
    CI_bc_u = 2 * bv[i] - qhi[i]))), cohs)
}

ss79 <- ss_build("Data/NLSY_estimation/nlsy79_income_unadj_baseline.rds",
                 "Data/NLSY_estimation/nlsy79_income_unadj_boot.rds",
                 c("Full", as.character(1961:1964)), rho = NULL)
ss97 <- ss_build("Data/NLSY_estimation/nlsy97_income_unadj_baseline.rds",
                 "Data/NLSY_estimation/nlsy97_income_unadj_boot.rds",
                 c("Full", as.character(1980:1984)), rho = 0.5)
saveRDS(list(ss79 = ss79, ss97 = ss97), "Data/NLSY_estimation/intermediate/ss_exchange_unadj.rds")

cat("Steady-state exchange mobility, UNADJUSTED (mean_bc [CI]):\n")
show <- function(ss, cohs) for (c in cohs) { d <- ss[[c]]$ss_df
  cat(sprintf("  %-5s %.3f [%.3f, %.3f]\n", c, d$mean_bc, d$CI_bc_l, d$CI_bc_u)) }
cat("NLSY79:\n"); show(ss79, c("Full", as.character(1961:1964)))
cat("NLSY97:\n"); show(ss97, c("Full", as.character(1980:1984)))


## ===================== [precompute_mte_delta] =====================
# ============================================================================
# NLSY/_precompute_mte_delta.R
# ----------------------------------------------------------------------------
# Validation of the bootstrap standard errors for the per-class MTE (mean time
# to exit) against the analytical DELTA-METHOD SE.
#
# The estimator defines  MTE_i = 1 / (1 - p_ii),  where p_ii = P(child in class
# i | parent in class i) is the diagonal of the row-conditional transition
# matrix (synthetic_dynasty_estimation.R, line 410). With g(p) = 1/(1-p),
# g'(p) = 1/(1-p)^2 = MTE_i^2, so by the delta method
#         SE(MTE_i) = MTE_i^2 * SE(p_ii),
# where SE(p_ii) is the DESIGN-BASED SE of the diagonal conditional proportion
# (survey mean of the indicator child==i over the domain parent==i), using the
# same design as the bootstrap: household clusters for NLSY79, VSTRAT/VPSU for
# NLSY97. We compare this to the bootstrap SD (MTE_df$se_adj), pooled "Full".
# Saves Data/NLSY_estimation/intermediate/mte_delta_check.rds (survey, class, mte, se_boot,
# se_delta, ratio).
# ============================================================================
suppressPackageStartupMessages({library(dplyr); library(survey)})
source("NLSY_fun/nlsy_estimation_functions.R")
options(survey.lonely.psu = "adjust")

# assign 5 income classes with the estimator's scheme: pooled parent+child
# weighted quintiles, then DESCENDING so Class 1 = highest income.
add_classes <- function(d) {
  cp <- pooled_weighted_cutpoints(d$parent, d$child, d$weight, 5)
  d$pc <- 6L - assign_tile(d$parent, cp)     # 6 - tile -> Class1 = highest
  d$cc <- 6L - assign_tile(d$child,  cp)
  d
}

# per-class delta-method SE(MTE) from a survey design
mte_delta <- function(des, label) {
  purrr::map_dfr(1:5, function(k) {
    sub <- subset(des, pc == k)
    m   <- svymean(~I(cc == k), sub)          # p_ii = P(child=k | parent=k)
    p   <- coef(m)[["I(cc == k)TRUE"]]; se_p <- SE(m)[["I(cc == k)TRUE"]]
    mte <- 1 / (1 - p)
    tibble::tibble(survey = label, class = paste0("Class", k),
                   p_ii = p, mte = mte, se_pii = se_p, se_delta = mte^2 * se_p)
  })
}

# NLSY79: the bootstrap now defines the pooled "Full" on the 1961-64 reporting
# domain (the 1957-60 cohorts are weighting scaffolding only), so the analytical
# MTE is computed on the same Zhou domain and reproduces the bootstrap baseline.
d79 <- add_classes(readRDS("Data/NLSY_estimation/inputs/nlsy79_income_input.rds") %>%
                     filter(zhou, parent > 0, child > 0))
d97 <- add_classes(readRDS("Data/NLSY_estimation/inputs/nlsy97_income_input.rds") %>%
                     filter(parent > 0, child > 0))
des79 <- svydesign(ids = ~hhid, weights = ~weight, data = d79)
des97 <- svydesign(ids = ~vpsu, strata = ~vstrat, weights = ~weight, data = d97, nest = TRUE)

delta <- bind_rows(mte_delta(des79, "NLSY79"), mte_delta(des97, "NLSY97"))

# bootstrap SE per class from the processed MTE_df (pooled "Full")
boot_se <- function(f, label) {
  vl <- readRDS(f)$validation$validation_lst
  vl[["Full"]]$MTE_df %>% filter(class != "AMTE") %>%
    transmute(survey = label, class, mte_boot = mean_baseline, se_boot = se_adj)
}
bse <- bind_rows(boot_se("Data/NLSY_estimation/processed/nlsy79_income_processed.rds", "NLSY79"),
                 boot_se("Data/NLSY_estimation/processed/nlsy97_income_processed.rds", "NLSY97"))

# individual (naive i.i.d.) bootstrap SE per class, same measure, for the
# three-way comparison: analytical delta vs design bootstrap vs naive bootstrap
ise <- bind_rows(boot_se("Data/NLSY_estimation/processed/nlsy79_income_ind_processed.rds", "NLSY79"),
                 boot_se("Data/NLSY_estimation/processed/nlsy97_income_ind_processed.rds", "NLSY97")) %>%
  rename(se_ind = se_boot) %>% select(survey, class, se_ind)

chk <- delta %>% left_join(bse, by = c("survey", "class")) %>%
  left_join(ise, by = c("survey", "class")) %>%
  transmute(survey, class, mte = round(mte, 3), se_boot = round(se_boot, 4),
            se_delta = round(se_delta, 4), se_ind = round(se_ind, 4),
            ratio = round(se_delta / se_boot, 2),
            ratio_ind = round(se_boot / se_ind, 2))
saveRDS(chk, "Data/NLSY_estimation/intermediate/mte_delta_check.rds")
print(as.data.frame(chk))


## ===================== [precompute_updown_decomp] =====================
# ============================================================================
# Up/down decomposition vs the PREVIOUS income construction (both surveys).
#
# VALIDATED replica: reproduces the previous version's spec from raw data --
#   child window : previous filter is on SURVEY year (birth+27..birth+39);
#                  TNFI/CV_INCOME report the PRIOR calendar year, so the
#                  effective income-age window is 26-38.
#   parent NLSY79: 4 survey rounds (1979/80/82/83) restricted to respondents
#                  living at home (HHI.24 == 1); survey 1981 skipped because
#                  HHI.24_1981 is absent from the extract.
#   parent NLSY97: 1997 gross household income primary, 1998-2001 mean fallback
#                  (same coalesce rule as ours).
#   CPI-U deflation, no size adjustment, previous sample + weights.
# Replica check: NLSY79 0.258/0.431 vs published 0.260/0.426;
#                NLSY97 matches to ~0.01.
# Ladder: replica -> our observation rules -> PCE -> equivalization.
# Saves Data/NLSY_estimation/intermediate/updown_decomposition.rds.
# ============================================================================
suppressPackageStartupMessages({library(dplyr)})
source("NLSY_fun/nlsy_estimation_functions.R")

price <- tibble::tribble(~y, ~pce, ~cpi,
  1978,22.826,65.2, 1979,24.653,72.6, 1980,27.288,82.4, 1981,29.694,90.9, 1982,31.339,96.5,
  1983,32.634,99.6, 1984,33.913,103.9,1985,35.140,107.6,1986,35.892,109.6,1987,37.010,113.6,
  1988,38.393,118.3,1989,40.028,124.0,1990,41.696,130.7,1991,43.129,136.2,1992,44.140,140.3,
  1993,45.207,144.5,1994,46.104,148.2,1995,47.045,152.4,1996,48.055,156.9,1997,48.937,160.5,
  1998,49.360,163.0,1999,50.093,166.6,2000,51.395,172.2,2001,52.545,177.1,2002,53.271,179.9,
  2003,54.416,184.0,2004,55.965,188.9,2005,57.892,195.3,2006,59.610,201.6,2007,61.395,207.3,
  2008,63.243,215.3,2009,63.170,214.5,2010,64.391,218.1,2011,66.087,224.9,2012,67.343,229.6,
  2013,68.263,233.0,2014,69.470,236.7,2015,69.711,237.0,2016,70.487,240.0,2017,71.916,245.1,
  2018,73.560,251.1,2019,74.777,255.7,2020,75.884,258.8,2021,79.017,271.0,2022,84.462,292.7,
  2023,88.017,304.7,2024,90.223,313.7,2025,92.400,321.9)
pidx <- function(y, i) approx(price$y, price[[i]], y, rule = 2)$y
BASE <- pidx(2023, "pce")

rap <- "Data/NLSY/previous_construction"
run_updown_decomp <- file.exists(file.path(rap, "NLSY79HHIncDyads.csv"))
if (!run_updown_decomp)
  message("[precompute] RA dyads absent - previous-construction decomposition skipped (optional)")
if (run_updown_decomp) {
ra79 <- read.csv(file.path(rap, "NLSY79HHIncDyads.csv")) %>%
  transmute(id = as.character(CASEID_1979), ra_par = as.numeric(Lifetime_ParInc_Real),
            ra_chd = as.numeric(Lifetime_ChildInc_Real), ra_w = Weights / 100)
ra97 <- read.csv(file.path(rap, "NLSY97HHIncDyads.csv")) %>%
  transmute(id = as.character(PUBID_1997), ra_par = as.numeric(Lifetime_ParInc_Real),
            ra_chd = as.numeric(Lifetime_ChildInc_Real), ra_w = Weights / 100)
l79 <- readRDS("Data/NLSY_estimation/intermediate/nlsy79_income_long.rds") %>% mutate(id = as.character(caseid))
l97 <- readRDS("Data/NLSY_estimation/intermediate/nlsy97_income_long.rds") %>% mutate(id = as.character(pubid))
# published-series boot inputs (comparable-civilian sample, validated custom weights)
bi79 <- readRDS("Data/NLSY_estimation/inputs/nlsy79_income_input.rds") %>% dplyr::filter(zhou)
bi97 <- readRDS("Data/NLSY_estimation/inputs/nlsy97_income_input.rds")

# previous-version NLSY79 parent: 4 screened lived-at-home rounds, from the extract
ext79 <- read.csv("Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_raw.csv")
prev_par79 <- function(defl) {
  vals <- sapply(c(1979, 1980, 1982, 1983), function(sy) {
    tn <- ext79[[paste0("TNFI_TRUNC_", sy)]]; hh <- ext79[[paste0("HHI.24_", sy)]]
    ifelse(!is.na(hh) & hh == 1 & !is.na(tn) & tn >= 0,
           tn * BASE / pidx(sy - 1, defl), NA) })
  tibble(id = as.character(ext79$CASEID_1979), par = rowMeans(vals, na.rm = TRUE))
}

chd_agg <- function(long, alo, ahi, defl, eq) long %>%
  filter(is_valid_income, age_in_income_year >= alo, age_in_income_year <= ahi) %>%
  mutate(r = income_nominal * BASE / pidx(income_year, defl),
         v = if (eq) r / sqrt(family_size) else r) %>%
  group_by(id) %>% summarise(chd = mean(v, na.rm = TRUE), .groups = "drop")
par_agg <- function(long, defl, eq) long %>% filter(is_parent_obs, is_valid_income) %>%
  mutate(r = income_nominal * BASE / pidx(income_year, defl),
         v = if (eq) r / sqrt(family_size) else r) %>%
  group_by(id) %>% summarise(par = mean(v, na.rm = TRUE), .groups = "drop")

ud <- function(d) { d <- d %>% filter(is.finite(par), is.finite(chd), is.finite(ra_w),
                                      par > 0, chd > 0)
  cp <- pooled_weighted_cutpoints(d$par, d$chd, d$ra_w, 5)
  pt <- assign_tile(d$par, cp); ct <- assign_tile(d$chd, cp); W <- sum(d$ra_w)
  c(up = sum(d$ra_w[ct > pt]) / W, down = sum(d$ra_w[ct < pt]) / W,
    imm = sum(d$ra_w[ct == pt]) / W,
    ratio = weighted.mean(d$chd, d$ra_w) / weighted.mean(d$par, d$ra_w)) }

ud_input <- function(d) {   # up/down on a boot-input file (parent/child/weight cols)
  d <- d %>% filter(is.finite(parent), is.finite(child), parent > 0, child > 0)
  cp <- pooled_weighted_cutpoints(d$parent, d$child, d$weight, 5)
  pt <- assign_tile(d$parent, cp); ct <- assign_tile(d$child, cp); W <- sum(d$weight)
  c(up = sum(d$weight[ct > pt]) / W, down = sum(d$weight[ct < pt]) / W,
    imm = sum(d$weight[ct == pt]) / W,
    ratio = weighted.mean(d$child, d$weight) / weighted.mean(d$parent, d$weight))
}

decomp_one <- function(long, ra, label, replica_par, boot_input) {
  # replica: previous spec (child income-age 26-38; previous parent rule; CPI)
  rep_par <- if (label == "NLSY79") replica_par("cpi") else par_agg(long, "cpi", FALSE)
  r_pub <- ud(ra %>% mutate(par = ra_par, chd = ra_chd))
  r_rep <- ud(ra %>% inner_join(rep_par, by = "id") %>%
                inner_join(chd_agg(long, 26, 38, "cpi", FALSE), by = "id"))
  mk <- function(defl, eq) ra %>% inner_join(par_agg(long, defl, eq), by = "id") %>%
    inner_join(chd_agg(long, 27, 38, defl, eq), by = "id")
  r_obs <- ud(mk("cpi", FALSE)); r_pce <- ud(mk("pce", FALSE)); r_eqv <- ud(mk("pce", TRUE))
  r_pub2 <- ud_input(boot_input)
  obs_lab <- "+ our child window (27-38; parent rule now identical)"
  tibble::tibble(survey = label,
    step = c("Previous version, as published",
             "Replica of previous version (validated)",
             obs_lab,
             "+ PCE deflator", "+ family-size equivalization",
             "+ comparable-civilian sample & custom weights (= published series)"),
    up = c(r_pub["up"], r_rep["up"], r_obs["up"], r_pce["up"], r_eqv["up"], r_pub2["up"]),
    down = c(r_pub["down"], r_rep["down"], r_obs["down"], r_pce["down"], r_eqv["down"], r_pub2["down"]),
    ratio = c(r_pub["ratio"], r_rep["ratio"], r_obs["ratio"], r_pce["ratio"], r_eqv["ratio"], r_pub2["ratio"]),
    immobility = c(r_pub["imm"], r_rep["imm"], r_obs["imm"], r_pce["imm"], r_eqv["imm"], r_pub2["imm"])) %>%
    mutate(delta_up = up - dplyr::lag(up))
}

decomp <- bind_rows(decomp_one(l79, ra79, "NLSY79", prev_par79, bi79),
                    decomp_one(l97, ra97, "NLSY97", NULL, bi97))
saveRDS(decomp, "Data/NLSY_estimation/intermediate/updown_decomposition.rds")
print(as.data.frame(decomp), digits = 3)
}


## ===================== [precompute_unadjusted_updown] =====================
# ============================================================================
# NLSY/_precompute_unadjusted_updown.R
# ----------------------------------------------------------------------------
# Baseline (point-estimate) upward vs downward mobility computed with the
# SAME machinery as the headline series (pooled weighted cut-points ->
# per-cohort weighted transition matrices -> single_simulation), but on the
# NO-family-size-adjustment income (real, size-UNadjusted) alongside the
# equivalized headline. Shows how much of the upward tilt is the equivalization.
#   NLSY79: mean_parent_income_real  / mean_adult_income_real   (zhou domain)
#   NLSY97: parent_income_real       / mean_adult_income_real
# Saves Data/NLSY_estimation/intermediate/unadjusted_updown.rds (tidy: survey,cohort,adj,up,down).
# ============================================================================
suppressPackageStartupMessages(library(dplyr))
source("NLSY_fun/nlsy_estimation_functions.R")
source("NLSY_fun/nlsy_data_functions.R")
source("NLSY_fun/nlsy_estimation_functions.R")
suppressWarnings(suppressMessages(load_analyst_estimators()))

# up/down (t = 1) per cohort + Full from a parent/child/weight/cohort frame
updown_baseline <- function(d, survey, adj) {
  cohorts <- sort(unique(d$cohort))
  bm <- baseline_measures(d, "parent", "child", "weight", "cohort", cohorts, type = "income")
  purrr::map_dfr(names(bm), function(nm) {
    mv <- tryCatch(bm[[nm]][["Movement Measure"]], error = function(e) NULL)
    if (is.null(mv)) return(NULL)
    tibble::tibble(survey = survey, cohort = nm, adj = adj,
                   up = mv[["Historical Upward Mobility"]][1],
                   down = mv[["Historical Downward Mobility"]][1])
  })
}

f79 <- readRDS("Data/NLSY_estimation/dyads/nlsy79_income_final.rds") %>%
  filter(nlsy79_zhou_income_sample)
n79_eq <- f79 %>% filter(!is.na(mean_parent_income_equiv), !is.na(mean_adult_income_equiv)) %>%
  transmute(cohort = birth_year, weight = custom_weight,
            parent = mean_parent_income_equiv, child = mean_adult_income_equiv)
n79_un <- f79 %>% filter(!is.na(mean_parent_income_real), !is.na(mean_adult_income_real)) %>%
  transmute(cohort = birth_year, weight = custom_weight,
            parent = mean_parent_income_real, child = mean_adult_income_real)

f97 <- readRDS("Data/NLSY_estimation/dyads/nlsy97_income_final.rds") %>% filter(has_custom_weight)
n97_eq <- f97 %>% filter(!is.na(parent_income_equiv), !is.na(mean_adult_income_equiv)) %>%
  transmute(cohort = birth_year, weight = custom_weight,
            parent = parent_income_equiv, child = mean_adult_income_equiv)
n97_un <- f97 %>% filter(!is.na(parent_income_real), !is.na(mean_adult_income_real)) %>%
  transmute(cohort = birth_year, weight = custom_weight,
            parent = parent_income_real, child = mean_adult_income_real)

updown <- bind_rows(
  updown_baseline(n79_eq, "NLSY79", "Equivalized (headline)"),
  updown_baseline(n79_un, "NLSY79", "Unadjusted (household)"),
  updown_baseline(n97_eq, "NLSY97", "Equivalized (headline)"),
  updown_baseline(n97_un, "NLSY97", "Unadjusted (household)"))

saveRDS(updown, "Data/NLSY_estimation/intermediate/unadjusted_updown.rds")
print(updown %>% filter(cohort == "Full") %>%
        mutate(across(c(up, down), ~ round(., 3))) %>% as.data.frame())


## ===================== [precompute_family_size] =====================
# ============================================================================
# NLSY/_precompute_family_size.R
# ----------------------------------------------------------------------------
# Household-size table behind the family-size sensitivity: NLSY79 responds more
# to √-size equivalization than NLSY97 because its PARENTAL households are larger
# (the 1961-64 cohorts grew up in bigger, baby-boom-era families), so dividing by
# √(size) gives the child a larger per-capita lift relative to the parent.
# Reports median parent/child household size, the equivalization boost
# √(parent/child), and the child/parent income ratio raw vs equivalized, per
# survey. All from the production data objects (no hard-coded values).
# Saves Data/NLSY_estimation/intermediate/family_size_table.rds.
# ============================================================================
suppressPackageStartupMessages(library(dplyr))

med_size <- function(long_f, id, cohs) {
  l <- readRDS(long_f) %>% mutate(id = as.character(.data[[id]]))
  p <- l %>% filter(is_parent_obs) %>% group_by(id) %>%
    summarise(fs = mean(family_size, na.rm = TRUE), by = first(birth_year)) %>% filter(by %in% cohs)
  a <- l %>% filter(is_adult_obs) %>% group_by(id) %>%
    summarise(fs = mean(family_size, na.rm = TRUE), by = first(birth_year)) %>% filter(by %in% cohs)
  c(parent = median(p$fs, na.rm = TRUE), child = median(a$fs, na.rm = TRUE))
}
ratios <- function(final_f, par_eq, chd_eq, w, filt) {
  f <- readRDS(final_f); f <- f[f[[filt]] %in% TRUE, ]
  wm <- function(x) weighted.mean(f[[x]], f[[w]], na.rm = TRUE)
  c(raw = wm(sub("equiv", "real", chd_eq)) / wm(sub("equiv", "real", par_eq)),
    eq  = wm(chd_eq) / wm(par_eq))
}

s79 <- med_size("Data/NLSY_estimation/intermediate/nlsy79_income_long.rds", "caseid", 1961:1964)
s97 <- med_size("Data/NLSY_estimation/intermediate/nlsy97_income_long.rds", "pubid", 1980:1984)
r79 <- ratios("Data/NLSY_estimation/dyads/nlsy79_income_final.rds",
              "mean_parent_income_equiv", "mean_adult_income_equiv",
              "custom_weight", "nlsy79_zhou_income_sample")
r97 <- ratios("Data/NLSY_estimation/dyads/nlsy97_income_final.rds",
              "parent_income_equiv", "mean_adult_income_equiv",
              "custom_weight", "has_custom_weight")

family_size_table <- tibble::tibble(
  Survey = c("NLSY79 (1961–64)", "NLSY97 (1980–84)"),
  `Parent household` = c(s79["parent"], s97["parent"]),
  `Child household`  = c(s79["child"],  s97["child"]),
  `Equiv. boost √(P/C)` = sqrt(c(s79["parent"] / s79["child"], s97["parent"] / s97["child"])),
  `Child/parent ratio, raw`         = c(r79["raw"], r97["raw"]),
  `Child/parent ratio, equivalized` = c(r79["eq"],  r97["eq"]))
saveRDS(family_size_table, "Data/NLSY_estimation/intermediate/family_size_table.rds")
print(as.data.frame(family_size_table), digits = 3)


## ===================== [precompute_income_windows] =====================
# ============================================================================
# NLSY/_precompute_income_windows.R
# ----------------------------------------------------------------------------
# Calendar-year capture windows behind the NLSY79 1961-64 upward-mobility rise.
# Parent income is drawn from a FIXED 1978-82 window (calendar-based), but later-
# born children have their adult income captured in later calendar years, and
# their parents are captured younger. Both raise the child/parent income ratio
# against the fixed pooled cut-points, producing the rising upward mobility.
# NLSY79 is biennial after 1994 (income-years become odd-only: 1995,1997,1999…),
# so an even-born cohort's age-38 income-year (even) is never collected and the
# window caps one year early -- documented in the `Child income yrs` column.
# All values from production objects. Saves Data/NLSY_estimation/intermediate/income_windows.rds.
# ============================================================================
suppressPackageStartupMessages(library(dplyr))

l <- readRDS("Data/NLSY_estimation/intermediate/nlsy79_income_long.rds")
win <- function(flag) l %>% filter(.data[[flag]]) %>% group_by(birth_year) %>%
  summarise(yrs = paste0(min(income_year), "–", max(income_year)), .groups = "drop") %>%
  filter(birth_year %in% 1961:1964) %>% arrange(birth_year)
pw <- win("is_parent_obs"); cw <- win("is_adult_obs")

f <- readRDS("Data/NLSY_estimation/dyads/nlsy79_income_final.rds") %>% filter(nlsy79_zhou_income_sample)
lev <- f %>% group_by(birth_year) %>% summarise(
  par = round(weighted.mean(mean_parent_income_equiv, custom_weight, na.rm = TRUE)),
  chd = round(weighted.mean(mean_adult_income_equiv,  custom_weight, na.rm = TRUE)),
  .groups = "drop") %>% arrange(birth_year)

vl <- readRDS("Data/NLSY_estimation/processed/nlsy79_income_processed.rds")$validation$validation_lst
up <- sapply(as.character(1961:1964), function(c) {
  d <- vl[[c]]$movement_df %>% filter(category == "Historical Upward Mobility")
  round(as.numeric(d$mean_bc[1]), 3) })

income_windows <- tibble::tibble(
  Cohort = 1961:1964,
  `Parent income yrs` = pw$yrs,
  `Child income yrs`  = cw$yrs,
  `Parent real $` = lev$par,
  `Child real $`  = lev$chd,
  `Child/parent`  = round(lev$chd / lev$par, 2),
  `Upward mobility` = up)
saveRDS(income_windows, "Data/NLSY_estimation/intermediate/income_windows.rds")
print(as.data.frame(income_windows), row.names = FALSE)


## ===================== [precompute_pooled_cutpoints] =====================
# ============================================================================
# NLSY/_precompute_pooled_cutpoints.R
# ----------------------------------------------------------------------------
# Baseline pooled weighted quintile cut-points (equivalized real 2023$) -- the
# income-class definition reported in the paper. For each survey, pool parent +
# child equivalized real income over the reporting domain (NLSY79 1961-64,
# NLSY97 1980-84) and take the four weighted quintile boundaries Q20/Q40/Q60/Q80.
# These are the SAME cut-points the baseline estimator computes inside tm_series
# for the pooled "Full" transition matrix (pooled_weighted_cutpoints(), no
# resampling), so the table matches the reported estimates exactly. Class 5 =
# below Q20 (lowest); Class 1 = above Q80 (highest).
# Saves Data/NLSY_estimation/intermediate/pooled_cutpoints.csv (+ .rds).
# ============================================================================
source("NLSY_fun/nlsy_estimation_functions.R")   # pooled_weighted_cutpoints()

cut_row <- function(f, cohs, label) {
  d  <- readRDS(f) %>%
    dplyr::filter(cohort %in% cohs, !is.na(parent), !is.na(child), !is.na(weight))
  cp <- pooled_weighted_cutpoints(d$parent, d$child, d$weight, 5L)
  tibble::tibble(Survey = label, Q20 = cp[1], Q40 = cp[2], Q60 = cp[3], Q80 = cp[4])
}
pooled_cutpoints <- dplyr::bind_rows(
  cut_row("Data/NLSY_estimation/inputs/nlsy79_income_input.rds", 1961:1964, "NLSY79 (1961-64)"),
  cut_row("Data/NLSY_estimation/inputs/nlsy97_income_input.rds", 1980:1984, "NLSY97 (1980-84)"))
saveRDS(pooled_cutpoints, "Data/NLSY_estimation/intermediate/pooled_cutpoints.rds")
readr::write_csv(pooled_cutpoints, "Data/NLSY_estimation/intermediate/pooled_cutpoints.csv")

cat("\nBaseline pooled weighted quintile cut-points (equivalized real 2023$).\n")
cat("Class 5 = below Q20; Class 1 = above Q80.\n")
usd <- function(x) paste0("$", formatC(round(x), format = "d", big.mark = ","))
print(data.frame(Survey = pooled_cutpoints$Survey,
                 Q20 = usd(pooled_cutpoints$Q20), Q40 = usd(pooled_cutpoints$Q40),
                 Q60 = usd(pooled_cutpoints$Q60), Q80 = usd(pooled_cutpoints$Q80)),
      row.names = FALSE)
