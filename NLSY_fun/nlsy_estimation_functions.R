# ============================================================================
# NLSY_fun/nlsy_estimation_functions.R
# ----------------------------------------------------------------------------
# Estimation-side machinery, consolidated:
#   1. weighted pooled cut-points, tile assignment, weighted transition
#      matrices                          (was nlsy_transition_helpers.R)
#   2. tm_series + the three resampling schemes (individual / household /
#      Fay-BRR survey), baseline_measures, and the analyst-estimator loader
#                                        (was nlsy_measure_bootstrap.R)
# The estimator itself is UNCHANGED in synthetic_dynasty_estimation.R.
# ============================================================================

## ===================== [nlsy_transition_helpers] =====================
# ============================================================================
# R/nlsy_transition_helpers.R
# ----------------------------------------------------------------------------
# Weighted parent->child transition matrices for the NLSY mobility analysis.
#
# Design (agreed with the analyst):
#   * Survey weights are the returned custom weights; the whole NLSY79 / NLSY97
#     cohort is the design, and each birth cohort is treated as a DOMAIN (subset).
#   * INCOME classes: pool parent AND child (equivalized real) income into ONE
#     weighted distribution over the whole analysis group, take 5 weighted-quantile
#     cut-points, and apply those SAME cut-points to classify parent and child
#     (project scheme D7). Cut-points are computed once per analysis group, then
#     reused across cohort domains so cohorts are on a common income scale.
#   * OCCUPATION classes are the 5-class EGP already built (no cut-points).
#
# Direction conventions (documented, not enforced):
#   income tile 1 = LOWEST .. 5 = HIGHEST (standard quintiles)
#   EGP class  1 = HIGHEST .. 5 = LOWEST  (Morgan 2017)
# The transition matrix P is row-stochastic: P[i,j] = weighted Pr(child = j | parent = i).
# ============================================================================

# ---- weighted quantiles -----------------------------------------------------
weighted_quantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

# pooled_weighted_cutpoints(): weighted quantile cut-points of the POOLED
# parent+child distribution (the income-class definition). Returns n_tiles-1 cuts.
pooled_weighted_cutpoints <- function(parent, child, weight, n_tiles = 5L) {
  probs <- seq_len(n_tiles - 1L) / n_tiles
  weighted_quantile(c(parent, child), c(weight, weight), probs)
}

# assign_tile(): bucket values into 1..n_tiles using cut-points (n_tiles-1 of them).
# tile 1 = below the first cut (lowest) .. tile n = at/above the last cut (highest).
assign_tile <- function(x, cutpoints) {
  cls <- findInterval(x, cutpoints) + 1L
  cls[is.na(x)] <- NA_integer_
  as.integer(cls)
}

# ---- the transition matrix core --------------------------------------------
# weighted_transition_matrix(): row-stochastic P from parent/child classes + weights.
weighted_transition_matrix <- function(parent_class, child_class, weight, n = 5L) {
  ok <- !is.na(parent_class) & !is.na(child_class) & !is.na(weight) & weight > 0
  pc <- factor(parent_class[ok], levels = seq_len(n))
  cc <- factor(child_class[ok],  levels = seq_len(n))
  counts <- tapply(weight[ok], list(parent = pc, child = cc), sum)
  counts[is.na(counts)] <- 0
  raw    <- table(parent = pc, child = cc)              # unweighted cell counts
  trans  <- sweep(counts, 1, rowSums(counts), "/")
  list(transition   = trans,                            # P(child|parent), row-stochastic
       counts       = counts,                           # weighted cell totals
       raw_counts   = raw,                              # unweighted n per cell
       parent_margin = rowSums(counts) / sum(counts),
       child_margin  = colSums(counts) / sum(counts),
       n            = sum(ok),
       n_weighted   = sum(weight[ok]))
}

# ---- income / occupation wrappers ------------------------------------------
income_transition_matrix <- function(data, parent_col, child_col, weight_col,
                                     cutpoints, n = 5L) {
  weighted_transition_matrix(
    assign_tile(data[[parent_col]], cutpoints),
    assign_tile(data[[child_col]],  cutpoints),
    data[[weight_col]], n)
}

occ_transition_matrix <- function(data, parent_col, child_col, weight_col, n = 5L) {
  weighted_transition_matrix(data[[parent_col]], data[[child_col]], data[[weight_col]], n)
}

# ---- cohort-domain driver ---------------------------------------------------
# cohort_transition_matrices(): the full-sample matrix plus one per cohort domain,
# all sharing the SAME income cut-points (computed once on the whole `data`).
cohort_transition_matrices <- function(data, parent_col, child_col, weight_col,
                                       cohort_col, type = c("income", "occ"),
                                       n = 5L, cutpoints = NULL) {
  type <- match.arg(type)
  if (type == "income" && is.null(cutpoints)) {
    cutpoints <- pooled_weighted_cutpoints(data[[parent_col]], data[[child_col]],
                                           data[[weight_col]], n)
  }
  tm_one <- function(d) if (type == "income")
    income_transition_matrix(d, parent_col, child_col, weight_col, cutpoints, n)
  else occ_transition_matrix(d, parent_col, child_col, weight_col, n)

  out <- list(Full = tm_one(data))
  for (co in sort(unique(data[[cohort_col]])))
    out[[as.character(co)]] <- tm_one(data[data[[cohort_col]] == co, , drop = FALSE])
  attr(out, "cutpoints") <- cutpoints
  out
}

# ---- small validity + summary utilities ------------------------------------
# check every row sums to 1 (within tolerance) and no empty parent rows.
validate_transition <- function(tm, tol = 1e-8) {
  rs <- rowSums(tm$transition)
  ok_rows <- all(abs(rs[rs > 0] - 1) < tol)
  list(rows_sum_to_1 = ok_rows,
       empty_parent_rows = sum(rowSums(tm$counts) == 0),
       n = tm$n)
}

# common scalar mobility summaries from a row-stochastic P
mobility_summaries <- function(P) {
  n <- nrow(P)
  immobility <- mean(diag(P), na.rm = TRUE)          # mean diagonal (staying)
  # mean absolute rank change of the child's class given parent (needs marginals);
  # here a simple trace-based immobility + off-diagonal mass:
  c(mean_diagonal = round(immobility, 3),
    off_diagonal  = round(1 - immobility, 3))
}

# pretty-print a transition matrix as a percentage table
format_transition <- function(tm, digits = 1) round(tm$transition * 100, digits)

## ===================== [nlsy_measure_bootstrap] =====================
# ============================================================================
# R/nlsy_measure_bootstrap.R
# ----------------------------------------------------------------------------
# Bootstrap adapter for the analyst's mobility measures. It does NOT touch the
# estimators: it only RESAMPLES the NLSY data and CONSTRUCTS weighted transition
# matrices in the exact object shape the analyst's `single_simulation()` expects
#   tm = list(num = cbind(M, validation = rowSums(M)),   # weighted counts + row totals
#             perc = M / rowSums(M))                      # row-stochastic P
# so we can `lapply(tm_series, single_simulation, size, tol, weighted = TRUE)`.
#
# Three resamplers differ only in HOW the data is resampled; each returns a
# weighted transition-matrix SERIES (one per cohort domain, plus "Full"):
#   replicate_individual()  -- individual bootstrap (rows with replacement)
#   replicate_household()   -- NLSY79 weighted household-cluster (siblings together)
#   replicate_survey()      -- NLSY97 Fay-BRR (swap in a replicate weight column)
# For INCOME, the pooled parent+child weighted quintile cut-points are RECOMPUTED
# inside each replicate (never frozen); for OCCUPATION the 5-class EGP is fixed.
#
# Requires: source the analyst's estimators (see load_analyst_estimators()).
# Requires helpers: pooled_weighted_cutpoints(), assign_tile() from
#   R/nlsy_transition_helpers.R.
# ============================================================================

# load_analyst_estimators(): source the analyst's UNCHANGED functions + packages.
load_analyst_estimators <- function(
    root = ".") {
  for (p in c("dplyr","tidyverse","markovchain","logmult","assertthat","pracma","ipfr","combinat"))
    suppressPackageStartupMessages(require(p, character.only = TRUE))
  suppressWarnings(suppressMessages({
    source(file.path(root, "Functions", "utils.R"))
    source(file.path(root, "Functions", "CopulaFunctions1.R"))
    source(file.path(root, "NLSY_fun", "synthetic_dynasty_estimation.R"))   # single_simulation (matrix, weighted)
  }))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 1. Build the weighted transition-matrix object in the analyst's shape
# ---------------------------------------------------------------------------
# d          : data with integer class columns pc_col (parent 1..size), cc_col (child)
# descending : reverse to "Class 1 = highest" (income tiles are 1=lowest; EGP already 1=highest)
# rescale    : scale weighted counts to the dyad count nrow(d) (matches the report's income build;
#              leaves P unchanged). Occupation uses raw weighted counts (rescale = FALSE).
build_weighted_tm <- function(d, pc_col, cc_col, weight_col, size = 5L,
                              descending = FALSE, rescale = FALSE) {
  f <- stats::as.formula(sprintf("%s ~ factor(%s, 1:%d) + factor(%s, 1:%d)",
                                 weight_col, pc_col, size, cc_col, size))
  M <- matrix(as.numeric(as.matrix(stats::xtabs(f, data = d))), size, size)
  if (descending) M <- M[size:1, size:1]                 # Class 1 = highest
  if (rescale && sum(M) > 0) M <- M * nrow(d) / sum(M)
  dimnames(M) <- list(paste0("Class", 1:size), paste0("Class", 1:size))
  list(num = cbind(M, validation = rowSums(M)), perc = M / rowSums(M))
}

# tm_series(): classify (income) or use EGP (occ), then build one weighted tm per
# cohort domain + a pooled "Full". Income cut-points are recomputed here from the
# (resampled) data and shared across the cohort domains.
#
# DOMAIN RESTRICTION. The analysis is restricted to the reporting cohorts
# (`cohorts`, e.g. 1961-64) BEFORE the income cut-points and the pooled "Full"
# are computed. The NLSY79 inputs carry a wider birth-cohort span (1957-64) so
# that the returned custom weights are stable at the 1960 boundary; those extra
# cohorts are weighting scaffolding only and must not enter the class cut-points
# or the pooled statistic. Because resampling (household / survey / individual)
# is done UPSTREAM on the full input, this is standard domain estimation:
# resample the whole design, then analyse the domain. Set `domain = FALSE` to
# recover the legacy behaviour (cut-points / Full on the whole input).
tm_series <- function(d, parent_col, child_col, weight_col, cohort_col, cohorts,
                      type = c("income", "occ"), size = 5L, cutpoints = NULL,
                      domain = TRUE) {
  type <- match.arg(type)
  if (isTRUE(domain))
    d <- d[d[[cohort_col]] %in% cohorts, , drop = FALSE]  # reporting domain only
  if (type == "income") {
    if (is.null(cutpoints))
      cutpoints <- pooled_weighted_cutpoints(d[[parent_col]], d[[child_col]], d[[weight_col]], size)
    d$.pc <- assign_tile(d[[parent_col]], cutpoints)
    d$.cc <- assign_tile(d[[child_col]], cutpoints)
    descending <- TRUE; rescale <- TRUE                  # tiles 1=lowest -> reverse; rescale counts
  } else {
    d$.pc <- as.integer(d[[parent_col]]); d$.cc <- as.integer(d[[child_col]])
    descending <- FALSE; rescale <- FALSE                # EGP already 1=highest
  }
  d <- d[!is.na(d$.pc) & !is.na(d$.cc) & !is.na(d[[weight_col]]), , drop = FALSE]
  bwt <- function(dd) build_weighted_tm(dd, ".pc", ".cc", weight_col, size, descending, rescale)
  out <- list(Full = bwt(d))
  for (co in cohorts) out[[as.character(co)]] <- bwt(d[d[[cohort_col]] == co, , drop = FALSE])
  out
}

# ---------------------------------------------------------------------------
# 2. The three resamplers -> each returns a weighted tm SERIES for ONE replicate
# ---------------------------------------------------------------------------

# individual bootstrap: respondents sampled with replacement (uniform, NOT weight-
# proportional); custom weight carried unchanged.
replicate_individual <- function(data, parent_col, child_col, weight_col, cohort_col,
                                 cohorts, type = "income", size = 5L) {
  d <- data[sample(nrow(data), replace = TRUE), , drop = FALSE]
  tm_series(d, parent_col, child_col, weight_col, cohort_col, cohorts, type, size)
}

# weighted household-cluster bootstrap (NLSY79): resample HHIDs with replacement,
# carry all sampled siblings together (preserves sibling dependence).
replicate_household <- function(data, hhid_col, parent_col, child_col, weight_col,
                                cohort_col, cohorts, type = "income", size = 5L) {
  hh_idx <- split(seq_len(nrow(data)), data[[hhid_col]])
  drawn <- sample(names(hh_idx), length(hh_idx), replace = TRUE)
  d <- data[unlist(hh_idx[drawn], use.names = FALSE), , drop = FALSE]
  tm_series(d, parent_col, child_col, weight_col, cohort_col, cohorts, type, size)
}

# NLSY97 Fay-BRR replicate: the "resample" is a reweighting -- swap the custom
# weight for the b-th Fay replicate weight (from as.svrepdesign). Data rows are
# fixed; only the weight column changes.
replicate_survey <- function(data, repweight_vec, parent_col, child_col,
                             cohort_col, cohorts, type = "income", size = 5L) {
  d <- data
  d$.repw <- repweight_vec
  tm_series(d, parent_col, child_col, ".repw", cohort_col, cohorts, type, size)
}

# ---------------------------------------------------------------------------
# 3. ONE bootstrap replicate -> named list (cohort -> single_simulation result).
#    Mirrors the analyst's single_simulation_lst() role (resample + build weighted
#    tm series + apply single_simulation). For scheme = "survey", pass the i-th Fay
#    replicate-weight column via `repweight_vec`.
# ---------------------------------------------------------------------------
single_bootstrap_replicate <- function(data, scheme = c("individual", "household", "survey"),
                                       parent_col, child_col, weight_col, cohort_col, cohorts,
                                       type = "income", size = 5L, tol = 0.01,
                                       hhid_col = NULL, repweight_vec = NULL) {
  scheme <- match.arg(scheme)
  tms <- switch(scheme,
    individual = replicate_individual(data, parent_col, child_col, weight_col, cohort_col, cohorts, type, size),
    household  = replicate_household(data, hhid_col, parent_col, child_col, weight_col, cohort_col, cohorts, type, size),
    survey     = replicate_survey(data, repweight_vec, parent_col, child_col, cohort_col, cohorts, type, size))
  lapply(tms, function(tm) tryCatch(single_simulation(tm, size = size, tol = tol, weighted = TRUE),
                                    error = function(e) NA))
}

# Convenience driver used by the notebook (mclapply over single_bootstrap_replicate).
bootstrap_measures <- function(data, scheme = c("individual", "household", "survey"),
                               parent_col, child_col, weight_col, cohort_col, cohorts,
                               type = "income", size = 5L, tol = 0.01, B = 500L,
                               hhid_col = NULL, repweights = NULL, seed = 20260722L,
                               mc_cores = 1L) {
  scheme <- match.arg(scheme)
  one <- function(b) {
    set.seed(seed + b)
    single_bootstrap_replicate(data, scheme, parent_col, child_col, weight_col, cohort_col,
                               cohorts, type, size, tol, hhid_col,
                               repweight_vec = if (scheme == "survey") repweights[, b] else NULL)
  }
  if (mc_cores > 1L && requireNamespace("parallel", quietly = TRUE))
    parallel::mclapply(seq_len(B), one, mc.cores = mc_cores)
  else lapply(seq_len(B), one)
}

# baseline_measures(): the POINT estimate -- build the weighted tm series from the
# (un-resampled) data and apply single_simulation. Mirrors calculate_baseline().
baseline_measures <- function(data, parent_col, child_col, weight_col, cohort_col,
                              cohorts, type = "income", size = 5L, tol = 0.01) {
  tms <- tm_series(data, parent_col, child_col, weight_col, cohort_col, cohorts, type, size)
  lapply(tms, function(tm) tryCatch(single_simulation(tm, size = size, tol = tol, weighted = TRUE),
                                    error = function(e) NA))
}

# boot_ci(): extract a scalar measure (via `accessor`) for a cohort across all
# replicates and summarise with a percentile CI. `accessor` maps a single_simulation
# result list -> a scalar (e.g. function(r) r$AMTE).
#
# fay_rho: set to the Fay coefficient (0.5 here) for SURVEY (Fay-BRR) replicates.
#   Fay replicate weights are perturbed by only +/-(1 - rho), so raw replicate
#   deviations understate the design variance by (1 - rho)^2. The fix rescales
#   deviations around the full-sample estimate (mse = TRUE convention):
#     v* = est + (v - est) / (1 - rho)
#   so that sd(v*) matches the Fay variance  Var = sum((v - est)^2) / (R (1-rho)^2)
#   and the percentile CI is scaled consistently. Leave NULL for the individual
#   and household bootstrap (raw replicate sd is already correct there).
boot_ci <- function(boot, accessor, cohort = "Full", point = NULL, level = 0.95,
                    fay_rho = NULL) {
  v <- vapply(boot, function(rep) {
    r <- rep[[as.character(cohort)]]
    if (is.list(r)) suppressWarnings(as.numeric(accessor(r))) else NA_real_
  }, numeric(1))
  est <- if (is.null(point)) NA_real_ else
    suppressWarnings(as.numeric(accessor(point[[as.character(cohort)]])))
  if (!is.null(fay_rho) && !is.na(est)) v <- est + (v - est) / (1 - fay_rho)
  a <- (1 - level) / 2
  tibble::tibble(
    cohort   = as.character(cohort),
    estimate = est,
    boot_mean = round(mean(v, na.rm = TRUE), 4),
    se       = round(stats::sd(v, na.rm = TRUE), 4),
    lower    = round(stats::quantile(v, a, na.rm = TRUE), 4),
    upper    = round(stats::quantile(v, 1 - a, na.rm = TRUE), 4),
    n_valid  = sum(!is.na(v))
  )
}
