# ============================================================================
# NLSY_fun/nlsy_data_functions.R
# ----------------------------------------------------------------------------
# ALL data-construction functions for the NLSY pipeline, consolidated:
#   1. raw-extract reading, ID standardization, audited merges, custom-weight
#      ID exports                       (was nlsy_data_helpers.R)
#   2. income maps + builders: TNFI/CV_INCOME reshaping, generation tagging,
#      lived-at-home screen, NLSY97 parent-income rules, PCE deflation +
#      sqrt-size equivalization         (was nlsy_income_helpers.R)
#   3. occupation EGP crosswalk chain   (was nlsy_occupation_helpers.R)
#   4. harmonized education             (was nlsy_education_helpers.R)
#   5. bootstrap/survey-design input builders (was nlsy_bootstrap_helpers.R)
# ============================================================================

## ===================== [nlsy_data_helpers] =====================
# ============================================================================
# R/nlsy_data_helpers.R
# ----------------------------------------------------------------------------
# Generic, reusable helpers for the NLSY79 / NLSY97 validation data build.
#
# These functions implement the audited-merge, ID-integrity, missing-code, and
# custom-weight-export utilities required by the pipeline specification
# (see nlsy_data_pipeline_prompt_regenerated.md, Section 8).
#
# Design principles:
#   * Every respondent-wide file must have a unique, non-missing character ID.
#   * NLS negative "codes" (-1 refusal, -2 don't know, -3 invalid skip,
#     -4 valid skip / not-in-universe, -5 non-interview) are NEVER collapsed to
#     a single NA before they have been audited. nls_missing_reason() preserves
#     the distinction; nls_clean_numeric() only nulls the negatives once the
#     reason has been captured.
#   * HHID (NLSY79) / SIDCODE (NLSY97) are household / sibling-cluster
#     identifiers and are NEVER used as a respondent merge key.
#
# Nothing here reads or mutates a raw NLS file in place; callers pass data in.
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Small internal utilities
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Ensure an output directory exists (used before every write_*).
.ensure_dir <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# ---------------------------------------------------------------------------
# 1. NLS missing-code handling
# ---------------------------------------------------------------------------

# Map a raw NLS numeric value to a human-readable missing reason.
# Positive values (and 0) are "valid"; each negative code has a distinct
# meaning that must be preserved for the diagnostics tables (spec Section 17.5).
nls_missing_reason <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    is.na(x)   ~ "na_source",
    x == -1    ~ "refusal",
    x == -2    ~ "dont_know",
    x == -3    ~ "invalid_skip",
    x == -4    ~ "valid_skip_not_in_universe",
    x == -5    ~ "noninterview",
    x  < 0     ~ "other_negative_code",
    TRUE       ~ "valid"
  )
}

# Convert a raw NLS numeric field to a clean numeric: every negative code
# becomes NA, positives/zero pass through. Always call nls_missing_reason()
# FIRST if you need to keep the reason (this function discards it by design).
nls_clean_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(x) & x < 0, NA_real_, x)
}

# ---------------------------------------------------------------------------
# 2. Respondent-ID integrity
# ---------------------------------------------------------------------------

# assert_unique_id(): guarantee `id` exists, is non-missing, and is unique.
# On failure, dump the offending rows to output/audit/duplicates_<label>.csv
# and stop() -- a duplicated respondent key would silently multiply rows on the
# next join (spec Section 8 / 17 stop-condition).
assert_unique_id <- function(data, id, label,
                             audit_dir = file.path("Data", "NLSY_estimation", "audit")) {
  stopifnot(is.data.frame(data))
  if (!id %in% names(data)) {
    stop(sprintf("[assert_unique_id:%s] id column '%s' not found.", label, id))
  }
  idv <- data[[id]]
  if (any(is.na(idv) | idv == "")) {
    stop(sprintf("[assert_unique_id:%s] '%s' has %d missing values.",
                 label, id, sum(is.na(idv) | idv == "")))
  }
  dup <- idv[duplicated(idv)]
  if (length(dup) > 0L) {
    out <- file.path(audit_dir, sprintf("duplicates_%s.csv", label))
    .ensure_dir(out)
    readr::write_csv(dplyr::filter(data, .data[[id]] %in% unique(dup)), out)
    stop(sprintf("[assert_unique_id:%s] %d duplicate '%s' values. See %s",
                 label, length(unique(dup)), id, out))
  }
  invisible(TRUE)
}

# standardize_id(): create a canonical character respondent key `new_id` while
# preserving the raw ID column. Character keys avoid float/precision surprises
# and make joins deterministic. Refuses to use household/cluster IDs as keys.
standardize_id <- function(data, raw_id, new_id) {
  forbidden <- c("HHID", "SIDCODE", "R0000149", "R1193000", "hhid", "sidcode")
  if (raw_id %in% forbidden) {
    stop(sprintf("[standardize_id] '%s' is a household/cluster ID; never use it as a respondent key.", raw_id))
  }
  if (!raw_id %in% names(data)) {
    stop(sprintf("[standardize_id] raw id '%s' not found.", raw_id))
  }
  data[[new_id]] <- as.character(data[[raw_id]])
  data
}

# ---------------------------------------------------------------------------
# 3. Audited left join
# ---------------------------------------------------------------------------

# left_join_with_audit(): join augmentation columns onto an RA base file with a
# full audit trail. Guarantees no row multiplication and no unique-ID growth,
# resolves overlapping non-key columns explicitly, and writes match-rate /
# unmatched-ID reports to output/audit (spec Section 8 / 10).
#
#   base    : RA file (left-hand side, authoritative respondent universe)
#   augment : augmentation file
#   by      : named or unnamed join key (respondent ID; NEVER a household ID)
#   label   : used for audit file names
#   overlap : how to resolve non-key columns present in BOTH inputs:
#             "keep_base" (default) or "prefer_augment"
left_join_with_audit <- function(base, augment, by, label,
                                 overlap = c("keep_base", "prefer_augment"),
                                 audit_dir = file.path("Data", "NLSY_estimation", "audit")) {
  overlap <- match.arg(overlap)
  key_base <- if (is.null(names(by))) by else names(by)
  key_aug  <- unname(by)

  # (a) Respondent-wide inputs must have unique keys before the join.
  assert_unique_id(base,    key_base, paste0(label, "_base"),    audit_dir)
  assert_unique_id(augment, key_aug,  paste0(label, "_augment"), audit_dir)

  n_base_rows <- nrow(base)
  n_base_ids  <- dplyr::n_distinct(base[[key_base]])

  # (b) Identify overlapping non-key columns and resolve them explicitly.
  overlap_cols <- setdiff(intersect(names(base), names(augment)), key_base)
  overlap_log <- tibble::tibble(column = character(), resolution = character())
  if (length(overlap_cols) > 0L) {
    if (overlap == "keep_base") {
      augment <- dplyr::select(augment, -dplyr::all_of(overlap_cols))
      overlap_log <- tibble::tibble(column = overlap_cols, resolution = "kept_base")
    } else {
      base <- dplyr::select(base, -dplyr::all_of(overlap_cols))
      overlap_log <- tibble::tibble(column = overlap_cols, resolution = "preferred_augment")
    }
  }

  # (c) The join itself.
  joined <- dplyr::left_join(base, augment, by = by)

  # (d) Hard invariants: rows and unique base IDs must not increase.
  if (nrow(joined) != n_base_rows) {
    stop(sprintf("[left_join_with_audit:%s] row multiplication: %d -> %d.",
                 label, n_base_rows, nrow(joined)))
  }
  if (dplyr::n_distinct(joined[[key_base]]) != n_base_ids) {
    stop(sprintf("[left_join_with_audit:%s] unique base IDs changed: %d -> %d.",
                 label, n_base_ids, dplyr::n_distinct(joined[[key_base]])))
  }

  # (e) Match-rate / unmatched-ID diagnostics.
  base_ids <- unique(as.character(base[[key_base]]))
  aug_ids  <- unique(as.character(augment[[key_aug]]))
  matched  <- intersect(base_ids, aug_ids)
  audit <- tibble::tibble(
    label            = label,
    n_base_rows      = n_base_rows,
    n_augment_rows   = nrow(augment),
    n_base_ids       = length(base_ids),
    n_augment_ids    = length(aug_ids),
    n_matched        = length(matched),
    match_rate_base  = round(length(matched) / length(base_ids), 4),
    n_base_unmatched = length(setdiff(base_ids, aug_ids)),
    n_augment_only   = length(setdiff(aug_ids, base_ids))
  )
  .ensure_dir(file.path(audit_dir, "x"))
  readr::write_csv(audit, file.path(audit_dir, sprintf("merge_audit_%s.csv", label)))
  readr::write_csv(tibble::tibble(base_unmatched_id = setdiff(base_ids, aug_ids)),
                   file.path(audit_dir, sprintf("unmatched_base_%s.csv", label)))
  readr::write_csv(tibble::tibble(augment_only_id = setdiff(aug_ids, base_ids)),
                   file.path(audit_dir, sprintf("augment_only_%s.csv", label)))
  if (nrow(overlap_log) > 0L) {
    readr::write_csv(overlap_log, file.path(audit_dir, sprintf("overlap_resolution_%s.csv", label)))
  }

  attr(joined, "merge_audit") <- audit
  joined
}

# ---------------------------------------------------------------------------
# 4. Parental education (kept as TWO separate variables, never collapsed)
# ---------------------------------------------------------------------------

# prepare_parent_education(): clean mother's and father's highest grade as two
# independent variables, each with its own observed-value flag and missing
# reason. The spec forbids any max-across-parents / "highest-parent" collapse
# in the primary demographic module (Section 8 / 11 / 19).
prepare_parent_education <- function(data, mother_raw, father_raw) {
  data %>%
    dplyr::mutate(
      mother_education_missing_reason = nls_missing_reason(.data[[mother_raw]]),
      father_education_missing_reason = nls_missing_reason(.data[[father_raw]]),
      mother_highest_education        = nls_clean_numeric(.data[[mother_raw]]),
      father_highest_education        = nls_clean_numeric(.data[[father_raw]]),
      mother_education_observed       = mother_education_missing_reason == "valid",
      father_education_observed       = father_education_missing_reason == "valid"
    )
}

# ---------------------------------------------------------------------------
# 5. Weighted rank / cutpoint / class helpers (used inside every bootstrap)
# ---------------------------------------------------------------------------

# weighted_midrank(): Zhou-style weighted mid-ranks with deterministic tie
# handling. Returns ranks on 0-1 (default) or 1-100 scale. Ties receive the
# average cumulative weight of the tie block (a "mid" rank).
weighted_midrank <- function(x, w = NULL, scale = c("unit", "percent")) {
  scale <- match.arg(scale)
  n <- length(x)
  if (is.null(w)) w <- rep(1, n)
  ok <- !is.na(x) & !is.na(w) & w > 0
  r <- rep(NA_real_, n)
  if (!any(ok)) return(r)
  xo <- x[ok]; wo <- w[ok]
  ord <- order(xo)
  xs <- xo[ord]; ws <- wo[ord]
  csum <- cumsum(ws)
  below <- csum - ws                 # total weight strictly before each element
  mid <- below + ws / 2              # mid-weight position within its own mass
  # Assign every tied value the mean mid-position of its tie block.
  mid_by_tie <- ave(mid, xs, FUN = mean)
  W <- sum(ws)
  rr <- mid_by_tie / W               # 0-1 weighted mid-rank
  out_ok <- numeric(length(xo))
  out_ok[ord] <- rr
  r[ok] <- if (scale == "percent") out_ok * 100 else out_ok
  r
}

# weighted_cutpoints(): weighted quantile cut-points for `probs`, with
# diagnostics for ties / zero-width classes (which would collapse a quintile).
weighted_cutpoints <- function(x, w = NULL, probs = seq(0.2, 0.8, by = 0.2)) {
  ok <- !is.na(x) & (is.null(w) | !is.na(w))
  x2 <- x[ok]; w2 <- if (is.null(w)) rep(1, length(x2)) else w[ok]
  ord <- order(x2); x2 <- x2[ord]; w2 <- w2[ord]
  cw <- cumsum(w2) / sum(w2)
  cuts <- vapply(probs, function(p) x2[which(cw >= p)[1]], numeric(1))
  zero_width <- any(diff(c(min(x2), cuts, max(x2))) <= 0)
  attr(cuts, "zero_width_class") <- zero_width
  attr(cuts, "n_ties") <- sum(duplicated(cuts))
  cuts
}

# assign_income_class(): bucket a value into classes defined by cutpoints, with
# explicit boundary behaviour (value <= cut_k -> class k). Missing input stays
# missing and is kept distinct from a valid zero / negative income.
assign_income_class <- function(x, cutpoints) {
  brks <- c(-Inf, cutpoints, Inf)
  cls <- cut(x, breaks = brks, labels = FALSE, right = TRUE, include.lowest = TRUE)
  cls[is.na(x)] <- NA_integer_
  as.integer(cls)
}

# ---------------------------------------------------------------------------
# 6. Custom-weight ID export (THE key deliverable for the NLS system)
# ---------------------------------------------------------------------------

# write_custom_weight_ids(): write the exact analytic-sample respondent IDs in
# the two forms the workflow needs:
#   * <path_stub>.txt : ONE numeric ID per line, headerless, sorted, unique.
#                       This is the file uploaded to the NLS "Custom Weighting"
#                       tool ("Upload a list of IDs").
#   * <path_stub>.csv : the same IDs with a header, for the analyst's records.
# It also writes a small composition summary (counts by cohort/sex/race/sample
# type) next to the files so the request is self-documenting (spec Section 8).
write_custom_weight_ids <- function(data, id, path_stub,
                                    cohort = NULL, sex = NULL,
                                    race = NULL, sample_type = NULL) {
  stopifnot(id %in% names(data))
  ids <- sort(unique(suppressWarnings(as.integer(as.character(data[[id]])))))
  ids <- ids[!is.na(ids)]
  .ensure_dir(paste0(path_stub, ".txt"))

  # (a) Headerless one-ID-per-line TXT -- the NLS upload format.
  writeLines(as.character(ids), paste0(path_stub, ".txt"))

  # (b) Headed CSV for the record.
  readr::write_csv(tibble::tibble(!!id := ids), paste0(path_stub, ".csv"))

  # (c) Composition summary for the request.
  comp <- tibble::tibble(n_ids = length(ids))
  summarise_by <- function(col, nm) {
    if (!is.null(col) && col %in% names(data)) {
      t <- data %>% dplyr::count(.data[[col]], name = "n") %>%
        dplyr::mutate(variable = nm) %>%
        dplyr::rename(level = 1)
      t %>% dplyr::mutate(level = as.character(level)) %>%
        dplyr::select(variable, level, n)
    } else NULL
  }
  parts <- Filter(Negate(is.null), list(
    summarise_by(cohort, "cohort"),
    summarise_by(sex, "sex"),
    summarise_by(race, "race"),
    summarise_by(sample_type, "sample_type")
  ))
  if (length(parts) > 0L) {
    readr::write_csv(dplyr::bind_rows(parts),
                     paste0(path_stub, "_composition.csv"))
  }
  message(sprintf("[write_custom_weight_ids] %d IDs -> %s.{txt,csv}",
                  length(ids), basename(path_stub)))
  invisible(ids)
}

# ---------------------------------------------------------------------------
# 7. Small assertion helper for inline QC
# ---------------------------------------------------------------------------
assert_that <- function(cond, msg) {
  if (!isTRUE(cond)) stop(paste0("[assertion failed] ", msg))
  invisible(TRUE)
}

## ===================== [nlsy_income_helpers] =====================
# ============================================================================
# R/nlsy_income_helpers.R
# ----------------------------------------------------------------------------
# Variable maps and Bloome-Dyer-Zhou (2018) income-timing helpers.
#
# The raw NLS extract CSVs carry columns named by *reference number*
# (e.g. "R0217900"). These maps translate reference numbers into
# (survey_year, income reference year, concept) so the wide extracts can be
# reshaped to one row per respondent-income-year (spec Section 13).
#
# Timing rules implemented (spec Section 1.1):
#   * NLSY79 parent window : survey 1979-1983  -> income years 1978-1982
#   * NLSY97 parent window : survey 1997-2001  -> income years 1996-2000
#   * Adult income         : respondent age 27-32 in the INCOME reference year
#   * >= 2 valid parent obs AND >= 2 valid adult obs per respondent
#
# NLS "total net family income" (NLSY79 TNFI_TRUNC) and NLSY97
# CV_INCOME_GROSS_YR / CV_INCOME_FAMILY all refer to the PAST calendar year, so
# income_reference_year = survey_year - 1 (verified against the codebooks; the
# maps below store the reference year explicitly rather than assuming it).
# ============================================================================

# ---------------------------------------------------------------------------
# 1. NLSY79 TNFI_TRUNC reference-number -> survey/income-year map
#    (source: NLSY79HHINC.sdf; income = "total net family income, past year")
# ---------------------------------------------------------------------------
nlsy79_tnfi_map <- tibble::tribble(
  ~refnum,     ~survey_year,
  "R0217900",  1979L,
  "R0406010",  1980L,
  "R0618410",  1981L,
  "R0898600",  1982L,
  "R1144500",  1983L,
  "R1519700",  1984L,
  "R1890400",  1985L,
  "R2257500",  1986L,
  "R2444700",  1987L,
  "R2870200",  1988L,
  "R3074000",  1989L,
  "R3400700",  1990L,
  "R3656100",  1991L,
  "R4006600",  1992L,
  "R4417700",  1993L,
  "R5080700",  1994L,
  "R5166000",  1996L,
  "R6478700",  1998L,
  "R7006500",  2000L,
  "R7703700",  2002L,
  "R8496100",  2004L,
  "T0987800",  2006L,
  "T2210000",  2008L,
  "T3107800",  2010L,
  "T4112300",  2012L,
  "T5022600",  2014L,
  "T5770800",  2016L,
  "T8218700",  2018L,
  "T8787900",  2020L,
  "T9299700",  2022L
) %>%
  dplyr::mutate(income_year = survey_year - 1L, concept = "TNFI_TRUNC")

# ---------------------------------------------------------------------------
# 2. NLSY97 income reference-number maps (source: NLSY97HHINC.sdf)
#    CV_INCOME_GROSS_YR  1997-2003  : gross HH income (parental HH in early yrs)
#    CV_INCOME_FAMILY    2004-2023  : gross family income of R's OWN family
# ---------------------------------------------------------------------------
nlsy97_income_map <- tibble::tribble(
  ~refnum,     ~survey_year, ~concept,
  "R1204500",  1997L, "CV_INCOME_GROSS_YR",
  "R2563300",  1998L, "CV_INCOME_GROSS_YR",
  "R3884900",  1999L, "CV_INCOME_GROSS_YR",
  "R5464100",  2000L, "CV_INCOME_GROSS_YR",
  "R7227800",  2001L, "CV_INCOME_GROSS_YR",
  "S1541700",  2002L, "CV_INCOME_GROSS_YR",
  "S2011500",  2003L, "CV_INCOME_GROSS_YR",
  "S3812400",  2004L, "CV_INCOME_FAMILY",
  "S5412800",  2005L, "CV_INCOME_FAMILY",
  "S7513700",  2006L, "CV_INCOME_FAMILY",
  "T0014100",  2007L, "CV_INCOME_FAMILY",
  "T2016200",  2008L, "CV_INCOME_FAMILY",
  "T3606500",  2009L, "CV_INCOME_FAMILY",
  "T5206900",  2010L, "CV_INCOME_FAMILY",
  "T6656700",  2011L, "CV_INCOME_FAMILY",
  "T8129100",  2013L, "CV_INCOME_FAMILY",
  "U0008900",  2015L, "CV_INCOME_FAMILY",
  "U1845500",  2017L, "CV_INCOME_FAMILY",
  "U3444000",  2019L, "CV_INCOME_FAMILY",
  "U4949700",  2021L, "CV_INCOME_FAMILY",
  "U6356400",  2023L, "CV_INCOME_FAMILY"
) %>%
  dplyr::mutate(income_year = survey_year - 1L)

# NLSY97 household-income source flag (1997): tells whether CV_INCOME_GROSS_YR
# for 1997 came from the parent questionnaire or the youth (spec Section 13.2).
nlsy97_income_source_1997 <- "R1204600"   # CV_HH_INCOME_SOURCE (1=Parent, 2=Youth)

# CV_INCOME_GROSS_YR reference numbers for the PARENT window (survey 1997-2001).
# The parent questionnaire was fielded ONLY in 1997, so 1997 is the sole
# genuinely parent-reported year; 1998-2001 are youth-reported and noisier
# (see analysis/investigate_nlsy97_parent_income.R).
nlsy97_parent_gross_map <- c(`1997` = "R1204500", `1998` = "R2563300",
                             `1999` = "R3884900", `2000` = "R5464100",
                             `2001` = "R7227800")

# ---------------------------------------------------------------------------
# NLSY97 parent income: 1997-primary coalesce with 1998-2001 imputation
# ---------------------------------------------------------------------------

# build_nlsy97_parent_income(): construct a SINGLE parental household-income
# value per respondent using the reliable 1997 parent-reported figure as the
# primary source and coalescing to 1998->1999->2000->2001 ONLY when 1997 is
# missing. Rationale (empirically confirmed): only Round 1 carried the parent
# questionnaire; later "household income" is youth-reported, has an 8% don't-know
# rate (vs 0% for parents), is mostly valid-skipped, and correlates only ~0.44
# with the parent figure. We therefore do NOT average 1996-2000 (the RA's
# original approach); we prefer 1997 and impute only to keep otherwise-complete
# dyads. Full provenance is recorded so the imputation is auditable.
#
#   hh        : wide NLSY97 HHINC data, already standardized to carry `id_col`
#   id_col    : respondent key ("pubid")
#   require_parent_source_1997 : if TRUE, the 1997 value counts as primary only
#               when CV_HH_INCOME_SOURCE == 1 (Parent); youth-reported 1997
#               values (source == 2) fall through to imputation.
#   fallback_priority : survey years tried, in order, when 1997 is unavailable.
build_nlsy97_parent_income <- function(hh, id_col = "pubid",
                                       gross_map = nlsy97_parent_gross_map,
                                       source_col = nlsy97_income_source_1997,
                                       require_parent_source_1997 = TRUE,
                                       fallback_priority = c("1998", "1999", "2000", "2001")) {
  stopifnot(id_col %in% names(hh))
  # clean each candidate year (negative NLS codes -> NA)
  vals <- lapply(names(gross_map), function(y) nls_clean_numeric(hh[[gross_map[[y]]]]))
  names(vals) <- names(gross_map)

  src97 <- suppressWarnings(as.integer(hh[[source_col]]))   # 1 = Parent, 2 = Youth

  # primary = 1997 value; require parent-reported when asked to
  primary <- vals[["1997"]]
  if (require_parent_source_1997) primary[!(src97 %in% 1L)] <- NA_real_

  chosen_val  <- primary
  chosen_iy   <- ifelse(!is.na(primary), 1996L, NA_integer_)   # 1997 survey -> income yr 1996
  chosen_svy  <- ifelse(!is.na(primary), 1997L, NA_integer_)
  # coalesce: fill remaining NAs from the fallback years, in priority order
  for (y in fallback_priority) {
    take <- is.na(chosen_val) & !is.na(vals[[y]])
    chosen_val[take] <- vals[[y]][take]
    chosen_iy[take]  <- as.integer(y) - 1L
    chosen_svy[take] <- as.integer(y)
  }

  tibble::tibble(
    !!id_col                     := as.character(hh[[id_col]]),
    parent_income_nominal        = chosen_val,
    parent_income_year           = chosen_iy,          # income reference year
    parent_income_source_round   = chosen_svy,         # survey year value came from
    parent_income_is_parent_reported = chosen_svy == 1997L & !is.na(chosen_svy),
    parent_income_is_imputed     = !is.na(chosen_svy) & chosen_svy != 1997L,
    n_parent_income_years_available =
      Reduce(`+`, lapply(vals, function(v) as.integer(!is.na(v)))),
    has_parent_income            = !is.na(chosen_val)
  )
}

# crossvalidate_nlsy97_parent_income(): assess whether each 1998-2001 youth
# report is a trustworthy stand-in for the 1997 parent-reported value, using
# only respondents who have BOTH. Reports pair count, log correlation, median
# ratio (fallback / 1997), and the share landing within +/-25% and +/-50% of
# the parent figure. High corr / ratio ~ 1 / high within-band share => the
# fallback year is a valid imputation source.
crossvalidate_nlsy97_parent_income <- function(hh,
                                               source_col = nlsy97_income_source_1997,
                                               gross_map = nlsy97_parent_gross_map) {
  v97 <- nls_clean_numeric(hh[[gross_map[["1997"]]]])
  src <- suppressWarnings(as.integer(hh[[source_col]]))
  base <- ifelse(src %in% 1L, v97, NA_real_)           # parent-reported 1997 only
  purrr::map_dfr(c("1998", "1999", "2000", "2001"), function(y) {
    vy <- nls_clean_numeric(hh[[gross_map[[y]]]])
    ok <- !is.na(base) & base > 0 & !is.na(vy) & vy > 0
    lr <- log(vy[ok] / base[ok])
    tibble::tibble(
      fallback_year = y,
      n_pairs       = sum(ok),
      corr_log      = round(stats::cor(log(base[ok]), log(vy[ok])), 3),
      median_ratio  = round(stats::median(vy[ok] / base[ok]), 3),
      pct_within_25 = round(mean(abs(lr) < log(1.25)) * 100, 1),
      pct_within_50 = round(mean(abs(lr) < log(1.50)) * 100, 1)
    )
  })
}

# build_nlsy97_parent_income_avg(): the equivalized/real PARENT income values
# under the "1997-primary, AVERAGED fallback" rule. Identical to
# build_nlsy97_parent_income for the primary case (1997 parent-reported when
# source == Parent), but when 1997 is missing it imputes with the MEAN of the
# available 1998-2001 youth-reported years -- each deflated to real$ and
# size-equivalized AT ITS OWN survey year BEFORE averaging -- rather than the
# first available single year. This uses more of the observed information for
# the ~1.4k imputed dyads and is the construction that matches the RA's intended
# design. Returns one row per respondent with the equivalized (primary) and
# real (robustness) parent income plus the fallback-year count for auditing.
# `size_by_year` must carry (id_col, survey_year, household_size).
build_nlsy97_parent_income_avg <- function(hh, size_by_year, pce, id_col = "pubid",
                                           gross_map = nlsy97_parent_gross_map,
                                           source_col = nlsy97_income_source_1997,
                                           require_parent_source_1997 = TRUE) {
  stopifnot(id_col %in% names(hh), id_col %in% names(size_by_year))
  yrs  <- names(gross_map)
  vals <- lapply(yrs, function(y) nls_clean_numeric(hh[[gross_map[[y]]]]))
  names(vals) <- yrs
  src97 <- suppressWarnings(as.integer(hh[[source_col]]))
  if (require_parent_source_1997) vals[["1997"]][!(src97 %in% 1L)] <- NA_real_

  # long: one row per (id, candidate survey year) with a non-missing gross value
  long <- purrr::map_dfr(yrs, function(y) tibble::tibble(
    !!id_col := as.character(hh[[id_col]]),
    survey_year = as.integer(y),
    income_year = as.integer(y) - 1L,       # 1997 survey -> income year 1996, etc.
    is_primary  = y == "1997",
    parent_income_nominal = vals[[y]]))
  long <- long[!is.na(long$parent_income_nominal), , drop = FALSE]

  sby <- size_by_year
  sby[[id_col]] <- as.character(sby[[id_col]])
  long <- dplyr::left_join(long, sby[, c(id_col, "survey_year", "household_size")],
                           by = c(id_col, "survey_year"))
  long <- add_real_equiv(long, "parent_income_nominal", "income_year", "household_size", pce,
                         out_real = "parent_income_real", out_equiv = "parent_income_equiv")

  # aggregate: 1997 primary when present; else MEAN of available 1998-2001 years
  long %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      parent_income_real  = if (any(is_primary)) parent_income_real [is_primary][1]
                            else mean(parent_income_real [!is_primary], na.rm = TRUE),
      parent_income_equiv = if (any(is_primary)) parent_income_equiv[is_primary][1]
                            else mean(parent_income_equiv[!is_primary], na.rm = TRUE),
      n_parent_fallback_yrs = sum(!is_primary),
      .groups = "drop") %>%
    dplyr::mutate(dplyr::across(c(parent_income_real, parent_income_equiv),
                                ~ ifelse(is.nan(.x), NA_real_, .x)))
}

# ---------------------------------------------------------------------------
# 3. Core identifier / demographic reference numbers
# ---------------------------------------------------------------------------
nlsy79_core_refnums <- c(
  caseid       = "R0000100",
  birth_month  = "R0000300",
  birth_year   = "R0000500",
  sample_id    = "R0173600",
  race         = "R0214700",
  sex          = "R0214800"
)

nlsy79_demo_refnums <- c(
  caseid         = "R0000100",
  hhid           = "R0000149",   # household ID -- sibling clustering ONLY
  hgc_mother     = "R0006500",
  hgc_father     = "R0007900",
  hgc_respondent = "T9900000"    # HGC_EVER (highest grade ever completed, XRND)
)

nlsy97_core_refnums <- c(
  pubid        = "R0000100",
  sex          = "R0536300",
  birth_month  = "R0536401",
  birth_year   = "R0536402",
  sample_type  = "R1235800",
  race         = "R1482600"
)

nlsy97_demo_refnums <- c(
  pubid          = "R0000100",
  sidcode        = "R1193000",   # household/sibling ID -- clustering ONLY
  hh_size_1997   = "R1205400",   # CV_HH_SIZE 1997
  vstrat         = "R1489700",   # VSTRAT
  vpsu           = "R1489800",   # VPSU
  hgc_father     = "R1302400",   # CV_HGC_BIO_DAD
  hgc_mother     = "R1302500",   # CV_HGC_BIO_MOM
  weight_cc_1997 = "R1236101"    # SAMPLING_WEIGHT_CC round 1
)

# ---------------------------------------------------------------------------
# 3b. Family-size / household-size maps and the PCE deflator
# ---------------------------------------------------------------------------

# Need-adjustment denominator: NLSY79 FAMSIZE and NLSY97 CV_HH_SIZE, one per
# survey year, from the demographic augmentation extracts. Survey years line up
# with the income survey years, so size can be joined to income by survey_year.
nlsy79_famsize_map <- tibble::tribble(
  ~refnum,~survey_year, "R0217502",1979L,"R0405210",1980L,"R0647103",1981L,
  "R0896710",1982L,"R1144410",1983L,"R1519610",1984L,"R1890210",1985L,
  "R2257410",1986L,"R2444610",1987L,"R2870110",1988L,"R3073910",1989L,
  "R3400600",1990L,"R3656000",1991L,"R4006500",1992L,"R4417600",1993L,
  "R5080600",1994L,"R5165900",1996L,"R6478600",1998L,"R7006400",2000L,
  "R7703600",2002L,"R8496000",2004L)

nlsy97_hhsize_map <- tibble::tribble(
  ~refnum,~survey_year, "R1205400",1997L,"R2563700",1998L,"R3885300",1999L,
  "R5464500",2000L,"R7228200",2001L,"S1542100",2002L,"S2011900",2003L,
  "S3813400",2004L,"S5413000",2005L,"S7513900",2006L,"T0014300",2007L,
  "T2016400",2008L,"T3606700",2009L,"T5207100",2010L,"T6656900",2011L,
  "T8129300",2013L,"U0009100",2015L,"U1845700",2017L,"U3444200",2019L,
  "U4949900",2021L,"U6356600",2023L)

# reshape_size_long(): pivot the mapped size reference numbers to one row per
# (id, survey_year), cleaned of NLS negative codes.
reshape_size_long <- function(data, id_col, map) {
  present <- map[map$refnum %in% names(data), , drop = FALSE]
  data %>%
    dplyr::select(dplyr::all_of(c(id_col, present$refnum))) %>%
    tidyr::pivot_longer(-dplyr::all_of(id_col),
                        names_to = "refnum", values_to = "size_raw") %>%
    dplyr::left_join(present, by = "refnum") %>%
    dplyr::transmute(!!id_col := .data[[id_col]], survey_year,
                     family_size = nls_clean_numeric(size_raw))
}

# load_pce_index(): read the BEA PCE chain-type price index (annual, 2017=100).
#   Source: U.S. Bureau of Economic Analysis, series DPCERG (NIPA Table 2.3.4),
#   downloaded from FRED series DPCERG3A086NBEA:
#   https://fred.stlouisfed.org/series/DPCERG3A086NBEA
# FRED CSV format is (observation_date, DPCERG3A086NBEA). Returns a tidy
# (year, pce_index) tibble; attaches the base-year index for `base_year`.
load_pce_index <- function(path, base_year) {
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  names(raw)[1:2] <- c("observation_date", "pce_index")
  out <- raw %>%
    dplyr::mutate(year = as.integer(substr(as.character(observation_date), 1, 4)),
                  pce_index = as.numeric(pce_index)) %>%
    dplyr::filter(!is.na(pce_index)) %>%
    dplyr::select(year, pce_index)
  base <- out$pce_index[out$year == base_year]
  if (length(base) != 1) stop(sprintf("[load_pce_index] base year %d not found.", base_year))
  attr(out, "base_year") <- base_year
  attr(out, "base_index") <- base
  attr(out, "source") <- "BEA DPCERG3A086NBEA (2017=100); https://fred.stlouisfed.org/series/DPCERG3A086NBEA"
  out
}

# add_real_equiv(): attach PCE-deflated real income and sqrt(size)-equivalized
# income to a frame. Deflation base = pce base_index / pce at the INCOME year.
# Stores BOTH the size-unadjusted real income (robustness) and the equivalized
# income (primary). Missing/zero size -> equiv is NA (real is still kept).
add_real_equiv <- function(df, nominal_col, income_year_col, size_col, pce,
                           out_real = "income_real", out_equiv = "income_equiv") {
  base_index <- attr(pce, "base_index")
  pce_lookup <- stats::setNames(pce$pce_index, pce$year)
  iy   <- df[[income_year_col]]
  pcey <- as.numeric(pce_lookup[as.character(iy)])
  real <- df[[nominal_col]] * base_index / pcey
  sz   <- df[[size_col]]
  df[[out_real]]  <- real
  df[[out_equiv]] <- ifelse(!is.na(sz) & sz > 0, real / sqrt(sz), NA_real_)
  df
}

# summarise_income_person_real(): like summarise_income_person() but averages
# BOTH the size-unadjusted real income and the equivalized income within person,
# for each generation. Requires income_real / income_equiv columns present.
summarise_income_person_real <- function(long, id_col,
                                         min_parent = 2L, min_adult = 2L) {
  long %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      n_parent_income_obs = sum(is_parent_obs, na.rm = TRUE),
      n_adult_income_obs  = sum(is_adult_obs,  na.rm = TRUE),
      mean_parent_income_real  = mean(income_real[is_parent_obs],  na.rm = TRUE),
      mean_parent_income_equiv = mean(income_equiv[is_parent_obs], na.rm = TRUE),
      mean_adult_income_real   = mean(income_real[is_adult_obs],   na.rm = TRUE),
      mean_adult_income_equiv  = mean(income_equiv[is_adult_obs],  na.rm = TRUE),
      parent_income_years = paste(sort(income_year[is_parent_obs]), collapse = ";"),
      adult_income_years  = paste(sort(income_year[is_adult_obs]),  collapse = ";"),
      .groups = "drop") %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("mean_"),
                                ~ ifelse(is.nan(.x), NA_real_, .x)),
                  complete_income_dyad = n_parent_income_obs >= min_parent &
                                         n_adult_income_obs  >= min_adult)
}

# ---------------------------------------------------------------------------
# 4. Reshape a wide extract to long income observations
# ---------------------------------------------------------------------------

# reshape_income_long(): given a wide data frame keyed by `id_col`, a map with
# columns (refnum, survey_year, income_year, concept), pivot the mapped columns
# to one row per (id, survey_year). Cleans NLS negatives and records the
# missing reason. Only reference numbers actually present are used.
reshape_income_long <- function(data, id_col, map) {
  present <- map[map$refnum %in% names(data), , drop = FALSE]
  long <- data %>%
    dplyr::select(dplyr::all_of(c(id_col, present$refnum))) %>%
    tidyr::pivot_longer(-dplyr::all_of(id_col),
                        names_to = "refnum", values_to = "income_raw") %>%
    dplyr::left_join(present, by = "refnum") %>%
    dplyr::mutate(
      missing_reason = nls_missing_reason(income_raw),
      income_nominal = nls_clean_numeric(income_raw)
    )
  long
}

# ---------------------------------------------------------------------------
# 5. Bloome-Dyer-Zhou generation tagging
# ---------------------------------------------------------------------------

# tag_income_generation(): flag each long income observation as a valid PARENT
# or ADULT observation under the BDZ rules.
#   parent_survey_years : e.g. 1979:1983 (NLSY79) or 1997:2001 (NLSY97)
#   birth_year          : respondent birth year (to derive age in income year)
#   adult_age_min/max   : 27 / 32 in the income reference year
# A `parent_source_ok` vector (already computed by the caller from
# co-residence / income-source evidence) gates parental observations.
tag_income_generation <- function(long, birth_year_lookup, id_col,
                                  parent_survey_years,
                                  adult_age_min = 27L, adult_age_max = 32L,
                                  parent_source_ok = NULL) {
  long <- long %>%
    dplyr::left_join(birth_year_lookup, by = id_col) %>%
    dplyr::mutate(
      age_in_income_year = income_year - birth_year,
      is_valid_income    = missing_reason == "valid" & !is.na(income_nominal),
      is_parent_window   = survey_year %in% parent_survey_years,
      is_adult_window    = age_in_income_year >= adult_age_min &
                           age_in_income_year <= adult_age_max
    )
  if (!is.null(parent_source_ok)) {
    long$parent_source_ok <- parent_source_ok
  } else {
    long$parent_source_ok <- TRUE
  }
  long %>%
    dplyr::mutate(
      is_parent_obs = is_valid_income & is_parent_window & parent_source_ok,
      is_adult_obs  = is_valid_income & is_adult_window
    )
}

# ---------------------------------------------------------------------------
# 6. Collapse long income obs to one row per respondent
# ---------------------------------------------------------------------------

# summarise_income_person(): count valid parent/adult obs and average the
# (already real, size-adjusted) income within respondent. Averaging is
# UNWEIGHTED across a person's own annual observations, per spec Section 13.5.
#   value_col : the real-equivalized income column to average
summarise_income_person <- function(long, id_col, value_col) {
  long %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      n_parent_income_obs = sum(is_parent_obs, na.rm = TRUE),
      n_adult_income_obs  = sum(is_adult_obs,  na.rm = TRUE),
      mean_parent_income  = mean(.data[[value_col]][is_parent_obs], na.rm = TRUE),
      mean_adult_income   = mean(.data[[value_col]][is_adult_obs],  na.rm = TRUE),
      parent_income_years = paste(sort(income_year[is_parent_obs]), collapse = ";"),
      adult_income_years  = paste(sort(income_year[is_adult_obs]),  collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_parent_income = ifelse(is.nan(mean_parent_income), NA_real_, mean_parent_income),
      mean_adult_income  = ifelse(is.nan(mean_adult_income),  NA_real_, mean_adult_income),
      complete_income_dyad = n_parent_income_obs >= 2L & n_adult_income_obs >= 2L
    )
}

## ===================== [nlsy_occupation_helpers] =====================
# ============================================================================
# R/nlsy_occupation_helpers.R
# ----------------------------------------------------------------------------
# Build NLSY79 EGP occupation measures FROM SCRATCH via the crosswalk chain,
# so we can create measures the RA dyad does not contain:
#   * child EGP closest to age 40
#   * child longest-duration EGP (exposure-weighted modal class)
# alongside the RA's mode/best/avg (for validation).
#
# Crosswalk chain (replicates TransitionMatrices.Rmd, chunk nlsy79occupationv2):
#   raw Census code --Dorn--> OCC1990dd --IPUMS--> OCC2010 --Occ10EGP--> EGP-11
#   --map5--> 5-class EGP (Class 1 = highest)
#
# Census system by source:
#   father FAMOCC-26  : Census 1970  (Dorn 1970->1990dd)
#   child  CPSOCC80   : Census 1980  (Dorn 1980->1990dd), survey years 1982-2000
#   child  OCCALL-EMP : Census 2000  (Dorn 2000->1990dd), survey years 2002-2022
#                       (2004-2022 are 4-digit -> divide by 10 to 3-digit)
# ============================================================================

# ---------------------------------------------------------------------------
# 1. Crosswalk lookups
# ---------------------------------------------------------------------------
load_occ_crosswalks <- function(root) {
  rc <- function(p) readr::read_csv(file.path(root, p), show_col_types = FALSE, progress = FALSE)
  d70 <- rc("Data/NLSY/NLSY79/NLSY79_occ/Crosswalks/Dorn/19701990Crosswalkdd.csv")
  d80 <- rc("Data/NLSY/NLSY79/NLSY79_occ/Crosswalks/Dorn/19801990Crosswalkdd.csv")
  d00 <- rc("Data/NLSY/NLSY79/NLSY79_occ/Crosswalks/Dorn/20001990Crosswalkdd.csv")
  ipums <- rc("Data/NLSY/NLSY79/NLSY79_occ/Crosswalks/IPUMS/Crosswalk19902010IPUMS.csv")
  names(ipums)[1] <- "OCC1990"                    # header has a BOM; force name
  egp <- rc("Data/NLSY/NLSY79/NLSY79_occ/Crosswalks/EGP/Occ10EGPCrosswalk.csv")
  mklu <- function(key, val) stats::setNames(as.integer(val), as.character(as.integer(key)))
  list(
    d70      = mklu(d70$occ_1970, d70$occ1990dd),
    d80      = mklu(d80$occ_1980, d80$occ1990dd),
    d00      = mklu(d00$occ_2000, d00$occ1990dd),
    to2010   = mklu(ipums$OCC1990, ipums$OCC2010),
    to_egp11 = mklu(egp$occ10,     egp$egp10_10)
  )
}

# EGP-11 -> 5-class (Morgan 2017; verbatim from the RA's EGPRescaled case_when)
egp11_to_5 <- function(e) {
  dplyr::case_when(
    e == 1 ~ 1L, e == 2 ~ 2L,
    e %in% c(3L, 8L)      ~ 3L,
    e %in% c(9L, 7L)      ~ 4L,
    e %in% c(4L, 10L, 11L) ~ 5L,
    TRUE ~ NA_integer_
  )
}

# Crosswalk a vector of raw census codes (one census `system`) to 5-class EGP.
# Negative NLS missing codes -> NA before crosswalking; unmatched codes -> NA.
census_to_egp5 <- function(code, system, cw) {
  code <- suppressWarnings(as.numeric(code))
  code[!is.na(code) & code < 0] <- NA            # NLS negative missing codes
  if (system == "2000_4digit") { code <- floor(code / 10); system <- "2000" }
  d <- switch(system, "1970" = cw$d70, "1980" = cw$d80, "2000" = cw$d00,
              stop("unknown census system: ", system))
  occ90   <- d[as.character(as.integer(code))]
  occ2010 <- cw$to2010[as.character(occ90)]
  egp11   <- cw$to_egp11[as.character(occ2010)]
  egp11_to_5(egp11)
}

# ---------------------------------------------------------------------------
# 2. Column -> (survey year, census system) map for the labelled extract
# ---------------------------------------------------------------------------
# Given the labelled occupation column names, classify each child-occupation
# column by survey year and census coding system (for the divide-by-10 rule).
nlsy79_child_occ_columns <- function(cols) {
  tibble::tibble(col = cols) %>%
    dplyr::filter(grepl("^CPSOCC80", col) | grepl("^OCCALL", col)) %>%
    dplyr::mutate(
      survey_year = as.integer(stringr::str_extract(col, "[0-9]{4}$")),
      system = dplyr::case_when(
        grepl("^CPSOCC80", col)            ~ "1980",
        grepl("^OCCALL", col) & survey_year == 2002 ~ "2000",       # 3-digit
        grepl("^OCCALL", col) & survey_year >= 2004 ~ "2000_4digit" # 4-digit /10
      )
    )
}

# ---------------------------------------------------------------------------
# 3. Child EGP life-trajectory (long: one row per respondent-survey-year)
# ---------------------------------------------------------------------------
# data : labelled NLSY79 occupation extract with an added character `caseid`
#        and integer `birth_year`.
build_child_egp_trajectory <- function(data, id_col, cw) {
  occ_cols <- nlsy79_child_occ_columns(names(data))
  purrr::pmap_dfr(occ_cols, function(col, survey_year, system) {
    tibble::tibble(
      !!id_col := data[[id_col]],
      birth_year  = data$birth_year,
      survey_year = survey_year,
      egp5 = census_to_egp5(data[[col]], system, cw)
    )
  }) %>%
    dplyr::filter(!is.na(egp5)) %>%
    dplyr::mutate(age = survey_year - birth_year) %>%
    dplyr::arrange(.data[[id_col]], survey_year)
}

# ---------------------------------------------------------------------------
# 4. Person-level child EGP measures from the trajectory
# ---------------------------------------------------------------------------
# For each respondent, within the adult career window [age_min, age_max]:
#   child_egp_age40          : EGP at the observation closest to age 40
#   child_egp_longest_dur    : EGP class with the greatest exposure (gap-weighted)
#   child_egp_mode/best/avg  : RA-style measures (for validation)
# Tie-breaks documented inline.
summarise_child_egp <- function(traj, id_col, age_min = 25L, age_max = 55L,
                                target_age = 40L) {
  traj <- dplyr::filter(traj, age >= age_min, age <= age_max)

  traj %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::group_modify(function(g, key) {
      g <- dplyr::arrange(g, survey_year)
      yrs <- g$survey_year
      # exposure = years until the next observation; last obs gets the prior gap
      gaps <- diff(yrs)
      exposure <- if (length(gaps) == 0) 1 else c(gaps, gaps[length(gaps)])

      # --- longest-duration: max total exposure by class; tie -> higher status ---
      dur <- tapply(exposure, g$egp5, sum)
      longest <- as.integer(names(dur)[which.max(dur + 1e-9 * (max(as.integer(names(dur))) -
                                                               as.integer(names(dur))))])
      # (the tiny term breaks ties toward the lower class number = higher status)

      # --- closest to age 40: min |age-40|; tie -> older side, then later year ---
      d40 <- abs(g$age - target_age)
      cand <- which(d40 == min(d40))
      pick <- cand[which.max(g$age[cand] + g$survey_year[cand] * 1e-6)]  # older, then later

      # --- RA-style measures ---
      egp <- g$egp5
      best <- min(egp)                              # lowest code = highest status
      avgv <- mean(egp)
      ux <- unique(egp); modev <- mean(ux[tabulate(match(egp, ux)) == max(tabulate(match(egp, ux)))])

      tibble::tibble(
        n_child_occ_obs        = nrow(g),
        child_egp_age40        = g$egp5[pick],
        child_occ_age_at_meas  = g$age[pick],
        child_occ_dist_from_40 = g$age[pick] - target_age,
        child_egp_longest_dur  = longest,
        child_egp_longest_share = as.numeric(max(dur) / sum(dur)),
        child_egp_mode         = as.integer(floor(modev + 0.5)),
        child_egp_best         = as.integer(best),
        child_egp_avg          = as.integer(floor(avgv + 0.5))
      )
    }) %>%
    dplyr::ungroup()
}

## ===================== [nlsy_education_helpers] =====================
# ============================================================================
# R/nlsy_education_helpers.R
# ----------------------------------------------------------------------------
# Harmonize NLSY education to a GSS-comparable scheme.
#
# The two cohorts code education differently:
#   * years of schooling ("highest grade completed", 0-20) for the NLSY79
#     respondent and ALL parents (NLSY79 + NLSY97);
#   * a DEGREE code (0 None ... 7 Professional) for the NLSY97 respondent
#     (CVC_HIGHEST_DEGREE_EVER).
# GSS carries both `educ` (years, 0-20) and `degree` (0 LT-HS .. 4 graduate).
#
# We map everything to a single 5-category scheme that (a) both years and degree
# sources can produce and (b) lines up with GSS `degree`/BDZ's education groups:
#
#   0 = Less than high school        (< 12 years,           or degree "None")
#   1 = High school (diploma/GED)    (12 years,             or GED / HS diploma)
#   2 = Some college / associate     (13-15 years,          or Associate/AA)
#   3 = Bachelor's                   (16 years,             or Bachelor's)
#   4 = Graduate / advanced          (>= 17 years,          or Master's/PhD/Prof)
#
# Category 2 folds "some college, no degree" and "associate" together, because
# years alone cannot distinguish them (this matches GSS's own educ-vs-degree
# ambiguity). Raw values and native years are preserved so either GSS variable
# (`educ` years or `degree`) can be reconstructed.
# ============================================================================

educ5_labels <- c(
  "0" = "Less than high school",
  "1" = "High school (diploma/GED)",
  "2" = "Some college / associate",
  "3" = "Bachelor's",
  "4" = "Graduate / advanced"
)

# clean a "highest grade completed" (years) field: negatives = NLS missing,
# 95 = ungraded/other -> NA; keep 0-20.
clean_hgc_years <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.na(y) & (y < 0 | y == 95)] <- NA_real_
  y
}

# years of schooling -> 5-category harmonized education
years_to_educ5 <- function(x) {
  y <- clean_hgc_years(x)
  dplyr::case_when(
    is.na(y)            ~ NA_integer_,
    y < 12              ~ 0L,
    y == 12             ~ 1L,
    y >= 13 & y <= 15   ~ 2L,
    y == 16             ~ 3L,
    y >= 17             ~ 4L
  )
}

# NLSY97 respondent degree code -> 5-category harmonized education
#   0 None -> 0 ; 1 GED, 2 HS -> 1 ; 3 AA -> 2 ; 4 BA -> 3 ; 5 MA, 6 PhD, 7 Prof -> 4
nlsy97_degree_to_educ5 <- function(code) {
  c <- suppressWarnings(as.integer(code))
  dplyr::case_when(
    is.na(c) | c < 0 ~ NA_integer_,
    c == 0           ~ 0L,
    c %in% c(1L, 2L) ~ 1L,
    c == 3L          ~ 2L,
    c == 4L          ~ 3L,
    c %in% c(5L, 6L, 7L) ~ 4L
  )
}

# ---------------------------------------------------------------------------
# 4-category scheme: fully consistent across years AND degree sources (no
# some-college/associate ambiguity), so it is comparable ACROSS cohorts and with
# GSS without the years-vs-degree asymmetry that affects category 2 of educ5.
#   0 = Less than high school
#   1 = High school or some college (no 4-year degree)   [12-15 yrs, or GED/HS/AA]
#   2 = Bachelor's                                        [16 yrs, or BA]
#   3 = Graduate / advanced                               [>=17 yrs, or MA/PhD/Prof]
# ---------------------------------------------------------------------------
educ4_labels <- c("0" = "Less than high school",
                  "1" = "High school / some college",
                  "2" = "Bachelor's",
                  "3" = "Graduate / advanced")

years_to_educ4 <- function(x) {
  y <- clean_hgc_years(x)
  dplyr::case_when(is.na(y) ~ NA_integer_, y < 12 ~ 0L, y >= 12 & y <= 15 ~ 1L,
                   y == 16 ~ 2L, y >= 17 ~ 3L)
}
nlsy97_degree_to_educ4 <- function(code) {
  c <- suppressWarnings(as.integer(code))
  dplyr::case_when(is.na(c) | c < 0 ~ NA_integer_, c == 0 ~ 0L,
                   c %in% c(1L, 2L, 3L) ~ 1L, c == 4L ~ 2L, c %in% c(5L, 6L, 7L) ~ 3L)
}

# generic dispatchers
to_educ5 <- function(x, type = c("years", "nlsy97_degree")) {
  type <- match.arg(type)
  if (type == "years") years_to_educ5(x) else nlsy97_degree_to_educ5(x)
}
to_educ4 <- function(x, type = c("years", "nlsy97_degree")) {
  type <- match.arg(type)
  if (type == "years") years_to_educ4(x) else nlsy97_degree_to_educ4(x)
}

# add_harmonized_education(): attach, for one raw education field, the harmonized
# 5-category value, the native years (NA when the source is a degree code), the
# raw value, an observed flag, and the NLS missing reason -- prefixed by `who`
# (e.g. "respondent", "mother", "father"). Keeps generations/parents SEPARATE.
add_harmonized_education <- function(data, raw_col, who,
                                     type = c("years", "nlsy97_degree")) {
  type <- match.arg(type)
  raw <- data[[raw_col]]
  data[[paste0(who, "_educ_raw")]]            <- raw
  data[[paste0(who, "_educ5")]]               <- to_educ5(raw, type)   # GSS degree-style (5-cat)
  data[[paste0(who, "_educ4")]]               <- to_educ4(raw, type)   # fully cross-consistent (4-cat)
  data[[paste0(who, "_educ_years")]]          <- if (type == "years") clean_hgc_years(raw) else NA_real_
  data[[paste0(who, "_educ_source")]]         <- type
  data[[paste0(who, "_educ_missing_reason")]] <- nls_missing_reason(raw)
  data[[paste0(who, "_educ_observed")]]       <- !is.na(data[[paste0(who, "_educ5")]])
  data
}

## ===================== [nlsy_bootstrap_helpers] =====================
# ============================================================================
# R/nlsy_bootstrap_helpers.R
# ----------------------------------------------------------------------------
# Resampling / survey-design preparation helpers (spec Section 16).
#
# Three inference machineries are prepared here (NOT run at scale in the build;
# only tiny smoke tests):
#   1. Individual bootstrap (NLSY79 & NLSY97): resample respondents uniformly
#      with replacement; replicate weight = custom_weight * multiplicity.
#      Respondents are NEVER sampled proportional to their custom weight.
#   2. NLSY79 household-cluster bootstrap: resample original HHIDs with
#      replacement, carrying all siblings together (preserves sibling
#      dependence). This is explicitly NOT a full survey-design bootstrap --
#      HHID does not reproduce the unavailable NLSY79 geographic PSU design.
#   3. NLSY97 Fay-BRR replicate design from VSTRAT / VPSU (two VPSUs per
#      VSTRAT -> balanced repeated replication is the natural design replicate).
#
# The transition estimator recomputes weighted cut-points and class labels
# INSIDE every replicate (spec Section 1.3 / 13.6): cut-points are never fixed
# and reused across replicates.
# ============================================================================

# ---------------------------------------------------------------------------
# 1. Individual bootstrap replicate weights
# ---------------------------------------------------------------------------

# make_individual_bootstrap_input(): attach the pieces the individual bootstrap
# needs. The replicate loop itself lives in the later inferential notebook;
# here we just validate and stamp metadata so the input object is self-describing.
make_individual_bootstrap_input <- function(data, id_col, weight_col,
                                            cohort_col = NULL, scheme = "uniform_individual") {
  stopifnot(id_col %in% names(data), weight_col %in% names(data))
  attr(data, "bootstrap_scheme") <- scheme
  attr(data, "id_col")           <- id_col
  attr(data, "weight_col")       <- weight_col
  attr(data, "cohort_col")       <- cohort_col
  data
}

# draw_individual_replicate(): one uniform-with-replacement respondent resample.
# Returns the data with a replicate weight = custom_weight * multiplicity m_i.
draw_individual_replicate <- function(data, id_col, weight_col) {
  ids <- data[[id_col]]
  drawn <- sample(ids, size = length(ids), replace = TRUE)      # uniform, NOT weight-proportional
  m <- table(factor(drawn, levels = ids))
  data$m_i <- as.integer(m[match(data[[id_col]], names(m))])
  data$m_i[is.na(data$m_i)] <- 0L
  data$weight_rep <- data[[weight_col]] * data$m_i
  data[data$m_i > 0L, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 2. NLSY79 household-cluster bootstrap
# ---------------------------------------------------------------------------

# make_household_bootstrap_input(): index respondents by original HHID so the
# replicate loop can resample whole households (siblings move together).
make_household_bootstrap_input <- function(data, id_col, hhid_col, weight_col) {
  stopifnot(hhid_col %in% names(data))
  attr(data, "bootstrap_scheme") <-
    "weighted original-household cluster bootstrap preserving sibling dependence"
  attr(data, "cluster_col") <- hhid_col
  attr(data, "id_col")      <- id_col
  attr(data, "weight_col")  <- weight_col
  data
}

# draw_household_replicate(): resample HHIDs with replacement; each respondent's
# replicate weight = custom_weight * (number of times its household was drawn).
draw_household_replicate <- function(data, hhid_col, weight_col) {
  hh <- unique(data[[hhid_col]])
  drawn <- sample(hh, size = length(hh), replace = TRUE)
  mult <- table(factor(drawn, levels = hh))
  data$hh_mult <- as.integer(mult[match(data[[hhid_col]], names(mult))])
  data$hh_mult[is.na(data$hh_mult)] <- 0L
  data$weight_rep <- data[[weight_col]] * data$hh_mult
  data[data$hh_mult > 0L, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 3. NLSY97 survey design + Fay-BRR replicate design
# ---------------------------------------------------------------------------

# build_nlsy97_design(): pooled svydesign on VPSU nested in VSTRAT, weighted by
# the returned custom weight. Requires the `survey` package.
build_nlsy97_design <- function(data, weight_col = "custom_weight",
                                psu = "vpsu", strata = "vstrat") {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("[build_nlsy97_design] package 'survey' is required.")
  }
  ok <- !is.na(data[[psu]]) & !is.na(data[[strata]]) & !is.na(data[[weight_col]])
  survey::svydesign(
    ids     = stats::as.formula(paste0("~", psu)),
    strata  = stats::as.formula(paste0("~", strata)),
    weights = stats::as.formula(paste0("~", weight_col)),
    data    = data[ok, , drop = FALSE],
    nest    = TRUE
  )
}

# build_nlsy97_fay_brr(): Fay-adjusted BRR replicate design from the svydesign.
build_nlsy97_fay_brr <- function(design, fay_rho = 0.5) {
  survey::as.svrepdesign(design, type = "Fay", fay.rho = fay_rho, mse = TRUE)
}

# ---------------------------------------------------------------------------
# 4. Replicate-safe transition estimator
# ---------------------------------------------------------------------------

# transition_matrix_reest(): compute a parent x adult income transition matrix,
# recomputing weighted quintile cut-points and class labels FROM SCRATCH on the
# supplied (possibly resampled / reweighted) data. This is the function every
# replicate calls, guaranteeing cut-points are never frozen across replicates.
transition_matrix_reest <- function(data, parent_col, adult_col, weight_col,
                                    n_classes = 5L, scheme = c("pooled_parent_child",
                                                               "generation_specific")) {
  scheme <- match.arg(scheme)
  probs <- seq_len(n_classes - 1L) / n_classes
  w <- data[[weight_col]]
  if (scheme == "pooled_parent_child") {
    pooled <- c(data[[parent_col]], data[[adult_col]])
    pooled_w <- c(w, w)
    cuts_p <- cuts_a <- weighted_cutpoints(pooled, pooled_w, probs)
  } else {
    cuts_p <- weighted_cutpoints(data[[parent_col]], w, probs)
    cuts_a <- weighted_cutpoints(data[[adult_col]],  w, probs)
  }
  pc <- assign_income_class(data[[parent_col]], cuts_p)
  ac <- assign_income_class(data[[adult_col]],  cuts_a)
  # Weighted joint distribution, row-normalised to transition probabilities.
  tab <- tapply(w, list(parent = factor(pc, 1:n_classes),
                        adult  = factor(ac, 1:n_classes)),
                FUN = sum)
  tab[is.na(tab)] <- 0
  trans <- sweep(tab, 1, rowSums(tab), "/")
  list(counts = tab, transition = trans, cut_parent = cuts_p, cut_adult = cuts_a)
}

# pooled_income_classes(): the PROJECT income-class definition (cfg$income_class_scheme
# = "pooled_parent_child_per_survey"). Pool father and child (equivalized) real
# income WITHIN one survey into a single weighted distribution, take `n_classes`
# weighted-quantile cutoffs, and assign BOTH generations' class from those SAME
# cutoffs. Call it once per survey (NLSY79, NLSY97) so pooling is per-survey,
# across cohorts -- NOT within-cohort. class 1 = lowest .. n = highest; set
# descending = TRUE for the "Class 1 = highest" labeling used in the RA's TM report.
pooled_income_classes <- function(data, parent_col, adult_col, weight_col,
                                  n_classes = 5L, descending = FALSE) {
  w <- data[[weight_col]]
  probs <- seq_len(n_classes - 1L) / n_classes
  cuts <- weighted_cutpoints(c(data[[parent_col]], data[[adult_col]]),
                             c(w, w), probs)
  pc <- assign_income_class(data[[parent_col]], cuts)
  ac <- assign_income_class(data[[adult_col]],  cuts)
  if (descending) { pc <- n_classes + 1L - pc; ac <- n_classes + 1L - ac }
  list(cutpoints = cuts, parent_class = pc, adult_class = ac)
}
