# nlsy97_clean_household_income.R
# ---------------------------------------------------------------------------
# Final NLSY97 income construction: RELEASED summary variables + a per-round
# family-only household eligibility screen. No component reconstruction.
#
#   parent income = CV_INCOME_GROSS_YR, Round-1 parent-reported PRIMARY +
#                   Rounds 2-5 averaged fallback (the baseline rule),
#                   with each household-YEAR excluded when that round's roster
#                   contains a confirmed nonfamily member (Roommate/Other
#                   non-relative).
#   adult income  = CV_INCOME_FAMILY, ages 27-38 (baseline rule, cfg).
#
# Relationship source: the consolidated NLSY97 input (round-specific HHI2_RELY in
# Round 1, HHI_RELY in Rounds 2-5). Occupancy = RELY present and not a valid
# -4/-5 skip (verified to match the Round-1 HHI2_ID occupancy to 99.97%).
#
# This file is the SOLE active NLSY97 robustness income construction; it depends
# only on shared baseline utilities (add_real_equiv, build_nlsy97_parent_income_avg,
# nlsy97_parent_gross_map, reshape_income_long, cfg, pce) and NEVER touches the
# NLSY79 pipeline or the shared bootstrap/estimator functions.
# ---------------------------------------------------------------------------

# roster relationship code classification (from the local value labels):
#   nonfamily: 68 Roommate, 85 Other non-relative
#   unknown  : 99 relationship missing, and negative refusal/DK codes (-1,-2,-3)
#   family   : every other substantive code (0=self, parents/step/adoptive/foster,
#              siblings, in-laws, grandparents, aunts/uncles/cousins/nieces/
#              nephews, 69 Lover/partner, 84 Other relative, ...)
REL_NONFAMILY_CODES <- c(68, 85)
REL_UNKNOWN_CODES   <- c(99, -1, -2, -3)   # occupied but relationship not ascertained
REL_SKIP_CODES      <- c(-4, -5)           # empty roster slot / non-interview => UNOCCUPIED

# survey year <-> round <-> income reference year map (Rounds 1-5), with the
# HHI2_RELY (R1) / HHI_RELY (R2-5) roster relationship column BASE R-number and
# position count. Roster column for position k = sprintf("R%05d00", base+k-1)
# (the CSV column-name form, e.g. R13158.00 -> "R1315800"). These come from the
# consolidated NLSY97 input, so no codebook is needed at run time.
NLSY97_ROUND_MAP <- tibble::tribble(
  ~round, ~survey_year, ~income_year, ~rely_base, ~n_pos,
  1L, 1997L, 1996L, 13158L, 17L,
  2L, 1998L, 1997L, 24163L, 14L,
  3L, 1999L, 1998L, 37269L, 14L,
  4L, 2000L, 1999L, 51918L, 14L,
  5L, 2001L, 2000L, 69197L, 16L
)
.rely_col <- function(base, k) sprintf("R%05d00", base + (k - 1L))

# ---------------------------------------------------------------------------
# build_nlsy97_round_clean_flags(): one row per respondent with, for each of
# Rounds 1-5, clean_family_household_r{1..5} (primary rule = no confirmed
# nonfamily member) plus has_nonfamily / has_unknown per round and a strict
# variant. Households with NO occupied roster in a round -> that round's flag is
# NA. `data` is the consolidated NLSY97 input (carries the RELY R-number columns
# + R0000100); a path is also accepted.
# ---------------------------------------------------------------------------
build_nlsy97_round_clean_flags <- function(data = NULL, id_col = "pubid") {
  if (is.null(data)) data <- "Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_raw.csv"
  d <- if (is.character(data)) as.data.frame(data.table::fread(data, showProgress = FALSE)) else as.data.frame(data)
  g  <- function(c) if (!c %in% names(d)) rep(NA_real_, nrow(d)) else suppressWarnings(as.numeric(d[[c]]))
  out <- tibble::tibble(!!id_col := as.character(g("R0000100")))

  for (rr in seq_len(nrow(NLSY97_ROUND_MAP))) {
    row <- NLSY97_ROUND_MAP[rr, ]
    rely <- lapply(seq_len(row$n_pos), function(i) g(.rely_col(row$rely_base, i)))
    occ  <- lapply(rely, function(v) !is.na(v) & !(v %in% REL_SKIP_CODES))          # occupied roster slot
    nonf <- Reduce(`|`, Map(function(o, v) o & !is.na(v) & v %in% REL_NONFAMILY_CODES, occ, rely))
    unk  <- Reduce(`|`, Map(function(o, v) o & !is.na(v) & v %in% REL_UNKNOWN_CODES,   occ, rely))
    n_occ <- Reduce(`+`, lapply(occ, as.integer))
    observed <- n_occ > 0                                    # a household roster exists this round
    clean       <- ifelse(observed, !nonf, NA)              # primary rule: no confirmed nonfamily member
    clean_strict<- ifelse(observed, !nonf & !unk, NA)
    out[[sprintf("clean_family_household_r%d", row$round)]]        <- clean
    out[[sprintf("clean_family_household_strict_r%d", row$round)]] <- clean_strict
    out[[sprintf("has_nonfamily_r%d", row$round)]] <- ifelse(observed, nonf, NA)
    out[[sprintf("has_unknown_r%d",   row$round)]] <- ifelse(observed, unk,  NA)
    out[[sprintf("roster_observed_r%d", row$round)]] <- observed
  }
  out
}

# round -> released CV_INCOME_GROSS_YR column (colname form), for the screen
NLSY97_ROUND_GROSS_COL <- c(`1` = "R1204500", `2` = "R2563300",
                            `3` = "R3884900", `4` = "R5464100", `5` = "R7227800")

# ---------------------------------------------------------------------------
# build_nlsy97_income_dyads_robust(): the SOLE active NLSY97 robustness income
# construction. Restores the baseline released-summary-variable pipeline
# and adds ONE screen: each early household-YEAR whose roster contains a
# CONFIRMED nonfamily member (Roommate/Other non-relative) is excluded from the
# parent-income candidate set BEFORE the existing 1997-primary / 1998-2001-mean
# fallback rule is applied.  Adult income (CV_INCOME_FAMILY, 27-38) is untouched
# (§7: the early screen is NOT applied to later CV_INCOME_FAMILY).
#
# It reuses the baseline builders line-for-line (build_nlsy97_parent_income,
# build_nlsy97_parent_income_avg) on the screened data; the ONLY change is the
# household-year screen.  Returns the baseline preweight person record
# schema exactly, with the parent side re-derived on screened data and
# complete_income_dyad recomputed.
#
#   preweight   : the baseline nlsy97_income_person_preweight.rds (adult income,
#                 demographics, sample flags, design vars — all preserved).
#   hh_raw      : raw NLSY97 extract carrying the released gross columns
#                 (nlsy97_parent_gross_map) + CV_HH_INCOME_SOURCE + R0000100.
#   round_flags : output of build_nlsy97_round_clean_flags().
# ---------------------------------------------------------------------------
build_nlsy97_income_dyads_robust <- function(preweight, hh_raw, round_flags,
                                             size_by_year, pce, cfg,
                                             id_col = "pubid",
                                             audit_dir = NULL) {
  hh <- standardize_id(hh_raw, "R0000100", id_col)
  fl <- round_flags; fl[[id_col]] <- as.character(fl[[id_col]])
  hh[[id_col]] <- as.character(hh[[id_col]])
  hh <- dplyr::left_join(hh, fl, by = id_col)

  # ---- THE SCREEN: null a round's gross value where that round is CONFIRMED
  #      nonfamily (has_nonfamily_r == TRUE). NA / not-observed rounds are NOT
  #      confirmed nonfamily and are left eligible (exclude only confirmed).
  n_screened <- integer(0)
  for (r in as.character(1:5)) {
    gcol <- NLSY97_ROUND_GROSS_COL[[r]]; fcol <- sprintf("has_nonfamily_r%s", r)
    if (!gcol %in% names(hh) || !fcol %in% names(hh)) next
    contaminated <- hh[[fcol]] %in% TRUE
    hh[[gcol]][contaminated] <- NA
    n_screened[r] <- sum(contaminated)
  }

  # ---- re-derive parent income on the SCREENED data with the UNCHANGED builders
  parent <- build_nlsy97_parent_income(hh, id_col = id_col) %>%
    dplyr::select(-dplyr::any_of(c("parent_income_real", "parent_income_equiv"))) %>%
    dplyr::left_join(build_nlsy97_parent_income_avg(hh, size_by_year, pce, id_col = id_col),
                     by = id_col)

  # ---- swap the screened parent side into the preweight record, recompute dyad
  keep_parent <- setdiff(names(parent), id_col)
  pw <- preweight
  pw[[id_col]] <- as.character(pw[[id_col]])
  pw <- pw %>% dplyr::select(-dplyr::any_of(keep_parent)) %>%
    dplyr::left_join(parent, by = id_col) %>%
    dplyr::mutate(
      n_adult_income_obs = dplyr::coalesce(n_adult_income_obs, 0L),
      parent_dyad_ok       = has_parent_income,   # default rule: coalesce_1997_primary
      complete_income_dyad = dplyr::coalesce(parent_dyad_ok &
                              n_adult_income_obs >= cfg$min_adult_income_obs, FALSE),
      nlsy97_zhou_income_sample = is_comparable_civilian_sample & complete_income_dyad)

  if (!is.null(audit_dir)) {
    dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(round_flags, file.path(audit_dir, "nlsy97_family_relationship_audit.rds"))
    saveRDS(data.frame(round = 1:5, gross_col = unname(NLSY97_ROUND_GROSS_COL),
                       n_household_years_screened = as.integer(n_screened[as.character(1:5)])),
            file.path(audit_dir, "nlsy97_income_year_screen_audit.rds"))
  }
  attr(pw, "n_screened_by_round") <- n_screened
  pw
}

# ---------------------------------------------------------------------------
# build_nlsy97_adult_released_income(): adult (child) income from RELEASED
# CV_INCOME_FAMILY for an arbitrary age window, reusing the exact baseline
# machinery (reshape_income_long + nlsy97_income_map + tag_income_generation +
# add_real_equiv). Used to build the BDZ-style 27-32 window; the 27-38 window it
# produces reproduces the baseline preweight adult income. Returns per-respondent
# n_adult_income_obs + mean real/equiv adult income for the window.
# ---------------------------------------------------------------------------
build_nlsy97_adult_released_income <- function(hh_raw, birth_year_df, size_by_year, pce, cfg,
                                               age_lo = 27L, age_hi = 38L, id_col = "pubid") {
  hh <- standardize_id(hh_raw, "R0000100", id_col)
  long <- hh %>%
    reshape_income_long(id_col, nlsy97_income_map) %>%
    tag_income_generation(birth_year_df, id_col,
                          parent_survey_years = cfg$parent_rounds_97,
                          adult_age_min = age_lo, adult_age_max = age_hi) %>%
    dplyr::left_join(dplyr::rename(size_by_year, family_size = household_size),
                     by = c(id_col, "survey_year")) %>%
    add_real_equiv("income_nominal", "income_year", "family_size", pce)
  long %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      n_adult_income_obs      = sum(is_adult_obs, na.rm = TRUE),
      mean_adult_income_real  = mean(income_real[is_adult_obs],  na.rm = TRUE),
      mean_adult_income_equiv = mean(income_equiv[is_adult_obs], na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("mean_"), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

# ---------------------------------------------------------------------------
# build_nlsy97_parent_measureB_screened(): BDZ's parental-income rule
# ("average all years in which parents reported their income", BDZ fn20) =
# the mean of parent-reported-year income across 1996-2000 (MEASURE B), with the
# family-only household-YEAR screen applied FIRST: any parent-year whose roster
# has a confirmed nonfamily member is dropped BEFORE averaging (§2.3). A
# respondent is retained if >= minobs clean parent-years remain.
#   income_long : the baseline nlsy97_income_long (is_parent_obs, survey_year,
#                 income_real, income_equiv).
#   round_flags : build_nlsy97_round_clean_flags() output.
# Returns per-respondent parent_income_real/equiv (measure B, screened) + np.
# ---------------------------------------------------------------------------
build_nlsy97_parent_measureB_screened <- function(income_long, round_flags,
                                                  id_col = "pubid", minobs = 1L, screen = TRUE) {
  rmap <- c(`1`=1997L,`2`=1998L,`3`=1999L,`4`=2000L,`5`=2001L)
  flong <- dplyr::bind_rows(lapply(1:5, function(r) tibble::tibble(
    !!id_col := as.character(round_flags[[id_col]]),
    survey_year = rmap[[as.character(r)]],
    nonfam = round_flags[[sprintf("has_nonfamily_r%d", r)]] %in% TRUE)))
  p <- income_long %>% dplyr::mutate(!!id_col := as.character(.data[[id_col]])) %>%
    dplyr::filter(is_parent_obs)
  if (screen) p <- p %>% dplyr::left_join(flong, by = c(id_col, "survey_year")) %>%
    dplyr::filter(!(nonfam %in% TRUE))                          # drop contaminated parent-YEARS
  p %>% dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(parent_income_real  = mean(income_real,  na.rm = TRUE),
                     parent_income_equiv = mean(income_equiv, na.rm = TRUE),
                     np_parent_measureB  = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(np_parent_measureB >= minobs)
}
