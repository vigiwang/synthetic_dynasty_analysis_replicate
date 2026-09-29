#!/usr/bin/env Rscript
# ===========================================================================
# 02b_apply_nlsy97_nonrelative_clean_household.R
# ---------------------------------------------------------------------------
# Applies the NLSY97 family-only household screen to the income dyads BEFORE any
# weight attachment (clean the dyads first, then weight in 03).
#
# Each early parent-income household-YEAR whose roster (HHI2_RELY / HHI_RELY,
# Rounds 1-5) contains a confirmed non-relative (code 68 Roommate / 85 Other
# non-relative) is excluded before the parent measure is formed. A respondent is
# dropped only if no valid clean parent-year remains. This OVERWRITES the NLSY97
# income preweight with the cleaned sample so 03 attaches the existing weights to
# the cleaned dyads (no weight recalibration for the ~90 removed obs).
# ===========================================================================
HPC      <- TRUE
dir_root <- if (HPC) "/home/weiqiw/synthetic_dynasty_analysis" else "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(purrr); library(tibble); library(readr)})
options(warn = 1)
source(file.path("NLSY_fun", "nlsy_config.R"))
source(file.path("NLSY_fun", "nlsy_data_functions.R"))
source(file.path("NLSY_fun", "nlsy97_clean_household_income.R"))

ANALY <- file.path(dir_root, "Data", "NLSY_estimation", "dyads")
V     <- file.path(dir_root, "Data", "NLSY", "_validation"); dir.create(V, recursive = TRUE, showWarnings = FALSE)
CANON <- file.path(dir_root, "Data", "NLSY", "NLSY97", "NLSY97_inc")

preweight <- readRDS(file.path(ANALY, "nlsy97_income_person_preweight.rds"))   # baseline (pre-screen)
hh_raw    <- as.data.frame(data.table::fread(p$n97_inc, showProgress = FALSE)) # consolidated NLSY97 input
size_by_year <- reshape_size_long(standardize_id(hh_raw, "R0000100", "pubid"), "pubid", nlsy97_hhsize_map) %>%
  rename(household_size = family_size)

message("[02b] building NLSY97 family-only household flags (R1-5) ...")
round_flags <- build_nlsy97_round_clean_flags(hh_raw)

message("[02b] applying non-relative screen to NLSY97 income dyads ...")
ch <- build_nlsy97_income_dyads_robust(preweight, hh_raw, round_flags, size_by_year, pce, cfg)

## exclusion audit: exactly which baseline-complete dyads are dropped, and why
base_ids  <- preweight %>% filter(complete_income_dyad) %>% pull(pubid) %>% as.character()
clean_ids <- ch        %>% filter(complete_income_dyad) %>% pull(pubid) %>% as.character()
excluded  <- setdiff(base_ids, clean_ids)
cons <- hh_raw; cons$pubid <- as.character(cons[["R0000100"]]); lab <- c("68"="Roommate","85"="Other non-relative")
aud <- bind_rows(lapply(excluded, function(id) {
  i <- which(cons$pubid == id); hits <- c()
  for (r in 1:5) { base <- NLSY97_ROUND_MAP$rely_base[r]; np <- NLSY97_ROUND_MAP$n_pos[r]
    codes <- vapply(seq_len(np), function(k) suppressWarnings(as.numeric(cons[[sprintf("R%05d00", base+k-1)]][i])), 0)
    h <- codes[codes %in% c(68,85)]; if (length(h)) hits <- c(hits, sprintf("R%d:%s", r, paste(h, collapse="/"))) }
  cc <- unique(unlist(regmatches(hits, gregexpr("68|85", hits))))
  tibble(pubid = id, affected_rounds = paste(hits, collapse=";"), relationship_codes = paste(cc, collapse="/"),
         derived_label = paste(lab[cc], collapse="/"), in_former_complete_dyads = TRUE, final_exclusion = TRUE,
         reason = "confirmed non-relative in early parent-income round(s); no clean parent-year remained") }))
write.csv(aud, file.path(V, "NLSY97_nonrelative_exclusion_audit.csv"), row.names = FALSE)

## OVERWRITE the NLSY97 income preweight with the cleaned sample (03 weights this)
saveRDS(ch, file.path(ANALY, "nlsy97_income_person_preweight.rds"))
## canonical preweight CSV (no weights)
dir.create(CANON, recursive = TRUE, showWarnings = FALSE)
ch %>% filter(complete_income_dyad) %>%
  transmute(id = as.character(pubid), vstrat, vpsu, cohort = birth_year,
            parent = parent_income_equiv, child = mean_adult_income_equiv) %>%
  write.csv(file.path(CANON, "NLSY97_inc_complete_dyads_preweight.csv"), row.names = FALSE)

message(sprintf("[02b] NLSY97 income: baseline complete=%d -> clean complete=%d (excluded %d non-relative; audit written)",
                length(base_ids), length(clean_ids), length(excluded)))
