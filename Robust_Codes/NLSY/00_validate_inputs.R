#!/usr/bin/env Rscript
# ===========================================================================
# 00_validate_inputs.R  — validate canonical inputs before any construction.
# Validates; it does NOT silently repair. Fails loudly on a missing/incomplete
# input. Run first (master --stage=validate).
# ===========================================================================
root <- Sys.getenv("NLSY_PROJECT_ROOT", unset = getwd())
suppressPackageStartupMessages({library(dplyr)})
V <- file.path(root, "Data", "NLSY", "_validation"); dir.create(file.path(V, "logs"), recursive = TRUE, showWarnings = FALSE)
fail <- function(msg) stop("[00_validate] ", msg, call. = FALSE)
ok   <- function(msg) cat("  OK:", msg, "\n")

cat("== 00_validate_inputs ==\nroot:", root, "\n")

## required canonical raw + supporting files (3 datasets; NLSY97 occ not in this pipeline)
req <- c(
  "Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_raw.csv",
  "Data/NLSY/NLSY79/NLSY79_inc/NLSY79_inc_complete_dyads_preweight.csv",
  "Data/NLSY/NLSY79/NLSY79_occ/NLSY79_occ_complete_dyads_preweight.csv",
  "Data/NLSY/NLSY97/NLSY97_inc/DPCERG3A086NBEA.csv",
  "Data/NLSY/NLSY79/NLSY79_occ/Crosswalks/crosswalk_manifest.csv",
  "Data/main_results/main_rst_bc.rds")
for (f in req) if (!file.exists(file.path(root, f))) fail(paste("missing required input:", f)) else ok(f)

## NLSY97 consolidated raw: key + income + relationship-roster columns present
d <- data.table::fread(file.path(root, "Data/NLSY/NLSY97/NLSY97_inc/NLSY97_inc_raw.csv"), nrows = 5)
needcol <- c("R0000100","R1204500","R1204600","R2563300","R7227800","T0014100","R1205400","R1315800","R2416300","R6919700")
miss <- setdiff(needcol, names(d)); if (length(miss)) fail(paste("NLSY97 raw missing columns:", paste(miss, collapse=","))) else ok("NLSY97 raw has income + RELY roster columns")

## PCE coverage for all income reference years 1978-2023
pce <- data.table::fread(file.path(root, "Data/NLSY/NLSY97/NLSY97_inc/DPCERG3A086NBEA.csv"))
yrs <- as.integer(substr(as.character(pce[[1]]), 1, 4))
if (!all(1978:2023 %in% yrs)) fail("PCE index does not cover 1978-2023") else ok("PCE covers 1978-2023")

## Functions/ present (frozen benchmark core)
if (length(list.files(file.path(root, "Functions"), pattern="\\.R$")) < 5) fail("Functions/ incomplete") else ok("Functions/ present")

cat("[00_validate] ALL INPUT CHECKS PASSED\n")
