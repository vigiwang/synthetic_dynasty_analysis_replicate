# ============================================================================
# 07_generate_results.R      (create_figX tradition: figures as fig*.rds)
# ----------------------------------------------------------------------------
# Builds EVERY Section-2 and Section-3 figure of the NLSY validation analysis
# directly from the processed bootstrap objects (NO Rmd on the pipeline), in
# the synthetic-dynasty paper cosmetics (Functions/plot.R: no grid, black axis
# lines, outside ticks, serif/Times), and saves each as a ggplot object:
#
#   Plot/robustness_check/NLSY/fig2_*.rds   S2.1 absolute mobility (adjusted):
#                                  UPDOWN OM UP DOWN ADJ SM EM EM_SS
#   Plot/robustness_check/NLSY/fig3_*.rds   S2.2 persistence/association:
#                                  AIM GEENENS RANKRANK LAMBDA2 AMTE HELLINGER
#   Plot/robustness_check/NLSY/fig4_*.rds   S2.3 income vs occupation vs GSS:
#                                  OM UP DOWN SM EM AIM AMTE LAMBDA2 GEENENS
#                                  HELLINGER MTE_CLS IM_CLS
#   Plot/robustness_check/NLSY/fig5_*.rds   S2.4 class-specific (adjusted):  MTE_CLS IM_CLS
#   Plot/robustness_check/NLSY/fig6_*.rds   S3.1 absolute mobility (unadjusted)
#   Plot/robustness_check/NLSY/fig7_*.rds   S3.2 persistence/association (unadjusted)
#   Plot/robustness_check/NLSY/fig8_*.rds   S3.3 class-specific (unadjusted)
#
# MAKE_PLOTLY: also saves ggplotly twins to Plot/robustness_check/NLSY/plotly/ (needs
# plotly >= 4.12 with ggplot2 >= 4.0; skipped gracefully otherwise).
# LOCAL VERIFICATION (HPC = FALSE, or NLSY_VERIFY env set): each figure's
# ggplot_build() layer data is md5-checked against
# NLSY_fun/figure_data_manifest.csv (the canonical values), and 300-dpi PNGs
# are rendered to Plot/robustness_check/NLSY/png/.
# ============================================================================

#----------Preliminaries----------#
HPC <- TRUE   # set FALSE for a local test run (enables manifest check + PNGs)
if (HPC) {
  dir_root <- "/home/weiqiw/synthetic_dynasty_analysis"
} else {
  dir_root <- "/Users/wangweiqi/Desktop/NLSY Replication Package Final"
}
if (nzchar(Sys.getenv("NLSY_PROJECT_ROOT"))) dir_root <- Sys.getenv("NLSY_PROJECT_ROOT"); setwd(dir_root)
MAKE_PLOTLY  <- TRUE
LOCAL_VERIFY <- !HPC   # manifest check + PNGs only on local runs

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(stringr); library(ggplot2)
})
root <- dir_root
source(file.path(root, "NLSY_fun", "nlsy_estimation_functions.R"))
source(file.path(root, "NLSY_fun", "nlsy_plot_functions.R"))
dir_plot <- file.path(root, "Plot", "robustness_check", "NLSY")
dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------#
#     LOAD PROCESSED OBJECTS (same as the report)       #
#-------------------------------------------------------#
.vl <- function(f) readRDS(file.path(root, "Data/NLSY_estimation", "processed", f))$validation$validation_lst
inc79  <- .vl("nlsy79_income_processed.rds");        inc97  <- .vl("nlsy97_income_processed.rds")
occ79  <- .vl("nlsy79_occ_processed.rds")
inc79u <- .vl("nlsy79_income_unadj_processed.rds");  inc97u <- .vl("nlsy97_income_unadj_processed.rds")
bench  <- readRDS(file.path(root, "Data", "main_results", "main_rst_bc.rds"))$validation_lst

INTERM <- file.path(root, "Data/NLSY_estimation/intermediate")
updown_adj        <- readRDS(file.path(INTERM, "unadjusted_updown.rds"))
ss_exchange       <- readRDS(file.path(INTERM, "ss_exchange.rds"))
ss_exchange_unadj <- readRDS(file.path(INTERM, "ss_exchange_unadj.rds"))
rr_cohort         <- readRDS(file.path(INTERM, "rankrank_by_cohort.rds"))

nlk79 <- setNames(as.character(1961:1964), 1961:1964)
nlk97 <- setNames(as.character(1980:1984), 1980:1984)
bkA   <- bench_key(c(1961:1964, 1980:1984))
bk5   <- bench_key(1961:1964)

sp_income  <- list(list(series = "NLSY income", vl = inc79, keys = nlk79),
                   list(series = "NLSY income", vl = inc97, keys = nlk97),
                   list(series = "Benchmark: GSS occupation", vl = bench, keys = bkA))
lev_income <- c("NLSY income", "Benchmark: GSS occupation")
sp_income_u  <- list(list(series = "NLSY income (unadjusted)", vl = inc79u, keys = nlk79),
                     list(series = "NLSY income (unadjusted)", vl = inc97u, keys = nlk97),
                     list(series = "Benchmark: GSS occupation", vl = bench, keys = bkA))
lev_income_u <- c("NLSY income (unadjusted)", "Benchmark: GSS occupation")
sp_incocc  <- list(list(series = "NLSY income", vl = inc79, keys = nlk79),
                   list(series = "NLSY income", vl = inc97, keys = nlk97),
                   list(series = "NLSY occupation", vl = occ79, keys = nlk79),
                   list(series = "Benchmark: GSS occupation", vl = bench, keys = bkA))
lev_incocc <- c("NLSY income", "NLSY occupation", "Benchmark: GSS occupation")
# paper (0_main.pdf) y-axis convention: measure name -> title; the y-axis shows
# the UNIT. Movement measures -> "Probability to move"; memory -> "Memory";
# mean-time-to-exit -> "Number of generations"; distinct indices keep their name.
Y_UNIT <- function(measure) {
  if (grepl("Overall|Upward|Downward|Structural|Exchange", measure)) "Probability to move"
  else if (grepl("memory|AIM|\\bIM\\b", measure, ignore.case = TRUE)) "Memory"
  else if (grepl("time to exit|MTE|AMTE", measure, ignore.case = TRUE)) "Number of generations"
  else measure
}
# measure -> y-axis UNIT only; NO on-plot title/subtitle (the measure and the
# income spec live in the figure caption + the fig*_{equalized,unequalized} filename).
P <- function(getr, specs, measure, levs)
  measure_plot(build_measure_df(getr, specs), Y_UNIT(measure), levs)

#-------------------------------------------------------#
#     BUILD THE FIGURES  (label, width, height, plot)   #
#-------------------------------------------------------#
figs <- list(
  # ---- fig2: S2.1 absolute mobility, adjusted ----
  fig2_UPDOWN = list("m-updown", 7.6, 4.8, updown_plot(inc79, inc97, bench, nlk79, nlk97, bkA)),
  fig2_OM     = list("m-overall",    7.6, 4.6, P(mg$overall,    sp_income, "Overall mobility  (t = 1)",    lev_income)),
  fig2_UP     = list("m-upward",     7.6, 4.6, P(mg$upward,     sp_income, "Upward mobility  (t = 1)",     lev_income)),
  fig2_DOWN   = list("m-downward",   7.6, 4.6, P(mg$downward,   sp_income, "Downward mobility  (t = 1)",   lev_income)),
  fig2_ADJ    = list("m-updown-adj", 7.6, 4.8, updown_adj_plot(updown_adj)),
  fig2_SM     = list("m-structural", 7.6, 4.6, P(mg$structural, sp_income, "Structural mobility  (t = 1)", lev_income)),
  fig2_EM     = list("m-exchange",   7.6, 4.6, P(mg$exchange,   sp_income, "Exchange mobility  (t = 1)",   lev_income)),
  fig2_EM_SS  = list("m-ss-exchange", 7.6, 4.6,
                     ss_exchange_plot(inc79, inc97, ss_exchange$ss79, ss_exchange$ss97, nlk79, nlk97)),
  # ---- fig3: S2.2 persistence and association, adjusted ----
  fig3_AIM       = list("m-aim",      7.6, 4.6, P(mg$aim,       sp_income, "Average individual memory (AIM, t = 1)", lev_income)),
  fig3_GEENENS   = list("m-geenens",  7.6, 4.6, P(mg$geenens,   sp_income, "Geenens D",                     lev_income)),
  fig3_RANKRANK  = list("m-rankrank", 9,   4.6, rankrank_plot(rr_cohort)),
  fig3_LAMBDA2   = list("m-lambda2",  7.6, 4.6, P(mg$lambda2,   sp_income, "Second eigenvalue",             lev_income)),
  fig3_AMTE      = list("m-amte",     7.6, 4.6, P(mg$amte,      sp_income, "Average mean time to exit (AMTE)", lev_income)),
  fig3_HELLINGER = list("m-hellinger", 7.6, 4.6, P(mg$hellinger, sp_income, "Hellinger dependence",         lev_income)),
  # ---- fig4: S2.3 NLSY income vs NLSY occupation vs GSS ----
  fig4_OM        = list("c-overall",   7.6, 4.6, P(mg$overall,    sp_incocc, "Overall mobility  (t = 1)",    lev_incocc)),
  fig4_UP        = list("c-upward",    7.6, 4.6, P(mg$upward,     sp_incocc, "Upward mobility  (t = 1)",     lev_incocc)),
  fig4_DOWN      = list("c-downward",  7.6, 4.6, P(mg$downward,   sp_incocc, "Downward mobility  (t = 1)",   lev_incocc)),
  fig4_SM        = list("c-structural", 7.6, 4.6, P(mg$structural, sp_incocc, "Structural mobility  (t = 1)", lev_incocc)),
  fig4_EM        = list("c-exchange",  7.6, 4.6, P(mg$exchange,   sp_incocc, "Exchange mobility  (t = 1)",   lev_incocc)),
  # steady-state (converged) mobility compared across the three measures
  fig4_OM_SS     = list("c-overall-ss",  7.6, 4.6, P(mg_ss("Historical Mobility"), sp_incocc, "Overall mobility  (steady state)",  lev_incocc)),
  fig4_EM_SS     = list("c-exchange-ss", 7.6, 4.6, P(mg_ss("Exchange Mobility"),   sp_incocc, "Exchange mobility  (steady state)", lev_incocc)),
  fig4_AIM       = list("c-aim",       7.6, 4.6, P(mg$aim,        sp_incocc, "Average individual memory (AIM, t = 1)", lev_incocc)),
  fig4_AMTE      = list("c-amte",      7.6, 4.6, P(mg$amte,       sp_incocc, "Average mean time to exit (AMTE)", lev_incocc)),
  fig4_LAMBDA2   = list("c-lambda2",   7.6, 4.6, P(mg$lambda2,    sp_incocc, "Second eigenvalue",            lev_incocc)),
  fig4_GEENENS   = list("c-geenens",   7.6, 4.6, P(mg$geenens,    sp_incocc, "Geenens D",                    lev_incocc)),
  fig4_HELLINGER = list("c-hellinger", 7.6, 4.6, P(mg$hellinger,  sp_incocc, "Hellinger dependence",         lev_incocc)),
  fig4_MTE_CLS   = list("c-cls-mte", 10, 4, class_compare_plot(mte_by_class, occ79, bench, nlk79, bk5,
                                                               y_title = "Number of generations")),
  fig4_IM_CLS    = list("c-cls-imm", 10, 4, class_compare_plot(imm_by_class, occ79, bench, nlk79, bk5,
                                                               y_title = "Memory")),
  # ---- fig5: S2.4 class-specific, adjusted ----
  fig5_MTE_CLS = list("cls-mte", 9, 4.6, class_by_survey_plot(mte_by_class, inc79, inc97, nlk79, nlk97,
                                                              "Number of generations")),
  fig5_IM_CLS  = list("cls-imm", 9, 4.6, class_by_survey_plot(imm_by_class, inc79, inc97, nlk79, nlk97,
                                                              "Memory")),
  # ---- fig6: S3.1 absolute mobility, unadjusted ----
  fig6_UPDOWN = list("m-updown-u", 7.6, 4.8, updown_plot(inc79u, inc97u, bench, nlk79, nlk97, bkA)),
  fig6_OM     = list("m-overall-u",    7.6, 4.6, P(mg$overall,    sp_income_u, "Overall mobility  (t = 1)",    lev_income_u)),
  fig6_UP     = list("m-upward-u",     7.6, 4.6, P(mg$upward,     sp_income_u, "Upward mobility  (t = 1)",     lev_income_u)),
  fig6_DOWN   = list("m-downward-u",   7.6, 4.6, P(mg$downward,   sp_income_u, "Downward mobility  (t = 1)",   lev_income_u)),
  fig6_SM     = list("m-structural-u", 7.6, 4.6, P(mg$structural, sp_income_u, "Structural mobility  (t = 1)", lev_income_u)),
  fig6_EM     = list("m-exchange-u",   7.6, 4.6, P(mg$exchange,   sp_income_u, "Exchange mobility  (t = 1)",   lev_income_u)),
  fig6_EM_SS  = list("m-ss-exchange-u", 7.6, 4.6,
                     ss_exchange_plot(inc79u, inc97u, ss_exchange_unadj$ss79, ss_exchange_unadj$ss97, nlk79, nlk97)),
  # ---- fig7: S3.2 persistence and association, unadjusted ----
  fig7_AIM       = list("m-aim-u",      7.6, 4.6, P(mg$aim,       sp_income_u, "Average individual memory (AIM, t = 1)", lev_income_u)),
  fig7_GEENENS   = list("m-geenens-u",  7.6, 4.6, P(mg$geenens,   sp_income_u, "Geenens D",             lev_income_u)),
  fig7_LAMBDA2   = list("m-lambda2-u",  7.6, 4.6, P(mg$lambda2,   sp_income_u, "Second eigenvalue",     lev_income_u)),
  fig7_AMTE      = list("m-amte-u",     7.6, 4.6, P(mg$amte,      sp_income_u, "Average mean time to exit (AMTE)", lev_income_u)),
  fig7_HELLINGER = list("m-hellinger-u", 7.6, 4.6, P(mg$hellinger, sp_income_u, "Hellinger dependence", lev_income_u)),
  # ---- fig8: S3.3 class-specific, unadjusted ----
  fig8_MTE_CLS = list("cls-mte-u", 9, 4.6, class_by_survey_plot(mte_by_class, inc79u, inc97u, nlk79, nlk97,
                                                                "Number of generations")),
  fig8_IM_CLS  = list("cls-imm-u", 9, 4.6, class_by_survey_plot(imm_by_class, inc79u, inc97u, nlk79, nlk97,
                                                                "Memory"))
)

#-----------------------------#
#  Style + save the figures:  #
#-----------------------------#
# income-specification tag: fig2-5 = size-equalized (adjusted); fig6-8 =
# unequalized (family-size-unadjusted). Tag the FILENAME only (no on-plot text).
spec_tag <- function(nm) if (grepl("^fig[678]", nm)) "unequalized" else "equalized"
for (nm in names(figs)) {
  p <- style_dynasty(figs[[nm]][[4]]) + labs(title = NULL, subtitle = NULL)
  figs[[nm]][[4]] <- p
  saveRDS(p, file.path(dir_plot, paste0(nm, "_", spec_tag(nm), ".rds")))
}
cat(sprintf("[07] %d figures saved to %s as fig*_{equalized,unequalized}.rds (styled ggplot objects)\n",
            length(figs), dir_plot))

# ---- optional plotly twins --------------------------------------------------
if (MAKE_PLOTLY) {
  if (dir.exists(file.path(root, "rlibs"))) .libPaths(c(file.path(root, "rlibs"), .libPaths()))
  if (requireNamespace("plotly", quietly = TRUE) && packageVersion("plotly") >= "4.12.0") {
    dir.create(file.path(dir_plot, "plotly"), showWarnings = FALSE)
    n_pl <- 0L
    for (nm in names(figs)) {
      pl <- tryCatch(suppressWarnings(suppressMessages(plotly::ggplotly(figs[[nm]][[4]]))),
                     error = function(e) NULL)
      if (!is.null(pl)) { saveRDS(pl, file.path(dir_plot, "plotly", paste0(nm, ".rds"))); n_pl <- n_pl + 1L }
    }
    cat(sprintf("[07] %d plotly twins saved to %s/plotly\n", n_pl, dir_plot))
  } else cat("[07] plotly >= 4.12 not available - twins skipped",
             "(install.packages('plotly', lib = 'rlibs'))\n")
}

#-----------------------------#
#  Local verification only:   #
#-----------------------------#
if (LOCAL_VERIFY) {
  manifest <- read.csv(file.path(root, "NLSY_fun", "figure_data_manifest.csv"))
  n_ok <- 0L
  dir.create(file.path(dir_plot, "png"), showWarnings = FALSE)
  for (nm in names(figs)) {
    lab <- figs[[nm]][[1]]; fw <- figs[[nm]][[2]]; fh <- figs[[nm]][[3]]
    h <- digest::digest(ggplot_build(figs[[nm]][[4]])$data, algo = "md5")
    ref <- manifest$data_md5[manifest$label == lab]
    if (length(ref) == 1 && h == ref) n_ok <- n_ok + 1L else message("  DATA MISMATCH: ", nm, " (", lab, ")")
    grDevices::png(file.path(dir_plot, "png", paste0(nm, "_", spec_tag(nm), ".png")),
                   width = fw, height = fh, units = "in", res = 300)
    print(figs[[nm]][[4]]); dev.off()
  }
  cat(sprintf("[07] verification: %d/%d figures data-identical to canonical; PNGs in %s/png\n",
              n_ok, length(figs), dir_plot))
}
