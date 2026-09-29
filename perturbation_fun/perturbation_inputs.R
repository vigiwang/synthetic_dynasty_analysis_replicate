#-------------------------------------------------------------------------------
# perturbation_fun/perturbation_inputs.R
#
# Configuration, measure dictionary, input loading/validation, and the fixed
# sparse-cell mask for the Appendix J perturbation pipeline (Tables J.1-J.3,
# Figures J.1-J.3). All empirical inputs are the stored Midway main-analysis
# objects under Data/main_results/ (read-only).
#-------------------------------------------------------------------------------

#----------------------------- configuration ----------------------------------#

get_perturbation_config <- function() {
  list(
    # verified production input filenames under Data/main_results/
    input_files = c(
      raw      = "gss_fc_occ10_5class.rds",
      baseline = "main_rst_baseline.rds",
      boot     = "main_rst_boot.rds",
      bc       = "main_rst_bc.rds"
    ),
    # paper sample contract (Appendix J is attached to the main GSS analysis)
    age_min = 25L, age_max = 55L,
    sample_n = 24265L,
    cohort_n_range = c(184L, 788L),
    cohorts = 1945:1990,
    n_classes = 5L,
    class_ids = paste0("C", 1:5),
    class_labels = c(
      "Higher-grade professionals, administrators, and managers (EGP I)",
      "Lower-grade professionals, administrators, and managers (EGP II)",
      "Higher-grade routine non-manual workers, technicians, and supervisors (EGP V + IIIa)",
      "Higher-grade manual workers, lower-grade technicians, and farmers (EGP VI + IVc)",
      "Lower-grade routine non-manual and manual workers (EGP VII + IIIb)"
    ),
    # eras (Table J.1 summaries) and trend windows (Tables J.2-J.3, figures)
    early_era   = c(1945L, 1974L),
    late_era    = c(1975L, 1990L),
    full_window = c(1945L, 1990L),
    late_window = c(1975L, 1990L),
    late_start  = 1975L,
    contribution_divider_full = 1975L,
    contribution_divider_late = 1983L,
    # perturbation design
    kappa = 1,
    sparse_threshold = 5L,
    n_boot = 2000L,
    # production steady-state stop rule (paper convention)
    steady_state_tol = 0.01,
    # single numerical tolerance for identity guards
    tol = 1e-10,
    # mechanism labels; Table J.3 lists Augmentation before Redistribution
    mechanisms = c("Redistribution", "Augmentation"),
    table_J3_mechanism_order = c("Augmentation", "Redistribution"),
    # paper artifact filenames. Figures are lettered K (current manuscript
    # appendix lettering) and split into four measure blocks per figure
    # family (a movement, b MTE, c memory, d Class-4) so each sub-figure
    # renders large enough at text width; the canonical rds keeps the
    # original figure_J*_data key names.
    artifacts = list(
      results_rds = "perturbation_results.rds",
      table_J1 = "table_J1_sparse_cell_design.tex",
      table_J2 = "table_J2_cohort_trend_results.tex",
      table_J3 = "table_J3_slope_bootstrap_calibration.tex",
      figures = c(
        K1a = "figure_K1a_trend_changes_movement",
        K1b = "figure_K1b_trend_changes_mte",
        K1c = "figure_K1c_trend_changes_memory",
        K1d = "figure_K1d_trend_changes_class4",
        K2a = "figure_K2a_slope_calibration_movement",
        K2b = "figure_K2b_slope_calibration_mte",
        K2c = "figure_K2c_slope_calibration_memory",
        K2d = "figure_K2d_slope_calibration_class4",
        K3a = "figure_K3a_cohort_contributions_movement",
        K3b = "figure_K3b_cohort_contributions_mte",
        K3c = "figure_K3c_cohort_contributions_memory",
        K3d = "figure_K3d_cohort_contributions_class4"
      )
    )
  )
}

#--------------------------- measure dictionary -------------------------------#
# The 20 paper measures in the exact Appendix J order, with paper-facing
# labels, groups, and per-decade display units. decade_factor converts the
# per-decade OLS slope (decade-centered regressor) into display units:
# proportions -> percentage points (x100); generations and index units -> x1.

get_perturbation_measure_dictionary <- function() {
  dplyr::tibble(
    measure_id = c("historical", "structural", "exchange", "upward", "downward",
                   "AMTE", paste0("MTE_C", 1:5),
                   "AIM", paste0("IM_C", 1:5),
                   "SSM", "ss_C4_share", "mfp_to_C4"),
    label = c("Overall mobility (OM$_1$)",
              "Structural mobility (SM$_1$)",
              "Exchange mobility (EM$_1$)",
              "Upward mobility", "Downward mobility",
              "Average mean time to exit (AMTE)",
              sprintf("MTE, Class %d", 1:5),
              "Aggregate intergenerational memory (AIM$_1$)",
              sprintf("IM, Class %d", 1:5),
              "Steady-state mobility (SSM)",
              "Steady-state share, Class 4",
              "Mean first passage time into Class 4"),
    label_plain = c("Overall mobility (OM1)", "Structural mobility (SM1)",
                    "Exchange mobility (EM1)", "Upward mobility", "Downward mobility",
                    "Average mean time to exit (AMTE)",
                    sprintf("MTE, Class %d", 1:5),
                    "Aggregate intergenerational memory (AIM1)",
                    sprintf("IM, Class %d", 1:5),
                    "Steady-state mobility (SSM)",
                    "Steady-state share, Class 4",
                    "Mean first passage time into Class 4"),
    group = c(rep("Measures of movement", 5),
              rep("Mean time to exit", 6),
              rep("Intergenerational memory", 6),
              rep("Steady-state quantities", 2),
              "Mean first passage time"),
    unit = c(rep("proportion", 5),
             rep("generations", 6),
             rep("index", 6),
             "proportion", "proportion", "generations"),
    unit_label = c(rep("pp per decade", 5),
                   rep("generations per decade", 6),
                   rep("index units per decade", 6),
                   "pp per decade", "pp per decade", "generations per decade"),
    decade_factor = c(rep(100, 5), rep(1, 6), rep(1, 6), 100, 100, 1),
    display_order = 1:20
  )
}

#------------------------------ input loading ---------------------------------#

load_perturbation_inputs <- function(cfg, dir_main_results) {
  paths <- file.path(dir_main_results, cfg$input_files)
  names(paths) <- names(cfg$input_files)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing production inputs under Data/main_results/: ",
         paste(basename(missing), collapse = ", "))
  }
  raw <- readRDS(paths[["raw"]])
  raw$cohort   <- as.integer(haven::zap_labels(raw$cohort))
  raw$age      <- as.integer(haven::zap_labels(raw$age))
  raw$status_p <- as.integer(haven::zap_labels(raw$status_p))
  raw$status_c <- as.integer(haven::zap_labels(raw$status_c))
  list(raw = raw,
       baseline = readRDS(paths[["baseline"]]),
       boot = readRDS(paths[["boot"]]),
       bc_path = paths[["bc"]],
       input_paths = paths)
}

#----------------------- standardization + validation -------------------------#
# Reconstruct the production smoothed matrices/marginals and the raw
# cohort-by-origin-by-destination dyad counts; confirm the paper sample.

standardize_perturbation_inputs <- function(inputs, cfg) {
  n_cls <- cfg$n_classes
  cls <- cfg$class_ids

  ## smoothed production matrices: P_c = row-normalized tm_num; mu0 = row shares
  smoothed <- lapply(seq_along(cfg$cohorts), function(i) {
    num <- matrix(inputs$baseline[[i]]$tm_num$mean, n_cls, n_cls,
                  byrow = TRUE, dimnames = list(cls, cls))
    rs <- rowSums(num)
    list(cohort = cfg$cohorts[i], tm_num = num, P = num / rs, mu0 = rs / sum(rs))
  })
  names(smoothed) <- as.character(cfg$cohorts)

  ## raw dyad-count matrices (paper analytic sample: ages 25-55, cohorts 1945-1990)
  gg <- inputs$raw[inputs$raw$age >= cfg$age_min & inputs$raw$age <= cfg$age_max &
                   inputs$raw$cohort >= min(cfg$cohorts) &
                   inputs$raw$cohort <= max(cfg$cohorts) &
                   inputs$raw$status_p %in% 1:n_cls &
                   inputs$raw$status_c %in% 1:n_cls, ]
  raw_mats <- lapply(cfg$cohorts, function(ch) {
    ci <- gg[gg$cohort == ch, ]
    tb <- table(factor(ci$status_p, 1:n_cls), factor(ci$status_c, 1:n_cls))
    m <- matrix(as.integer(tb), n_cls, n_cls, dimnames = list(cls, cls))
    list(cohort = ch, N_raw = m, n_ic = rowSums(m), N_c = sum(m))
  })
  names(raw_mats) <- as.character(cfg$cohorts)

  ## ---- paper sample contract checks (stop on mismatch, never silent) ----
  n_total <- sum(vapply(raw_mats, function(r) r$N_c, numeric(1)))
  if (n_total != cfg$sample_n) {
    stop("Production sample mismatch: reconstructed n = ", n_total,
         " but the paper reports n = ", cfg$sample_n,
         " (ages ", cfg$age_min, "-", cfg$age_max, ", cohorts ",
         min(cfg$cohorts), "-", max(cfg$cohorts), ").")
  }
  n_by_cohort <- vapply(raw_mats, function(r) r$N_c, numeric(1))
  if (min(n_by_cohort) != cfg$cohort_n_range[1] ||
      max(n_by_cohort) != cfg$cohort_n_range[2]) {
    stop("Cohort sample-size range mismatch: observed [", min(n_by_cohort), ", ",
         max(n_by_cohort), "], paper reports [",
         cfg$cohort_n_range[1], ", ", cfg$cohort_n_range[2], "].")
  }
  ## five-class coding check: Class 4 must be EGP VI + IVc in the raw labels
  lab4 <- sort(unique(as.character(inputs$raw$label_p[inputs$raw$status_p == 4])))
  if (!all(grepl("IVc|VI$|^VI", lab4))) {
    stop("Class 4 coding mismatch: expected EGP VI + IVc, found labels: ",
         paste(lab4, collapse = " + "))
  }
  ## bootstrap contract: 2,000 replicates, identical cohort labels
  if (length(inputs$boot) != cfg$n_boot) {
    stop("Bootstrap replicate count is ", length(inputs$boot),
         "; the paper analysis uses ", cfg$n_boot, ".")
  }
  if (!identical(as.integer(names(inputs$boot[[1]])), cfg$cohorts)) {
    stop("Bootstrap cohort labels do not match 1945-1990.")
  }

  list(smoothed = smoothed, raw_mats = raw_mats,
       n_by_cohort = n_by_cohort, n_total = n_total)
}

#--------------------------- fixed sparse-cell mask ---------------------------#
# A_ijc = 1{N_raw_ijc <= threshold}, constructed ONCE from the deterministic
# original raw counts and held fixed across mechanisms, kappa values, and
# all bootstrap replicates.

build_fixed_sparse_mask <- function(raw_mats, cfg) {
  masks <- lapply(raw_mats, function(r) r$N_raw <= cfg$sparse_threshold)
  inventory <- dplyr::bind_rows(lapply(names(raw_mats), function(nm) {
    r <- raw_mats[[nm]]; A <- masks[[nm]]
    m_ic <- rowSums(A)
    gr <- expand.grid(j = seq_len(cfg$n_classes), i = seq_len(cfg$n_classes))
    dplyr::tibble(
      cohort = r$cohort,
      origin = cfg$class_ids[gr$i], destination = cfg$class_ids[gr$j],
      raw_count = r$N_raw[cbind(gr$i, gr$j)],
      sparse_cell = A[cbind(gr$i, gr$j)],
      affected_row = m_ic[gr$i] > 0,
      m_ic = m_ic[gr$i], K_c = sum(A), n_ic = r$n_ic[gr$i])
  }))
  list(masks = masks, inventory = inventory)
}
