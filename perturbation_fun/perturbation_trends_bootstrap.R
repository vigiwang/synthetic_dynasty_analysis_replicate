#-------------------------------------------------------------------------------
# perturbation_fun/perturbation_trends_bootstrap.R
#
# Trend estimation, paired bootstrap, stationary-TV design calibration
# (Table J.1), slope calibration (Table J.3 / Figure J.2), and cohort
# contributions (Figure J.3) for the Appendix J pipeline. Reuses the
# 2,000 stored GAM-smoothed bootstrap replicates: no resampling, no GAM
# refits, replicate identity preserved, paired within replicate.
#-------------------------------------------------------------------------------

#------------------------------ trend fitting ---------------------------------#
# Paper trend: theta_c = alpha + beta * x_c + e, with the centered
# per-decade regressor x_c = (c - mean(c in W)) / 10; beta is per decade.

fit_paper_trend <- function(theta, cohorts) {
  x <- (cohorts - mean(cohorts)) / 10
  sum(x * (theta - mean(theta))) / sum(x^2)
}

# OLS slope weight h_c implied by the centered per-decade regressor
# (beta = sum_c h_c * theta_c)
trend_slope_weights <- function(cohorts) {
  x <- (cohorts - mean(cohorts)) / 10
  x / sum(x^2)
}

window_cohorts <- function(cfg, window) {
  if (window == "full") cfg$full_window[1]:cfg$full_window[2]
  else cfg$late_window[1]:cfg$late_window[2]
}

#-------------------------- deterministic trends ------------------------------#

calculate_deterministic_trends <- function(measures_long, dictionary, cfg) {
  out <- list()
  for (w in c("full", "late")) {
    cw <- window_cohorts(cfg, w)
    for (m in dictionary$measure_id) {
      base <- measures_long[measures_long$mechanism == "Baseline" &
                            measures_long$measure_id == m &
                            measures_long$cohort %in% cw, ]
      base <- base[order(base$cohort), ]
      b0 <- fit_paper_trend(base$value, base$cohort)
      for (mech in cfg$mechanisms) {
        pert <- measures_long[measures_long$mechanism == mech &
                              measures_long$measure_id == m &
                              measures_long$cohort %in% cw, ]
        pert <- pert[order(pert$cohort), ]
        stopifnot(identical(pert$cohort, base$cohort))
        b1 <- fit_paper_trend(pert$value, pert$cohort)
        ## exact identity: difference of slopes = slope of differences
        b_diff <- fit_paper_trend(pert$value - base$value, base$cohort)
        stopifnot(abs((b1 - b0) - b_diff) < cfg$tol)
        out[[length(out) + 1]] <- dplyr::tibble(
          window = w, measure_id = m, mechanism = mech,
          beta_base = b0, beta_pert = b1, delta_beta = b1 - b0)
      }
    }
  }
  dplyr::bind_rows(out)
}

#--------------------------- paired bootstrap ---------------------------------#
# Within each stored replicate: baseline and perturbed slopes from the SAME
# replicate; paired delta beta; u_total = replicate baseline slope minus the
# deterministic baseline slope; direction_retained relative to the sign of
# the DETERMINISTIC baseline slope.

calculate_paired_bootstrap_trends <- function(boot, std, masks, dictionary, cfg,
                                              deterministic_trends,
                                              log_every = 200) {
  cls_cols <- paste0("Class", 1:cfg$n_classes)
  cohorts <- cfg$cohorts; n_coh <- length(cohorts)
  late_idx <- which(cohorts >= cfg$late_start)
  mids <- dictionary$measure_id; M <- length(mids)
  n_rep <- length(boot)
  n_ic_list <- lapply(std$raw_mats, function(r) r$n_ic)
  N_c_vec <- vapply(std$raw_mats, function(r) r$N_c, numeric(1))

  det_base <- deterministic_trends[!duplicated(deterministic_trends[,
    c("window", "measure_id")]), c("window", "measure_id", "beta_base")]
  det_base_full <- stats::setNames(
    det_base$beta_base[det_base$window == "full"],
    det_base$measure_id[det_base$window == "full"])[mids]
  det_base_late <- stats::setNames(
    det_base$beta_base[det_base$window == "late"],
    det_base$measure_id[det_base$window == "late"])[mids]

  variants <- c("Baseline", cfg$mechanisms)
  mk <- function() matrix(NA_real_, n_rep, M, dimnames = list(NULL, mids))
  sf <- lapply(variants, function(v) mk()); names(sf) <- variants
  sl <- lapply(variants, function(v) mk()); names(sl) <- variants

  for (b in seq_len(n_rep)) {
    rb <- boot[[b]]
    vals <- lapply(variants, function(v)
      matrix(NA_real_, n_coh, M, dimnames = list(NULL, mids)))
    names(vals) <- variants
    for (ci in seq_len(n_coh)) {
      MM <- as.matrix(rb[[ci]]$tm_num[, cls_cols]); rs <- rowSums(MM)
      P_b <- MM / rs; mu0_b <- rs / sum(rs)
      vals$Baseline[ci, ] <-
        calculate_reported_paper_measures(P_b, mu0_b, cfg$steady_state_tol)[mids]
      for (mech in cfg$mechanisms) {
        pb <- apply_perturbation_to_bootstrap_replicate(
          mech, P_b, mu0_b, n_ic_list[[ci]], N_c_vec[ci], masks[[ci]], cfg$kappa)
        vals[[mech]][ci, ] <- calculate_reported_paper_measures(
          pb$transition_matrix, pb$father_origin_marginal,
          cfg$steady_state_tol)[mids]
      }
    }
    for (v in variants) for (m in mids) {
      sf[[v]][b, m] <- fit_paper_trend(vals[[v]][, m], cohorts)
      sl[[v]][b, m] <- fit_paper_trend(vals[[v]][late_idx, m], cohorts[late_idx])
    }
    if (b %% log_every == 0) message("    paired bootstrap replicate ", b, "/", n_rep)
  }

  out <- list()
  for (w in c("full", "late")) {
    S <- if (w == "full") sf else sl
    d0 <- if (w == "full") det_base_full else det_base_late
    for (mech in cfg$mechanisms) for (m in mids) {
      bb <- S$Baseline[, m]; bp <- S[[mech]][, m]
      out[[length(out) + 1]] <- dplyr::tibble(
        replicate_id = seq_len(n_rep), window = w, measure_id = m,
        mechanism = mech, beta_base = bb, beta_pert = bp,
        delta_beta = bp - bb,
        direction_retained = sign(bp) == sign(d0[[m]]),
        u_total = bb - d0[[m]])
    }
  }
  dplyr::bind_rows(out)
}

#----------------------- targeted slope displacement --------------------------#
# Mechanism-specific targeted benchmarks (validated production designs):
#   Redistribution — flagged cells take replicate transition-law values,
#     donors rescale within the row, father-origin marginal held fixed.
#   Augmentation  — flagged cells take replicate joint pseudo-count values
#     on the raw total-dyad scale; the father-origin marginal and the
#     transition law derive from that targeted joint table.

calculate_targeted_slope_bootstrap <- function(boot, std, masks, dictionary, cfg,
                                               deterministic_trends,
                                               log_every = 200) {
  cls_cols <- paste0("Class", 1:cfg$n_classes)
  cohorts <- cfg$cohorts; n_coh <- length(cohorts)
  late_idx <- which(cohorts >= cfg$late_start)
  mids <- dictionary$measure_id; M <- length(mids)
  n_rep <- length(boot)
  N_c_vec <- vapply(std$raw_mats, function(r) r$N_c, numeric(1))
  C0 <- lapply(names(std$smoothed), function(nm)
    (N_c_vec[[nm]] * std$smoothed[[nm]]$mu0) * std$smoothed[[nm]]$P)
  names(C0) <- names(std$smoothed)

  det_base <- deterministic_trends[!duplicated(deterministic_trends[,
    c("window", "measure_id")]), c("window", "measure_id", "beta_base")]

  mk <- function() matrix(NA_real_, n_rep, M, dimnames = list(NULL, mids))
  store <- list(Redistribution_full = mk(), Redistribution_late = mk(),
                Augmentation_full = mk(), Augmentation_late = mk())

  for (b in seq_len(n_rep)) {
    rb <- boot[[b]]
    vr <- matrix(NA_real_, n_coh, M, dimnames = list(NULL, mids))
    va <- matrix(NA_real_, n_coh, M, dimnames = list(NULL, mids))
    for (ci in seq_len(n_coh)) {
      nm <- names(std$smoothed)[ci]
      s <- std$smoothed[[nm]]; A <- masks[[ci]]
      MM <- as.matrix(rb[[ci]]$tm_num[, cls_cols]); rs <- rowSums(MM)
      P_b <- MM / rs; mu0_b <- rs / sum(rs)
      if (!any(A)) {  # K_c = 0: targeted projections equal the baseline
        v0 <- calculate_reported_paper_measures(s$P, s$mu0, cfg$steady_state_tol)[mids]
        vr[ci, ] <- v0; va[ci, ] <- v0
        next
      }
      ## redistribution: transition-law-only targeted projection
      Pc <- s$P; ok <- TRUE
      for (i in seq_len(cfg$n_classes)) {
        T_i <- which(A[i, ]); if (!length(T_i)) next
        D_i <- setdiff(seq_len(cfg$n_classes), T_i)
        if (!length(D_i)) { ok <- FALSE; break }
        pv <- P_b[i, T_i]
        Pc[i, T_i] <- pv
        Pc[i, D_i] <- s$P[i, D_i] * (1 - sum(pv)) / (1 - sum(s$P[i, T_i]))
      }
      if (ok && !any(Pc < -1e-12))
        vr[ci, ] <- calculate_reported_paper_measures(Pc, s$mu0,
                                                      cfg$steady_state_tol)[mids]
      ## augmentation: targeted joint-table projection
      C_b <- (N_c_vec[ci] * mu0_b) * P_b
      C_t <- C0[[nm]]; flag <- which(A)
      C_t[flag] <- C_b[flag]
      n_t <- rowSums(C_t)
      if (all(n_t > 0)) {
        va[ci, ] <- calculate_reported_paper_measures(
          C_t / n_t, n_t / sum(n_t), cfg$steady_state_tol)[mids]
      }
    }
    for (m in mids) {
      store$Redistribution_full[b, m] <- fit_paper_trend(vr[, m], cohorts)
      store$Redistribution_late[b, m] <- fit_paper_trend(vr[late_idx, m], cohorts[late_idx])
      store$Augmentation_full[b, m] <- fit_paper_trend(va[, m], cohorts)
      store$Augmentation_late[b, m] <- fit_paper_trend(va[late_idx, m], cohorts[late_idx])
    }
    if (b %% log_every == 0) message("    targeted slope bootstrap replicate ", b, "/", n_rep)
  }

  out <- list()
  for (w in c("full", "late")) for (mech in cfg$mechanisms) {
    S <- store[[paste0(mech, "_", w)]]
    for (m in mids) {
      d0 <- det_base$beta_base[det_base$window == w & det_base$measure_id == m]
      out[[length(out) + 1]] <- dplyr::tibble(
        replicate_id = seq_len(nrow(S)), window = w, measure_id = m,
        mechanism = mech, u_target = S[, m] - d0)
    }
  }
  dplyr::bind_rows(out)
}

#------------------------ paired delta-beta summary ---------------------------#
# Unadjusted paired percentile interval (2.5th/97.5th) and P(direction
# retained) as a sign-stability frequency (not a p-value).

summarize_paired_delta_beta <- function(bootstrap_trends) {
  dplyr::summarise(
    dplyr::group_by(bootstrap_trends, window, measure_id, mechanism),
    delta_beta_lo = stats::quantile(delta_beta, 0.025, names = FALSE),
    delta_beta_hi = stats::quantile(delta_beta, 0.975, names = FALSE),
    p_direction_retained = mean(direction_retained),
    n_replicates = dplyr::n(),
    .groups = "drop")
}

#--------------------- Table J.1 design calibration ---------------------------#

tv_distance <- function(a, b) 0.5 * sum(abs(a - b))

# stationary TV of the redistribution direction at magnitude kappa
.redistribution_direction <- function(P, n_ic, mask) {
  D <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P))
  for (i in seq_len(nrow(P))) {
    T_i <- which(mask[i, ]); m <- length(T_i)
    if (m == 0 || m == ncol(P)) next
    S_T <- sum(P[i, T_i]); D_i <- setdiff(seq_len(ncol(P)), T_i)
    D[i, T_i] <- 1 / n_ic[i]
    D[i, D_i] <- -(m / n_ic[i]) * P[i, D_i] / (1 - S_T)
  }
  D
}

# mechanism-specific kappa at which stationary TV reaches `target`
# (validated interpolation/root-finding conventions)
.kappa_at_stationary_tv <- function(mechanism, s, r, mask, target, n_grid = 80,
                                    aug_search_max = 4096) {
  pi_base <- stationary_distribution_exact(s$P)
  if (mechanism == "Redistribution") {
    D <- .redistribution_direction(s$P, r$n_ic, mask)
    km <- redistribution_feasibility_ceiling(s$P, r$n_ic, mask)$ceiling
    grid <- seq(0, km, length.out = n_grid)
    tvv <- vapply(grid, function(k)
      tv_distance(stationary_distribution_exact(s$P + k * D), pi_base), numeric(1))
    hit <- which(tvv >= target)
    if (!length(hit)) return(NA_real_)
    i <- hit[1]
    if (i == 1) return(0)
    tryCatch(stats::uniroot(function(k)
      tv_distance(stationary_distribution_exact(s$P + k * D), pi_base) - target,
      c(grid[i - 1], grid[i]))$root, error = function(e) grid[i])
  } else {
    f <- function(k) tv_distance(stationary_distribution_exact(
      apply_augmentation_perturbation(s$P, s$mu0, r$N_c * s$mu0, mask,
                                      k)$transition_matrix), pi_base) - target
    if (f(0) >= 0) return(0)
    k_hi <- 8
    while (f(k_hi) < 0 && k_hi < aug_search_max) k_hi <- k_hi * 2
    if (f(k_hi) < 0) return(NA_real_)
    stats::uniroot(f, c(0, k_hi))$root
  }
}

calculate_stationary_design_calibration <- function(std, masks, mask_inventory,
                                                    boot, cfg, log_every = 500) {
  cls_cols <- paste0("Class", 1:cfg$n_classes)
  cohorts <- cfg$cohorts; n_coh <- length(cohorts); n_rep <- length(boot)

  ## perturbation-free total-system stationary-TV bootstrap per cohort
  pi_base <- lapply(std$smoothed, function(s) stationary_distribution_exact(s$P))
  tv_draws <- matrix(NA_real_, n_rep, n_coh)
  for (b in seq_len(n_rep)) {
    rb <- boot[[b]]
    for (ci in seq_len(n_coh)) {
      MM <- as.matrix(rb[[ci]]$tm_num[, cls_cols])
      tv_draws[b, ci] <- tv_distance(
        stationary_distribution_exact(MM / rowSums(MM)), pi_base[[ci]])
    }
    if (b %% log_every == 0) message("    total-system TV bootstrap replicate ", b, "/", n_rep)
  }
  total_q50 <- apply(tv_draws, 2, stats::median)
  total_q95 <- apply(tv_draws, 2, stats::quantile, probs = 0.95, names = FALSE)

  out <- list()
  for (ci in seq_len(n_coh)) {
    nm <- as.character(cohorts[ci])
    s <- std$smoothed[[nm]]; r <- std$raw_mats[[nm]]; A <- masks[[ci]]
    K_c <- sum(A); m_rows <- sum(rowSums(A) > 0)
    mu1_base <- as.numeric(s$mu0 %*% s$P)
    # hoist scalars: inside dplyr::tibble() the freshly created columns
    # would mask the outer total_q50/total_q95 vectors
    q50_ci <- total_q50[ci]; q95_ci <- total_q95[ci]
    for (mech in cfg$mechanisms) {
      p <- apply_perturbation_to_cohort(mech, s$P, s$mu0, r$n_ic, r$N_c, A, cfg$kappa)
      tv1 <- tv_distance(stationary_distribution_exact(p$transition_matrix),
                         pi_base[[ci]])
      kb <- if (K_c == 0) NA_real_ else
        .kappa_at_stationary_tv(mech, s, r, A, q95_ci)
      out[[length(out) + 1]] <- dplyr::tibble(
        cohort = cohorts[ci], mechanism = mech,
        era = ifelse(cohorts[ci] < cfg$late_start, "early", "late"),
        n_sparse_cells = K_c, n_affected_rows = m_rows,
        budget = cfg$kappa * K_c,
        stationary_tv_k1 = tv1,
        tv_mu0 = tv_distance(p$father_origin_marginal, s$mu0),
        tv_mu1 = tv_distance(p$child_destination_marginal, mu1_base),
        total_q50 = q50_ci, total_q95 = q95_ci,
        R_total_stationary = tv1 / q95_ci,
        kappa_boot95_total = kb,
        capped = p$capped)
    }
  }
  dplyr::bind_rows(out)
}

#----------------- Table J.3 / Figure J.2 slope calibration -------------------#

calculate_slope_calibration <- function(bootstrap_trends, deterministic_trends) {
  agg <- dplyr::summarise(
    dplyr::group_by(bootstrap_trends, window, measure_id, mechanism),
    total_beta_lo = stats::quantile(beta_base, 0.025, names = FALSE),
    total_beta_hi = stats::quantile(beta_base, 0.975, names = FALSE),
    q95_total = stats::quantile(abs(u_total), 0.95, names = FALSE),
    q95_target = stats::quantile(abs(u_target), 0.95, names = FALSE),
    u_total_lo = stats::quantile(u_total, 0.025, names = FALSE),
    u_total_hi = stats::quantile(u_total, 0.975, names = FALSE),
    u_target_lo = stats::quantile(u_target, 0.025, names = FALSE),
    u_target_hi = stats::quantile(u_target, 0.975, names = FALSE),
    .groups = "drop")
  j <- dplyr::left_join(deterministic_trends, agg,
                        by = c("window", "measure_id", "mechanism"))
  dplyr::mutate(j,
    R_total = abs(delta_beta) / q95_total,
    R_target = abs(delta_beta) / q95_target,
    outside_total_95 = beta_pert < total_beta_lo | beta_pert > total_beta_hi)
}

#--------------------- Figure J.3 cohort contributions ------------------------#
# q_delta,c = h_W,c * (theta_pert,c - theta_base,c); sum_c q = delta_beta.
# ALL 20 measures, both mechanisms, both windows; no flag-based filtering.

calculate_cohort_contributions <- function(measures_long, dictionary, cfg,
                                           deterministic_trends) {
  out <- list()
  for (w in c("full", "late")) {
    cw <- window_cohorts(cfg, w)
    h <- trend_slope_weights(cw)
    for (mech in cfg$mechanisms) for (m in dictionary$measure_id) {
      base <- measures_long[measures_long$mechanism == "Baseline" &
                            measures_long$measure_id == m &
                            measures_long$cohort %in% cw, ]
      pert <- measures_long[measures_long$mechanism == mech &
                            measures_long$measure_id == m &
                            measures_long$cohort %in% cw, ]
      base <- base[order(base$cohort), ]; pert <- pert[order(pert$cohort), ]
      d <- pert$value - base$value
      q <- h * d
      db <- deterministic_trends$delta_beta[
        deterministic_trends$window == w &
        deterministic_trends$measure_id == m &
        deterministic_trends$mechanism == mech]
      stopifnot(abs(sum(q) - db) < 1e-10)
      out[[length(out) + 1]] <- dplyr::tibble(
        window = w, mechanism = mech, measure_id = m, cohort = cw,
        ols_slope_weight = h, cohort_difference = d, q_delta = q)
    }
  }
  dplyr::bind_rows(out)
}
