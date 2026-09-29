#-------------------------------------------------------------------------------
# perturbation_fun/perturbation_measures.R
#
# Thin wrappers around the paper's production estimators in
# Functions/synthetic_dynasty_estimation.R. No mobility formula is
# rederived here: OM/SM/EM/upward/downward come from
# compute_movement_mobility(); the production steady-state stop rule
# (tol = 0.01) and stationary shares come from
# calculate_the_steady_state_tm(); AIM/IM from calculate_memory_curve();
# MTE_i = 1/(1 - P_ii) as in the paper. The exact invariant distribution
# (left Perron eigenvector) and the fundamental-matrix mean first passage
# time are the linear-algebra utilities behind SSM and MFPT, matching the
# validated Appendix J production values.
#-------------------------------------------------------------------------------

# production functions expect a data.frame with columns Class1..Class5
construct_production_measure_input <- function(P) {
  d <- as.data.frame(P)
  colnames(d) <- paste0("Class", seq_len(ncol(P)))
  d
}

# exact invariant distribution: left eigenvector of P for eigenvalue 1
stationary_distribution_exact <- function(P) {
  ev <- eigen(t(P))
  v <- Re(ev$vectors[, which.max(Re(ev$values))])
  v / sum(v)
}

# mean first passage time into class `target` (fundamental matrix)
mean_first_passage_into <- function(P, target) {
  idx <- setdiff(seq_len(nrow(P)), target)
  Q <- P[idx, idx, drop = FALSE]
  mean(solve(diag(nrow(Q)) - Q, rep(1, nrow(Q))))
}

# all 20 Appendix J measures for one (P, mu0) pair
calculate_reported_paper_measures <- function(P, mu0, steady_state_tol = 0.01) {
  n_cls <- nrow(P)
  Pdf <- construct_production_measure_input(P)
  ss  <- calculate_the_steady_state_tm(Pdf, n_cls, steady_state_tol)
  pinf <- ss$p_infty
  mv  <- compute_movement_mobility(Pdf, mu0, 1, n_cls)
  mem <- calculate_memory_curve(Pdf, 1, pinf[, 1:n_cls], n_cls, mu0)
  pi_inf <- stationary_distribution_exact(P)
  mte <- unname(1 / (1 - diag(P)))
  im  <- unname(as.numeric(mem$IM_vec))
  out <- c(
    historical  = mv$historical_mobility,
    structural  = mv$structural_mobility,
    exchange    = mv$exchange_mobility,
    upward      = mv$historical_upward_mobility,
    downward    = mv$historical_downward_mobility,
    AMTE        = mean(mte),
    stats::setNames(mte, paste0("MTE_C", 1:n_cls)),
    AIM         = as.numeric(mem$AIM),
    stats::setNames(im, paste0("IM_C", 1:n_cls)),
    SSM         = 1 - sum(pi_inf * diag(P)),
    ss_C4_share = unname(pinf[n_cls, 4]),
    mfp_to_C4   = mean_first_passage_into(P, 4)
  )
  stopifnot(all(is.finite(out)))
  out
}

# tidy cohort-level measure results for baseline + both mechanisms
calculate_cohort_measure_results <- function(cohort_objects, dictionary, kappa) {
  dplyr::bind_rows(lapply(cohort_objects, function(ch) {
    dplyr::bind_rows(lapply(names(ch$results), function(mech) {
      v <- ch$results[[mech]]$measures
      dplyr::tibble(cohort = ch$cohort, mechanism = mech,
                    kappa = ifelse(mech == "Baseline", 0, kappa),
                    measure_id = names(v), value = as.numeric(v))
    }))
  }))
}
