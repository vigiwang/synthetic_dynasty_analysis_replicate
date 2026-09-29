#-------------------------------------------------------------------------------
# perturbation_fun/perturbation_mechanisms.R
#
# The two Appendix J perturbation mechanisms, exactly as in the paper:
#   Redistribution — within-father-origin-row reallocation toward sparse cells
#                    (+kappa/n_ic per target; proportional donors; mu0 fixed)
#   Augmentation   — kappa pseudo-dyads added to every sparse cell of the
#                    fitted joint pseudo-count table C0 = w * P with
#                    w_ic = N_c * mu0_hat_ic (P, mu0, mu1 all respond)
# One active implementation of each mechanism; no trend, plotting, table,
# or file-writing code here.
#-------------------------------------------------------------------------------

#---------------------------- redistribution ----------------------------------#

# feasibility ceiling: first father-origin row to exhaust donor mass
redistribution_feasibility_ceiling <- function(P, n_ic, mask) {
  vals <- c(); rows <- c()
  for (i in seq_len(nrow(P))) {
    T_i <- which(mask[i, ]); m <- length(T_i)
    if (m == 0 || m == ncol(P)) next
    vals <- c(vals, n_ic[i] * (1 - sum(P[i, T_i])) / m)
    rows <- c(rows, i)
  }
  if (!length(vals)) return(list(ceiling = Inf, binding_row = NA_integer_))
  list(ceiling = min(vals), binding_row = rows[which.min(vals)])
}

apply_redistribution_perturbation <- function(P, mu0, n_ic, mask, kappa, tol = 1e-10) {
  feas <- redistribution_feasibility_ceiling(P, n_ic, mask)
  if (kappa > feas$ceiling + tol) {
    stop("Redistribution infeasible: kappa = ", kappa,
         " exceeds the row feasibility ceiling ", signif(feas$ceiling, 4),
         " (binding father-origin row ", feas$binding_row, ").")
  }
  P_R <- P
  for (i in seq_len(nrow(P))) {
    T_i <- which(mask[i, ]); m <- length(T_i)
    if (m == 0) next
    if (m == ncol(P)) {
      stop("Redistribution infeasible: every destination in father-origin row ",
           i, " is flagged (no donor cells).")
    }
    S_T <- sum(P[i, T_i]); D_i <- setdiff(seq_len(ncol(P)), T_i)
    P_R[i, T_i] <- P[i, T_i] + kappa / n_ic[i]
    P_R[i, D_i] <- P[i, D_i] * (1 - S_T - kappa * m / n_ic[i]) / (1 - S_T)
  }
  stopifnot(all(is.finite(P_R)), max(abs(rowSums(P_R) - 1)) < tol,
            min(P_R) > -tol)
  list(transition_matrix = P_R,
       father_origin_marginal = mu0,                       # unchanged by design
       child_destination_marginal = as.numeric(mu0 %*% P_R),
       perturbation_budget = kappa * sum(mask),
       feasibility_ceiling = feas$ceiling,
       binding_row = feas$binding_row,
       capped = FALSE)
}

#------------------------------ augmentation ----------------------------------#

apply_augmentation_perturbation <- function(P, mu0, w_ic, mask, kappa, tol = 1e-10) {
  C0 <- w_ic * P                                # fitted joint pseudo-count table
  CA <- C0 + kappa * (mask * 1)
  n_A <- rowSums(CA)                            # = w_ic + kappa * m_ic
  P_A <- CA / n_A
  N_A <- sum(n_A)                               # = sum(w) + kappa * K_c
  mu0_A <- n_A / N_A
  stopifnot(all(is.finite(P_A)), max(abs(rowSums(P_A) - 1)) < tol,
            abs(sum(mu0_A) - 1) < tol, min(P_A) >= 0)
  list(transition_matrix = P_A,
       father_origin_marginal = mu0_A,
       child_destination_marginal = as.numeric(mu0_A %*% P_A),
       fitted_joint_table = CA,
       perturbation_budget = kappa * sum(mask),
       feasibility_ceiling = Inf,               # no algebraic ceiling
       binding_row = NA_integer_,
       capped = FALSE)
}

#------------------------------ dispatchers -----------------------------------#

apply_perturbation_to_cohort <- function(mechanism, P, mu0, n_ic, N_c, mask, kappa) {
  if (mechanism == "Redistribution") {
    apply_redistribution_perturbation(P, mu0, n_ic, mask, kappa)
  } else if (mechanism == "Augmentation") {
    apply_augmentation_perturbation(P, mu0, N_c * mu0, mask, kappa)
  } else stop("Unknown mechanism: ", mechanism)
}

# Replicate-level application (paired bootstrap). The redistribution rule
# follows the validated production behavior: the same fixed mask and raw
# n_ic ruler, with the row-specific feasibility bound applied within the
# replicate law (replicates can sit closer to their row ceiling than the
# deterministic law does). Augmentation uses w_b = N_c * mu0_b so replicate
# variation in both the transition law and the father-origin marginal
# propagates.
apply_perturbation_to_bootstrap_replicate <- function(mechanism, P_b, mu0_b,
                                                      n_ic, N_c, mask, kappa) {
  if (mechanism == "Redistribution") {
    P_R <- P_b
    for (i in seq_len(nrow(P_b))) {
      T_i <- which(mask[i, ]); m <- length(T_i)
      if (m == 0) next
      D_i <- setdiff(seq_len(ncol(P_b)), T_i)
      if (!length(D_i)) next
      S_b <- sum(P_R[i, T_i])
      k_row <- min(kappa, n_ic[i] * (1 - S_b) / m)
      if (k_row <= 0) next
      add <- k_row / n_ic[i]
      P_R[i, D_i] <- P_R[i, D_i] * (1 - S_b - m * add) / (1 - S_b)
      P_R[i, T_i] <- P_R[i, T_i] + add
    }
    list(transition_matrix = P_R, father_origin_marginal = mu0_b)
  } else if (mechanism == "Augmentation") {
    a <- apply_augmentation_perturbation(P_b, mu0_b, N_c * mu0_b, mask, kappa)
    list(transition_matrix = a$transition_matrix,
         father_origin_marginal = a$father_origin_marginal)
  } else stop("Unknown mechanism: ", mechanism)
}
