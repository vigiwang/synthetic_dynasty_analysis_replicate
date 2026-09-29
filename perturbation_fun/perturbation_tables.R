#-------------------------------------------------------------------------------
# perturbation_fun/perturbation_tables.R
#
# Table preparation and LaTeX export for Appendix J (Tables J.1-J.3).
# All numbers derive from the unrounded objects saved by script 01; no
# value is typed manually and no perturbation is recomputed here.
#-------------------------------------------------------------------------------

#------------------------------- formatting -----------------------------------#

.tex_num <- function(x, d = 3, signed = TRUE) {
  ifelse(is.na(x), "--", {
    s <- formatC(x, format = "f", digits = d, flag = if (signed) "+" else "")
    sub("^-(?=0(\\.0+)?$)", if (signed) "+" else "", s, perl = TRUE)  # no -0.000
  })
}
.tex_ci <- function(lo, hi, d = 3)
  ifelse(is.na(lo), "--",
         paste0("[", .tex_num(lo, d, FALSE), ", ", .tex_num(hi, d, FALSE), "]"))
.tex_lines <- function(lines, path) writeLines(enc2utf8(lines), path, useBytes = TRUE)

#------------------------------- Table J.1 ------------------------------------#

prepare_table_J1_sparse_cell_design <- function(table_J1_data) {
  agg <- dplyr::summarise(
    dplyr::group_by(table_J1_data, mechanism, era),
    mu0_changes = ifelse(dplyr::first(mechanism) == "Augmentation", "yes", "no"),
    med_sparse_cells = stats::median(n_sparse_cells),
    med_affected_rows = stats::median(n_affected_rows),
    med_budget = stats::median(budget),
    med_stationary_tv = stats::median(stationary_tv_k1),
    total_q50 = stats::median(total_q50),
    total_q95 = stats::median(total_q95),
    med_R_total = stats::median(R_total_stationary),
    kappa_boot95_total = stats::median(kappa_boot95_total, na.rm = TRUE),
    med_tv_mu0 = stats::median(tv_mu0),
    med_tv_mu1 = stats::median(tv_mu1),
    n_capped = sum(capped),
    .groups = "drop")
  agg[order(agg$mechanism == "Augmentation", agg$era), ]   # Redistribution first
}

write_table_J1_tex <- function(table_J1, path) {
  rows <- sprintf("%s & %s & %s & %d & %d & %d & %s & %s & %s & %s & %s & %s & %s & %d \\\\",
    table_J1$mechanism, table_J1$era, table_J1$mu0_changes,
    table_J1$med_sparse_cells, table_J1$med_affected_rows, table_J1$med_budget,
    .tex_num(table_J1$med_stationary_tv, 4, FALSE),
    .tex_num(table_J1$total_q50, 4, FALSE), .tex_num(table_J1$total_q95, 4, FALSE),
    .tex_num(table_J1$med_R_total, 3, FALSE),
    .tex_num(table_J1$kappa_boot95_total, 2, FALSE),
    .tex_num(table_J1$med_tv_mu0, 4, FALSE), .tex_num(table_J1$med_tv_mu1, 4, FALSE),
    table_J1$n_capped)
  .tex_lines(c(
    "\\begin{landscape}",
    "\\begin{table}[!htbp]",
    "\\centering\\small",
    "\\caption{Sparse-cell perturbation design and bootstrap calibration at $\\kappa=1$.}",
    "\\label{tab:sparse-design-calibration}",
    "\\begin{tabular}{llcccccccccccc}",
    "\\toprule",
    paste0("Mechanism & Era & $\\mu_0$ changes & \\shortstack{med.\\\\sparse cells} & ",
           "\\shortstack{med.\\\\affected rows} & \\shortstack{med.\\\\budget} & ",
           "\\shortstack{med.\\\\stationary TV} & \\shortstack{total\\\\q50} & ",
           "\\shortstack{total\\\\q95} & \\shortstack{med.\\\\$R_{\\mathrm{total}}(1)$} & ",
           "$\\kappa_{\\mathrm{boot95}}$ & \\shortstack{med.\\\\TV$(\\mu_0)$} & ",
           "\\shortstack{med.\\\\TV$(\\mu_1)$} & \\shortstack{capped\\\\cohorts} \\\\"),
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    "\\begin{minipage}{\\linewidth}\\vspace{4pt}\\footnotesize",
    paste0("\\emph{Notes:} A cell is sparse when its raw dyad count is at most five; ",
      "the mask is fixed across mechanisms, $\\kappa$ values, and bootstrap replicates. ",
      "At $\\kappa=1$ each cohort receives $K_c$ dyad-equivalents (one per sparse cell), ",
      "so the budget column is the median $K_c$ within the era. Redistribution moves ",
      "probability toward sparse cells within each father-origin row, holding the ",
      "smoothed marginal distribution of fathers across occupational origin classes ",
      "fixed; augmentation adds $\\kappa$ pseudo-dyads per sparse cell to the fitted ",
      "pseudo-count table, so the transition matrix and both marginals respond. ",
      "Stationary-distribution total variation (TV) compares the perturbed invariant ",
      "distribution with the baseline; total q50/q95 are the perturbation-free ",
      "bootstrap TV quantiles; $R_{\\mathrm{total}}(1)$ divides the TV at $\\kappa=1$ ",
      "by the total q95; $\\kappa_{\\mathrm{boot95}}$ is the $\\kappa$ whose TV equals ",
      "the total q95. Medians are taken across cohorts within era ",
      "(early = 1945--1974, late = 1975--1990)."),
    "\\end{minipage}", "\\end{table}", "\\end{landscape}"), path)
  invisible(path)
}

#------------------------------- Table J.2 ------------------------------------#

prepare_table_J2_cohort_trend_results <- function(table_J2_data, dictionary) {
  d <- dplyr::left_join(table_J2_data,
    dictionary[, c("measure_id", "label", "group", "decade_factor", "display_order")],
    by = "measure_id")
  d <- dplyr::mutate(d,
    dplyr::across(c(beta_base, beta_pert, delta_beta, delta_beta_lo, delta_beta_hi),
                  ~ .x * decade_factor, .names = "{.col}_display"))
  d[order(match(d$window, c("full", "late")), d$display_order), ]
}

write_table_J2_tex <- function(table_J2, path) {
  wide <- function(w) {
    d <- table_J2[table_J2$window == w, ]
    red <- d[d$mechanism == "Redistribution", ]
    aug <- d[d$mechanism == "Augmentation", ]
    stopifnot(identical(red$measure_id, aug$measure_id))
    out <- c()
    for (g in unique(red$group)) {
      out <- c(out, sprintf("\\addlinespace\\multicolumn{10}{l}{\\emph{%s}} \\\\", g))
      i <- red$group == g
      out <- c(out, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        red$label[i],
        .tex_num(red$beta_base_display[i]),
        .tex_num(red$beta_pert_display[i]), .tex_num(red$delta_beta_display[i]),
        .tex_ci(red$delta_beta_lo_display[i], red$delta_beta_hi_display[i]),
        .tex_num(red$p_direction_retained[i], 2, FALSE),
        .tex_num(aug$beta_pert_display[i], 3),
        .tex_num(aug$delta_beta_display[i], 3),
        .tex_ci(aug$delta_beta_lo_display[i], aug$delta_beta_hi_display[i]),
        .tex_num(aug$p_direction_retained[i], 2, FALSE)))
    }
    out
  }
  hdr <- paste0("Measure & \\shortstack{Baseline\\\\$\\beta$} & ",
    "\\shortstack{Redist.\\\\$\\beta$} & \\shortstack{Redist.\\\\$\\Delta\\beta$} & ",
    "\\shortstack{Redist.\\\\$\\Delta\\beta$ 95\\% CI} & \\shortstack{Redist.\\\\$P(\\mathrm{dir.})$} & ",
    "\\shortstack{Augm.\\\\$\\beta$} & \\shortstack{Augm.\\\\$\\Delta\\beta$} & ",
    "\\shortstack{Augm.\\\\$\\Delta\\beta$ 95\\% CI} & \\shortstack{Augm.\\\\$P(\\mathrm{dir.})$} \\\\")
  .tex_lines(c(
    "\\begin{landscape}",
    "\\begin{table}[!htbp]",
    "\\centering\\footnotesize",
    "\\caption{Cohort-trend results under both sparse-cell perturbation mechanisms ($\\kappa=1$).}",
    "\\label{tab:sparse-trend-results}",
    "\\begin{tabular}{lccccccccc}",
    "\\toprule",
    "\\multicolumn{10}{c}{\\emph{Panel A: full period, 1945--1990}} \\\\",
    "\\midrule", hdr, "\\midrule",
    wide("full"),
    "\\midrule",
    "\\multicolumn{10}{c}{\\emph{Panel B: late period, 1975--1990}} \\\\",
    "\\midrule", hdr, "\\midrule",
    wide("late"),
    "\\bottomrule", "\\end{tabular}",
    "\\begin{minipage}{\\linewidth}\\vspace{4pt}\\footnotesize",
    paste0("\\emph{Notes:} Slopes are per decade in panel-specific units: percentage ",
      "points for proportion measures, generations for mean-time-to-exit and ",
      "first-passage measures, and index units for intergenerational-memory measures. ",
      "Redist.\\ = within-row redistribution; Augm.\\ = augmented pseudo-dyad table; ",
      "both mechanisms use the same fixed sparse-cell mask and the same cohort-specific ",
      "input budget at $\\kappa=1$. $\\Delta\\beta=\\beta_{\\mathrm{pert}}-\\beta_{\\mathrm{base}}$. ",
      "The 95\\% interval is an unadjusted paired percentile interval for $\\Delta\\beta$ ",
      "(baseline and perturbed slopes estimated within the same stored replicate and ",
      "differenced); it is \\emph{not} an interval for the final perturbed slope. ",
      "$P(\\mathrm{dir.})$ = $P(\\text{direction retained})$, the share of replicates ",
      "whose perturbed slope keeps the sign of the deterministic baseline slope --- ",
      "a sign-stability frequency, not a $p$-value."),
    "\\end{minipage}", "\\end{table}", "\\end{landscape}"), path)
  invisible(path)
}

#------------------------------- Table J.3 ------------------------------------#

prepare_table_J3_slope_bootstrap_calibration <- function(table_J3_data, dictionary,
                                                         mechanism_order) {
  d <- dplyr::left_join(table_J3_data,
    dictionary[, c("measure_id", "label", "decade_factor", "display_order")],
    by = "measure_id")
  d <- dplyr::mutate(d,
    dplyr::across(c(beta_base, beta_pert, delta_beta, delta_beta_lo, delta_beta_hi,
                    total_beta_lo, total_beta_hi, q95_total, q95_target),
                  ~ .x * decade_factor, .names = "{.col}_display"))
  d[order(match(d$window, c("full", "late")), d$display_order,
          match(d$mechanism, mechanism_order)), ]
}

write_table_J3_tex <- function(table_J3, path) {
  row_of <- function(d) sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    d$label, ifelse(d$mechanism == "Redistribution", "Redist.", "Augm."),
    .tex_num(d$beta_base_display),
    .tex_ci(d$total_beta_lo_display, d$total_beta_hi_display),
    .tex_num(d$beta_pert_display), .tex_num(d$delta_beta_display),
    .tex_ci(d$delta_beta_lo_display, d$delta_beta_hi_display),
    .tex_num(d$q95_total_display, 3, FALSE), .tex_num(d$R_total, 2, FALSE),
    .tex_num(d$q95_target_display, 3, FALSE), .tex_num(d$R_target, 2, FALSE),
    ifelse(d$outside_total_95, "Yes", "--"))
  hdr <- paste0("Measure & Mechanism & \\shortstack{Baseline\\\\$\\beta$} & ",
    "\\shortstack{Total-bootstrap\\\\$\\beta$ 95\\% range} & \\shortstack{Perturbed\\\\$\\beta$} & ",
    "$\\Delta\\beta$ & \\shortstack{Paired $\\Delta\\beta$\\\\95\\% CI} & ",
    "\\shortstack{Total\\\\$|u|$ q95} & $R_{\\mathrm{total}}$ & ",
    "\\shortstack{Targeted\\\\$|u|$ q95} & $R_{\\mathrm{target}}$ & ",
    "\\shortstack{Outside\\\\total 95\\%} \\\\")
  .tex_lines(c(
    "\\begin{landscape}",
    "\\footnotesize",
    "\\begin{longtable}{llcccccccccc}",
    "\\caption{Perturbation-induced trend changes relative to bootstrap slope variation ($\\kappa=1$).}",
    "\\label{tab:sparse-slope-calibration} \\\\",
    "\\toprule",
    "\\multicolumn{12}{c}{\\emph{Panel A: full period, 1945--1990}} \\\\",
    "\\midrule", hdr, "\\midrule", "\\endfirsthead",
    "\\toprule", hdr, "\\midrule", "\\endhead",
    "\\midrule \\multicolumn{12}{r}{\\emph{continued}} \\\\ \\endfoot",
    "\\bottomrule \\endlastfoot",
    row_of(table_J3[table_J3$window == "full", ]),
    "\\midrule",
    "\\multicolumn{12}{c}{\\emph{Panel B: late period, 1975--1990}} \\\\",
    "\\midrule",
    row_of(table_J3[table_J3$window == "late", ]),
    "\\end{longtable}",
    "\\begin{minipage}{\\linewidth}\\footnotesize",
    paste0("\\emph{Notes:} Per-decade units as in the trend table. The total-bootstrap ",
      "$\\beta$ range is the unadjusted central 95\\% percentile range of the baseline ",
      "bootstrap slopes (not the bias-adjusted interval of the main analysis). ",
      "$u_{\\mathrm{total}}$ is the displacement of each replicate's baseline slope ",
      "from the deterministic baseline slope; $u_{\\mathrm{target}}$ is the ",
      "mechanism-specific targeted displacement (redistribution: targeted ",
      "transition-law variation with the father-origin marginal held fixed; ",
      "augmentation: targeted joint-table variation with replicate-specific ",
      "father-origin marginals). $R_{\\mathrm{total}}=|\\Delta\\beta|/q_{.95}(|u_{\\mathrm{total}}|)$ ",
      "and $R_{\\mathrm{target}}=|\\Delta\\beta|/q_{.95}(|u_{\\mathrm{target}}|)$ are ",
      "effect-size calibrations, not hypothesis tests. ``Outside total 95\\%'' marks ",
      "deterministic perturbed slopes outside the total-bootstrap range. The benchmark ",
      "quantities are not confidence intervals for $\\Delta\\beta$; the paired ",
      "$\\Delta\\beta$ interval is the uncertainty statement about the ",
      "perturbation-induced change."),
    "\\end{minipage}",
    "\\end{landscape}"), path)
  invisible(path)
}
