#-------------------------------------------------------------------------------
# perturbation_fun/perturbation_figures.R
#
# Figure constructors for the perturbation appendix (manuscript lettering
# K). Each function consumes already computed tidy data from
# perturbation_results.rds (whose keys keep the original figure_J*_data
# names) and returns ONE final ggplot object. Each figure family (K.1
# trend changes, K.2 slope calibration, K.3 cohort contributions) is
# rendered as three sub-figures via the measure groups below, so every
# sub-figure is large enough to read at text width.
# Style follows Functions/plot.R (the paper's visual authority, sourced
# by script 02): Times New Roman, blank grid and background, black
# x-axis line, outward black ticks, leading zeros stripped from numeric
# tick labels, grayscale only. Mechanism marks are fixed across all
# figures:
#   Redistribution — black (#1F1F1F) filled circle, solid interval
#   Augmentation   — gray (#7A7A7A) filled square, dashed interval
#-------------------------------------------------------------------------------

# sub-figure measure groups (user-chosen split, four blocks):
#   movement = all movement measures incl. SSM (pp per decade)
#   mte      = class-specific mean time to exit (generations per decade)
#   memory   = AIM + class-specific IM (index units per decade)
#   class4   = the Class-4 targets (mixed units)
# AMTE appears in the tables and the canonical rds but is excluded from
# the figures (user decision). Vector order = row order in the figures.
perturbation_figure_groups <- function() {
  list(
    movement = c("historical", "structural", "exchange", "upward", "downward",
                 "SSM"),
    mte      = paste0("MTE_C", 1:5),
    memory   = c("AIM", paste0("IM_C", 1:5)),
    class4   = c("ss_C4_share", "mfp_to_C4")
  )
}

perturbation_mech_shapes <- c(Redistribution = 16, Augmentation = 15)
perturbation_mech_lines  <- c(Redistribution = "solid", Augmentation = "dashed")
perturbation_mech_cols   <- c(Redistribution = "#1F1F1F", Augmentation = "#7A7A7A")

# numeric tick labels, paper rounding rule: every break on an axis uses
# the SAME number of decimals — two, or more only where the panel scale
# needs them to distinguish breaks (e.g. the .004-scale panels). The
# integer zero is omitted (plot.R to_text convention, ".65") unless the
# axis reaches |value| >= 1, in which case all labels keep it ("0.25",
# "1.25") for within-axis consistency.
perturbation_tick_text <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_character_, length(x))
  if (!any(ok)) return(out)
  decimals_needed <- function(v) {
    for (d in 0:6) if (isTRUE(all.equal(v, round(v, d), tolerance = 1e-8)))
      return(d)
    6L
  }
  d <- max(2L, vapply(x[ok], decimals_needed, integer(1)))
  out[ok] <- formatC(x[ok], format = "f", digits = d)
  out[ok] <- sub("^-(0?\\.?0*)$", "\\1", out[ok])   # no negative zero
  if (max(abs(x[ok])) < 1) out[ok] <- sub("^(-?)0\\.", "\\1.", out[ok])
  out
}

# serif family with a Linux/HPC fallback. Midway nodes ship no serif
# font, so PDFs are never rendered there (script 02 saves the ggplot rds
# only under HPC); the rds is downloaded and rendered on the local Mac,
# where "Times New Roman" (Darwin branch) or the "serif" family baked in
# by an HPC run both resolve to Times at render time.
perturbation_serif_family <- function() {
  if (Sys.info()[["sysname"]] == "Darwin") "Times New Roman" else "serif"
}

# figure-facing measure labels: true Unicode subscripts on the acronym
# indices (OM1 -> OM<sub>1</sub> etc.), matching the tables' TeX
# subscripts (OM$_1$). Class numbers ("MTE, Class 1") are not indices
# and stay as-is. AIM uses the bare acronym in figures (user request);
# tables keep the full "Aggregate intergenerational memory (AIM$_1$)".
perturbation_measure_label <- function(x) {
  x <- sub("Aggregate intergenerational memory (AIM1)", "AIM1", x, fixed = TRUE)
  x <- sub("Mean first passage time into Class 4", "MFPT, Class 4", x,
           fixed = TRUE)
  gsub("(OM|SM|EM|AIM)1", "\\1₁", x)
}

# theme wrapper mirroring the Functions/plot.R conventions used by the
# paper cross-cohort figures: Times family; bold size-26 axis titles
# (titlefont size 26, bold); PLAIN tick labels (tickfont: x size 15,
# y size 17); legend text size 20; panel annotations size 20; outward
# black ticks; no grid; black axis line
perturbation_theme_paper <- function() {
  ggplot2::theme_minimal(base_family = perturbation_serif_family()) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(hjust = 0.5, size = 15, color = "black"),
      # y tick labels (the measure names) enlarged to 20 for readability
      axis.text.y = ggplot2::element_text(size = 20, color = "black"),
      axis.ticks.x = ggplot2::element_line(color = "black"),
      axis.ticks.y = ggplot2::element_line(color = "black"),
      axis.ticks.length = ggplot2::unit(0.2, "cm"),
      # both axis lines drawn, matching the plotly tradition in plot.R
      # (showline = TRUE on the x AND y axes)
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.title.x = ggplot2::element_text(size = 26, face = "bold",
                                           margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(size = 26, face = "bold",
                                           margin = ggplot2::margin(r = 10)),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(3, "lines"),
      legend.text = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_blank(),
      # 20 bold = plot.R plotly subplot annotation font. The strip title
      # is stacked (group / unit / window lines) so each line fits a
      # single panel width on the split sub-figure canvases; strip.clip
      # lets a borderline line overhang instead of being cut
      strip.text = ggplot2::element_text(face = "bold", size = 20),
      strip.clip = "off",
      # edge tick labels are centered on the outermost panel boundary and
      # would otherwise be half-clipped at the device edge (e.g. "1.25")
      plot.margin = ggplot2::margin(t = 12, r = 34, b = 12, l = 14))
}

# shared panel/measure factors from the dictionary (single label source);
# measure_ids restricts a figure to one sub-figure group AND fixes the
# row order (vector order governs, e.g. SSM before AIM in the summary)
.figure_panel_data <- function(d, dictionary, measure_ids = NULL) {
  dict <- dictionary[order(dictionary$display_order), ]
  if (!is.null(measure_ids)) {
    dict <- dict[match(measure_ids, dict$measure_id), ]
    d <- d[d$measure_id %in% measure_ids, ]
  }
  # explicit factor levels: without this, ggplot orders the character
  # column alphabetically (Augmentation first), flipping the legend and
  # the within-row dodge arrangement relative to the printed paper
  d$mechanism <- factor(d$mechanism,
                        levels = c("Redistribution", "Augmentation"))
  # panels split by DISPLAY UNIT only (measures sharing a unit share one
  # panel and one free x axis): pp measures merge into a single block,
  # while e.g. AIM (index units) gets its own panel in the summary figure
  panel_of <- stats::setNames(dict$unit_label, dict$measure_id)
  d$panel <- factor(panel_of[d$measure_id], levels = unique(unname(panel_of)))
  labs <- perturbation_measure_label(dict$label_plain)
  d$measure_lab <- factor(
    stats::setNames(labs, dict$measure_id)[d$measure_id],
    levels = rev(labs))
  d$window_lab <- factor(ifelse(d$window == "full",
    "Full period (1945–1990)", "Late period (1975–1990)"),
    levels = c("Full period (1945–1990)", "Late period (1975–1990)"))
  # facet id keeps the panel-major grid (separate free axes per unit
  # family) but the displayed strip shows ONLY the window line — the
  # group/unit headings are dropped from the figures (user request; units
  # are documented in the LaTeX captions)
  d$facet_id <- droplevels(interaction(d$panel, d$window_lab,
                                       sep = "|||", lex.order = TRUE))
  attr(d, "facet_labels") <-
    stats::setNames(sub("^.*\\|\\|\\|", "", levels(d$facet_id)),
                    levels(d$facet_id))
  d
}

# fix every y-axis label column of a built K1/K2 figure to ONE width:
# with the common canvas width this makes the panel positions and x-axis
# lengths identical across all K1/K2 sub-figures BY CONSTRUCTION (no
# font metrics involved), so the stacked sub-figures align in the
# manuscript. Returns a gtable; render it with grid::grid.draw() on a
# cairo_pdf device (the canonical ggplot object in the rds is untouched).
perturbation_align_axis_grob <- function(p, label_col_in = 3.4) {
  g <- ggplot2::ggplotGrob(p)
  ax_cols <- unique(g$layout$l[grepl("^axis-l", g$layout$name)])
  for (cl in ax_cols) g$widths[cl] <- grid::unit(label_col_in, "in")
  g
}

# x-axis title: states the unit directly when the figure is
# unit-homogeneous (one panel family), generic otherwise
.figure_x_title <- function(d) {
  if (nlevels(d$panel) == 1) {
    bquote(Delta * beta * "  (" * .(levels(d$panel)[1]) * ")")
  } else {
    expression(Delta * beta * "  (per decade, panel-specific units)")
  }
}

#------------------------------- Figure K.1 -----------------------------------#
# Deterministic delta beta with paired-bootstrap 95% intervals; both
# windows, both mechanisms; per-decade display units. measure_ids selects
# the sub-figure group (K.1a/b/c); x_title carries the group-specific
# axis title (e.g. a single unit when the group is unit-homogeneous).

make_figure_K1_trend_changes <- function(figure_J1_data, dictionary,
                                         measure_ids = NULL,
                                         x_title = NULL) {
  d <- .figure_panel_data(figure_J1_data, dictionary, measure_ids)
  facet_labels <- attr(d, "facet_labels")
  if (is.null(x_title)) x_title <- .figure_x_title(d)
  dg <- ggplot2::position_dodge(width = 0.78)
  ggplot2::ggplot(d, ggplot2::aes(x = delta_beta_display, y = measure_lab,
                                  group = mechanism)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = delta_beta_lo_display, xmax = delta_beta_hi_display,
                   linetype = mechanism, color = mechanism),
      position = dg, linewidth = 0.8) +
    # end caps drawn separately and ALWAYS solid: a dashed cap segment
    # renders only a fragment of the dash pattern, which made the bar
    # appear off-center at its ends
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = delta_beta_lo_display, xmax = delta_beta_lo_display,
                   color = mechanism),
      width = 0.2, orientation = "y", position = dg, linewidth = 0.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = delta_beta_hi_display, xmax = delta_beta_hi_display,
                   color = mechanism),
      width = 0.2, orientation = "y", position = dg, linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = mechanism, color = mechanism),
                        size = 2.5, stroke = 0.6, position = dg) +
    ggplot2::scale_shape_manual(values = perturbation_mech_shapes) +
    ggplot2::scale_linetype_manual(values = perturbation_mech_lines) +
    ggplot2::scale_color_manual(values = perturbation_mech_cols) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 1.4)),
      shape = ggplot2::guide_legend(override.aes = list(size = 3.5))) +
    ggplot2::scale_x_continuous(labels = perturbation_tick_text) +
    ggplot2::facet_wrap(~ facet_id, ncol = 2, scales = "free",
                        labeller = ggplot2::as_labeller(facet_labels)) +
    ggplot2::labs(x = x_title, y = NULL) +
    perturbation_theme_paper()
}

#------------------------------- Figure K.2 -----------------------------------#
# Figure J.1 marks over the bootstrap benchmark bands, exactly as printed
# in the paper:
#   light gray: signed total-system range   [q.025(u_total),  q.975(u_total)]
#   dark gray:  signed targeted range       [q.025(u_target), q.975(u_target)]
# The bands calibrate effect size; they are not confidence intervals for
# the plotted deterministic delta beta.

make_figure_K2_slope_calibration <- function(figure_J2_data, dictionary,
                                             measure_ids = NULL,
                                             x_title = NULL) {
  d <- .figure_panel_data(figure_J2_data, dictionary, measure_ids)
  facet_labels <- attr(d, "facet_labels")
  if (is.null(x_title)) x_title <- .figure_x_title(d)
  dg <- ggplot2::position_dodge(width = 0.78)
  ggplot2::ggplot(d, ggplot2::aes(y = measure_lab, group = mechanism)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = u_total_lo_display, xmax = u_total_hi_display),
      color = "grey85", linewidth = 7, position = dg) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = u_target_lo_display, xmax = u_target_hi_display),
      color = "grey62", linewidth = 3, position = dg) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = delta_beta_lo_display, xmax = delta_beta_hi_display,
                   linetype = mechanism, color = mechanism),
      position = dg, linewidth = 0.8) +
    # solid end caps (see make_figure_K1_trend_changes)
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = delta_beta_lo_display, xmax = delta_beta_lo_display,
                   color = mechanism),
      width = 0.2, orientation = "y", position = dg, linewidth = 0.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = delta_beta_hi_display, xmax = delta_beta_hi_display,
                   color = mechanism),
      width = 0.2, orientation = "y", position = dg, linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(x = delta_beta_display, shape = mechanism, color = mechanism),
      size = 2.5, stroke = 0.6, position = dg) +
    ggplot2::scale_shape_manual(values = perturbation_mech_shapes) +
    ggplot2::scale_linetype_manual(values = perturbation_mech_lines) +
    ggplot2::scale_color_manual(values = perturbation_mech_cols) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 1.4)),
      shape = ggplot2::guide_legend(override.aes = list(size = 3.5))) +
    ggplot2::scale_x_continuous(labels = perturbation_tick_text) +
    ggplot2::facet_wrap(~ facet_id, ncol = 2, scales = "free",
                        labeller = ggplot2::as_labeller(facet_labels)) +
    ggplot2::labs(x = x_title, y = NULL) +
    perturbation_theme_paper()
}

#------------------------------- Figure K.3 -----------------------------------#
# Cohort-specific OLS slope contributions q_delta, four columns:
# full/Redistribution, full/Augmentation, late/Redistribution,
# late/Augmentation. Signed bars sum exactly to the deterministic delta
# beta. measure_ids selects the sub-figure group; no flag-based filtering.

make_figure_K3_cohort_contributions <- function(figure_J3_data, dictionary,
                                                measure_ids = NULL,
                                                divider_full = 1975,
                                                divider_late = 1983) {
  dict <- dictionary[order(dictionary$display_order), ]
  d <- figure_J3_data
  if (!is.null(measure_ids)) {
    dict <- dict[match(measure_ids, dict$measure_id), ]
    d <- d[d$measure_id %in% measure_ids, ]
  }
  labs <- perturbation_measure_label(dict$label_plain)
  d$measure_lab <- factor(
    stats::setNames(labs, dict$measure_id)[d$measure_id],
    levels = labs)
  col_levels <- c("Full window (1945–1990)\nRedistribution",
                  "Full window (1945–1990)\nAugmentation",
                  "Late window (1975–1990)\nRedistribution",
                  "Late window (1975–1990)\nAugmentation")
  d$col <- factor(paste0(ifelse(d$window == "full",
    "Full window (1945–1990)\n", "Late window (1975–1990)\n"),
    d$mechanism), levels = col_levels)
  vlines <- data.frame(
    col = factor(col_levels, levels = col_levels),
    x = c(divider_full, divider_full, divider_late, divider_late) - 0.5)
  ggplot2::ggplot(d, ggplot2::aes(cohort, q_delta_display)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
    ggplot2::geom_col(fill = "grey40", width = 0.85) +
    ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = x),
                        linetype = "dashed", color = "grey45", linewidth = 0.35) +
    ggplot2::scale_y_continuous(n.breaks = 3,
                                labels = perturbation_tick_text) +
    ggplot2::facet_grid(measure_lab ~ col, scales = "free",
      labeller = ggplot2::labeller(measure_lab = ggplot2::label_wrap_gen(20))) +
    ggplot2::labs(x = "Birth cohorts",
                  y = expression(q[Delta][","*c] * "  (per decade)")) +
    perturbation_theme_paper() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(hjust = 0.5, size = 15,
                                          angle = 90, vjust = 0.5,
                                          color = "black"),
      axis.text.y = ggplot2::element_text(size = 15, color = "black"),
      strip.text.y = ggplot2::element_text(face = "bold", size = 13, angle = 0),
      strip.text.x = ggplot2::element_text(face = "bold", size = 16))
}
