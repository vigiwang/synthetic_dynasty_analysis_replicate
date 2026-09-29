# ============================================================================
# R/nlsy_measure_plot.R
# ----------------------------------------------------------------------------
# Two-segment (NLSY79 + NLSY97) cross-cohort measure plots with bias-corrected
# confidence bands, in the house style of Functions/plot.R (serif, black primary
# series, grey benchmark, shaded CI ribbons). Every NLSY panel can be overlaid
# with the GSS occupation benchmark (data/main_rst_bc.rds) at the corresponding
# cohorts.
#
# Input everywhere is the shared `validation_lst[[cohort]]$<measure>_df` object
# (columns mean_bc / CI_bc_l / CI_bc_u), produced by generate_validation_data()
# for both the NLSY processed results and the GSS benchmark.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tibble); library(ggplot2)
})

# benchmark validation_lst is UNNAMED, ordered by cohort 1945..1990 (index = year - 1944)
bench_key <- function(years) setNames(as.integer(years) - 1944L, as.character(years))

# Even-spacing x positions: the two NLSY segments are far apart in calendar time
# (1961-64 vs 1980-84) with nothing between, which crowds each cluster. Map each
# cohort to an evenly-spaced ordinal position instead, with a single-slot break
# (position 5) between the segments:
#   1961..1964 -> 1..4      |  gap @ 5  |  1980..1984 -> 6..10
POS_GAP    <- 5
POS_BREAKS <- c(1:4, 6:10)
POS_LABELS <- c(as.character(1961:1964), as.character(1980:1984))
cohort_pos <- function(coh) ifelse(coh <= 1970, coh - 1960, coh - 1980 + 6)

# extract one bias-corrected scalar series (mean_bc + CI) across a set of cohorts.
#   vl   : validation_lst (named "1960".. for NLSY, or integer-indexed for benchmark)
#   keys : named vector mapping cohort-year -> element key (name or index)
#   getr : function(v) -> one row with mean_bc / CI_bc_l / CI_bc_u
extract_series <- function(vl, keys, getr) {
  imap_dfr(keys, function(key, coh) {
    v <- vl[[key]]; if (is.null(v)) return(NULL)
    r <- getr(v)
    tibble(cohort = as.numeric(coh), pos = cohort_pos(as.numeric(coh)),
           mean_bc = r$mean_bc[1], CI_bc_l = r$CI_bc_l[1], CI_bc_u = r$CI_bc_u[1])
  })
}

# assemble a multi-series data frame; specs is a list of list(series=, vl=, keys=).
# A series that spans BOTH cohort segments gets an NA row inserted at the gap slot
# so the line/ribbon breaks instead of bridging the two segments.
build_measure_df <- function(getr, specs) {
  df <- bind_rows(lapply(specs, function(s)
    extract_series(s$vl, s$keys, getr) %>% mutate(series = s$series)))
  brk <- df %>% group_by(series) %>%
    summarise(span = any(cohort < 1970) & any(cohort > 1975), .groups = "drop") %>%
    filter(span) %>%
    transmute(series, cohort = NA_real_, pos = POS_GAP,
              mean_bc = NA_real_, CI_bc_l = NA_real_, CI_bc_u = NA_real_)
  bind_rows(df, brk) %>% arrange(series, pos)
}

# house-style palettes (subset to the series actually present)
.mp_col <- c(`NLSY income` = "#000000", `NLSY income (unadjusted)` = "#000000",
             `NLSY occupation` = "#3D3D3D",
             `Benchmark: GSS occupation` = "#8a8a8a",
             `Design-based` = "#000000", `Individual bootstrap` = "#8a8a8a",
             `Income: upward` = "#000000", `Income: downward` = "#000000",
             `Benchmark: upward` = "#8a8a8a", `Benchmark: downward` = "#8a8a8a",
             `Exchange (t = 1)` = "#c2c2c2", `Exchange (t = 2)` = "#8a8a8a",
             `Exchange (t = 3)` = "#565656", `Exchange (steady state)` = "#000000")
.mp_lty <- c(`NLSY income` = "solid", `NLSY income (unadjusted)` = "solid",
             `NLSY occupation` = "dotdash",
             `Benchmark: GSS occupation` = "dashed",
             `Design-based` = "solid", `Individual bootstrap` = "dashed",
             `Income: upward` = "solid", `Income: downward` = "longdash",
             `Benchmark: upward` = "solid", `Benchmark: downward` = "longdash",
             `Exchange (t = 1)` = "dotted", `Exchange (t = 2)` = "dashed",
             `Exchange (t = 3)` = "longdash", `Exchange (steady state)` = "solid")
.mp_shp <- c(`NLSY income` = 16, `NLSY income (unadjusted)` = 16,
             `NLSY occupation` = 18,
             `Benchmark: GSS occupation` = 15,
             `Design-based` = 16, `Individual bootstrap` = 17,
             `Income: upward` = 17, `Income: downward` = 15,
             `Benchmark: upward` = 2,  `Benchmark: downward` = 0,
             `Exchange (t = 1)` = 1, `Exchange (t = 2)` = 2,
             `Exchange (t = 3)` = 5, `Exchange (steady state)` = 16)

# ---------------------------------------------------------------------------
# Paper typography (matches Functions/plot.R + the GSS2024 figures):
#   Times New Roman; BOLD axis titles; plain (non-bold) tick labels; plain
#   legend text with NO legend title; bold facet-strip labels; no panel grid;
#   black L-shaped axis lines with outside ticks. Sizes are tuned for the
#   GSS-matching 5000x3000 px render (10 x 6 in @ 500 dpi); if you render at a
#   different physical size, scale FONT_SCALE accordingly.
# ---------------------------------------------------------------------------
PAPER_FAMILY <- "Times New Roman"   # falls back to a serif face if unavailable
FONT_SCALE   <- 1.0                 # global multiplier for all text sizes
.ps <- function(x) x * FONT_SCALE   # pt-size helper
theme_paper <- function(legend_position = "top", x_angle = 0) {
  x_hjust <- if (x_angle == 0) 0.5 else 1
  theme_minimal(base_family = PAPER_FAMILY, base_size = .ps(18)) +
    theme(
      panel.grid        = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.6),
      axis.ticks        = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title.x = element_text(face = "bold", size = .ps(24), margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = .ps(24), margin = margin(r = 8)),
      axis.text.x  = element_text(size = .ps(17), colour = "black",
                                  angle = x_angle, hjust = x_hjust,
                                  vjust = if (x_angle == 0) NA_real_ else 0.5),
      axis.text.y  = element_text(size = .ps(17), colour = "black"),
      legend.position = legend_position,
      legend.title    = element_blank(),
      legend.text     = element_text(size = .ps(17)),
      legend.key.width = unit(2.2, "lines"),
      strip.text   = element_text(face = "bold", size = .ps(19)),
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption  = element_text(size = .ps(13), colour = "grey40"),
      plot.margin   = margin(10, 14, 10, 10)
    )
}

# the plot. df: cohort, pos, mean_bc, CI_bc_l, CI_bc_u, series.
#   y_limits      : optional c(lo, hi) to zoom the y-axis (coord_cartesian; does not
#                   drop data from stats) -- used to keep a degenerate off-scale point
#                   (e.g. sparse-cell Altham) from wrecking the readable range.
#   ribbon_series : optional subset of series that get a CI band (default: all).
measure_plot <- function(df, y_title, series_levels = NULL, seg_labels = TRUE,
                         y_zero = FALSE, y_limits = NULL, ribbon_series = NULL,
                         title = NULL, subtitle = NULL) {
  if (is.null(series_levels)) series_levels <- unique(df$series)
  if (!is.null(title))    title    <- gsub("\\s{2,}", " ", title)     # collapse double spaces
  if (!is.null(subtitle)) subtitle <- gsub("\\s{2,}", " ", subtitle)
  df$series <- factor(df$series, levels = series_levels)
  col <- .mp_col[series_levels]; lty <- .mp_lty[series_levels]; shp <- .mp_shp[series_levels]
  rib_df <- if (is.null(ribbon_series)) df else df[df$series %in% ribbon_series, ]
  p <- ggplot(df, aes(pos, mean_bc, color = series, fill = series,
                      linetype = series, shape = series)) +
    geom_vline(xintercept = POS_GAP, linetype = "dotted", color = "grey75", linewidth = 0.5) +
    geom_ribbon(data = rib_df, aes(ymin = CI_bc_l, ymax = CI_bc_u),
                alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 3.8, na.rm = TRUE) +
    scale_color_manual(values = col, drop = FALSE) +
    scale_fill_manual(values = col, drop = FALSE) +
    scale_linetype_manual(values = lty, drop = FALSE) +
    scale_shape_manual(values = shp, drop = FALSE) +
    scale_x_continuous(breaks = POS_BREAKS, labels = POS_LABELS,
                       limits = c(0.5, 10.5), expand = c(0.01, 0)) +
    # paper cosmetics (Functions/plot.R): x = "Birth cohorts"; unit on the
    # y-axis; bold serif axis titles, plain ticks (title/subtitle stay off-plot).
    labs(x = "Birth cohorts", y = y_title, title = title, subtitle = subtitle,
         color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
    theme_paper(legend_position = "top")
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  if (y_zero) p <- p + expand_limits(y = 0)
  if (seg_labels)
    p <- p +
      annotate("text", x = 2.5, y = Inf, label = "NLSY79", vjust = 1.4,
               family = PAPER_FAMILY, fontface = "bold", size = .ps(6), color = "grey35") +
      annotate("text", x = 8, y = Inf, label = "NLSY97", vjust = 1.4,
               family = PAPER_FAMILY, fontface = "bold", size = .ps(6), color = "grey35")
  p
}

# like build_measure_df() but each spec carries its OWN getter (series may pull
# different quantities, e.g. upward vs downward). specs: list(series, vl, keys, getr).
build_measure_df_multi <- function(specs) {
  df <- bind_rows(lapply(specs, function(s)
    extract_series(s$vl, s$keys, s$getr) %>% mutate(series = s$series)))
  brk <- df %>% group_by(series) %>%
    summarise(span = any(cohort < 1970) & any(cohort > 1975), .groups = "drop") %>%
    filter(span) %>%
    transmute(series, cohort = NA_real_, pos = POS_GAP,
              mean_bc = NA_real_, CI_bc_l = NA_real_, CI_bc_u = NA_real_)
  bind_rows(df, brk) %>% arrange(series, pos)
}

# upward-vs-downward panel: NLSY income (both segments) up & down, plus the GSS
# benchmark up & down; CI bands only on the two NLSY series to keep it readable.
updown_plot <- function(inc79, inc97, bench, nlk79, nlk97, bkA) {
  gU <- function(v) v$movement_df %>% filter(category == "Historical Upward Mobility",   t == 1)
  gD <- function(v) v$movement_df %>% filter(category == "Historical Downward Mobility", t == 1)
  specs <- list(
    list(series = "Income: upward",     vl = inc79, keys = nlk79, getr = gU),
    list(series = "Income: upward",     vl = inc97, keys = nlk97, getr = gU),
    list(series = "Income: downward",   vl = inc79, keys = nlk79, getr = gD),
    list(series = "Income: downward",   vl = inc97, keys = nlk97, getr = gD),
    list(series = "Benchmark: upward",  vl = bench, keys = bkA,  getr = gU),
    list(series = "Benchmark: downward",vl = bench, keys = bkA,  getr = gD))
  lv <- c("Income: upward", "Income: downward",
          "Benchmark: upward", "Benchmark: downward")
  measure_plot(build_measure_df_multi(specs),
               "Probability to move  (t = 1)",
               series_levels = lv,
               ribbon_series = lv)   # CI bands on all series, incl. GSS occ benchmark
}

# class-specific small multiples: one facet per origin/destination class, each a
# two-segment NLSY-vs-benchmark panel. getr_factory(k) -> getter for class k.
class_facet_plot <- function(specs, classes, getr_factory, y_title, series_levels) {
  df <- bind_rows(lapply(classes, function(k)
    build_measure_df_multi(lapply(specs, function(s) c(s, list(getr = getr_factory(k))))) %>%
      mutate(class = k)))
  df$series <- factor(df$series, levels = series_levels)
  df$class  <- factor(df$class,  levels = classes)
  col <- .mp_col[series_levels]; lty <- .mp_lty[series_levels]; shp <- .mp_shp[series_levels]
  ggplot(df, aes(pos, mean_bc, color = series, fill = series,
                 linetype = series, shape = series)) +
    geom_vline(xintercept = POS_GAP, linetype = "dotted", color = "grey80", linewidth = 0.4) +
    geom_ribbon(aes(ymin = CI_bc_l, ymax = CI_bc_u), alpha = 0.13, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_point(size = 2.6, na.rm = TRUE) +
    scale_color_manual(values = col, drop = FALSE) +
    scale_fill_manual(values = col, drop = FALSE) +
    scale_linetype_manual(values = lty, drop = FALSE) +
    scale_shape_manual(values = shp, drop = FALSE) +
    scale_x_continuous(breaks = POS_BREAKS, labels = POS_LABELS, limits = c(0.5, 10.5)) +
    facet_wrap(~class, scales = "free_y", ncol = 3) +
    labs(x = "Birth cohorts", y = y_title,
         color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
    theme_paper(legend_position = "top", x_angle = 90) +
    theme(axis.text.x = element_text(size = .ps(12)))   # denser panels: smaller ticks
}

# weighted pooled father+child income distribution with the quintile cut-points,
# one panel per segment (mirrors the NLSY Final Report's distribution figure).
# inputs: named list; each element a data.frame with parent, child, weight.
income_dist_plot <- function(inputs,
                             xlab  = "Equivalized real income  (log scale, 2023$)",
                             title = "Pooled father + child income distribution and quintile cut-points") {
  dens <- imap_dfr(inputs, function(d, nm) {
    d <- d %>% filter(is.finite(parent), is.finite(child), is.finite(weight),
                      parent > 0, child > 0)
    bind_rows(tibble(segment = nm, gen = "Parent (father)", income = d$parent, w = d$weight),
              tibble(segment = nm, gen = "Child (adult)",   income = d$child,  w = d$weight))
  })
  cuts <- imap_dfr(inputs, function(d, nm) {
    d  <- d %>% filter(is.finite(parent), is.finite(child), is.finite(weight),
                       parent > 0, child > 0)
    cp <- pooled_weighted_cutpoints(d$parent, d$child, d$weight, 5)
    tibble(segment = nm, cut = as.numeric(cp), lab = paste0("Q", c(20, 40, 60, 80)),
           vj = c(1.4, 2.9, 1.4, 2.9))          # stagger to avoid label overlap
  })
  dens$segment <- factor(dens$segment, levels = names(inputs))
  cuts$segment <- factor(cuts$segment, levels = names(inputs))
  ggplot(dens, aes(income, weight = w, color = gen, fill = gen)) +
    geom_vline(data = cuts, aes(xintercept = cut), linetype = "dashed",
               color = "grey55", linewidth = 0.4, inherit.aes = FALSE) +
    geom_density(alpha = 0.18, adjust = 1.2, linewidth = 0.8, na.rm = TRUE) +
    geom_text(data = cuts, aes(x = cut, y = Inf, label = lab, vjust = vj), inherit.aes = FALSE,
              size = 2.9, color = "grey45", family = "serif") +
    scale_x_log10(labels = function(x) paste0("$", formatC(x, format = "d", big.mark = ","))) +
    scale_color_manual(values = c(`Parent (father)` = "#8a8a8a", `Child (adult)` = "#000000")) +
    scale_fill_manual(values  = c(`Parent (father)` = "#8a8a8a", `Child (adult)` = "#000000")) +
    facet_wrap(~segment, scales = "free", ncol = 1) +
    labs(x = xlab, y = "Weighted density", color = NULL, fill = NULL, title = title) +
    theme_paper(legend_position = "top") +
    theme(plot.title = element_text(face = "bold", size = .ps(20), hjust = 0.5))
}

# class-by-survey panel (NLSY Final Report "By origin class" style): two panels
# side by side (NLSY79 | NLSY97), each showing the five class-specific series
# across that segment's cohorts. getr_factory(k) -> getter for class k.
class_by_survey_plot <- function(getr_factory, vl79, vl97, keys79, keys97,
                                 y_title, classes = paste0("Class", 1:5)) {
  one <- function(vl, keys, lab) bind_rows(lapply(classes, function(k)
    extract_series(vl, keys, getr_factory(k)) %>% mutate(class = k))) %>%
    mutate(survey = lab)
  df <- bind_rows(one(vl79, keys79, "NLSY79"), one(vl97, keys97, "NLSY97"))
  df$class  <- factor(df$class,  levels = classes)
  df$survey <- factor(df$survey, levels = c("NLSY79", "NLSY97"))
  greys  <- setNames(c("#9e9e9e", "#767676", "#525252", "#2e2e2e", "#000000"), classes)
  shapes <- setNames(c(16, 15, 18, 17, 25), classes)
  ltys   <- setNames(c("solid", "dashed", "dotdash", "dotted", "longdash"), classes)
  ggplot(df, aes(cohort, mean_bc, color = class, shape = class,
                 linetype = class, group = class)) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_point(size = 3.2, na.rm = TRUE) +
    scale_color_manual(values = greys, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    scale_linetype_manual(values = ltys, name = NULL) +
    scale_x_continuous(breaks = c(1961:1964, 1980:1984)) +
    facet_wrap(~survey, scales = "free_x") +
    labs(x = "Birth cohorts", y = y_title) +
    theme_paper(legend_position = "right")
}

# class-specific two-series comparison (e.g. NLSY79 occupation vs GSS occupation):
# one facet per origin class, both series overlaid with bias-corrected CI ribbons.
# getr_factory(k) -> getter for class k (same factories as class_by_survey_plot).
class_compare_plot <- function(getr_factory, vlA, vlB, keysA, keysB,
                               labels = c("NLSY79 occupation", "Benchmark: GSS occupation"),
                               y_title, classes = paste0("Class", 1:5)) {
  one <- function(vl, keys, lab) bind_rows(lapply(classes, function(k)
    extract_series(vl, keys, getr_factory(k)) %>% mutate(class = k))) %>%
    mutate(series = lab)
  df <- bind_rows(one(vlA, keysA, labels[1]), one(vlB, keysB, labels[2]))
  df$class  <- factor(df$class, levels = classes)
  df$series <- factor(df$series, levels = labels)
  pal <- setNames(c("#000000", "#8a8a8a"), labels)
  lty <- setNames(c("solid", "dashed"), labels)
  shp <- setNames(c(16, 15), labels)
  ggplot(df, aes(cohort, mean_bc, color = series, fill = series,
                 linetype = series, shape = series, group = series)) +
    geom_ribbon(aes(ymin = CI_bc_l, ymax = CI_bc_u), alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_point(size = 3.0, na.rm = TRUE) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = lty, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    facet_wrap(~class, nrow = 1) +
    scale_x_continuous(breaks = 1961:1964) +
    labs(x = "Birth cohorts", y = y_title) +
    theme_paper(legend_position = "top", x_angle = 45)
}

# ---- standard measure getters (work on both NLSY and benchmark validation_lst) ----
mg <- list(
  overall   = function(v) v$movement_df %>% filter(category == "Historical Mobility",          t == 1),
  upward    = function(v) v$movement_df %>% filter(category == "Historical Upward Mobility",    t == 1),
  downward  = function(v) v$movement_df %>% filter(category == "Historical Downward Mobility",  t == 1),
  structural= function(v) v$movement_df %>% filter(category == "Structural Mobility",           t == 1),
  exchange  = function(v) v$movement_df %>% filter(category == "Exchange Mobility",             t == 1),
  aim       = function(v) v$memory_df   %>% filter(class == "AIM",  t == 1),
  lambda2   = function(v) v$lambda2_df,
  amte      = function(v) v$MTE_df      %>% filter(class == "AMTE"),
  altham    = function(v) v$altham,
  geenens   = function(v) v$GeenensD,
  hellinger = function(v) v$HellingersDep
)

# steady-state getter factory: the converged value of a movement measure = the row
# at the LARGEST generation t that still has a non-NA bias-corrected estimate. At
# the steady state the marginal distribution has converged (structural mobility
# -> 0), so overall mobility equals the long-run "circulation" (exchange) mobility.
# Works uniformly on the NLSY income, NLSY occupation and GSS benchmark
# validation_lst objects (income runs to t = 7, occupation / GSS to t = 5), so it
# lets fig4 compare the STEADY-STATE mobility of all three measures on one panel.
mg_ss <- function(cat_name) function(v) {
  d <- v$movement_df %>%
    dplyr::filter(category == cat_name, !is.na(mean_bc)) %>%
    dplyr::mutate(t = as.integer(t))
  d %>% dplyr::filter(t == max(t)) %>% dplyr::slice(1)
}

# class-specific getter factories (return a getter for a given "Class1".."Class5")
mte_by_class <- function(k) function(v) v$MTE_df    %>% filter(class == k)              # marginal transition effect
imm_by_class <- function(k) function(v) v$memory_df %>% filter(class == k, t == 1)      # first-generation immobility
CLASSES <- paste0("Class", 1:5)

# ---------------------------------------------------------------------------
# updown_adj_plot(): baseline upward vs downward mobility, unadjusted (household)
# vs equivalized (headline), one panel per survey. Input is the tidy frame from
# NLSY/_precompute_unadjusted_updown.R (survey, cohort, adj, up, down);
# the pooled "Full" row is dropped (per-cohort trend only).
# ---------------------------------------------------------------------------
updown_adj_plot <- function(df) {
  d <- df %>%
    dplyr::filter(cohort != "Full") %>%
    dplyr::mutate(cohort = as.integer(cohort)) %>%
    tidyr::pivot_longer(c(up, down), names_to = "dir", values_to = "p") %>%
    dplyr::mutate(dir = factor(ifelse(dir == "up", "Upward", "Downward"),
                               c("Upward", "Downward")),
                  adj = factor(adj, c("Unadjusted (household)", "Equivalized (headline)")),
                  survey = factor(survey, c("NLSY79", "NLSY97")))
  ggplot(d, aes(cohort, p, colour = dir, linetype = adj, shape = adj, group = interaction(dir, adj))) +
    geom_line(linewidth = 1.1, na.rm = TRUE) +
    geom_point(size = 3.2, fill = "white", na.rm = TRUE) +
    scale_colour_manual(values = c(Upward = "#1a1a1a", Downward = "#c0392b"), name = NULL) +
    scale_linetype_manual(values = c(`Unadjusted (household)` = "solid",
                                     `Equivalized (headline)` = "22"), name = NULL) +
    scale_shape_manual(values = c(`Unadjusted (household)` = 16,
                                  `Equivalized (headline)` = 21), name = NULL) +
    scale_x_continuous(breaks = c(1961:1964, 1980:1984)) +
    facet_wrap(~survey, scales = "free_x") +
    labs(x = "Birth cohorts", y = "Probability to move") +
    theme_paper(legend_position = "top") +
    theme(legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
}

# ---------------------------------------------------------------------------
# rankrank_plot(): BDZ 2018 rank-rank slope per cohort, two panels (NLSY79 ·
# NLSY97), with a 95% design-based CI band and a dotted line at BDZ's own
# per-survey Table-1 value. Input: Data/NLSY_estimation/intermediate/rankrank_by_cohort.rds
# (survey, cohort, rr, rr_se); the pooled "Full" row is dropped.
# ---------------------------------------------------------------------------
rankrank_plot <- function(df) {
  d <- df %>% dplyr::filter(cohort != "Full") %>%
    dplyr::mutate(cohort = as.integer(cohort), lo = rr - 1.96 * rr_se,
                  hi = rr + 1.96 * rr_se, survey = factor(survey, c("NLSY79", "NLSY97")))
  bdz <- tibble::tibble(survey = factor(c("NLSY79", "NLSY97"), c("NLSY79", "NLSY97")),
                        b = c(0.426, 0.424))
  ggplot(d, aes(cohort, rr)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#1a1a1a") +
    geom_hline(data = bdz, aes(yintercept = b), linetype = "dotted",
               color = "#c0392b", linewidth = 0.7) +
    geom_line(linewidth = 1.1, color = "#1a1a1a", na.rm = TRUE) +
    geom_point(size = 3.4, color = "#1a1a1a", na.rm = TRUE) +
    scale_x_continuous(breaks = c(1961:1964, 1980:1984)) +
    facet_wrap(~survey, scales = "free_x") +
    labs(x = "Birth cohorts", y = "Rank-rank slope",
         caption = "Dotted red = BDZ 2018 Table 1 per-survey rank-rank slope.") +
    theme_paper(legend_position = "none")
}

# ---------------------------------------------------------------------------
# exchange-mobility-across-synthetic-generations panel
# (report chunks m-ss-exchange / m-ss-exchange-u): exchange mobility at each
# synthetic generation t = 1, 2, 3 overlaid with the steady-state value, both
# surveys -- the same layout as the Overall-mobility-across-generations figure
# (OM_t, paper Fig. 3), just for the exchange component. t=1,2,3 are read from
# movement_df; the steady-state line is the ss_exchange precompute (v$ss_df).
# vl79/vl97: processed validation lists; ss79/ss97: steady-state lists
# (ss_exchange rds); keys: cohort keys.
ss_exchange_plot <- function(vl79, vl97, ss79, ss97, keys79, keys97) {
  ex_t <- function(k) function(v)             # exchange mobility at generation t = k
    v$movement_df %>% filter(category == "Exchange Mobility", t == k)
  gen_lvls <- c("Exchange (t = 1)", "Exchange (t = 2)",
                "Exchange (t = 3)", "Exchange (steady state)")
  df_ss <- dplyr::bind_rows(
    build_measure_df(ex_t(1),
      list(list(series = "Exchange (t = 1)", vl = vl79, keys = keys79),
           list(series = "Exchange (t = 1)", vl = vl97, keys = keys97))),
    build_measure_df(ex_t(2),
      list(list(series = "Exchange (t = 2)", vl = vl79, keys = keys79),
           list(series = "Exchange (t = 2)", vl = vl97, keys = keys97))),
    build_measure_df(ex_t(3),
      list(list(series = "Exchange (t = 3)", vl = vl79, keys = keys79),
           list(series = "Exchange (t = 3)", vl = vl97, keys = keys97))),
    build_measure_df(function(v) v$ss_df,
      list(list(series = "Exchange (steady state)", vl = ss79, keys = keys79),
           list(series = "Exchange (steady state)", vl = ss97, keys = keys97)))) %>%
    dplyr::mutate(series = factor(series, gen_lvls))
  # route through the shared house-style plotter so cosmetics (axis-title weight,
  # tick-label size, segment labels) are identical to every other figure; the
  # four exchange series draw from the .mp_col/.mp_lty/.mp_shp palettes.
  measure_plot(df_ss, "Probability to move", series_levels = gen_lvls)
}

# synthetic-dynasty paper cosmetics (Functions/plot.R): no background grid,
# black axis lines, outside tick marks (5pt), serif/Times family
style_dynasty <- function(p) p +
  theme(panel.grid        = element_blank(),
        axis.line         = element_line(color = "black", linewidth = 0.6),
        axis.ticks        = element_line(color = "black", linewidth = 0.6),
        axis.ticks.length = unit(0.15, "cm"))
