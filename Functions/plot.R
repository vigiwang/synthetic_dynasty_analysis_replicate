## Plot for Main results:
## @copyright: Weiqi Wang



#=============================================================================#
#                         Helper Functions
#=============================================================================#

#-----------------------------------------------------------------------------#
#           Function 1: Fetch the Cross Cohort Data from Bootstrap List
#-----------------------------------------------------------------------------#


cross_cohort_plot_data<- function(cohort_lsts, measure, attribute_col, 
                                  metric_name,t_num){
  # Step 1: Generate the Cross Cohort Dataframe:
  attribute_dfs <-
      mclapply(
        cohort_lsts,
        function(cohort_lst){
          attribute_sub_df <- 
            cohort_lst[[measure]] %>%
            rename(
              measure = !!sym(attribute_col)
            ) %>%
            group_by(measure) %>%
            filter(
              if_all(any_of("t"), ~ .x %in% t_num),
              measure %in% metric_name
            )
        }
      )
  attribute_df <- 
    bind_rows(attribute_dfs, .id = "cohort_index") %>%
    mutate(
      cohort_index = as.numeric(cohort_index)
    ) 
  return(attribute_df)
}

#-----------------------------------------------------------------------------#
#           Function 2: Tick Values To text
#-----------------------------------------------------------------------------#

to_text <- function(tick_vals){
  return(
    sub("^0(?=\\.)", "", sprintf("%.2f", tick_vals), perl = TRUE)
  )
}

#=============================================================================#
#                         Figure 1: Marginal Distribution
#=============================================================================#

plot_marginal <- function(data, group_col, label_vec, group_num){
  # Process the data:
  group_col_str <- rlang::as_string(rlang::ensym(group_col))
  obs <- 
    data %>%
    dplyr::group_by(cohort, {{ group_col }}) %>%
    dplyr::summarise(count = n(), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(count_all = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pct = count / count_all) %>%
    dplyr::group_by(cohort) %>%
    dplyr::arrange(cohort, pct) %>%
    dplyr::mutate(
      status_stack = factor(.data[[group_col_str]],
                            levels = .data[[group_col_str]][order(-pct)])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cohort = as.character(cohort))
  
  obs$status_legend <- factor(
    obs[[group_col_str]],
    levels = label_vec,
    labels = label_vec
  )
  
  # ---- Your exact style mapping (≤ 7 classes) ----
  idx <- seq_len(length(label_vec))
  uniq_measures <- label_vec[idx]
  
  # colors (reversed grayscale, dark -> light)
  gray_vec <- c("#B7B7B7","#999999","#7A7A7A","#5C5C5C","#3D3D3D","#1F1F1F","#2E2E2E")
  class_colors <- setNames(gray_vec[idx], uniq_measures)
  
  # linetypes: plotly -> ggplot equivalents
  lty_plotly <- c("dash","dot","dashdot","solid","dash","longdash","longdashdot")
  lty_gg     <- c("dashed","dotted","dotdash","solid","dashed","longdash","twodash")
  class_linetypes <- setNames(lty_gg[idx], uniq_measures)
  
  # shapes to mirror your symbols
  shp_vec <- c(16, 15, 18, 17, 25, 16, 15)
  class_shapes <- setNames(shp_vec[idx], uniq_measures)
  
  max_obs  <- (round(max(obs$count)/50) + 1) * 50
  #tick_vals <- seq(1945, 1985, by = 2)
  #year_vec  <- as.character(tick_vals)
  #tick_vals <- seq(1,34,2)
  #year_vec <- cohort_bins_vec[tick_vals]
  tick_vals <- seq(1945, 1990, by = 3)
  year_vec  <- as.character(tick_vals)
  # ---- Top panel: counts (lines + markers) ----
  p_count <- ggplot(
    obs,
    aes(x = as.numeric(cohort),
        y = count,
        color = status_legend,
        linetype = status_legend,
        shape = status_legend,            
        fill = status_legend,             
        group = status_legend)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4, stroke = 0.6) +
    scale_color_manual(values = class_colors, name = "Class") +
    scale_fill_manual(values = class_colors,  name = "Class") +
    scale_linetype_manual(values = class_linetypes, name = "Class") +
    scale_shape_manual(values = class_shapes,   name = "Class") +
    guides(
      color    = guide_legend(override.aes = list(linewidth = 1.4)),
      shape    = guide_legend(override.aes = list(size = 5))
    ) +
    scale_y_continuous(
      #breaks = seq(0, max_obs, by = 50),
      breaks = seq(0, max_obs, by = 100),
      limits = c(5, max_obs)
    ) +
    labs(y = "Number of observations", x = NULL) +
    scale_x_continuous(
      breaks = tick_vals,
      labels = year_vec,
      limits = c(min(tick_vals) - 0.5, max(tick_vals) + 0.5),
      expand = c(0, 0)
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(face = "bold", size = 20),
      axis.ticks.y = element_line(color = "black"),
      panel.grid   = element_blank(),
      axis.title.y = element_text(size = 35, margin = margin(r = 10)),
      legend.position = c(0.9, 0.8),
      legend.key.size = unit(3, "lines"),
      legend.text  = element_text(size = 20),
      legend.title = element_blank(),
      axis.ticks.length = unit(0.2, "cm"),
      panel.background = element_blank()
    )
  
  # ---- Bottom panel: composition (stacked bars) ----
  p_comp <- ggplot(
    obs, aes(x = as.numeric(cohort), y = pct, fill = status_legend)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = class_colors, name = "Class") +
    scale_x_continuous(
      breaks = tick_vals,
      labels = year_vec,
      limits = c(min(tick_vals) - 0.5, max(tick_vals) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(y = "Composition", x = "Birth cohorts") +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      axis.text.x  = element_text(hjust = 0.5, face = "bold", size = 20), # previously 23
      axis.text.y  = element_text(face = "bold", size = 25),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.title.y = element_text(size = 35, margin = margin(r = 10)),
      axis.line.x  = element_line(color = "black", linewidth = 0.6),
      axis.title.x = element_text(size = 26, margin = margin(t = 10)),
      axis.ticks.length = unit(0.2, "cm")
    )
  p_comb <- p_count / p_comp + patchwork::plot_layout(heights = c(0.5, 0.1))
  return(
    list(
      p = p_comb,
      data = obs
    )
  )
}

plot_marginal_weight <- function(data, group_col, label_vec, group_num){
  # Process the data:
  group_col_str <- rlang::as_string(rlang::ensym(group_col))
  obs <- 
    data %>%
    dplyr::group_by(cohort, {{ group_col }}) %>%
    dplyr::summarise(count = round(sum(wn)), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(count_all = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pct = count / count_all) %>%
    dplyr::group_by(cohort) %>%
    dplyr::arrange(cohort, pct) %>%
    dplyr::mutate(
      status_stack = factor(.data[[group_col_str]],
                            levels = .data[[group_col_str]][order(-pct)])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cohort = as.character(cohort))
  
  obs$status_legend <- factor(
    obs[[group_col_str]],
    levels = label_vec,
    labels = label_vec
  )
  
  # ---- Your exact style mapping (≤ 7 classes) ----
  idx <- seq_len(length(label_vec))
  uniq_measures <- label_vec[idx]
  
  # colors (reversed grayscale, dark -> light)
  gray_vec <- c("#B7B7B7","#999999","#7A7A7A","#5C5C5C","#3D3D3D","#1F1F1F","#2E2E2E")
  class_colors <- setNames(gray_vec[idx], uniq_measures)
  
  # linetypes: plotly -> ggplot equivalents
  lty_plotly <- c("dash","dot","dashdot","solid","dash","longdash","longdashdot")
  lty_gg     <- c("dashed","dotted","dotdash","solid","dashed","longdash","twodash")
  class_linetypes <- setNames(lty_gg[idx], uniq_measures)
  
  # shapes to mirror your symbols
  shp_vec <- c(16, 15, 18, 17, 25, 16, 15)
  class_shapes <- setNames(shp_vec[idx], uniq_measures)
  
  max_obs  <- (round(max(obs$count)/50) + 1) * 50
  #tick_vals <- seq(1945, 1985, by = 2)
  #year_vec  <- as.character(tick_vals)
  #tick_vals <- seq(1,34,2)
  #year_vec <- cohort_bins_vec[tick_vals]
  tick_vals <- seq(1945, 1990, by = 3)
  year_vec  <- as.character(tick_vals)
  # ---- Top panel: counts (lines + markers) ----
  p_count <- ggplot(
    obs,
    aes(x = as.numeric(cohort),
        y = count,
        color = status_legend,
        linetype = status_legend,
        shape = status_legend,            
        fill = status_legend,             
        group = status_legend)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4, stroke = 0.6) +
    scale_color_manual(values = class_colors, name = "Class") +
    scale_fill_manual(values = class_colors,  name = "Class") +
    scale_linetype_manual(values = class_linetypes, name = "Class") +
    scale_shape_manual(values = class_shapes,   name = "Class") +
    guides(
      color    = guide_legend(override.aes = list(linewidth = 1.4)),
      shape    = guide_legend(override.aes = list(size = 5))
    ) +
    scale_y_continuous(
      breaks = seq(0, max_obs, by = 50),
      #breaks = seq(0, max_obs, by = 100),
      limits = c(5, max_obs)
    ) +
    labs(y = "Number of observations", x = NULL) +
    scale_x_continuous(
      breaks = tick_vals,
      labels = year_vec,
      limits = c(min(tick_vals) - 0.5, max(tick_vals) + 0.5),
      expand = c(0, 0)
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(face = "bold", size = 20),
      axis.ticks.y = element_line(color = "black"),
      panel.grid   = element_blank(),
      axis.title.y = element_text(size = 35, margin = margin(r = 10)),
      legend.position = c(0.9, 0.8),
      legend.key.size = unit(3, "lines"),
      legend.text  = element_text(size = 20),
      legend.title = element_blank(),
      axis.ticks.length = unit(0.2, "cm"),
      panel.background = element_blank()
    )
  
  # ---- Bottom panel: composition (stacked bars) ----
  p_comp <- ggplot(
    obs, aes(x = as.numeric(cohort), y = pct, fill = status_legend)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = class_colors, name = "Class") +
    scale_x_continuous(
      breaks = tick_vals,
      labels = year_vec,
      limits = c(min(tick_vals) - 0.5, max(tick_vals) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(y = "Composition", x = "Birth cohorts") +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      axis.text.x  = element_text(hjust = 0.5, face = "bold", size = 20), # previously 23
      axis.text.y  = element_text(face = "bold", size = 25),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.title.y = element_text(size = 35, margin = margin(r = 10)),
      axis.line.x  = element_line(color = "black", linewidth = 0.6),
      axis.title.x = element_text(size = 26, margin = margin(t = 10)),
      axis.ticks.length = unit(0.2, "cm")
    )
  p_comb <- p_count / p_comp + patchwork::plot_layout(heights = c(0.5, 0.1))
  return(
    list(
      p = p_comb,
      data = obs
    )
  )
}

#=============================================================================#
# Figure 2: Overall, Structural, Exchange, Upward and Downward Mobility
#=============================================================================#


#-----------------------------------------------------------------------------#
#           Function 1: Overall Mobility
#-----------------------------------------------------------------------------#


plot_overall_mobility <- function(
    data, initial, measure, attribute_col, 
    metric_name,label_name, t_num,x_tick_vals, y_tick_vals, x_range, y_range, y_tick_text,
    legend_x, legend_y, legend_size){
  # Fetch the relevant data:
  plot_data <-
    cross_cohort_plot_data(
      data,
      measure,
      attribute_col,
      metric_name,
      t_num
    ) %>%
    group_by(
      cohort_index
    ) %>%
    mutate(
      t_ss = as.integer(t == max(t[!is.na(mean_bc)], na.rm = TRUE))
    ) %>%
    left_join(
      cohort_bins,
      by = "cohort_index"
    ) %>%
    ungroup()
  if(initial == TRUE){
    fig <-
      plot_ly() %>%
      add_ribbons(
        data = plot_data[plot_data$t == 1, ],
        x = ~cohort_bins_vec,
        ymin = ~CI_bc_l,  
        ymax = ~CI_bc_u,   
        fillcolor = "rgba(0,0,0,0.1)", 
        line = list(color = "transparent"),
        showlegend = FALSE
      ) %>%
      add_trace(
        data = plot_data[plot_data$t == 1, ],
        x = ~cohort_bins_vec,
        y = ~mean_bc,
        type = "scatter",
        mode = "lines+markers",
        line = list(color = "black", dash = "solid"), 
        marker = list(symbol = "circle", size = 10, color = "black"),
        name = paste(label_name, "(t = 1)")
      )
  }else{
    # ---- steady-state: one black line + markers by period shape ----
    ss <- plot_data %>% filter(t_ss == 1) %>% arrange(cohort_index)
    if (nrow(ss) > 0) {
      # Section 1: Plot 3 lines:
      extra_ts <- intersect(c(1, 2, 3), sort(unique(plot_data$t)))
      
      make_greys <- function(n) {
        # from #CCCCCC to #333333
        hex <- function(x) paste0("#", paste(rep(x, 3), collapse = ""))
        shades <- grDevices::colorRampPalette(c("#CCCCCC", "#333333"))(n)
        shades
      }
      greys <- make_greys(3)
      fig <- plot_ly()
      for (tt in extra_ts) {
        dat_full <- plot_data %>% dplyr::filter(t == tt) %>% dplyr::arrange(cohort_index)
        if (nrow(dat_full) == 0) next
        
        col_tt <- greys[as.numeric(tt)]
        
        fig <- fig %>%
          add_trace(
            data = dat_full,
            x = ~cohort_bins_vec, y = ~mean_bc,
            type = "scatter", mode = "lines+markers",
            line = list(
              color = col_tt,
              dash  = "dash",
              width = 2
            ),
            marker = list(
              size   = 8,
              color  = col_tt,
              line   = list(width = 0)
            ),
            name = paste0("Overall mobility(t = ", tt, ")"),
            legendgroup = paste0("t_full_", tt),
            showlegend = TRUE,
            hoverinfo = "x+y+name"
          )
      }
      
      # Section 2: Plot the steady state:
      # map period t -> marker shape (extend as you like)
      shape_map <- c(
        "3" = "triangle-down", "4" = "circle", "5" = "triangle-up",
        "6" = "diamond", "7" = "square", "8" = "x"
      )
      
      # --- CI ribbon under SS
      fig <- 
        fig %>%
        add_ribbons(
          data = ss,
          x = ~cohort_bins_vec, ymin = ~CI_bc_l, ymax = ~CI_bc_u,
          fillcolor = "rgba(0,0,0,0.10)",
          line = list(color = "rgba(0,0,0,0)"),
          showlegend = FALSE
        )
      
      # --- single black line connecting all SS points (no markers here)
      fig <- fig %>%
        add_trace(
          data = ss,
          x = ~cohort_bins_vec, y = ~mean_bc,
          type = "scatter", mode = "lines+markers",
          line = list(color = "#000000"),
          marker = list(color = "#000000", size = 8),
          name = "Steady state mobility",
          legendgroup = "steady",
          showlegend = TRUE
        )
      
  }
}
  fig <-
    fig %>%
    layout(
      title = list(text = "",
                   font = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
                   y = 0.97, x = 0.5),
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside",
        ticklen = 5,
        tickwidth = 1,
        showline = TRUE,
        titlefont = list(
          size = 26,
          family = "Times New Roman",
          color = "black",
          bold = TRUE),
        tickfont = list(size = 15),
        mirror = FALSE,
        showgrid = FALSE,
        zeroline = TRUE,
        constrain = "range",
        tickangle = 0, tickmode = "array",
        tickvals = x_tick_vals,
        range = x_range,
        cliponaxis = TRUE,
        ticktext = year_vec , tickangle = 0),
      yaxis = list(
        title = "Probability to move",
        mirror = FALSE,
        tickvals = y_tick_vals,
        range = y_range,
        cliponaxis = FALSE,
        ticktext = y_tick_text,
        tickes = "outside", ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
        tickfont = list(size = 17), showgrid = FALSE, zeroline = TRUE),
      legend = list(
        orientation = "v",
        x = legend_x,
        y = legend_y,
        constrain = "range",
        font = list(size = legend_size, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      ),
      margin = list(b = 30)
    )
  return(fig)
}

#=============================================================================#
#                      Figure 2.3: Structural Mobility
#=============================================================================#


plot_single_measure <- function(data, t_num, measure, attribute_col, 
                                metric_name, y_title, steady_state = FALSE,
                                x_tick_vals, y_tick_vals,
                                x_range, y_range, y_tick_text) {
  # Fetch the relevant data:
  if(steady_state == FALSE){
    plot_data <-
      cross_cohort_plot_data(
        data,
        measure,
        attribute_col,
        metric_name,
        t_num
      )
  }else{
    plot_data <-
      cross_cohort_plot_data(
        data,
        measure,
        attribute_col,
        metric_name,
        t_num
      ) %>%
      group_by(
        cohort_index
      ) %>%
      mutate(
        t_ss = as.integer(t == max(t[!is.na(mean_bc)], na.rm = TRUE))
      ) %>%
      filter(
        t_ss == 1
      ) %>%
      left_join(
        cohort_bins,
        by = "cohort_index"
      ) %>%
      ungroup()
  }
  
  # Generate the plot:
  fig <-
    plot_ly() %>%
    add_ribbons(
      data = plot_data,
      x = cohort_bins_vec,
      ymin = ~CI_bc_l,
      ymax = ~CI_bc_u,
      fillcolor = "rgba(0,0,0,0.1)",
      line = list(color = "transparent"), 
      name = "Confidence Interval",
      showlegend = FALSE
    ) %>%
    add_trace(
      data = plot_data,
      x = cohort_bins_vec,
      y = ~mean_bc,
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "black", dash = "solid", width = 2), 
      marker = list(
        symbol = "circle", size = 10, color = "black"         
      ),
      name = "Mean",
      showlegend = FALSE
    )  %>%
    layout(
      title = list(
        text = "",
        font = list(
          size = 20, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        y = 0.97, x = 0.5),
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", 
        ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 15), 
        showgrid = FALSE, 
        tickangle = 0, 
        tickmode = "array", 
        tickvals = x_tick_vals, 
        range = x_range,
        tickangle = 0),
      yaxis = list(title = y_title, 
                   tickvals = y_tick_vals,
                   range = y_range,
                   ticktext = y_tick_text,
                   tickes = "outside", 
                   ticklen = 5, tickwidth = 1, showline = TRUE,
                   titlefont = list(
                     size = 26, 
                     family = "Times New Roman", 
                     color = "black", bold = TRUE),
                   tickfont = list(size = 17), 
                   showgrid = FALSE, zeroline = TRUE),
      margin = list(b = 30)
    ) 
  return(fig)
}


#=============================================================================#
#                      Figure 2.4: Upward Downward Mobility
#=============================================================================#

plot_updown_mobility <- function(data, t_num,
                                 x_tick_vals, y_tick_vals,
                                 x_range, y_range, y_tick_text,
                                 legend_x, legend_y) {
  # Fetch the relevant data:
  plot_data <-
    cross_cohort_plot_data(
      data,
      "movement_df",
      "category",
      c("Historical Upward Mobility","Historical Downward Mobility"),
      c(1)
    )
  
  fig <- 
    plot_ly() %>%
    add_ribbons(
      data = 
        plot_data[plot_data$measure == 
                    "Historical Upward Mobility", ],
      x = cohort_bins_vec,
      ymin = ~CI_bc_l,  
      ymax = ~CI_bc_u,  
      fillcolor = "rgba(0,0,0,0.1)",  
      line = list(color = "transparent"),
      showlegend = FALSE
    ) %>%
    add_ribbons(
      data = plot_data[plot_data$measure == 
                         "Historical Downward Mobility", ],
      x = cohort_bins_vec,
      ymin = ~CI_bc_l,  
      ymax = ~CI_bc_u, 
      fillcolor = "rgba(0,0,0,0.1)",  
      line = list(color = "transparent"),
      showlegend = FALSE
    ) %>%
    add_trace(
      data = plot_data[plot_data$measure == 
                         "Historical Upward Mobility", ],
      x = cohort_bins_vec,
      y = ~mean_bc,
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "black", dash = "solid"), 
      marker = list(symbol = "triangle-up", size = 10, color = "black"),
      name = paste0("Upward mobility")
    ) %>%
    add_trace(
      data = plot_data[plot_data$measure == 
                         "Historical Downward Mobility", ],
      x = cohort_bins_vec,
      y = ~mean_bc,
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "black", dash = "solid"), 
      marker = list(symbol = "square", size = 10, color = "black"),
      name = paste0("Downward mobility")
    ) %>%
    layout(
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 15), 
        showgrid = FALSE, tickangle = 0,tickmode = "array", 
        tickvals = x_tick_vals,  
        range = x_range,
        tickangle = 0),
      yaxis = list(
        title = list(text = "Probability to move", 
                     standoff = 10),automargin = TRUE,
        tickvals = y_tick_vals,
        ticktext = y_tick_text,
        range = y_range,
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 17), 
        showgrid = FALSE, zeroline = TRUE),
      margin = list( l = 30, b = 30),
      legend = list(
        orientation = "v",
        x = legend_x,  
        y = legend_y,     
        font = list(size = 20, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      )
    )
  return(fig)
}


#=============================================================================#
#                      Figure 4 & 6: Class specific MTE & IM:
#=============================================================================#


plot_measure_by_class <- function(
    data, measure, attribute_col,metric_name,y_title, 
    t_num,x_tick_vals, y_tick_vals, x_range, y_range,
    y_tick_text, legend_x, legend_y){
  # Fetch the data:
  plot_data <-
    cross_cohort_plot_data(
      data,
      measure,
      attribute_col,
      metric_name,
      t_num
    )
  # Generate the plot:
  # Measures order: put AMTE last (to draw on top)
  unique_measures <- c(unique(plot_data$measure))
  
  # ---- MANUAL STYLE BLOCKS ---------------------------------------------------

  other_colors  <- rev(c("#000000", "#1F1F1F", "#3D3D3D", "#5C5C5C", "#7A7A7A", "#999999", "#B7B7B7"))[seq_along(unique_measures)]    
  other_dashes  <- c("dash", "dot", "dashdot", "solid", "dash", "longdash", "longdashdot")[seq_along(unique_measures)]              
  other_symbols <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "circle", "square")[seq_along(unique_measures)] 
  
  
  fig <- plot_ly()
  
  # Add each non-AMTE series separately so line type / marker can differ per measure
  for(i in seq_along(unique_measures)){
    m_i <- unique_measures[i]
    dat_i <- plot_data %>% dplyr::filter(measure == m_i)
    
    fig <- fig %>%
      add_trace(
        data = dat_i,
        x = cohort_bins_vec,
        y = ~mean_bc,
        type = "scatter",
        mode = "lines+markers",
        name = m_i,
        line   = list(color = other_colors[i], dash = other_dashes[i], width = 2),   
        marker = list(color = other_colors[i], symbol = other_symbols[i], size = 7,
                      line = list(color = "#000000", width = 0.5))  
      )
  }
  
  # ---- Layout
  fig <- 
    fig %>%
    layout(
      title = list(
        text = "", font = list(
          size = 23, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        y = 0.97, y = 0.5),
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(size = 26,
                         family = "Times New Roman",
                         color = "black", bold = TRUE),
        tickfont = list(size = 14, face = "bold"), 
        showgrid = FALSE, tickangle = 0, 
        tickmode = "array", 
        tickvals = x_tick_vals, 
        range = x_range,
        tickangle = 0),
      yaxis = list(title = y_title, 
                   tickes = "inside", ticklen = 5, 
                   tickwidth = 1, showline = TRUE,
                   range = y_range,
                   tickvals = y_tick_vals,
                   ticktext = y_tick_text,
                   titlefont = list(
                     size = 26, 
                     family = "Times New Roman", 
                     color = "black", bold = TRUE),
                   tickfont = list(size = 17), 
                   showgrid = FALSE, zeroline = TRUE),
      margin = list( b = 30),
      legend = list(
        orientation = "v",
        x = legend_x,  
        y = legend_y,     
        font = list(
          size = 20, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      )
    )
  return(fig)
}

#=============================================================================#
#                      Figure 7: IM & AIM comparison
#=============================================================================#

plot_memory_contrast <- function(data, initial, t_num, y_tick_vals, y_range,
                                 y_tick_text, legend_x, legend_y, individual, 
                                 x_tick_vals,x_range){
  # Process data:
  plot_data <-
    data %>%
    filter(
      (class == "AIM") | str_starts(as.character(class), "Class")
    ) %>%
    filter(
      t %in% t_num
    ) %>%
    mutate(
      across(where(is.numeric), ~ replace_na(.x, 0))
    )
  # Measures order: put AMTE last (to draw on top)
  unique_measures <- c(unique(plot_data$class))
  other_measures  <- setdiff(unique_measures, "AIM")
  
  # ---- MANUAL STYLE BLOCKS ---------------------------------------------------
  
  other_colors  <- rev(c("#000000", "#1F1F1F", "#3D3D3D", "#5C5C5C", "#7A7A7A", "#999999", "#B7B7B7"))[seq_along(other_measures)]     
  other_dashes  <- c("dash", "dot", "dashdot", "solid", "dash", "longdash", "longdashdot")[seq_along(other_measures)]              
  other_symbols <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "circle", "square")[seq_along(other_measures)] 
  
  # AIM styles
  AIM_fill_rgba <- "rgba(0,0,0,0.12)"   
  AIM_line_col  <- "#000000"           
  AIM_line_dash <- "solid"             
  AIM_symbol    <- "square"           
  
  if(individual == TRUE){
    # ---- Build plot ------------------------------------------------------------
    fig <- plot_ly()
    
    # Add each non-AMTE series separately so line type / marker can differ per measure
    for(i in seq_along(other_measures)){
      m_i <- other_measures[i]
      dat_i <- plot_data %>% dplyr::filter(class == m_i)
      
      fig <- fig %>%
        add_trace(
          data = dat_i,
          x = ~t,
          y = ~mean_bc,
          type = "scatter",
          mode = "lines+markers",
          name = m_i,
          line   = list(color = other_colors[i], dash = other_dashes[i], width = 2),   
          marker = list(color = other_colors[i], symbol = other_symbols[i], size = 7,
                        line = list(color = "#000000", width = 0.5))                    # thin border for print clarity
        )
    }
  }else{
    # AMTE ribbon (CI)
    fig <- plot_ly()
    dat_AIM <- plot_data %>% dplyr::filter(class == "AIM")
    if(nrow(dat_AIM) > 0){
      fig <- 
        fig %>%
        add_ribbons(
          data = dat_AIM,
          x = ~t,
          ymin = ~CI_bc_l,
          ymax = ~CI_bc_u,
          fillcolor = AIM_fill_rgba,               # <<< EDIT HERE
          line = list(color = "transparent"),
          showlegend = FALSE,
          name = "AIM 95% CI"
        ) %>%
        # AMTE line + markers on top
        add_trace(
          data = dat_AIM,
          x = ~t,
          y = ~mean_bc,
          type = "scatter",
          mode = "lines+markers",
          name = "AIM",
          line   = list(color = AIM_line_col, dash = AIM_line_dash, width = 3), 
          marker = list(color = AIM_line_col, symbol = AIM_symbol, size = 9,
                        line = list(color = "#000000", width = 0.6))                 # <<< EDIT HERE (marker)
        )
    }
  }
  # ---- Layout 
  fig <- 
    fig %>%
    layout(
      title = list(text = "",
                   font = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
                   y = 0.97, x = 0.5),
      xaxis = list(
        title = "t",
        tickes = "outside", 
        ticklen = 5, 
        tickwidth = 1, 
        showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", 
          bold = TRUE),
        tickfont = list(size = 20, bold = TRUE), 
        tickvals = x_tick_vals,
        mirror = FALSE,
        showgrid = FALSE, 
        zeroline = TRUE,
        constrain = "range",
        tickangle = 0, tickmode = "array", 
        cliponaxis = TRUE,
        tickangle = 0),
      yaxis = list(
        title = "Memory", 
        mirror = FALSE,
        tickvals = y_tick_vals, # y_tick_vals
        range = y_range, # y_range
        cliponaxis = FALSE,
        ticktext = y_tick_text,
        tickes = "outside", ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
        tickfont = list(size = 17), showgrid = FALSE, zeroline = FALSE),
      legend = list(
        orientation = "v",
        x = legend_x,  
        y = legend_y, 
        constrain = "range",
        font = list(size = 25, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      ),
      margin = list( b = 30)
    )
  return(fig)
}



#=============================================================================#
#                      Figure 8: Outflow Comparison
#=============================================================================#

plot_outflow_contrast <- function(data, y_tick_vals, y_range,y_tick_text, legend_x, legend_y, x_tick_vals, x_range){
  # Process data:
  plot_data <- 
    data %>%
    mutate(
      cohort_index = as.numeric(cohort_index)
    ) %>%
    left_join(
      cohort_bins,
      by = "cohort_index"
    )
  # Measures order: put AMTE last (to draw on top)
  unique_measures <- c(unique(plot_data$class))
  
  # ---- MANUAL STYLE BLOCKS ---------------------------------------------------

  
  other_colors  <- rev(c("#000000", "#1F1F1F", "#3D3D3D", "#5C5C5C", "#7A7A7A", "#999999", "#B7B7B7"))[seq_along(unique_measures)]      
  other_dashes  <- c("dash", "dot", "dashdot", "solid", "dash", "longdash", "longdashdot")[seq_along(unique_measures)]              
  other_symbols <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "circle", "square")[seq_along(unique_measures)] 
  
  
  fill_colors <- c(
    Class1 = "rgba(183, 183, 183, 0.2)",  
    Class2 = "rgba(153, 153, 153, 0.2)",   
    Class3 = "rgba(122, 122, 122, 0.2)",    
    Class4 = "rgba(92, 92, 92, 0.2)",   
    Class5 = "rgba(61, 61, 61, 0.2)", 
    Class6 = "rgba(31, 31, 31, 0.2)",
    AIM = "rgba(225, 87, 89, 0.2)"
  )

    # ---- Build plot ------------------------------------------------------------
    fig <- plot_ly()
    
    # Add each non-AMTE series separately so line type / marker can differ per measure
    for(i in seq_along(unique_measures)){
      m_i <- unique_measures[i]
      dat_i <- plot_data %>% dplyr::filter(class == m_i)
      
      fig <- fig %>%
        add_trace(
          data = dat_i,
          x = ~cohort_bins_vec,
          y = ~mean_bc,
          type = "scatter",
          mode = "lines+markers",
          name = m_i,
          line   = list(color = other_colors[i], dash = other_dashes[i], width = 2),   
          marker = list(color = other_colors[i], symbol = other_symbols[i], size = 7,
                        line = list(color = "#000000", width = 0.5))                    # thin border for print clarity
        )
    }
   
   for (i in seq_along(unique_measures)) {
     m_i <- unique_measures[i]
     dat_i <- plot_data %>% dplyr::filter(class == m_i)
     fig <- 
      fig %>% add_ribbons(
      data = dat_i,
      x = ~cohort_bins_vec,
      ymin = ~CI_bc_l,
      ymax = ~CI_bc_u,
      fillcolor = fill_colors[i],
      line = list(color = "transparent"),
      name = m_i,
      showlegend = FALSE
    )
  }
  # ---- Layout 
  fig <- 
    fig %>%
    layout(
      title = list(text = "",
                   font = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
                   y = 0.97, x = 0.5),
      xaxis = list(
        title = "t",
        tickes = "outside", 
        ticklen = 5, 
        tickwidth = 1, 
        showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", 
          bold = TRUE),
        tickfont = list(size = 14), 
        tickvals = x_tick_vals,
        mirror = FALSE,
        showgrid = FALSE, 
        zeroline = TRUE,
        constrain = "range",
        tickangle = 0, tickmode = "array", 
        cliponaxis = TRUE,
        tickangle = 0),
      yaxis = list(
        title = "Outflow", 
        mirror = FALSE,
        tickvals = y_tick_vals, # y_tick_vals
        range = y_range, # y_range
        cliponaxis = FALSE,
        ticktext = y_tick_text,
        tickes = "outside", ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
        tickfont = list(size = 17), showgrid = FALSE, zeroline = FALSE),
      legend = list(
        orientation = "v",
        x = legend_x,  
        y = legend_y, 
        constrain = "range",
        font = list(size = 16, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      ),
      margin = list( b = 30)
    )
  return(fig)
}



#==============================================================================#
#                      Figure 8: IM & AIM comparison
#==============================================================================#


marg_class_comparison_evolution <- function(data, size){
  # Measures order: put AMTE last (to draw on top)
  unique_measures <- sort(unique(c(data$status)))
  # ---- MANUAL STYLE BLOCKS ---------------------------------------------------
  
  other_colors  <- rev(c("#000000", "#1F1F1F", "#3D3D3D", "#5C5C5C", "#7A7A7A", "#999999", "#B7B7B7"))[seq_along(unique_measures)]    
  other_dashes  <- c("dash", "dot", "dashdot", "solid", "dash", "longdash", "longdashdot")[seq_along(unique_measures)]            
  other_symbols <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "circle", "square")[seq_along(unique_measures)]
  
  # ---- Build plot ------------------------------------------------------------
  fig <- plot_ly()
  
  # Add each non-AMTE series separately so line type / marker can differ per measure
  for(i in seq_along(unique_measures)){
    m_i <- unique_measures[i]
    dat_i <- data %>% dplyr::filter(status == m_i)
    
    fig <- 
      fig %>%
      add_trace(
        data = dat_i,
        x = ~cohort,
        y = ~pct_change,
        type = "scatter",
        mode = "lines+markers",
        name = m_i,
        line   = list(color = other_colors[i], dash = other_dashes[i], width = 2),  
        marker = list(color = other_colors[i], symbol = other_symbols[i], size = 7,
                      line = list(color = "#000000", width = 0.5))                  
      )
  }
  
  # ---- Layout  --------------------------------
  fig <- 
    fig %>%
    layout(
      title = list(text = "",
                   font = list(size = 23, family = "Times New Roman", color = "black", bold = TRUE),
                   y = 0.97, y = 0.5),
      xaxis = list(title = "Birth cohorts",
                   tickes = "outside", ticklen = 5, tickwidth = 1, showline = TRUE,
                   titlefont = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
                   tickfont = list(size = 13, face = "bold"), showgrid = FALSE, tickangle = 0, 
                   tickmode = "array", tickvals = seq(1945, 1992, by = 3), 
                   range = c(1944, 1993),
                   ticktext = seq(1945, 1992, by = 3), tickangle = 0),
      yaxis = list(title = "Marginal Change", 
                   tickes = "inside", 
                   ticklen = 5, 
                   tickwidth = 1, 
                   showline = TRUE,
                   titlefont = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
                   tickfont = list(size = 17), showgrid = FALSE, zeroline = TRUE),
      margin = list( b = 30),
      legend = list(
        orientation = "v",
        x = 0.985,  
        y = 0.95,     
        font = list(size = 20, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      )
    )
  
  return(fig)
}

plot_contraction_speed <- function(
    data, y_title, t_index, individual, x_tick_vals, y_tick_vals, x_range,
    y_range, y_tick_text, legend_x, legend_y){
  # Fetch the relevant data:
  plot_data <-
    cross_cohort_plot_data(
      cohort_lsts,
      "memory_df",
      "class",
      c(paste0("Class",seq(1,5)),"AIM"),
      c(1,2,3,4)
    ) %>%
    mutate(
      across(where(is.numeric), ~ replace_na(.x, 0))
    ) %>%
    group_by(
      cohort_index,
      measure
    ) %>%
    arrange(
      t,
      .by_group = TRUE
    ) %>%
    mutate(
      contraction = lag(mean_bc) - mean_bc
    ) %>%
    dplyr::select(
      contraction,
      everything()
    ) %>%
    filter(
      t != 1
    ) %>%
    ungroup()
  unique_measures <- c(c(unique(plot_data$measure))[-1],"AIM")
  other_measures  <- setdiff(unique_measures, "AIM")
  plot_data_lst <- split(plot_data,plot_data$t)
  plot_data <- plot_data_lst[[t_index]]
  other_colors  <- rev(c("#000000", "#1F1F1F", "#3D3D3D", "#5C5C5C", "#7A7A7A", "#999999", "#B7B7B7"))[seq_along(other_measures)]     
  other_dashes  <- c("dash", "dot", "dashdot", "solid", "dash", "longdash", "longdashdot")[seq_along(other_measures)]              
  other_symbols <- c("circle", "square", "diamond", "triangle-up", "triangle-down", "circle", "square")[seq_along(other_measures)] 
  
  # AIM styles
  AIM_fill_rgba <- "rgba(0,0,0,0.12)"   
  AIM_line_col  <- "#000000"           
  AIM_line_dash <- "solid"             
  AIM_symbol    <- "square" 
  
  if(individual == TRUE){
    fig <- plot_ly()
    
    for(i in seq_along(other_measures)){
      m_i <- other_measures[i]
      dat_i <- plot_data %>% dplyr::filter(measure == m_i)
      
      fig <- fig %>%
        add_trace(
          data = dat_i,
          x = cohort_bins_vec,
          y = ~contraction,
          type = "scatter",
          mode = "lines+markers",
          name = m_i,
          line   = list(
            color = other_colors[i], dash = other_dashes[i], width = 2),
          marker = list(color = other_colors[i], 
                        symbol = other_symbols[i], size = 7,
                        line = list(color = "#000000", width = 0.5))                  
        )
    }
    fig
  }else{
    fig <- plot_ly()
    dat_AIM <- plot_data %>% dplyr::filter(measure == "AIM")
    if(nrow(dat_AIM) > 0){
      fig <- 
        fig %>%
        add_trace(
          data = dat_AIM,
          x = cohort_bins_vec,
          y = ~contraction,
          type = "scatter",
          mode = "lines+markers",
          name = "AIM",
          line   = 
            list(color = AIM_line_col, dash = AIM_line_dash, width = 3), 
          marker = list(color = AIM_line_col, symbol = AIM_symbol, 
                        size = 9,
                        line = list(color = "#000000", width = 0.6))
        )
    }
  }
  
  fig <- 
    fig %>%
    layout(
      title = list(
        text = "",
        font = list(size = 26, family = "Times New Roman", 
                    color = "black", bold = TRUE),
        y = 0.97, x = 0.5),
      xaxis = list(
        title = "t",
        tickes = "outside", 
        ticklen = 5, 
        tickwidth = 1, 
        showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", 
          bold = TRUE),
        tickfont = list(size = 15), 
        tickvals = x_tick_vals,
        mirror = FALSE,
        showgrid = FALSE, 
        zeroline = TRUE,
        constrain = "range",
        tickangle = 0, tickmode = "array", 
        cliponaxis = TRUE,
        tickangle = 0),
      yaxis = list(
        title = y_title, 
        mirror = FALSE,
        tickvals = y_tick_vals, 
        range = y_range, 
        cliponaxis = FALSE,
        ticktext = y_tick_text,
        tickes = "outside", ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(size = 26, family = "Times New Roman", color = "black", bold = TRUE),
        tickfont = list(size = 17), showgrid = FALSE, zeroline = FALSE),
      legend = list(
        orientation = "v",
        x = legend_x,  
        y = legend_y, 
        constrain = "range",
        font = list(size = 16, family = "Times New Roman", color = "black"),
        xanchor = "left",
        yanchor = "top"
      ),
      margin = list( b = 30)
    )
  return(fig)
}



#==============================================================================#
#                      Figure 9: Log multiplicative model:
#==============================================================================#

plot_log_multiplicative <- function(data, bin_width, year_l, year_u,
                                    age_l, age_u, cohort_bins,
                                    tick_vals, range){
  # Generate the regression data:
  if(bin_width == 5){
    gss_data_bin <-
      bind_rows(
        lapply(
          cross_tm_df(
            data %>%
              filter(
                cohort >= year_l & cohort <= year_u
              ),age_l,age_u, bin_width),
          function(df){
            df <-
              df %>%
              group_by(
                bins,
                status_c,
                status_p
              ) %>%
              summarise(
                Freq = n(),
                .groups = "drop"
              )
          }
        ),
        .id = "cohort_id"
      ) %>%
      mutate(
        bins = as.factor(bins),
        Freq = as.numeric(Freq)
      ) %>%
      mutate(
        bins = case_when(
          bins %in% c("[1983, 1988)","[1988, 1991]") ~ "[1983,1990]",
          bins %in% c("[1970, 1975)","[1975, 1978]") ~ "[1970,1977]",
          .default = bins
        ))
  }else{
    gss_data_bin <-
      data %>%
      filter(
        age >= 25 & age <= 55
      ) %>%
      filter(
        cohort >= 1945  & cohort <= 1990
      ) %>%
      mutate(
        bins = case_when(
          cohort >= 1945 & cohort <= 1954 ~ "[1945,1954]",
          cohort >= 1955 & cohort <= 1964 ~ "[1955,1964]",
          cohort >= 1965 & cohort <= 1974 ~ "[1965,1974]",
          cohort >= 1975 & cohort <= 1990 ~ "[1975,1990]",
          .default = NA
        )
      ) %>%
      group_by(
        bins,
        status_c,
        status_p
      )  %>%
      summarise(
        Freq = n(),
        .groups = "drop"
      ) %>%
      mutate(
        bins = as.factor(bins),
        Freq = as.numeric(Freq)
      ) 
  }
  
  gss_occ10_tab <- xtabs(Freq ~ status_c + status_p + bins, data = gss_data_bin)
  # Fit the log-multiplicative model:
  set.seed(98625)
  full_interaction2 <- unidiff(gss_occ10_tab, diagonal = "included")
  # Generate the plot data:
  phi_df <-
    full_interaction2$unidiff$layer$qvframe %>%
    mutate(
      mean_phi = exp(estimate),
      se_phi = mean_phi * quasiSE,
      se_phi = quasiSE
    ) %>%
    mutate(
      phi_CI_l = mean_phi - 1.96 * se_phi,
      phi_CI_u = mean_phi + 1.96 * se_phi,
      err_up = phi_CI_u - mean_phi,
      err_dn = mean_phi - phi_CI_l,
      cohort = cohort_bins
    )
  # Generate the plot:
  log_multiplicative_fig <- 
    plot_ly(
      phi_df,
      x = ~cohort,
      y = ~mean_phi,
      type = "scatter",
      mode = "lines+markers",
      line = list(width = 3, color = "black"),
      marker = list(size = 9, color = "black" ),
      error_y = list(
        type = "data",
        symmetric = FALSE,
        array = ~err_up,
        arrayminus = ~err_dn,
        thickness = 1.5,   
        width = 6,
        color = "grey40"
      ),
      customdata = ~cbind(phi_CI_l, phi_CI_u),
      hovertemplate = paste(
        "Cohort: %{x}",
        "<br>Estimate: %{y:.3f}",
        "<br>95% CI: [%{customdata[0]:.3f}, %{customdata[1]:.3f}]",
        "<extra></extra>"
      )
    ) %>%
    layout(
      xaxis = list(
        title = "Birth cohort",
        categoryorder = "array",
        categoryarray = levels(phi_df$cohort),
        tickfont = list(size = 14),
        titlefont = list(size = 16)
      ),
      yaxis = list(
        title = "Estimate",
        tickfont = list(size = 14),
        titlefont = list(size = 16),
        zeroline = FALSE
      ),
      margin = list(l = 80, r = 30, t = 30, b = 70)
    ) %>%
    layout(
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 14), 
        showgrid = FALSE, 
        tickangle = 0,
        tickmode = "array", 
        tickangle = 0),
      yaxis = list(
        title = list(text = "log-multiplicative  layer effect", 
                     standoff = 10),automargin = TRUE,
        tickvals = tick_vals,
        range = range,
        #tickvals = seq(0.45,1.35,0.1),
        #range = c(0.45,1.35),
        #tickvals = seq(0.65,1.75,0.1),
        #range = c(0.65,1.75),
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 23, family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 17), 
        showgrid = FALSE, zeroline = TRUE),
      margin = list( l = 30, b = 30)
    ) %>% layout(width = 1000, height = 600)
  return(log_multiplicative_fig)
}


plot_log_multiplicative_weight <- function(data, bin_width, year_l, year_u,
                                    age_l, age_u, cohort_bins){
  # Generate the regression data:
  if(bin_width == 5){
  gss_data_bin <-
    bind_rows(
      lapply(
        cross_tm_df(
          data %>%
            filter(
              cohort >= year_l & cohort <= year_u
            ),age_l,age_u, bin_width),
        function(df){
          df <-
            df %>%
            group_by(
              bins,
              status_c,
              status_p
            ) %>%
            summarise(
              Freq = round(sum(wn)),
              .groups = "drop"
            )
        }
      ),
      .id = "cohort_id"
    ) %>%
    mutate(
      bins = as.factor(bins),
      Freq = as.numeric(Freq)
    ) %>%
    mutate(
      bins = case_when(
        bins %in% c("[1984, 1989)","[1989, 1991]") ~ "[1984,1990]",
        .default = bins
      ))}else{
        gss_data_bin <-
          data %>%
          filter(
            age >= age_l & age <= age_u
          ) %>%
          filter(
            cohort >= year_l  & cohort <= year_u
          ) %>%
          mutate(
            bins = case_when(
              cohort >= 1945 & cohort <= 1954 ~ "[1945,1954]",
              cohort >= 1955 & cohort <= 1964 ~ "[1955,1964]",
              cohort >= 1965 & cohort <= 1974 ~ "[1965,1974]",
              cohort >= 1975 & cohort <= 1990 ~ "[1975,1990]",
              .default = NA
            )
          ) %>%
          group_by(
            bins,
            status_c,
            status_p
          )  %>%
          summarise(
            Freq = round(sum(wn)),
            .groups = "drop"
          ) %>%
          mutate(
            bins = as.factor(bins),
            Freq = as.numeric(Freq)
          ) 
      }
  gss_occ10_tab <- xtabs(Freq ~ status_c + status_p + bins, data = gss_data_bin)
  # Fit the log-multiplicative model:
  set.seed(98625)
  full_interaction2 <- unidiff(gss_occ10_tab, diagonal = "included")
  # Generate the plot data:
  phi_df <-
    full_interaction2$unidiff$layer$qvframe %>%
    mutate(
      mean_phi = exp(estimate),
      se_phi = mean_phi * quasiSE,
      se_phi = quasiSE
    ) %>%
    mutate(
      phi_CI_l = mean_phi - 1.96 * se_phi,
      phi_CI_u = mean_phi + 1.96 * se_phi,
      err_up = phi_CI_u - mean_phi,
      err_dn = mean_phi - phi_CI_l,
      cohort = cohort_bins
    )
  # Generate the plot:
  log_multiplicative_fig <- 
    plot_ly(
      phi_df,
      x = ~cohort,
      y = ~mean_phi,
      type = "scatter",
      mode = "lines+markers",
      line = list(width = 3, color = "black"),
      marker = list(size = 9, color = "black" ),
      error_y = list(
        type = "data",
        symmetric = FALSE,
        array = ~err_up,
        arrayminus = ~err_dn,
        thickness = 1.5,   
        width = 6,
        color = "grey40"
      ),
      customdata = ~cbind(phi_CI_l, phi_CI_u),
      hovertemplate = paste(
        "Cohort: %{x}",
        "<br>Estimate: %{y:.3f}",
        "<br>95% CI: [%{customdata[0]:.3f}, %{customdata[1]:.3f}]",
        "<extra></extra>"
      )
    ) %>%
    layout(
      xaxis = list(
        title = "Birth cohort",
        categoryorder = "array",
        categoryarray = levels(phi_df$cohort),
        tickfont = list(size = 14),
        titlefont = list(size = 16)
      ),
      yaxis = list(
        title = "Estimate",
        tickfont = list(size = 14),
        titlefont = list(size = 16),
        zeroline = FALSE
      ),
      margin = list(l = 80, r = 30, t = 30, b = 70)
    ) %>%
    layout(
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 14), 
        showgrid = FALSE, 
        tickangle = 0,
        tickmode = "array", 
        tickangle = 0),
      yaxis = list(
        title = list(text = "log-multiplicative  layer effect", 
                     standoff = 10),automargin = TRUE,
        tickvals = seq(0.35,1.35,0.1),
        range = c(0.35,1.35),
        #tickvals = seq(0.65,1.75,0.1),
        #range = c(0.65,1.75),
        tickes = "outside", ticklen = 5, 
        tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 23, family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 17), 
        showgrid = FALSE, zeroline = TRUE),
      margin = list( l = 30, b = 30)
    ) %>% layout(width = 1000, height = 600)
  return(log_multiplicative_fig)
}


#==============================================================================#
#                      Figure 10: Mean First Passage Time:
#==============================================================================#

MFP_class_comparison_new <- function(cohort_lsts, attribute,q1, q2, q3, class_specific, class){
  # Prepare the Data:
  modify_lsts <- 
    mclapply(
      cohort_lsts,
      function(cohort_lst) {
        # Modify the "AEO" element of the cohort_lst
        cohort_lst[["MFP"]] <- 
          cohort_lst[["MFP"]] %>%
          mutate(
            parent_num = str_extract(parent_status, "\\d+"),
            children_num = str_extract(children_status, "\\d+")
          ) %>%
          mutate(
            measure =
              paste0(parent_num,"-",children_num)
          ) %>%
          mutate(
            category = attribute,
            mean_bc = round(mean_bc)
          )
      }
    )
  attribute_df <- 
    bind_rows(modify_lsts, .id = "cohort_index") %>%
    mutate(
      cohort_index = as.numeric(cohort_index)
    )
  attribute_MFP <- 
    left_join(
      attribute_df,
      cohort_bins,
      by = "cohort_index"
    ) %>%
    mutate(
      cohort_bins_vec = as.character(cohort_bins_vec)
    )
  if(class_specific == TRUE){
    attribute_MFP <-
      attribute_MFP %>%
      filter(
        parent_status == class
      )
  }else{
    attribute_MFP
  }
  
  # Generate the label:
  very_low_label    <- paste0("Very Fast transit: [1, ", q1, ")")
  low_label <-  paste0("Fast transit: [", q1, ", ", q2, ")")
  medium_label <- paste0("Moderate transit: [", q2, ", ", q3, ")")
  high_label   <- paste0("Slow transit: [", q3, ",", max(attribute_MFP$mean_bc), ")")
  
  # Keep the Key Information:
  attribute_MFP_key <- 
    attribute_MFP %>%
    select(
      measure,
      parent_status,
      children_status,
      cohort_index,
      cohort_bins_vec,
      mean_bc
    ) %>%
    mutate(
      cohort_bins_vec = 
        factor(
          as.character(cohort_bins_vec), 
          levels = as.character(sort(unique(cohort_bins_vec)))),
      direction = ifelse(
        parent_status < children_status, "Downward transition", "Upward transition")
    )  %>%
    mutate(
      mean_bc_group_manual = case_when(
        mean_bc < q1            ~  very_low_label,
        mean_bc >= q1 & mean_bc < q2 ~ low_label,
        mean_bc >= q2 & mean_bc <= q3 ~ medium_label,
        mean_bc > q3            ~ high_label
      )
    ) %>%
    mutate(
      mean_bc_group_manual = factor(
        mean_bc_group_manual,
        levels = c(
          very_low_label,
          low_label,
          medium_label,
          high_label
        )
      )
    )
  # Generate the plot:
  fig_all <-
    ggplot(
      attribute_MFP_key, 
      aes(
        x = cohort_bins_vec, 
        y = measure, 
        fill = mean_bc_group_manual)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(
      values = setNames(
        c("#F0F0F0", "#BDBDBD", "#737373", "#252525"),
        c(very_low_label, low_label, medium_label, high_label)
      ),
      name = "MFP interval (# of generations):"
    ) +
    scale_x_discrete(
      breaks = 
        levels(attribute_MFP_key$cohort_bins_vec)[
          seq(1, length(levels(attribute_MFP_key$cohort_bins_vec)), by = 2)],
      expand = c(0, 0)
    )+
    labs(
      title = "Mean First Passage Time Cross-Cohort Comparison by Group",
      x = "Birth Cohorts",
      y = "Father–Child class transitions"
    ) +
    theme_minimal(base_family = "Times") +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      legend.key.size = unit(0.3, "cm"),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.2, "cm")  # Controls length of tick marks
    )
  
  fig_updown <-
    ggplot(
      attribute_MFP_key,
      aes(
        x = cohort_bins_vec, y = measure, fill = mean_bc_group_manual)) +
    geom_tile(color = "white") +
    facet_wrap(~ direction, scales = "free_y", ncol = 1, strip.position = "top") +
    scale_fill_manual(
      values = setNames(
        c("#F0F0F0", "#BDBDBD", "#737373", "#252525"),
        c(very_low_label, low_label, medium_label, high_label)
      ),
      name = "MFP interval (# of generations):"
    ) +
    scale_x_discrete(
      breaks = 
        levels(attribute_MFP_key$cohort_bins_vec)[
          seq(1, 
              length(levels(attribute_MFP_key$cohort_bins_vec)), 
              by = 2)],
      expand = c(0, 0)
    ) +
    labs(
      x = "Birth cohorts",
      y = "Origin–Destination class transitions",
      title = "Mobility Patterns by Direction Across Birth Cohorts"
    ) +
    labs(
      x = "Birth cohorts",
      y = "Origin–Destination class transitions",
      title = "Mobility Patterns by Direction Across Birth Cohorts"
    ) +
    theme_minimal(base_size = 10, base_family = "Times") +
    theme(
      legend.position = "top",                     
      legend.direction = "horizontal",                
      legend.box = "horizontal",      
      axis.text.x = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(size = 20, face = "bold"),
      axis.ticks.length = unit(0.2, "cm"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key.size = unit(0.3, "cm"),
      axis.title.x = element_text(face = "bold", size = 22,margin = margin(t = 10)),
      axis.title.y = element_text(face = "bold", size = 22, margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
    )
  return(
    list(
      all = fig_all,
      updown = fig_updown
    )
  )
}


#==============================================================================#
#                      Figure 11: Plot the d'2/d'1
#==============================================================================#

plot_d_prime_contraction <- function(data,  y_tick_vals, y_range,
                                     y_tick_text,  x_tick_vals,x_range,
                                     y_title){
  ## Step 1: Fetch a new data:
  d_prime_lsts <-
    lapply(
      data,
      function(lst){
        lst$d_prime_df
      }
    )
  
  ## Step 2: Generate the new data:
  plot_data <-
    bind_rows(
      lapply(
        d_prime_lsts,
        function(lst){
          df <-
            lst %>%
            mutate(
              d_prime_ratio_bc = 
                mean_bc[match(2, t)]/ mean_bc[match(1, t)],
              d_prime_ratio_CI_u = 
                CI_bc_u[match(2, t)]/ CI_bc_u[match(1, t)],
              d_prime_ratio_CI_l = 
                CI_bc_l[match(2, t)]/ CI_bc_l[match(1, t)]
            ) %>%
            dplyr::select(
              d_prime_ratio_bc,
              d_prime_ratio_CI_u,
              d_prime_ratio_CI_l
            ) %>%
            unique()
        }
      ),
      .id = "cohort_index"
    ) %>%
    mutate(
      cohort_index = as.numeric(cohort_index)
    ) %>%
    left_join(
      cohort_bins,
      by = "cohort_index"
    )
  
  ## Step 3: Generate the plot:
  fig <-
    plot_ly() %>%
    add_ribbons(
      data = plot_data,
      x = cohort_bins_vec,
      ymin = ~d_prime_ratio_CI_l,
      ymax = ~d_prime_ratio_CI_u,
      fillcolor = "rgba(0,0,0,0.1)",
      line = list(color = "transparent"), 
      name = "Confidence Interval",
      showlegend = FALSE
    ) %>%
    add_trace(
      data = plot_data,
      x = cohort_bins_vec,
      y = ~d_prime_ratio_bc,
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "black", dash = "solid", width = 2), 
      marker = list(
        symbol = "circle", size = 10, color = "black"         
      ),
      name = "Mean",
      showlegend = FALSE
    )  %>%
    layout(
      title = list(
        text = "",
        font = list(
          size = 20, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        y = 0.97, x = 0.5),
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", 
        ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 15), 
        showgrid = FALSE, 
        tickangle = 0, 
        tickmode = "array", 
        tickvals = x_tick_vals, 
        range = x_range,
        tickangle = 0),
      yaxis = list(
        title = y_title, 
        tickvals = y_tick_vals,
        range = y_range,
        ticktext = y_tick_text,
        tickes = "outside", 
        ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 17), 
        showgrid = FALSE, zeroline = TRUE),
      margin = list(b = 30)
    ) 
  return(fig)
}


#==============================================================================#
#                      Figure 11: Plot the d'2/d'1
#==============================================================================#

plot_d <- function(data, t_num, y_title, size, 
                   x_tick_vals, y_tick_vals,
                   x_range, y_range, y_tick_text) {
  # Fetch the relevant data:
  # fetch the data that is relevant:
  plot_data <-
    bind_rows(
      lapply(
        data,
        function(lst){
          df <-
            lst$memory_df %>%
            filter(
              t == t_num
            ) %>%
            filter(
              class %in% paste0("Class",seq(1,size))
            ) %>%
            slice_max(mean_bc, n = 1, with_ties = FALSE)
        }
      ),
      .id = "cohort_index"
    ) %>%
    mutate(
      cohort_index = as.numeric(cohort_index)
    ) %>%
    left_join(
      cohort_bins,
      by = "cohort_index"
    )
  
  # Generate the plot:
  fig <-
    plot_ly() %>%
    add_ribbons(
      data = plot_data,
      x = cohort_bins_vec,
      ymin = ~CI_bc_l,
      ymax = ~CI_bc_u,
      fillcolor = "rgba(0,0,0,0.1)",
      line = list(color = "transparent"), 
      name = "Confidence Interval",
      showlegend = FALSE
    ) %>%
    add_trace(
      data = plot_data,
      x = cohort_bins_vec,
      y = ~mean_bc,
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "black", dash = "solid", width = 2), 
      marker = list(
        symbol = "circle", size = 10, color = "black"         
      ),
      name = "Mean",
      showlegend = FALSE
    )  %>%
    layout(
      title = list(
        text = "",
        font = list(
          size = 20, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        y = 0.97, x = 0.5),
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", 
        ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 15), 
        showgrid = FALSE, 
        tickangle = 0, 
        tickmode = "array", 
        tickvals = x_tick_vals, 
        range = x_range,
        tickangle = 0),
      yaxis = list(title = y_title, 
                   tickvals = y_tick_vals,
                   range = y_range,
                   ticktext = y_tick_text,
                   tickes = "outside", 
                   ticklen = 5, tickwidth = 1, showline = TRUE,
                   titlefont = list(
                     size = 26, 
                     family = "Times New Roman", 
                     color = "black", bold = TRUE),
                   tickfont = list(size = 17), 
                   showgrid = FALSE, zeroline = TRUE),
      margin = list(b = 30)
    ) 
  return(fig)
}


#==============================================================================#
#                      Figure 12: Plot the Delta AIM
#==============================================================================#

plot_AIM_change <- function(
    plot_data, t, y_title, x_tick_vals, 
    mean_col, CI_u_col, CI_l_col,
    y_tick_vals, x_range, y_range, y_tick_text){ 
  
  fig <-
    plot_ly() %>%
    add_ribbons(
      data = plot_data[[t]],
      x = cohort_bins_vec,
      ymin = ~get(CI_l_col),
      ymax = ~get(CI_u_col),
      fillcolor = "rgba(0,0,0,0.1)",
      line = list(color = "transparent"), 
      name = "Confidence Interval",
      showlegend = FALSE
    ) %>%
    add_trace(
      data = plot_data[[t]],
      x = cohort_bins_vec,
      y = ~get(mean_col),
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "black", dash = "solid", width = 2),
      marker = list(
        symbol = "circle", size = 10, color = "black"         
      ),
      name = "Mean",
      showlegend = FALSE
    )  %>%
    layout(
      title = list(
        text = "",
        font = list(
          size = 20, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        y = 0.97, x = 0.5),
      xaxis = list(
        title = "Birth cohorts",
        tickes = "outside", 
        ticklen = 5, tickwidth = 1, showline = TRUE,
        titlefont = list(
          size = 26, 
          family = "Times New Roman", 
          color = "black", bold = TRUE),
        tickfont = list(size = 15), 
        showgrid = FALSE, 
        tickangle = 0, 
        tickmode = "array", 
        tickvals = x_tick_vals, 
        range = x_range,
        tickangle = 0),
      yaxis = list(title = y_title, 
                   tickvals = y_tick_vals,
                   range = y_range,
                   ticktext = y_tick_text,
                   tickes = "outside", 
                   ticklen = 5, tickwidth = 1, 
                   showline = TRUE,
                   titlefont = list(
                     size = 26, 
                     family = "Times New Roman", 
                     color = "black", bold = TRUE),
                   tickfont = list(size = 17), 
                   showgrid = FALSE, zeroline = TRUE),
      margin = list(b = 30)
    ) 
  return(fig)
}



















