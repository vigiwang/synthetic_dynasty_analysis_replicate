####### Mobility Measure Functions ############
####### @Copyright: Weiqi Wang

####================ Function 1: Calculate Transition Matrices  ============####

###### Section 1: Convert Raw Counts to Probabilities:

calculate_percentage <- function(row) {
  row_percentage <- row / sum(row)
  return(row_percentage)
}


trans_calc <- function(data){
  trans_matrix <- t(apply(data[, 1:ncol(data)], 1, calculate_percentage))
  return(trans_matrix)
}


####================ Function 2: Sample Dataframe  ============####

sample_data <- function(gss_data) {
  data <-
    gss_data %>%
    dplyr::select(
      status_p,
      status_c
    )
  sampled_data <- 
    data[sample(nrow(data), size = nrow(data), replace = TRUE), ]
  return(sampled_data)
}

####================ Function 3: Calculate Marginal Distribution  ==========####

marginal_dist <- function(tm_num, size){
  # Calculate the marginal distribution of Children:
  marginal_dist_c <-
    tm_num %>%
    colSums()
  marginal_dist_c_perc <- 
    marginal_dist_c[1:size] / marginal_dist_c["rowsum"]
  assert_that(
    abs(sum(marginal_dist_c_perc) - 1) < 1e-15
  )
  
  # Calculate the marginal distribution of Father:
  marginal_dist_f <-
    tm_num %>%
    pull(validation)
  marginal_dist_f_perc <- 
    marginal_dist_f / sum(marginal_dist_f)
  assert_that(
    abs(sum(marginal_dist_f_perc) - 1) < 1e-15
  )
  reslt_lst <- list()
  reslt_lst[["marginal_dist_c"]] <- marginal_dist_c_perc
  reslt_lst[["marginal_dist_f"]] <- marginal_dist_f_perc
  return(reslt_lst)
}

####========= Function 4: Calculate the Movement Mobility Measures  ==========####

compute_movement_mobility <- function(tm_perc, marginal_dist_f0,t,size){
  assert_that(
    length(marginal_dist_f0) == nrow(tm_perc)
  )
  # Section 1: Calculate Historical Mobility:
  # Iterate the transition matrix to get the marginal distribution of father 
  #up to time t:
  for (i in seq(1,t)){
    marginal_dist_ft <- marginal_dist_f0 %*% as.matrix(tm_perc[,1:size])
    marginal_dist_ft_pre <- marginal_dist_f0
    marginal_dist_f0 <- marginal_dist_ft
  }
  # Calculate the historical mobility:
  tm_temp <- t(marginal_dist_ft_pre) * tm_perc[,1:size]
  off_diags <- sum(tm_temp[row(tm_temp) == col(tm_temp)])
  historical_mobility <- 1 - off_diags
  ## calculate historical upward mobility:
  matrix_up_ind <- matrix(0, ncol = size, nrow = size)
  matrix_up_ind[col(matrix_up_ind) > row(matrix_up_ind)] <- 1
  matrix_up <- matrix_up_ind %*% as.matrix(tm_temp)
  historical_upward_mobility <- sum(matrix_up[row(matrix_up) == col(matrix_up)])
  ## calculate historical downward mobility:
  matrix_down_ind <- matrix(0, ncol = size, nrow = size)
  matrix_down_ind[col(matrix_down_ind) < row(matrix_down_ind)] <- 1
  matrix_down <- matrix_down_ind %*% as.matrix(tm_temp)
  historical_downward_mobility <- sum(matrix_down[row(matrix_down) == col(matrix_down)])
  
  # Section 2: Calculate Structural Mobility:
  marginal_dist_ft_next <- marginal_dist_ft %*% as.matrix(tm_perc[,1:size])
  structural_mobility <- 0.5 * (sum(abs(marginal_dist_ft - marginal_dist_ft_pre)))
  
  # Section 3: Calculate the Exchange Mobility:
  exchange_mobility <- historical_mobility - structural_mobility
  
  return(
    list(historical_mobility = historical_mobility, 
         historical_upward_mobility = historical_upward_mobility,
         historical_downward_mobility = historical_downward_mobility,
         structural_mobility = structural_mobility,
         exchange_mobility = exchange_mobility))
}



####================= Function 5: Matrix Power Calculation  ===============####

matrix_power = \(x, n) Reduce(`%*%`, replicate(n, x, simplify = FALSE))


####===============  Function 6: Calculate TV between 2 rows ===============####

row_difference <- function(row_index1, row_index2, tm_perc){
  return(0.5 * sum(abs(tm_perc[row_index1,] - tm_perc[row_index2,]), na.rm = TRUE))
}


#====================== Convert tm_num to Probability Matrix ==================#
convert_obs_to_prob <- function(tm_num){
  all_obs <-
    tm_num %>%
    summarise(
      all_obs = sum(rowsum)
    ) %>%
    pull()
  tm_num_prob <-
    tm_num/all_obs 
  
  tm_num_prob <-
    tm_num_prob %>%
    dplyr::select(
      -rowsum,
      -validation
    )
  return(tm_num_prob)
}

####===============  Function 7: Calculate d prime   ======================####

calculate_d_prime <- function(marginal_dist_f0, tm_perc, size, steady_state, t){
  # Step 1: Decide whether to use d_prime to measure pure exchange mobility or 
  #         steady_state:
  if(steady_state == FALSE){
    # Raise the transition matrix to a power of t:
    Pt <- matrix_power(as.matrix(tm_perc[,1:size]),t)
    mat <- t(marginal_dist_f0) * as.data.frame(Pt)
    print(mat)
  }else{
    mat <- matrix_power(as.matrix(tm_perc[,1:size]),t)
  }
  
  # Step 2: Define a function to calculate row difference:
  row_difference <- function(row_index1, row_index2, tm_perc){
    total_variation <- 
      0.5 * sum(abs(tm_perc[row_index1,] - tm_perc[row_index2,]), na.rm = TRUE)
    return(total_variation)
  }
  
  # Step 3: Generate the list of combination between rows:
  combinations <- t(combn(1:size, 2))
  combinations_list <- split(combinations, seq(nrow(combinations)))
  combination_index_list <- 
    lapply(
      combinations_list,
      function(sub_lst){
        index <- paste0(sub_lst[1],"-",sub_lst[2])
      }
    )
  difference_lst <-
    lapply(
      combinations_list, 
      function(x){
        row_difference(x[1],x[2],mat)
      }
    )
  d_prime <- max(unlist(difference_lst))
  d_prime_index <- which.max(unlist(difference_lst))
  d_prime_combination <- combination_index_list[[as.numeric(d_prime_index)]]
  return(list(
    d_prime = d_prime, 
    mat = mat,
    d_prime_combination = d_prime_combination
  ))
}


####================  Function 8: Calculate steady state   ===================####

calculate_the_steady_state_tm <- function(tm_perc, size, tol){
  row_difference <- 100
  t <- 1
  while(row_difference > tol & t < 50){
    d_prime_result <-  calculate_d_prime(NULL, tm_perc, size, TRUE, t)
    row_difference <- d_prime_result$d_prime
    p_infty <- d_prime_result$mat
    t <- t + 1
  }
  return(list(t = t-1, row_difference = row_difference, p_infty = p_infty)) 
}


####================  Function 9: Calculate the Memory curve:   ===============####

calculate_memory_curve <- function(tm_perc, t, steady_state_tm, size, marginal_dist_f0){
  # Step 1: Raise current transition matrix to power t:
  tm_perc_t <- matrix_power(as.matrix(tm_perc[,1:size]),t)
  # Step 2: Calculate the total variation between current distribution of class i and the steady state:
  IM_vec <- 0.5 * rowSums(abs(tm_perc_t[,1:size] - steady_state_tm))
  # Step 3: Calculate the AIM:
  AIM <- t(as.matrix(IM_vec)) %*% as.matrix(marginal_dist_f0)
  return(list(IM_vec = IM_vec, AIM = AIM ))
}

####=============  Function 10: Calculate the expected outcome :   ==========####

calculate_expected_outcome <- function(tm_perc, t, size, marginal_dist_f0){
  # Step 1: Calculate the Expected outcomes:
  outcomes <- seq(1, size)
  expected_outcomes <- as.matrix(tm_perc[,1:size]) %*% as.matrix(outcomes)
  # Step 2: Calculate the Average Expected Outcome:
  marginal_dist_ft <- t(marginal_dist_f0) %*% matrix_power(as.matrix(tm_perc[,1:size]),t)
  average_expected_outcomes <- as.matrix(marginal_dist_ft) %*% as.matrix(expected_outcomes)
  return(list(average_expected_outcomes = average_expected_outcomes,
              expected_outcomes = expected_outcomes))
}

####=============  Function 11: Calculate the MTE :   =======================####

MTE = \(p_ii,n) Reduce(`+`, seq(1,n) * p_ii ^(seq(1,n) - 1)) * (1 - p_ii)

calculate_MTE <- function(tm_perc, step, size){
  # Step 1: Define a function to calculate MTE for a given n:
  diag <- as.list(diag(as.matrix(tm_perc[,1:size])))
  # Step 2: Store the MTE results up to n:
  MTE_lst <- lapply(diag,function(x){MTE(x, step)})
  return(MTE_lst)
}

####=============  Function 11: Calculate the R-R \rho :   ===================####

calculate_rank_correlation <- function(tm_df){
  tm_df <- 
    tm_df %>%
    dplyr::select(
      status_c,
      status_p
    ) %>%
    mutate(
      status_c = rank(status_c, ties.method = "average"),
      status_p = rank(status_p, ties.method = "average")
    )
  rho <- 
    cov(tm_df$status_c, tm_df$status_p)/(sd(tm_df$status_c)*sd(tm_df$status_p))
  return(rho)
}


####========  Function 12: Calculate the 2nd Largest Eigenvalues :   ==========####

second_eigenvalues <- function(tm_perc,size){
  data <- 
    as.matrix(
      tm_perc[,1:size])
  data[is.na(data)] <- 0
  eigen_result <- eigen(data)$values
#  eigen_2nd <- sort(eigen_result, decreasing = TRUE)[2]
  eigen_2nd <- eigen_result[order(Mod(eigen_result), decreasing = TRUE)][2]
  eigen_2nd_Mod <- Mod(eigen_2nd)
  return(eigen_2nd_Mod)
}

second_eigenvec <- function(tm_perc,size){
  data <- 
    as.matrix(
      tm_perc[,1:size])
  data[is.na(data)] <- 0
  eigen_val <- eigen(data)$values
  eigen_vec <- eigen(data)$vectors
  eigen_vec_2nd <- eigen_vec[,order(Mod(eigen_val), decreasing = TRUE)[2]]
  return(eigen_vec_2nd)
}

#### ================= Final Simulation Methods ==================== ####

single_simulation <- function(tm_df_org, size, tol){
    tm_df <- tm_df_org
    if(size == 5){
      tm_results <- transition_matrix_5class_fun(tm_df)
      tm_perc <- tm_results$trans_matrix_perc
      tm_num <- tm_results$trans_matrix_num
    }else if(size == 6){
      tm_results <- transition_matrix_6class_fun(tm_df)
      tm_perc <- tm_results$trans_matrix_perc
      tm_num <- tm_results$trans_matrix_num
    }else{
      tm_results <- transition_matrix_7class_fun(tm_df)
      tm_perc <- tm_results$trans_matrix_perc
      tm_num <- tm_results$trans_matrix_num
    }
    # --- early check ---
    if (is.null(tm_perc) || (anyNA(tm_perc))) {
      return(NA)
    }
    # Step 3: Calculate marginal_dist_f0:  
    marginal_dist_f0 <- marginal_dist(tm_num, size)$marginal_dist_f
    # Step 4: Calculate the Steady State Results:
    steady_state_rst <- calculate_the_steady_state_tm(tm_perc, size, tol)
    t_max <- steady_state_rst$t
    steady_state_tm <- steady_state_rst$p_infty
    if (any(steady_state_tm == 0)) {
      return(NA)
    }
    #===========================Our Measures==============================#
    #----------------------t-vary measure(curve measure)-------------------#
    # For each t up to the t_max:
    historical_mobility_vec <- c()
    historical_upward_mobility_vec <- c()
    historical_downward_mobility_vec <- c()
    structural_mobility_vec <- c()
    exchange_mobility_vec <- c()
    IM_curve_lst <- list()
    AIM_vec <- c()
    AEO_vec <- c()
    d_prime_vec <- c()
    d_prime_combination_vec <- c()
    for(t in seq(t_max)){ # I will still iterate until t_max
      # Step 5: Calculate the movement mobility for current t;
      movement_mobility <- 
        compute_movement_mobility(tm_perc,marginal_dist_f0, t, size)
      historical_mobility_vec <- 
        c(historical_mobility_vec,movement_mobility$historical_mobility)
      historical_upward_mobility_vec <- 
        c(historical_upward_mobility_vec,movement_mobility$historical_upward_mobility)
      historical_downward_mobility_vec <- 
        c(historical_downward_mobility_vec,movement_mobility$historical_downward_mobility)
      structural_mobility_vec <- 
        c(structural_mobility_vec,movement_mobility$structural_mobility)
      exchange_mobility_vec <- 
        c(exchange_mobility_vec,movement_mobility$exchange_mobility)
      # Step 6: Calculate the memory curve for current t:
      memory_curve <-
        calculate_memory_curve(tm_perc, t, steady_state_tm, size, marginal_dist_f0)
      AIM_vec <- c(AIM_vec, as.numeric(memory_curve$AIM))
      IM_curve_lst[[t]] <- memory_curve$IM_vec
      # Step 8: Calculate the expected outcome for current t:
      expected_outcome <-
        calculate_expected_outcome(tm_perc, t, size, marginal_dist_f0)
      AEO_vec <-
        c(AEO_vec, as.numeric(expected_outcome$average_expected_outcomes))
      # Step 9: Record D prime:
      d_prime_rst <- calculate_d_prime(NULL, tm_perc, size, TRUE, t)
      d_prime <- d_prime_rst$d_prime
      d_prime_combination <- d_prime_rst$d_prime_combination
      d_prime_vec <- c(d_prime_vec, as.numeric(d_prime))
      d_prime_combination_vec <- c(d_prime_combination_vec, d_prime_combination)
    }
    # Organize the results:
    ## Movement Measure:
    movement_measure <-  
      tibble(
        `Historical Mobility` = historical_mobility_vec,
        `Historical Upward Mobility` = historical_upward_mobility_vec,
        `Historical Downward Mobility` = historical_downward_mobility_vec,
        `Structural Mobility` = structural_mobility_vec,
        `Exchange Mobility` = exchange_mobility_vec
      )
    ## D prime:
    d_prime_df <-
      tibble(
        `d prime` = d_prime_vec,
         d_prime_combination <- d_prime_combination_vec
      )
    ## Memory Curve:
    memory_curve <-
      cbind(
        do.call(rbind, IM_curve_lst),
        data_frame(
          AIM = AIM_vec)
      ) %>%
      mutate(
        t = seq(t_max)
      ) %>%
      mutate(
        flag = ifelse(t <= t_max, 0, 1)
      ) %>%
      mutate(
        across(
          -c(t,flag),
          ~ ifelse(flag == 1, NA, .x) 
        )
      ) %>%
      dplyr::select(
        -flag
      )
    
    ## Average Outcome:
    expected_outcome_tss <- 
      calculate_expected_outcome(steady_state_tm, 1, size, marginal_dist_f0)
    average_expected_outcomes <-
      list(
        AEO = AEO_vec,
        EO = t(expected_outcome$expected_outcomes),
        EO_tss = t(expected_outcome_tss$expected_outcomes)
      )
    #----------------------t-invariant measure(gloabl measure)-----------------#
    # Step 9: Calculate MTE, given current tm_perc
    # It will converge to 1/(1 - p_ii) till the end
    MTE_vec <- 1/(1 - diag(as.matrix(tm_perc[,1:size])))
    AMTE <- mean(MTE_vec)
    # Step 10: Calculate second_eigenvalues
    lambda2 <- second_eigenvalues(tm_perc, size)
    vec_2nd <- second_eigenvec(tm_perc, size)
    # Step 11: Calculate Rank-Rank Correlation
    rho <- calculate_rank_correlation(tm_df)
    # Step 12: Calculate MFP:
    
    MFP <- tryCatch({
      mc <- new("markovchain", transitionMatrix = as.matrix(tm_perc[, 1:size]))
      markovchain::meanFirstPassageTime(mc)
    }, error = function(e) NA)
    
    # Step 13: Calculate Altham index:
    P <- as.matrix(tm_num[,1:size])
    altham_index <- iac(P, weighting = "none") * 2 * sqrt(nrow(P) * ncol(P))
    
    #===========================Prof Coleman's Measures========================#
    # Step 1: Convert the tm_num to tm_prob
    tm_prob <- as.matrix(convert_obs_to_prob(tm_num))
    # Step 2: Calculate the Geenen's D:
    GeenensD <- round(calcGeenensStructZero(tm_prob)$GeenensD,3)
    # Step 3: Calculate Hellinger's Index:
    HellingerDep <- round(calcDiscreteHellinger(tm_prob)$HellingGeen,3)
    
    #==========================Store Results===================================#
    
    final_results <- list()
    final_results[["Movement Measure"]] <- movement_measure
    final_results[["Memory Curve"]] <- memory_curve
    final_results[["Expected Outcomes"]] <- average_expected_outcomes
    final_results[["MTE_vec"]] <- MTE_vec
    final_results[["AMTE"]] <- AMTE
    final_results[["lambda2"]] <- lambda2
    final_results[["marginal_dist_f0"]] <- marginal_dist_f0
    final_results[["d_prime"]] <- d_prime_df
    final_results[["rho"]] <- rho
    final_results[["MFP"]] <- MFP
    final_results[["GeenensD"]] <- GeenensD
    final_results[["HellingerDep"]] <- HellingerDep
    final_results[["t_max"]] <- t_max
    final_results[["tm_num"]] <- tm_num
    final_results[["vec_2nd"]] <- vec_2nd
    final_results[["altham_index"]] <- altham_index
    final_results[["Steady State TM"]] <- steady_state_tm
    return(final_results)
}



#### ======================== Average Steady States =========================####

average_steady_state_time <- function(tm_df_org, times, size, tol){
  single_steady_state <- function(tm_df_org, size, tol){
    # Step 1: Sample the Data:
    tm_df <- sample_data(tm_df_org)
    # Step 2: Calculate the tm_perc:
    if(size == 5){
      tm_perc <- transition_matrix_5class_fun(tm_df)$trans_matrix_perc
    }
    steady_state_rst <- calculate_the_steady_state_tm(tm_perc,size,tol)$t
  }
  t_vec <- 
    replicate(times, single_steady_state(tm_df_org, size, tol), simplify = FALSE)
  return(unlist(t_vec))
}


####### ============ Generate the Fitted TM DF lst ================= #######
generate_tm_df_fit <- function(fitted_tm_df){
  fitted_tm_df_lst_raw <- split(fitted_tm_df, fitted_tm_df$cohort)
  fitted_tm_df_lst <- 
    lapply(
      fitted_tm_df_lst_raw,
      function(fitted_tm_df_sub){
        data <-
          fitted_tm_df_sub %>%
          dplyr::select(
            round_predictions,
            parent_status,
            children_status,
            cohort
          ) %>%
          uncount(round_predictions) %>%
          rename(
            status_c = children_status,
            status_p = parent_status
          )
      }
    )
  return(fitted_tm_df_lst)
}


####### ============ Chunk Processing ================= #######

chunk_list <- function(data_list, n_cores) {
  # Total length of the list
  N <- length(data_list)
  # Calculate chunk size
  chunk_size <- n_cores
  # Split the list into chunks
  chunks <- split(data_list, ceiling(seq_along(data_list) / chunk_size))
  return(chunks)
}

####### ==================== Calculate MFP ============================ #######

#### Step 1: Sample Decedents:
sample_descedents <- function(initial_state, tm_perc, col_index, size){
  mat_t <- matrix_power(as.matrix(tm_perc), col_index)
  # Generate one decedent based on the column index
  if (is.null(initial_state) || is.na(initial_state)) {
    stop(paste(
      "Error: `initial_state` is NULL or NA.\n",
      "tm_perc:\n", paste(capture.output(print(tm_perc)), collapse = "\n"),
      "mat_t:\n", paste(capture.output(print(mat_t)), collapse = "\n"),
      "col_index:", col_index, "\n",
      "size:", size
    ))
  }
  
  # Assert the validity of the matrix indexing
  assert_that(
    length(mat_t[initial_state, ]) == size,
    msg = paste(
      "tm_perc:\n", paste(capture.output(print(tm_perc)), collapse = "\n"),
      "Length of mat_t[initial_state, ] is", 
      length(mat_t[initial_state, ]),
      "but size is", size, 
      "\nmat_t is:\n", paste(capture.output(print(mat_t)), collapse = "\n"),
      "\nThe initial state is", initial_state
    )
  )
  prob_vec <- mat_t[initial_state, ]
  if (all(prob_vec <= 0)) {
    message("Warning: All probabilities are zero for the given initial state. Returning NA.")
    return(NA)  
  }else{
    cell <- sample(1:size, 1, prob = prob_vec)
    return(cell)
  }
}

#### Step 2: Calculate the MFP for a single individual:
calculate_ind_MFP <- function(row_list, tm_perc, size, t){
  # Check `row_list`
  if (is.null(row_list$status_p) || is.null(row_list$status_c)) {
    stop("Error: `row_list` is missing required fields.")
  }
  
  if (is.na(row_list$status_p) || is.na(row_list$status_c)) {
    stop(paste(
      "Error: `row_list` contains NA values.\n",
      "row_list:\n", paste(capture.output(print(row_list)), collapse = "\n")
    ))
  }
  # Step 1: Generate the vector of the statuses of the decedents:
  des_vec <- 
    as.matrix(c(
      row_list$status_c,
      vapply(1:t, function(col_index) {
        sample_descedents(row_list$status_p, tm_perc, col_index, size)
      }, numeric(1))))
  
  # Step 2: track the first passage time:
  target_states <- 
    as.list(setdiff(seq(1, size), row_list$status_p))
  
  first_index_list <-
    setNames(
      lapply(
        target_states,
        function(target_state){
          first_index <-
            which(des_vec == target_state)[1]
          return(first_index)
        }
      ),
      target_states
    )
  return(first_index_list)
}

#### Step 3: Calculate the stats of the MFP:
calculate_MFP_stats <- function(tm_df_group, tm_perc, size, t){
  # Break it into row_lists:
  if (nrow(tm_df_group) == 0) {
    stop("Error: `tm_df_group` is empty.")
  }
  
  row_lists <- split(tm_df_group, seq(nrow(tm_df_group)))
  # Simulate_decedents:
  # Check each `row_list`
  lapply(row_lists, function(row_list) {
    if (is.null(row_list$status_p) || is.null(row_list$status_c)) {
      stop("Error: A `row_list` is missing required fields.")
    }
    if (any(is.na(row_list$status_p)) || any(is.na(row_list$status_c))) {
      stop("Error: A `row_list` contains NA values.")
    }
  })
  ind_MFP_lists <-
    lapply(
      row_lists,
      function(row_list){
        calculate_ind_MFP(row_list, tm_perc, size, t) 
      }
    )
  # Calculate the stats:
  statistics <- map(transpose(ind_MFP_lists), function(values) {
    values <- unlist(values)  # Flatten the values into a vector
    values <- values[!is.na(values)]  # Remove NA values
    list(
      mean = mean(values),
      median = median(values),
      percentile_90 = quantile(values, 0.9),
      percentile_95 = quantile(values, 0.95)
    )
  })
  return(statistics)
}

#### Step 4: Calculate the final step:
calculate_MFP_final <- function(tm_df_org, tm_perc, size, t){
  composition <- 
    tm_df_org %>%
    group_by(
      status_p, 
      status_c
    ) %>%
    summarise(
      count = n(), 
      .groups = "drop") %>%
    mutate(
      proportion = count / sum(count))
  
  sampled_data <- 
    composition %>%
    slice_sample(
      n = 200,
      weight_by = proportion,
      replace = TRUE
    ) %>%
    select(
      status_p, 
      status_c) %>%
    mutate(
      status_p = as.numeric(status_p),
      status_c = as.numeric(status_c))  
  
  tm_df_groups <-
    split(
      sampled_data %>%
        dplyr::select(
          status_c,
          status_p
        ),
      sampled_data %>%
        dplyr::select(
          status_p
        )
    )
  
  assert_that(
    all(sapply(tm_df_groups, nrow) > 0),
    msg = paste(
      "Error: One or more groups in `tm_df_groups` are empty.\n",
      "tm_df_groups structure:\n",
      paste(capture.output(print(tm_df_groups)), collapse = "\n")
    )
  )
  
  assert_that(
    all(sapply(tm_df_groups, function(group) all(complete.cases(group)))),
    msg = "Error: One or more groups in tm_df_groups contain NA values."
  )
  #print(unique(sampled_data$status_p))
  tm_perc <-
    tm_perc %>%
    dplyr::select(
      -rowsum
    )
  tm_perc[is.na(tm_perc)] <- 0
  assert_that(
    !any(is.na(tm_perc)),
    msg = paste(
      "Error: `tm_perc` contains NA values.\n",
      "tm_perc with NA values:\n", paste(capture.output(print(tm_perc)), collapse = "\n")
    )
  )
  
  MFP_stats_list <- 
    lapply(
      tm_df_groups,
      function(tm_df_group){
        calculate_MFP_stats(tm_df_group,tm_perc,size,t)
      }
    )
  return(MFP_stats_list)
}

#============================ 5 Class Transition Matrices =====================#

transition_matrix_5class_fun <- function(data){
  pc_status <- 
    data %>% 
    mutate(
      # Parent is in class 1:
      pc_11 = ifelse(status_p == 1 & status_c == 1, 1, 0),
      pc_12 = ifelse(status_p == 1 & status_c == 2, 1, 0),
      pc_13 = ifelse(status_p == 1 & status_c == 3, 1, 0),
      pc_14 = ifelse(status_p == 1 & status_c == 4, 1, 0),
      pc_15 = ifelse(status_p == 1 & status_c == 5, 1, 0),
      
      # Parent is in class 2:
      pc_21 = ifelse(status_p == 2 & status_c == 1, 1, 0),
      pc_22 = ifelse(status_p == 2 & status_c == 2, 1, 0),
      pc_23 = ifelse(status_p == 2 & status_c == 3, 1, 0),
      pc_24 = ifelse(status_p == 2 & status_c == 4, 1, 0),
      pc_25 = ifelse(status_p == 2 & status_c == 5, 1, 0),
      
      # Parent is in class 3:
      pc_31 = ifelse(status_p == 3 & status_c == 1, 1, 0),
      pc_32 = ifelse(status_p == 3 & status_c == 2, 1, 0),
      pc_33 = ifelse(status_p == 3 & status_c == 3, 1, 0),
      pc_34 = ifelse(status_p == 3 & status_c == 4, 1, 0),
      pc_35 = ifelse(status_p == 3 & status_c == 5, 1, 0),
      
      # Parent is in class 4:
      pc_41 = ifelse(status_p == 4 & status_c == 1, 1, 0),
      pc_42 = ifelse(status_p == 4 & status_c == 2, 1, 0),
      pc_43 = ifelse(status_p == 4 & status_c == 3, 1, 0),
      pc_44 = ifelse(status_p == 4 & status_c == 4, 1, 0),
      pc_45 = ifelse(status_p == 4 & status_c == 5, 1, 0),
      
      # Parent is in class 5:
      pc_51 = ifelse(status_p == 5 & status_c == 1, 1, 0),
      pc_52 = ifelse(status_p == 5 & status_c == 2, 1, 0),
      pc_53 = ifelse(status_p == 5 & status_c == 3, 1, 0),
      pc_54 = ifelse(status_p == 5 & status_c == 4, 1, 0),
      pc_55 = ifelse(status_p == 5 & status_c == 5, 1, 0)
    )
  
  
  trans_matrix_01_num <- 
    tibble(
      Class1 = c(
        sum(pc_status$pc_11, na.rm = TRUE),
        sum(pc_status$pc_21, na.rm = TRUE),
        sum(pc_status$pc_31, na.rm = TRUE),
        sum(pc_status$pc_41, na.rm = TRUE),
        sum(pc_status$pc_51, na.rm = TRUE)
      ),
      Class2 = c(
        sum(pc_status$pc_12, na.rm = TRUE),
        sum(pc_status$pc_22, na.rm = TRUE),
        sum(pc_status$pc_32, na.rm = TRUE),
        sum(pc_status$pc_42, na.rm = TRUE),
        sum(pc_status$pc_52, na.rm = TRUE)
      ),
      Class3 = c(
        sum(pc_status$pc_13, na.rm = TRUE),
        sum(pc_status$pc_23, na.rm = TRUE),
        sum(pc_status$pc_33, na.rm = TRUE),
        sum(pc_status$pc_43, na.rm = TRUE),
        sum(pc_status$pc_53, na.rm = TRUE)
      ),
      Class4 = c(
        sum(pc_status$pc_14, na.rm = TRUE),
        sum(pc_status$pc_24, na.rm = TRUE),
        sum(pc_status$pc_34, na.rm = TRUE),
        sum(pc_status$pc_44, na.rm = TRUE),
        sum(pc_status$pc_54, na.rm = TRUE)
      ),
      Class5 = c(
        sum(pc_status$pc_15, na.rm = TRUE),
        sum(pc_status$pc_25, na.rm = TRUE),
        sum(pc_status$pc_35, na.rm = TRUE),
        sum(pc_status$pc_45, na.rm = TRUE),
        sum(pc_status$pc_55, na.rm = TRUE)
      )
    ) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5,
           validation = 
             c(data %>% filter(status_p == 1) %>% nrow(),
               data %>% filter(status_p == 2) %>% nrow(), 
               data %>% filter(status_p == 3) %>% nrow(), 
               data %>% filter(status_p == 4) %>% nrow(),
               data %>% filter(status_p == 5) %>% nrow())
    )
  
  rownames(trans_matrix_01_num) <- c("Class1",
                                     "Class2",
                                     "Class3",
                                     "Class4",
                                     "Class5")
  
  
  trans_matrix_01_perc <- 
    as.data.frame(trans_calc(
      trans_matrix_01_num %>%
        dplyr::select(
          Class1,
          Class2,
          Class3,
          Class4,
          Class5))) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5) 
  
  rownames(trans_matrix_01_perc) <- 
    c("Class1",
      "Class2",
      "Class3",
      "Class4",
      "Class5")
  return(list(trans_matrix_num = trans_matrix_01_num,
              trans_matrix_perc = trans_matrix_01_perc))
}


#============================ 5 Class Transition Matrices: With weight =====================#

transition_matrix_5class_fun_weight <- function(data){
  pc_status <- 
    data %>% 
    mutate(
      # Parent is in class 1:
      pc_11 = ifelse(status_p == 1 & status_c == 1, wn, 0),
      pc_12 = ifelse(status_p == 1 & status_c == 2, wn, 0),
      pc_13 = ifelse(status_p == 1 & status_c == 3, wn, 0),
      pc_14 = ifelse(status_p == 1 & status_c == 4, wn, 0),
      pc_15 = ifelse(status_p == 1 & status_c == 5, wn, 0),
      
      # Parent is in class 2:
      pc_21 = ifelse(status_p == 2 & status_c == 1, wn, 0),
      pc_22 = ifelse(status_p == 2 & status_c == 2, wn, 0),
      pc_23 = ifelse(status_p == 2 & status_c == 3, wn, 0),
      pc_24 = ifelse(status_p == 2 & status_c == 4, wn, 0),
      pc_25 = ifelse(status_p == 2 & status_c == 5, wn, 0),
      
      # Parent is in class 3:
      pc_31 = ifelse(status_p == 3 & status_c == 1, wn, 0),
      pc_32 = ifelse(status_p == 3 & status_c == 2, wn, 0),
      pc_33 = ifelse(status_p == 3 & status_c == 3, wn, 0),
      pc_34 = ifelse(status_p == 3 & status_c == 4, wn, 0),
      pc_35 = ifelse(status_p == 3 & status_c == 5, wn, 0),
      
      # Parent is in class 4:
      pc_41 = ifelse(status_p == 4 & status_c == 1, wn, 0),
      pc_42 = ifelse(status_p == 4 & status_c == 2, wn, 0),
      pc_43 = ifelse(status_p == 4 & status_c == 3, wn, 0),
      pc_44 = ifelse(status_p == 4 & status_c == 4, wn, 0),
      pc_45 = ifelse(status_p == 4 & status_c == 5, wn, 0),
      
      # Parent is in class 5:
      pc_51 = ifelse(status_p == 5 & status_c == 1, wn, 0),
      pc_52 = ifelse(status_p == 5 & status_c == 2, wn, 0),
      pc_53 = ifelse(status_p == 5 & status_c == 3, wn, 0),
      pc_54 = ifelse(status_p == 5 & status_c == 4, wn, 0),
      pc_55 = ifelse(status_p == 5 & status_c == 5, wn, 0)
    )
  
  
  trans_matrix_01_num <- 
    tibble(
      Class1 = c(
        sum(pc_status$pc_11, na.rm = TRUE),
        sum(pc_status$pc_21, na.rm = TRUE),
        sum(pc_status$pc_31, na.rm = TRUE),
        sum(pc_status$pc_41, na.rm = TRUE),
        sum(pc_status$pc_51, na.rm = TRUE)
      ),
      Class2 = c(
        sum(pc_status$pc_12, na.rm = TRUE),
        sum(pc_status$pc_22, na.rm = TRUE),
        sum(pc_status$pc_32, na.rm = TRUE),
        sum(pc_status$pc_42, na.rm = TRUE),
        sum(pc_status$pc_52, na.rm = TRUE)
      ),
      Class3 = c(
        sum(pc_status$pc_13, na.rm = TRUE),
        sum(pc_status$pc_23, na.rm = TRUE),
        sum(pc_status$pc_33, na.rm = TRUE),
        sum(pc_status$pc_43, na.rm = TRUE),
        sum(pc_status$pc_53, na.rm = TRUE)
      ),
      Class4 = c(
        sum(pc_status$pc_14, na.rm = TRUE),
        sum(pc_status$pc_24, na.rm = TRUE),
        sum(pc_status$pc_34, na.rm = TRUE),
        sum(pc_status$pc_44, na.rm = TRUE),
        sum(pc_status$pc_54, na.rm = TRUE)
      ),
      Class5 = c(
        sum(pc_status$pc_15, na.rm = TRUE),
        sum(pc_status$pc_25, na.rm = TRUE),
        sum(pc_status$pc_35, na.rm = TRUE),
        sum(pc_status$pc_45, na.rm = TRUE),
        sum(pc_status$pc_55, na.rm = TRUE)
      )
    ) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5,
           validation = 
             c(sum(data %>% filter(status_p == 1) %>% pull(wn)),
               sum(data %>% filter(status_p == 2) %>% pull(wn)), 
               sum(data %>% filter(status_p == 3) %>% pull(wn)), 
               sum(data %>% filter(status_p == 4) %>% pull(wn)),
               sum(data %>% filter(status_p == 5) %>% pull(wn))
    ))
  
  rownames(trans_matrix_01_num) <- c("Class1",
                                     "Class2",
                                     "Class3",
                                     "Class4",
                                     "Class5")
  
  
  trans_matrix_01_perc <- 
    as.data.frame(trans_calc(
      trans_matrix_01_num %>%
        dplyr::select(
          Class1,
          Class2,
          Class3,
          Class4,
          Class5))) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5) 
  
  rownames(trans_matrix_01_perc) <- 
    c("Class1",
      "Class2",
      "Class3",
      "Class4",
      "Class5")
  return(list(trans_matrix_num = trans_matrix_01_num,
              trans_matrix_perc = trans_matrix_01_perc))
}


#=================== Generate the Bootstrap Fitted Data ==================#
generate_bootstrap_fitted_data <- function(
    gss_data,size,age_l,age_u, year_l, year_u,gamma_par){
  # Step 1: Cut the data by survey year, and shuffled within survey year:
  gss_sy_lst <- split(gss_data, gss_data$year)
  gss_shuffled_sy_lst <-
    lapply(
      gss_sy_lst,
      function(gss_sy){
        sampled_gss_sy <- 
          gss_sy[sample(nrow(gss_sy), 
                        size = nrow(gss_sy), 
                        replace = TRUE), ]
        return(sampled_gss_sy)
      }
    )
  gss_shuffled_sy <- bind_rows(gss_shuffled_sy_lst)
  gss_shuffled_sy <-
    gss_shuffled_sy %>%
    filter(
      cohort >= year_l & cohort <= (year_u + 1)
    )
  # Step 2: Generate the tm_num_lst for the regression:
  cs_tm_df_5_2555_1 <- cross_tm_df(gss_shuffled_sy, age_l, age_u, 1)
  cs_tm_5_2555_1_num_lst <- list()
  if(size == 5){
    for(i in seq(length(cs_tm_df_5_2555_1))){
      data <- cs_tm_df_5_2555_1[[i]]
      tms <- transition_matrix_5class_fun(data)
      cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
    }
  }else if(size == 6){
    for(i in seq(length(cs_tm_df_5_2555_1))){
      data <- cs_tm_df_5_2555_1[[i]]
      tms <- transition_matrix_6class_fun(data)
      cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
    }
  }else{
    for(i in seq(length(cs_tm_df_5_2555_1))){
      data <- cs_tm_df_5_2555_1[[i]]
      tms <- transition_matrix_7class_fun(data)
      cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
    }
  }
  # Step 3: Generate the regression data:
  regression_data_5_2555 <- 
    generate_model_data(FALSE, NULL, cs_tm_5_2555_1_num_lst, size, year_l, year_u)
  # Step 4: Fit the model:
  control_params <- gam.control(trace = FALSE)
  model_baseline <-
    gam(
      cell_count ~
        parent_children +
        s(cohort,
          by = interaction(parent_status, children_status),
          bs = "tp"),
      family = poisson(link="log"),
      gamma = gamma_par,
      select=FALSE, # Keep this as FALSE, avoiding oversmoothing
      method = "REML",
      data = regression_data_5_2555,
      control = control_params,
      gc.level = 2,
      nthreads = 3)
  # Step 5: Generate the model fitted data:
  fitted_data_5_2555 <- 
    generate_model_data(TRUE, model_baseline, cs_tm_5_2555_1_num_lst, size, year_l, year_u)
  return(fitted_data_5_2555)
}

#======================= Define the new bootstrap function =====================#
single_simulation_lst <- function(
    gss_data,size,age_l,age_u, year_l, year_u, est_l, est_u, tol, gamma_par){
  # Generate the fitted data:
  fitted_data_5_2555 <- generate_bootstrap_fitted_data(
    gss_data, size, age_l, age_u, year_l, year_u, gamma_par)
  # Convert the fitted_data to a list:
  cs_tm_df_5_2555_1_fitted_lst <- generate_tm_df_fit(fitted_data_5_2555)
  # Only focus on the target range:
  cs_tm_df_5_2555_1_fitted_lst <-
  cs_tm_df_5_2555_1_fitted_lst[names(cs_tm_df_5_2555_1_fitted_lst) >= est_l & 
                               names(cs_tm_df_5_2555_1_fitted_lst) <= est_u]
  # Generate the simulation result for each of the cohort:
  simulation_rst_lst <- 
    lapply(
      cs_tm_df_5_2555_1_fitted_lst,
      function(tm_df_org){
        rst <-
          single_simulation(tm_df_org, size, tol)
        print("1")
        return(rst)
      }
    )
  return(simulation_rst_lst)
}

# ========================= New Methods =========================== #

##=================##
#    Raw Data
##=================##

generate_raw_data <- function(gss_data, age_l, 
                              age_u, year_l, year_u, boot = TRUE){
  # Step 1: If bootstrap, will shuffle the data within each survey year:
  if(boot == TRUE){
    gss_sy_lst <- split(gss_data, gss_data$year)
    gss_shuffled_sy_lst <-
      lapply(
        gss_sy_lst,
        function(gss_sy){
          sampled_gss_sy <- 
            gss_sy[sample(nrow(gss_sy), 
                          size = nrow(gss_sy), 
                          replace = TRUE), ]
          return(sampled_gss_sy)
        }
      )
    # Generate the data frame:
    gss_df <-
      bind_rows(gss_shuffled_sy_lst) %>%
      filter(
        cohort >= year_l & cohort <= (year_u + 1)
      )
  }else{
    # Return a dataframe within desired year range without shuffling:
    gss_df <-
      gss_data %>%
      filter(
        cohort >= year_l & cohort <= (year_u + 1)
      )
  }
  
  # Step 2: Generate the tm_num_lst without smoothing:
  cs_tm_df_5_2555_1 <- cross_tm_df(gss_df, age_l, age_u, 1)
  cs_tm_5_2555_1_num_lst <- list()
  for(i in seq(length(cs_tm_df_5_2555_1))){
    data <- cs_tm_df_5_2555_1[[i]]
    tms <- transition_matrix_5class_fun(data)
    cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
  }
  regression_data_5_2555 <- 
    generate_model_data(FALSE, NULL, cs_tm_5_2555_1_num_lst, 5, year_l, year_u) %>%
    rename(
      round_predictions = cell_count
    )
  cs_tm_df_5_2555_1_fitted_lst <- generate_tm_df_fit(regression_data_5_2555)
  return(cs_tm_df_5_2555_1_fitted_lst)
}


##=================##
#    Survey Wave.  ##   
##=================##


generate_fitted_surveywave <- function(gss_data, age_l,
                                       age_u,
                                       boot = TRUE){
  # Step 1: process the data:
  if(boot == TRUE){
    gss_sy_lst <- split(gss_data, gss_data$year)
    gss_final_sy_lst <-
      lapply(
        seq(length(gss_sy_lst)),
        function(i){
          sampled_gss_sy <- 
            gss_sy_lst[[i]][
              sample(
                nrow(gss_sy_lst[[i]]), 
                size = nrow(gss_sy_lst[[i]]), 
                replace = TRUE), ] %>%
            filter(
              age >= age_l & age <= age_u
            ) %>%
            dplyr::select(
              -cohort
            ) %>%
            mutate(
              cohort = i 
            )
          return(sampled_gss_sy)
        }
      )
  }else{
    gss_sy_lst <- split(gss_data, gss_data$year)
    gss_final_sy_lst <-
      lapply(
        seq(length(gss_sy_lst)),
        function(i){
          gss_sy <- 
            gss_sy_lst[[i]] %>%
            filter(
              age >= age_l & age <= age_u
            ) %>%
            dplyr::select(
              -cohort
            ) %>%
            mutate(
              cohort = i 
            )
          return(gss_sy)
        }
      )
  }
  # Step 2: Generate the tm_num_lst for the regression:
  cs_tm_df_5_2555_1 <- gss_final_sy_lst 
  cs_tm_5_2555_1_num_lst <- list()
  for(i in seq(length(cs_tm_df_5_2555_1))){
    data <- cs_tm_df_5_2555_1[[i]]
    tms <- transition_matrix_5class_fun(data)
    cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
  }
  # Step 3: Generate the regression data:
  regression_data_5_2555 <- 
    generate_model_data(
      FALSE, NULL, cs_tm_5_2555_1_num_lst, 5, 1, length(unique(gss_data$year)))
  # Step 4: Fit the model:
  control_params <- gam.control(trace = FALSE)
  model_baseline <-
    gam(
      cell_count ~
        parent_children +
        s(cohort,
          by = interaction(parent_status, children_status),
          bs = "tp"),
      family = poisson(link="log"),
      gamma = 1,
      select=FALSE, 
      method = "REML",
      data = regression_data_5_2555,
      control = control_params,
      gc.level = 2,
      nthreads = 3)
  # Step 5: Generate the model fitted data:
  fitted_data_5_2555 <- 
    generate_model_data(
      TRUE, model_baseline, cs_tm_5_2555_1_num_lst, 5, 1, length(unique(gss_data$year)))
  cs_tm_df_5_2555_1_fitted_lst <- generate_tm_df_fit(fitted_data_5_2555)
  return(cs_tm_df_5_2555_1_fitted_lst)
}


##=================##
#    Bin.          ##   
##=================##



generate_bin <- function(gss_data, age_l,age_u, year_l, year_u,
                         boot = TRUE, bin_width){
  if(boot == TRUE){
    gss_sy_lst <- split(gss_data, gss_data$year)
    gss_shuffled_sy_lst <-
      lapply(
        gss_sy_lst,
        function(gss_sy){
          sampled_gss_sy <- 
            gss_sy[sample(nrow(gss_sy), 
                          size = nrow(gss_sy), 
                          replace = TRUE), ]
          return(sampled_gss_sy)
        }
      )
    gss_final_data <-
      bind_rows(gss_shuffled_sy_lst) %>%
      filter(
        cohort >= year_l & cohort <= year_u
      )
  }else{
    gss_final_data <-
      gss_data %>%
      filter(
        cohort >= year_l & cohort <= year_u
      )
  }
  
  # Generate transition matrices with assigned bin_width:
  if(bin_width == 5){
    cs_tm_df_5_2555 <- cross_tm_df(gss_final_data, age_l, age_u, 5)
  }else{
    cs_tm_df_5_2555_bin10 <-
      gss_final_data %>%
      filter(
        age >= age_l & age <= age_u
      ) %>%
      mutate(
        bins = case_when(
          cohort >= 1945 & cohort <= 1954 ~ "[1945,1954]",
          cohort >= 1955 & cohort <= 1964 ~ "[1955,1964]",
          cohort >= 1965 & cohort <= 1974 ~ "[1965,1974]",
          cohort >= 1975 & cohort <= 1990 ~ "[1975,1990]",
          .default = NA
        )) %>%
      filter(
        !is.na(bins)
      )
    cs_tm_df_5_2555 <- split(
      cs_tm_df_5_2555_bin10, cs_tm_df_5_2555_bin10$bins)
  }
  # Generate the transition matrices:
  cs_tm_5_2555_1_num_lst <- list()
  for(i in seq(length(cs_tm_df_5_2555))){
    data <- cs_tm_df_5_2555[[i]]
    tms <- transition_matrix_5class_fun(data)
    cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
  }
  regression_data_5_2555 <- 
    generate_model_data(FALSE, NULL, cs_tm_5_2555_1_num_lst, 5, 1, length(cs_tm_df_5_2555)) %>%
    rename(
      round_predictions = cell_count
    )
  cs_tm_df_5_2555_1_fitted_lst <- generate_tm_df_fit(regression_data_5_2555)
  return(cs_tm_df_5_2555_1_fitted_lst)
}



##=================##
#   Survey Weight. ##   
##=================##



generate_weight <- function(gss_data, age_l, age_u,
                            year_l, year_u, boot = TRUE){
  ## Step 1: Process the data:
  if(boot == FALSE){
    # Normalize the weight by survey year:
    gss_df_lst <- 
      lapply(
        split(gss_data,gss_data$year),
        function(df){
          df <-
            df %>%
            mutate(
              wn = wtssps / mean(wtssps)
            ) %>%
            mutate(
              flag = mean(wn)
            )
          assert_that(
            abs(unique(df$flag) -  1) < 1e-8
          )
          return(df)
        }
      )
    # Select the cohort with desired range:
    gss_df <-
      bind_rows(gss_df_lst) %>%
      filter(
        cohort >= year_l & cohort <= (year_u + 1)
      )
  }else{
    # Sample with replacement, normalize the weight:
    gss_shuffled_sy_lst <-
      lapply(
        split(gss_data, gss_data$year),
        function(gss_sy){
          sampled_gss_sy <- 
            gss_sy[sample(nrow(gss_sy), 
                          size = nrow(gss_sy), 
                          replace = TRUE), ] %>%
            mutate(
              wn = wtssps / mean(wtssps)
            ) %>%
            mutate(
              flag = mean(wn)
            )
          assert_that(
            abs(unique(sampled_gss_sy$flag) -  1) < 1e-8
          )
          
          return(sampled_gss_sy)
        }
      )
    # Select the desired cohort range:
    gss_df <-
      bind_rows(gss_shuffled_sy_lst) %>%
      filter(
        cohort >= year_l & cohort <= (year_u + 1)
      )
  }
  
  # Step 2: Generate the transition matrix:
  cs_tm_df_5_2555_1 <- cross_tm_df(gss_df, age_l, age_u, 1)
  cs_tm_5_2555_1_num_lst <- list()
  for(i in seq(length(cs_tm_df_5_2555_1))){
    data <- cs_tm_df_5_2555_1[[i]]
    tms <- transition_matrix_5class_fun_weight(data)
    cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
  }
  
  # Step 3: Generate the regression data:
  regression_data_5_2555 <- 
    generate_model_data(FALSE, NULL, cs_tm_5_2555_1_num_lst, 
                        size, year_l, year_u) %>%
    mutate(
      cell_count = round(cell_count)
    )
  # Step 4: Fit the model:
  control_params <- gam.control(trace = FALSE)
  model_baseline <-
    gam(
      cell_count ~
        parent_children +
        s(cohort,
          by = interaction(parent_status, children_status),
          bs = "tp"),
      family = poisson(link="log"),
      gamma = gamma_par,
      select = FALSE, # Keep this as FALSE, avoiding oversmoothing
      method = "REML",
      data = regression_data_5_2555,
      control = control_params,
      gc.level = 2,
      nthreads = 3)
  # Step 5: Generate the model fitted data:
  fitted_data_5_2555 <- 
    generate_model_data(
      TRUE, model_baseline, cs_tm_5_2555_1_num_lst, 
      size, year_l, year_u)
  cs_tm_df_5_2555_1_fitted_lst <-
    generate_tm_df_fit(fitted_data_5_2555)
  return(cs_tm_df_5_2555_1_fitted_lst)
}


single_simulation_lst_newform <- function(
    gss_data, age_l,age_u, year_l, year_u,
    boot = TRUE, bin_width, method){
  # Form 1: No Smooth:
  # Generate the fitted data:
  if(method == "raw"){
    cs_tm_df_5_2555_1_fitted_lst <- 
      generate_raw_data(
        gss_data = gss_data, 
        age_l = age_l, 
        age_u = age_u, 
        year_l = year_l,
        year_u = year_u,
        boot = boot)
  }else if(method == "sy"){
    cs_tm_df_5_2555_1_fitted_lst <- 
      generate_fitted_surveywave(
        gss_data = gss_data, 
        age_l = age_l,
        age_u = age_u, 
        boot = boot) 
  }else if(method == "weight"){
    cs_tm_df_5_2555_1_fitted_lst <- 
      generate_weight(
        gss_data = gss_data, 
        age_l = age_l,
        age_u = age_u, 
        year_l = year_l,
        year_u = year_u,
        boot = boot) 
    cs_tm_df_5_2555_1_fitted_lst <-
      cs_tm_df_5_2555_1_fitted_lst[names(cs_tm_df_5_2555_1_fitted_lst) >= "1945" & 
                                   names(cs_tm_df_5_2555_1_fitted_lst) <= "1990"]
  }else{
    cs_tm_df_5_2555_1_fitted_lst <- 
      generate_bin(
        gss_data = gss_data, 
        age_l = age_l,
        age_u = age_u, 
        year_l = year_l, 
        year_u = year_u,
        boot = boot, 
        bin_width = bin_width
      )
  }
  
  # Generate the simulation result for each of the cohort:
  simulation_rst_lst <- 
    lapply(
      cs_tm_df_5_2555_1_fitted_lst,
      function(tm_df_org){
        rst <-
          single_simulation(tm_df_org, 5, tol)
        return(rst)
      }
    )
  return(simulation_rst_lst)
}



