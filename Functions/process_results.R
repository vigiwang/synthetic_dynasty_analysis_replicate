### Clean Full Boot Stats
### @copyright weiqiw

##======================= Modify the Data Structure ==========================##
modify_cohort_lsts <- function(boot_rst, qt = 0.95) {
  all_years <- names(boot_rst[[1]])
  qt <- 0.95
  # Rearrange data: Outer list -> years, Inner list -> locations
  cohort_lsts <- map(set_names(all_years), function(year) {
    map(boot_rst, ~ .x[[year]])  # Extract each location's data for that year
  })
  cohort_lsts <- lapply(cohort_lsts, function(cohort_lst) {
    # Step 1: Extract "t_max" values from each sublist
    t_vec <- map(cohort_lst, "t_max") %>% unlist() %>% as.numeric()
    # Step 2: Compute the quantile threshold
    t_quantile <- round(quantile(t_vec, probs = qt, na.rm = TRUE))
    # Step 3: Identify sublists shorter than `t_quantile`
    index <- which(t_vec < t_quantile)
    # Step 4: Modify only those sublists
    cohort_lst[index] <- lapply(index, function(i) {
      # Modify the movement measure:
      cohort_lst[[i]]$`Movement Measure` <- bind_rows(
        cohort_lst[[i]]$`Movement Measure`,
        as.data.frame(matrix(
          NA,
          nrow = as.numeric(t_quantile - t_vec[i]),
          ncol = ncol(cohort_lst[[i]]$`Movement Measure`),
          dimnames = list(NULL, names(cohort_lst[[i]]$`Movement Measure`)))))
      
      # Modify the Memory Curve:
      cohort_lst[[i]]$`Memory Curve` <- bind_rows(
        cohort_lst[[i]]$`Memory Curve`,
        as.data.frame(matrix(
          NA,
          nrow = as.numeric(t_quantile - t_vec[i]),
          ncol = ncol(cohort_lst[[i]]$`Memory Curve`),
          dimnames = list(NULL, names(cohort_lst[[i]]$`Memory Curve`))))) %>%
        mutate(
          t = row_number()
        )
      # Modify the AEO vec:
      cohort_lst[[i]]$`Expected Outcomes`$AEO <-
        c(cohort_lst[[i]]$`Expected Outcomes`$AEO, 
          rep(tail(cohort_lst[[i]]$`Expected Outcomes`$AEO, 1),
              as.numeric(t_quantile - t_vec[i])))
      # Modify the d prime:
      cohort_lst[[i]]$d_prime <- 
        rbind(
          cohort_lst[[i]]$d_prime,
          as.data.frame(matrix(
            NA,
            nrow = as.numeric(t_quantile - t_vec[i]),
            ncol = ncol(cohort_lst[[i]]$d_prime),
            dimnames = list(NULL, names(cohort_lst[[i]]$d_prime))))
        )
      return(cohort_lst[[i]])  # Return the modified sublist
    })
    return(cohort_lst)  # Return the updated cohort list
  })
  return(cohort_lsts)  # Ensure function returns modified cohort_lsts
}

##=============================== Fetch the Results ==========================##
fetch_results <- function(bootstrap_result){
  result <- list()
  result[["Movement Measure"]] <- 
    lapply(bootstrap_result, function(x) x[["Movement Measure"]])
  result[["Memory Curve"]] <- 
    lapply(bootstrap_result, function(x) x[["Memory Curve"]])
  result[["Expected Outcomes"]] <- 
    lapply(bootstrap_result, function(x) x[["Expected Outcomes"]])
  result[["MTE_vec"]] <- 
    lapply(bootstrap_result, function(x) x[["MTE_vec"]])
  result[["AMTE"]] <- 
    lapply(bootstrap_result, function(x) x[["AMTE"]])
  result[["d_prime"]] <- 
    lapply(bootstrap_result, function(x) x[["d_prime"]])
  result[["marginal_dist_f0"]] <- 
    lapply(bootstrap_result, function(x) x[["marginal_dist_f0"]])
  result[["lambda2"]] <- 
    lapply(bootstrap_result, function(x) x[["lambda2"]])
  result[["rho"]] <- 
    lapply(bootstrap_result, function(x) x[["rho"]])
  result[["t_max"]] <- 
    lapply(bootstrap_result, function(x) x[["t_max"]] )
  result[["Steady State TM"]] <- 
    lapply(bootstrap_result, function(x) x[["Steady State TM"]] )
  result[["MFP"]] <-
    lapply(bootstrap_result, function(x) x[["MFP"]] )
  result[["GeenensD"]] <-
    lapply(bootstrap_result, function(x) x[["GeenensD"]] )
  result[["HellingerDep"]] <-
    lapply(bootstrap_result, function(x) x[["HellingerDep"]] )
   result[["altham_index"]] <-
     lapply(bootstrap_result, function(x) x[["altham_index"]] )
  result[["tm_num"]] <-
    lapply(bootstrap_result, function(x) x[["tm_num"]] )
  return(result)
}
##=============================== Clean the Results ==========================##
clean_results <- function(final_result, times, size){
  # 1. Extract the t_max df:
  t_loc_df <- 
    tibble(
      t_loc = unlist(final_result$t_max)
    )
  t_max <- max(unlist(final_result$t_max))
  # 2. Extract the Movement Measure:
  movement_boot <- final_result$`Movement Measure`[c(1:times)]
  movement_df<- 
    lapply(as.list(seq(1,t_max)), function(t){map_dfr(movement_boot, ~ .x[t, ])})
  # 3. Extract the Memory Curve:
  memory_boot <- final_result$`Memory Curve`[c(1:times)]
  memory_df <-
    lapply(as.list(seq(1,t_max)), function(t){map_dfr(memory_boot, ~ .x[t, ])})
  # 4. Extract AEO(vary by time):
  AEO_boot <- 
    lapply(final_result$`Expected Outcomes`[c(1:times)],function(x){x[["AEO"]]})
  AEO_df <- do.call(rbind, AEO_boot)
  colnames(AEO_df) <- paste0("t", seq(t_max))
  # 5. Extract EO(vary by class):
  EO_t1_boot <- 
    lapply(final_result$`Expected Outcomes`[c(1:times)],function(x){x[["EO"]]})
  EO_t1_df <- do.call(rbind, EO_t1_boot)
  colnames(EO_t1_df) <- paste("Class", seq(size))
  # 6. Extract EO_tss(Vary by class):
  EO_tss_boot <- 
    lapply(final_result$`Expected Outcomes`[c(1:times)],function(x){x[["EO_tss"]]})
  EO_tss_df <- do.call(rbind, EO_tss_boot)
  colnames(EO_tss_df) <- paste("Class", seq(size))
  # 7. Extract MTE(Vary by Class) & AMTE, put them into the same dataframe:
  MTE_boot <- final_result$`MTE_vec`[c(1:times)]
  MTE_df <- bind_rows(MTE_boot)
  AMTE <- unlist(final_result$AMTE[c(1:times)])
  MTE_df <- cbind(MTE_df,AMTE)
  # 8. Extract d prime:
  d_prime_boot <- final_result$d_prime[c(1:times)]
  d_prime_df_lst <-
    lapply(as.list(seq(1,t_max)), function(t){map_dfr(d_prime_boot, ~ .x[t, ])})
  d_prime_df <-
    lapply(
      seq_len(length(d_prime_df_lst)), 
      function(t){
        colnames(d_prime_df_lst[[t]]) <- c("d_prime","combination")
        d_prime_df_lst[[t]] %>% mutate(t = t)
      })
  
  # 9.  Lambda 2:
  lambda2_df <- tibble(lambda2 = unlist(final_result$lambda2[1:times]))
  # 10. rho:
  rho_df <- tibble(rho = unlist(final_result$rho[1:times]))
 
  # 11.  Extract the MTP result:
  MFP_lst <- final_result$MFP
  MFP_df <-
    map2_dfr(
      MFP_lst,
      .y = seq_along(MFP_lst),  # bootstrap index
      ~ as.data.frame(.x) %>%
        mutate(parent_status = rownames(.x)) %>%
        pivot_longer(
          cols = -parent_status,
          names_to = "children_status",
          values_to = "value"
        ) %>%
        mutate(bootstrap_index = .y)
    )
  # 12. Steady State TM:
  steady_state_TM <- final_result$`Steady State TM`
  # 13. tm_num:
  tm_num <- final_result$tm_num
  # 14. GeenensD:
  GeenensD_df <- tibble(GeenensD = unlist(final_result$GeenensD[1:times]))
  # 15. HellingerDep:
  HellingerDep_df <- tibble(HellingerDep = unlist(final_result$HellingerDep[1:times]))
  # Altham Index:
  AI_df <- tibble(altham = unlist(final_result$altham_index[1:times]))
  
  cleaned_result <- list()
  cleaned_result[["Movement Measure"]] = movement_df
  cleaned_result[["Memory Measure"]] = memory_df
  cleaned_result[["AEO"]] = AEO_df
  cleaned_result[["EO_t1"]] = EO_t1_df
  cleaned_result[["EO_tss"]] = EO_tss_df
  cleaned_result[["MTE"]] = MTE_df
  cleaned_result[["d_prime"]] = d_prime_df
  cleaned_result[["lambda2"]] = lambda2_df
  cleaned_result[["rho"]] = rho_df
  cleaned_result[["t_loc"]] = t_loc_df
  cleaned_result[["MFP"]] = MFP_df
  cleaned_result[["tm_num"]] = tm_num
  cleaned_result[["steady_state_TM"]] = steady_state_TM
  cleaned_result[["HellingerDep"]] = HellingerDep_df
  cleaned_result[["altham_index"]] <- AI_df
  cleaned_result[["GeenensD"]] = GeenensD_df
  return(cleaned_result)
}
##=============================== Calculate the Stats ========================##
calculate_stats <- function(cleaned_result,size){
  #-----------------Clean the Movement Measure--------------------#
  movement_lst <- cleaned_result$`Movement Measure`
  movement_df_lst <-
    lapply(
      seq(length(movement_lst)),
      function(i){
        movement_lst[[i]] <-
          movement_lst[[i]] %>%
          mutate(
            id = row_number(),
            t = i
          ) %>%
          drop_na()
      }
    )
  # Generate movement stat:
  movement_stats_df <-
    bind_rows(movement_df_lst, .id = "t") %>%
    group_by(
      t
    ) %>%
    summarise(
      across(
        -id,
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      cols = -t,  
      names_to = c("category", ".value"),
      names_pattern = "^(.*? Mobility)_(.*)$"
    ) %>%
    ungroup()
    
  

  
  #----------------------------Clean the Memory Measure-------------------------#
  #---------------------- Memory Across Different Generation-------------------#
  memory_lst <- cleaned_result$`Memory Measure`
  memory_df_lst <-
    lapply(
      memory_lst,
      function(memory_df){
        data <-
          memory_df %>%
          mutate(across(
            -t,                                 
            ~ log(.x),                         
            .names = "l{.col}"            
          )) %>%
          summarise(
            across(
              -t,
              list(
                mean = ~ mean(.x[is.finite(.x)], na.rm = TRUE),
                sd = ~ sd(.x[is.finite(.x)], na.rm = TRUE),
                CI_np_l = ~ quantile(.x[is.finite(.x)], 0.025, na.rm = TRUE),
                CI_np_u = ~ quantile(.x[is.finite(.x)], 0.975, na.rm = TRUE)
              ),
              .names = "{.col}_{.fn}"
            )) %>%
          pivot_longer(
            cols = everything(),
            names_to = c("class", ".value"),
            names_pattern = "(^[^_]+)_(.*)",
            values_to = "value"
          ) %>%
          mutate(
            CI_np = paste("(",CI_np_l,",",CI_np_u,")")
          )
      }
    )
  memory_df <-
    do.call(
      rbind,
      memory_df_lst 
    ) %>%
    mutate(
      t = rep(seq(length(memory_df_lst)),each = (size +1) * 2))

  #------------------------- Memory Change Speed ---------------------------#
  
  memory_df_lst <-
    lapply(
      memory_lst,
      function(memory_df){
        data <-
          memory_df %>%
          mutate(across(
            -t,                                 
            ~ log(.x),                         
            .names = "l{.col}"            
          )) 
      }
    )
  
  #memory_df_lst <-  cleaned_result$`Memory Measure`
  
  memory_slope_df <- reduce(
    list(
      memory_df_lst[[1]] %>% dplyr::select(-t),
      memory_df_lst[[2]] %>% dplyr::select(-t)
    ),
    `-`
    )  %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x[is.finite(.x)], na.rm = TRUE),
          sd = ~ sd(.x[is.finite(.x)], na.rm = TRUE),
          CI_np_l = ~ quantile(.x[is.finite(.x)], 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x[is.finite(.x)], 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("class", ".value"),
      names_pattern = "(^[^_]+)_(.*)",
      values_to = "value"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  
  #---------------------------Clean the AEO Measure----------------------------#
  AEO <- as.data.frame(cleaned_result$AEO)
  AEO_df <-
    AEO %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      cols = starts_with("t"),  
      names_to = c("time", ".value"),
      names_pattern = "t(\\d+)_(.*)"
    ) %>%
    mutate(time = as.integer(time)) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  #-----------------------------Clean the EO t1--------------------------------#
  EO_t1 <- as.data.frame(cleaned_result$EO_t1)
  EO_t1_df <-
    EO_t1 %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      cols = everything(),  
      names_to = c("Class", ".value"),
      names_pattern = "^(Class \\d+)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  #-----------------------------Clean the EO tss--------------------------------#
  EO_tss <- as.data.frame(cleaned_result$EO_tss)
  EO_tss_df <-
    EO_tss %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      cols = everything(),  
      names_to = c("Class", ".value"),
      names_pattern = "^(Class \\d+)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  
  #-----------------------------Clean the MTE--------------------------------#
  MTE_df <- cleaned_result$MTE
  clean_data <- function(data) {
    data <- lapply(data, function(col) {
      col <- ifelse(is.infinite(col), NA, col) # Replace Inf/-Inf
      return(col)
    })
    return(as.data.frame(data))
  }
  
  MTE_df <- clean_data(MTE_df)
  MTE_stat_df <-
    MTE_df %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("class", ".value"),
      names_pattern = "(^[^_]+)_(.*)",
      values_to = "value"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  
  #------------------------------Clean the lambda2-----------------------------#
  lambda2_df <- 
    cleaned_result$lambda2 %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        )
      )
    ) %>%
    pivot_longer(
      cols = everything(),  
      names_to = c("Attribute", ".value"),
      names_pattern = "^(lambda2)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  #-----------------------------Clean the rho--------------------------------#
  rho_df <-
    cleaned_result$rho %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        )
      )
    ) %>%
    pivot_longer(
      cols = everything(),  
      names_to = c("Attribute", ".value"),
      names_pattern = "^(rho)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  #-----------------------------Clean the t_loc--------------------------------#
  t_loc_df <- cleaned_result$t_loc
  t_loc_stat_df <-
    t_loc_df %>%
    tabyl(t_loc)
  
  #-----------------------------Clean the MFP--------------------------------#
 MFP_clean <- cleaned_result$MFP
 MFP_stat_df <-
   MFP_clean %>%
   group_by(
     parent_status,
     children_status
   ) %>%
   summarise(
     mean = mean(value, na.rm = TRUE),
     sd = sd(value, na.rm = TRUE),
     CI_np_l = quantile(value, 0.025, na.rm = TRUE),
     CI_np_u = quantile(value, 0.975, na.rm = TRUE)
   ) %>%
   filter(
     parent_status != children_status
   )

  #-----------------------------Clean the d prime--------------------------------#
  # Calculate stats for d prime:
  d_prime_lst <- cleaned_result$d_prime
  d_prime_stat_lst <-
    lapply(
      d_prime_lst,
      function(d_prime_sublst){
        d_prime_val <-
          d_prime_sublst %>%
          summarise(
            across(
              -c("t","combination"),
              list(
                mean = ~ mean(.x, na.rm = TRUE),
                sd = ~ sd(.x, na.rm = TRUE),
                CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
                CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
              ),
              .names = "{.col}_{.fn}"
            )
          )
        d_prime_comb <-
          d_prime_sublst %>%
          tabyl(
            combination
          ) %>%
          dplyr::select(
            -n
          )
        
        return(
          list(
            val = d_prime_val,
            comb = d_prime_comb
          )
        )
      }
    )
  
  d_prime_val_df <-
    bind_rows(
      lapply(
        seq_len(length(d_prime_stat_lst)),
        function(i){
          d_prime_stat_lst[[i]]$val %>%
            mutate(t = i)
        }
      ))
  d_prime_comb_df <-
    bind_rows(
      lapply(
        seq_len(length(d_prime_stat_lst)),
        function(i){
          d_prime_stat_lst[[i]]$comb %>%
            mutate(t = i)
        }
      ))
  
  #-----------------------------Clean the Geenen's D----------------------------# 
  GeenensD_df <- 
    cleaned_result$GeenensD %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        )
      )
    ) %>%
    pivot_longer(
      cols = everything(),  
      names_to = c("Attribute", ".value"),
      names_pattern = "^(GeenensD)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  #-----------------------Clean the Hellinger's Dependence---------------------# 
  HellingersDep_df <- 
    cleaned_result$HellingerDep %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        )
      )
    ) %>%
    pivot_longer(
      cols = everything(),  
      names_to = c("Attribute", ".value"),
      names_pattern = "^(HellingerDep)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )
  
  #-----------------------Clean the Altham Index---------------------# 
  AI_df <-
    cleaned_result$altham_index %>%
    summarise(
      across(
        everything(),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          CI_np_l = ~ quantile(.x, 0.025, na.rm = TRUE),
          CI_np_u = ~ quantile(.x, 0.975, na.rm = TRUE)
        )
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("Attribute", ".value"),
      names_pattern = "^(altham)_(.*)$"
    ) %>%
    mutate(
      CI_np = paste("(",CI_np_l,",",CI_np_u,")")
    )

  
  #-----------------------Clean the Steady State TM----------------------------# 
  steady_state_TM_lst <- 
    lapply(cleaned_result$steady_state_TM, function(tm){as.matrix(tm)})
  steady_state_TM_array <- simplify2array(steady_state_TM_lst)
  steady_state_TM_stats_lst <-
    list(
      apply(steady_state_TM_array, c(1, 2), mean),
      apply(steady_state_TM_array, c(1, 2), sd),
      apply(steady_state_TM_array, c(1, 2), function(x) quantile(x, probs = 0.025)),
      apply(steady_state_TM_array, c(1, 2), function(x) quantile(x, probs = 0.975))
    )
  attribute_lst <- c("mean","sd","CI_np_l","CI_np_u")
  steady_state_TM_stats_df_lst <-
    lapply(
      seq_len(length(steady_state_TM_stats_lst)),
      function(i){
        as.data.frame(
          steady_state_TM_stats_lst[[i]]) %>%
          mutate(
            Row = row_number()
          ) %>%
          pivot_longer(-Row, names_to = "Column", values_to = "Value") %>%
          mutate(Column = as.numeric(gsub("Class", "", Column)),   
                 Pair = paste(Row, Column, sep = "-")) %>%         
          dplyr::select(Pair, Value) %>%
          rename(
            !!attribute_lst[i] := Value
          )
      }
    )
  steady_state_TM_df <- reduce(steady_state_TM_stats_df_lst, left_join, by = "Pair")
  
  #-----------------------------Clean the tm_num-------------------------------# 
  tm_num_lst <- 
    lapply(cleaned_result$tm_num, function(tm){as.matrix(tm %>% dplyr::select(-rowsum,-validation))})
  tm_num_array <- simplify2array(tm_num_lst)
  tm_num_stats_lst <-
    list(
      apply(tm_num_array, c(1, 2), mean),
      apply(tm_num_array, c(1, 2), sd),
      apply(tm_num_array, c(1, 2), function(x) quantile(x, probs = 0.025)),
      apply(tm_num_array, c(1, 2), function(x) quantile(x, probs = 0.975))
    )
  attribute_lst <- c("mean","sd","CI_np_l","CI_np_u")
  tm_num_stats_df_lst <-
    lapply(
      seq_len(length(tm_num_stats_lst)),
      function(i){
        as.data.frame(
          tm_num_stats_lst[[i]]) %>%
          mutate(
            Row = row_number()
          ) %>%
          pivot_longer(-Row, names_to = "Column", values_to = "Value") %>%
          mutate(Column = as.numeric(gsub("Class", "", Column)),   
                 Pair = paste(Row, Column, sep = "-")) %>%         
          dplyr::select(Pair, Value) %>%
          rename(
            !!attribute_lst[i] := Value
          )
      }
    )
  tm_num_df <- reduce(tm_num_stats_df_lst, left_join, by = "Pair")
  stats_result <- list()
  stats_result[["Movement Measure"]] <- movement_stats_df
  stats_result[["t"]] <- t_loc_stat_df 
  stats_result[["Memory Measure"]] <- memory_df
  stats_result[["AEO"]] <- AEO_df
  stats_result[["EO_t1"]] <- EO_t1_df
  stats_result[["EO_tss"]] <- EO_tss_df
  stats_result[["MTE"]] <- MTE_stat_df
  stats_result[["lambda2"]] <- lambda2_df
  stats_result[["rho"]] <- rho_df
  stats_result[["d_prime_val"]] <- d_prime_val_df
  stats_result[["d_prime_comb"]] <- d_prime_comb_df
  stats_result[["GeenensD"]] <- GeenensD_df
  stats_result[["HellingersDep"]] <- HellingersDep_df
  stats_result[["steady_state_TM"]] <- steady_state_TM_df
  stats_result[["tm_num"]] <- tm_num_df
  stats_result[["MFP"]] <- MFP_stat_df
  stats_result[["altham_index"]] <- AI_df
  stats_result[["memory_slope_df"]] <- memory_slope_df
  return(stats_result)
}
##=============================== Final Function  ============================##
chunk_list <- function(data_list, n_cores) {
  # Total length of the list
  N <- length(data_list)
  # Calculate chunk size
  chunk_size <- n_cores
  # Split the list into chunks
  chunks <- split(data_list, ceiling(seq_along(data_list) / chunk_size))
  return(chunks)
}


process_results <- function(boot_rst, times, size, n_cores, reshape){
  # Step 0: Convert the input form:
  if(reshape == TRUE){
    bootstrap_result_lst <- modify_cohort_lsts(boot_rst, qt = 0.95)}else{
      bootstrap_result_lst <- boot_rst
    }
  bootstrap_result_lst <- bootstrap_result_lst[
    !sapply(bootstrap_result_lst, function(x) all(is.na(x[[1]])))
  ]
  cohort_vec <- names(bootstrap_result_lst)
  assert_that(
    all(diff(as.numeric(cohort_vec)) == 1)
  )
  print(paste("The current result starts from", min(as.numeric(cohort_vec)),
              "to", max(as.numeric(cohort_vec))))
  # Step 1: Fetch the results, remove the data to release memory:
  final_result_lst <- 
    mclapply(bootstrap_result_lst,fetch_results, mc.cores = n_cores)
  rm(bootstrap_result_lst)
  gc()
  # Step 2: Break the results to chunks for cleaning:
  chunk_list_for_clean <- chunk_list(final_result_lst, n_cores)
  cleaned_lst <- 
    lapply(
      chunk_list_for_clean, 
      function(chunk){
        mclapply(
          chunk, 
          function(final_result){
            rst <-
              clean_results(final_result,times,size)
            return(rst)
          }, 
          mc.cores = n_cores)
      }
    )
  cleaned_results <- unname(unlist(cleaned_lst, recursive = FALSE))
  rm(cleaned_lst)
  gc()
  # Step 3: Calculate the stats:
  chunk_list_for_stat <- chunk_list(cleaned_results, n_cores)
  stat_lst <- 
    lapply(
      chunk_list_for_stat, 
      function(chunk){
        mclapply(
          chunk,   
          function(cleaned_result){
            calculate_stats(cleaned_result,size)}, mc.cores = n_cores)
      }
    )
  stat_results <- unname(unlist(stat_lst, recursive = FALSE))
  return(stat_results)
}



##===========================Validation the Result============================##

##----------------- Section 1: Calculate the Point Estimate ------------------##
generate_baseline_fitted_data <- function(
    gss_data, size, age_l, age_u, year_l, year_u, gamma_par){
  # Filter the year:
  gss_data <-
    gss_data %>%
    filter(
      cohort >= year_l & cohort <= (year_u + 1)
    )
  # Step 2: Generate the tm_num_lst for the regression:
  cs_tm_df_5_2555_1 <- cross_tm_df(gss_data, age_l, age_u, 1)
  cs_tm_5_2555_1_num_lst <- list()
  for(i in seq(length(cs_tm_df_5_2555_1))){
    data <- cs_tm_df_5_2555_1[[i]]
    if(size == 5){
    tms <- transition_matrix_5class_fun(data)
    }else if(size == 6){
    tms <- transition_matrix_6class_fun(data)
    }else{
    tms <- transition_matrix_7class_fun(data)  
    }
    cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
  }
  # Step 3: Generate the regression data:
  regression_data_5_2555 <-
    generate_model_data(FALSE, NULL, cs_tm_5_2555_1_num_lst, size, year_l, year_u) %>%
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
      gamma = gamma_par,
      family = poisson(link="log"),
      select=FALSE, # Keep this as FALSE, avoiding oversmoothing
      method = "REML",
      data = regression_data_5_2555,
      control = control_params,
      gc.level = 2,
      nthreads = 3)
  # Step 5: Generate the model fitted data:
  fitted_data_5_2555 <-
    generate_model_data(
      TRUE, model_baseline, cs_tm_5_2555_1_num_lst, size, year_l, year_u)
  return(fitted_data_5_2555)
}


##----------------- Section 2: Calculate the Point Estimate ------------------##
single_simulation_lst_baseline <- function(
    gss_data,size,age_l,age_u, year_l, year_u, est_l, est_u, tol, gamma_par){
  # Generate the fitted data:
  fitted_data_5_2555 <- generate_baseline_fitted_data(
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
        return(rst)
      }
    )
  
  return(simulation_rst_lst)
}


calculate_baseline <- function(
  gss_data, size, age_l, age_u, year_l, year_u, est_l, est_u, tol, gamma_par){
  # Calculate baseline results:
  baseline_rst <- 
    mclapply(
      1,
      function(i){
        rst <-
          single_simulation_lst_baseline(
            gss_data,size,age_l,age_u, year_l, year_u, est_l, est_u, tol, gamma_par)
      },
      mc.cores = n_cores
    )
  # Process Baseline Rst:
  process_baseline_lst <- 
    process_results(baseline_rst, 1, size, n_cores, TRUE)
  
  return(process_baseline_lst)
}

calculate_baseline_newform <- function(
    gss_data, size, age_l, age_u,  year_l, year_u, bin_width, method){
  # Calculate baseline results:
  baseline_rst <- 
    mclapply(
      1,
      function(i){
        rst <-
          single_simulation_lst_newform(
            gss_data = gss_data,
            age_l = age_l,
            age_u = age_u, 
            year_l = year_l, 
            year_u = year_u, 
            boot = FALSE, 
            bin_width = bin_width, 
            method = method
            )
      },
      mc.cores = n_cores
    )
  # Process Baseline Rst:
  process_baseline_lst <- 
    process_results(baseline_rst, 1, size, n_cores, TRUE)
  
  return(process_baseline_lst)
}

##----------------- Section 2: Calculate the Raw Estimate --------------------##

single_simulation_lst_raw <- function(
    gss_data,size,age_l,age_u, year_l, year_u, tol){
  # Filter the year:
  gss_data <-
    gss_data %>%
    filter(
      cohort >= year_l & cohort <= (year_u + 1)
    )
  # Step 2: Generate the tm_num_lst for the regression:
  cs_tm_df_5_2555_1 <- cross_tm_df(gss_data, age_l, age_u, 1)
  cs_tm_5_2555_1_num_lst <- list()
  for(i in seq(length(cs_tm_df_5_2555_1))){
    data <- cs_tm_df_5_2555_1[[i]]
    tms <- transition_matrix_5class_fun(data)
    cs_tm_5_2555_1_num_lst[[i]] <- tms$trans_matrix_num
  }
  # Step 3: Generate the regression data:
  raw_data_5_2555 <- 
    generate_model_data(FALSE, NULL, cs_tm_5_2555_1_num_lst, size, year_l, year_u) %>%
    rename(
      round_predictions = cell_count
    )
  # Convert the fitted_data to a list:
  cs_tm_df_5_2555_1_raw_lst <- generate_tm_df_fit(raw_data_5_2555)
  # Only focus on the target range:
  cs_tm_df_5_2555_1_raw_lst <-
    cs_tm_df_5_2555_1_raw_lst[names(cs_tm_df_5_2555_1_raw_lst) >= "1945" & 
                                   names(cs_tm_df_5_2555_1_raw_lst) <= "1992"]
  # Generate the simulation result for each of the cohort:
  simulation_rst_lst <- 
    lapply(
      cs_tm_df_5_2555_1_raw_lst,
      function(tm_df_org){
        rst <-
          single_simulation(tm_df_org,size, tol)
        return(rst)
      }
    )
  return(simulation_rst_lst)
}


calculate_raw <- function(
    gss_data, size, age_l, age_u,  year_l, year_u, tol){
  # Calculate baseline results:
  raw_rst <- 
    mclapply(
      1,
      function(i){
        rst <-
          single_simulation_lst_raw(
            gss_data,size,age_l,age_u, year_l, year_u, tol)
      },
      mc.cores = n_cores
    )
  # Process Baseline Rst:
  process_raw_lst <- 
    process_results(raw_rst, 1, size, n_cores, TRUE)
  
  return(process_raw_lst)
}


##----------------------- Section 3: Validation the result -------------------##

valid_bootstrap_results <- function(boot_rst_df, baseline_rst_df, merge_vec){
  # Step 0: Merge two dataframe together:
  merged_df <-
    left_join(
      boot_rst_df,
      baseline_rst_df,
      by = merge_vec,
      suffix = c("_boot", "_baseline")
    ) %>%
    mutate(
      mean_bc = 2 * mean_baseline - mean_boot,
      CI_bc_l = 2 * mean_baseline - CI_np_u_boot,
      CI_bc_u = 2 * mean_baseline - CI_np_l_boot
    ) %>%
    mutate(
      valid =((mean_bc >= CI_bc_l - 1e-8) & (mean_bc <= CI_bc_u + 1e-8) | mean_baseline == 1)
    ) %>%
    dplyr::select(
      mean_baseline,
      mean_bc,
      CI_bc_u,
      CI_bc_l,
      everything()
    )
  return(merged_df)
}

valid_rst_per_cohort <- function(boot_rst, process_baseline){
  validation_lst <- list()
  record <- 
    tibble(
      Movement = NA,
      Memory = NA,
      AEO = NA,
      EO_t1 = NA,
      EO_tss = NA,
      MTE = NA,
      lambda2 = NA,
      rho = NA,
      d_prime = NA,
      MFP = NA,
      GeenensD = NA,
      HellingerDep = NA
    )
  # 1. Modify the Movement Measure:
  movement_df <- 
    valid_bootstrap_results(
      boot_rst$`Movement Measure`, 
      process_baseline$`Movement Measure`, 
      c("category","t"))
  # Record:
  record$Movement <- sum(!movement_df$valid, na.rm = TRUE)
  
  # 2. Modify the Memory Measure:
  memory_df <-
    valid_bootstrap_results(
      boot_rst$`Memory Measure`, 
      process_baseline$`Memory Measure` %>%
        mutate(
          mean = ifelse(mean == 0, NA, mean)
        ), 
      c("class","t")) %>%
    mutate(
      flag = case_when(
        class  == "Class1" & mean_bc < 0 & t %in% c(1) ~ 1,
        class  == "Class2" & mean_bc < 0 & t %in% c(1) ~ 1,
        class  == "Class3" & mean_bc < 0 & t %in% c(1) ~ 1,
        class  == "Class4" & mean_bc < 0 & t %in% c(1) ~ 1,
        class  == "Class5" & mean_bc < 0 & t %in% c(1) ~ 1,
        .default = 0
      )
    ) %>%
    dplyr::select(
      class,
      t,
      flag,
      everything()
    )
  # Record:
  record$Memory <- sum(memory_df$flag, na.rm = TRUE)
  # 3. Modify the AEO Measure:
  AEO_df <-
    valid_bootstrap_results(
      boot_rst$AEO, 
      process_baseline$AEO, 
      c("time"))
  # Record:
  record$AEO <- sum(!AEO_df$valid, na.rm = TRUE)
  # 4. EO_t1:
  EO_t1_df <-
    valid_bootstrap_results(
      boot_rst$EO_t1, 
      process_baseline$EO_t1, 
      c("Class"))
  # Record:
  record$EO_t1 <- sum(!EO_t1_df$valid,na.rm = TRUE)
  # 5. EO_t1:
  EO_tss_df <-
    valid_bootstrap_results(
      boot_rst$EO_tss, 
      process_baseline$EO_tss, 
      c("Class"))
  # Record:
  record$EO_tss <- sum(!EO_tss_df$valid, na.rm = TRUE)
  # 6. MTE:
  MTE_df <-
    valid_bootstrap_results(
      boot_rst$MTE, 
      process_baseline$MTE, 
      c("class"))
  # Record:
  record$MTE <- sum(!MTE_df$valid, na.rm = TRUE)
  # 7. lambda2:
  lambda2_df <-
    valid_bootstrap_results(
      boot_rst$lambda2, 
      process_baseline$lambda2, 
      c("Attribute"))
  # Record:
  record$lambda2 <- sum(!lambda2_df$valid, na.rm = TRUE)
  # 8. rho:
  rho_df <-
    valid_bootstrap_results(
      boot_rst$rho, 
      process_baseline$rho, 
      c("Attribute"))
  record$rho <- sum(!rho_df$valid, na.rm = TRUE)
  # 9. d prime:
  d_prime_df <-
    valid_bootstrap_results(
      boot_rst$d_prime_val %>%
        rename_with(~ str_replace_all(.x, "d_prime_", "")) %>%
        mutate(Attribute = "d_prime"), 
      process_baseline$d_prime_val %>%
        rename_with(~ str_replace_all(.x, "d_prime_", "")) %>%
        mutate(Attribute = "d_prime"), 
      c("t","Attribute")) %>%
    mutate(
      flag = ifelse(t == 1 & mean_baseline < mean_boot, 1, 0)
    )
  record$d_prime <- 
    sum(!d_prime_df$valid, na.rm = TRUE)
  # 10 MFP:
  MFP_df <-
    valid_bootstrap_results(
      boot_rst$MFP,
      process_baseline$MFP,
      c("children_status","parent_status")) %>%
    mutate(
      flag = ifelse(mean_bc < 0, 1, 0)
    )

  record$MFP <- sum(!MFP_df$valid, na.rm = TRUE)
  # 11. GeenensD:
  GeenensD_df <-
    valid_bootstrap_results(
      boot_rst$GeenensD, 
      process_baseline$GeenensD, 
      c("Attribute"))
  # Record:
  record$GeenensD <- sum(!GeenensD_df$valid, na.rm = TRUE)
  # 12. HellingerDep
  HellingersDep_df <-
    valid_bootstrap_results(
      boot_rst$HellingersDep, 
      process_baseline$HellingersDep, 
      c("Attribute"))
  
  # 13. Altham Index:
  AI_df <-
    valid_bootstrap_results(
      boot_rst$altham_index,
      process_baseline$altham_index,
      c("Attribute"))
  
  # Memory Slope Df:
  memory_slope_df <-
    valid_bootstrap_results(
      boot_rst$memory_slope_df,
      process_baseline$memory_slope_df,
      c("class")
    )
  
  # Record:
  record$HellingerDep <- sum(!HellingersDep_df$valid, na.rm = TRUE)
  
  validation_lst[["movement_df"]] <- movement_df
  validation_lst[["memory_df"]] <- memory_df
  validation_lst[["AEO_df"]] <- AEO_df
  validation_lst[["EO_t1_df"]] <- EO_t1_df
  validation_lst[["EO_tss_df"]] <- EO_tss_df
  validation_lst[["MTE_df"]] <- MTE_df
  validation_lst[["lambda2_df"]] <- lambda2_df
  validation_lst[["rho_df"]] <- rho_df
  validation_lst[["d_prime_df"]] <- d_prime_df
  validation_lst[["MFP"]] <- MFP_df
  validation_lst[["GeenensD"]] <- GeenensD_df
  validation_lst[["HellingersDep"]] <- HellingersDep_df
  validation_lst[["altham"]] <- AI_df
  validation_lst[["memory_slope"]] <- memory_slope_df
  validation_lst[["t_stat_df"]] <- boot_rst$t
  validation_lst[["t_baseline"]] <- process_baseline$t
  validation_lst[["record"]] <- record
  return(validation_lst)
}


##----------------------- Section 3: Validation the result -------------------##

generate_validation_data <- function(process_boot_rst_lst, process_baseline_lst){
  validation_lst <-
    lapply(
      seq(length(process_boot_rst_lst)),
      function(i){
        boot_rst <- process_boot_rst_lst[[i]]
        baseline_rst <- process_baseline_lst[[i]]
        validation_rst <- valid_rst_per_cohort(boot_rst, baseline_rst)
        return(validation_rst)
      }
    )
  
  record_lst <-
    lapply(
      validation_lst,
      function(validation_cohort){
        record <-
          validation_cohort$record
      }
    )
  record_df <- bind_rows(record_lst)
  return(
    list(
      validation_lst = validation_lst,
      record_df = record_df
    )
  )
}


##===================== Function : Generate Obs data =======================##

generate_obs_df <- function(boot_lst, baseline_lst, raw_data, size){
  # Boot # obs:
  boot_pairs <-
    bind_rows(
      lapply(
        boot_lst,
        function(x){
          x$tm_num %>%
            rename_with(~ paste0("boot_", .), -Pair)
        }
      )
    ) %>%
    mutate(
      cohort = rep(seq(1945,1992,1), each = size * size)
    )
  # Baseline # obs:
  baseline_pairs <-
    bind_rows(
      lapply(
        baseline_lst,
        function(x){
          x$tm_num %>%
            rename_with(~ paste0("baseline_", .), -Pair)
        }
      )
    ) %>%
    mutate(
      cohort = rep(seq(1945,1992,1), each = size * size)
    )
  # Raw Data # Obs:
  raw_pairs <-
    gss_data %>%
    filter(
      cohort >= 1945 & cohort <= 1992
    ) %>%
    filter(
      age >= 25 & age <= 55
    ) %>%
    mutate(
      Pair = paste0(status_p,"-",status_c)
    ) %>%
    group_by(
      cohort,
      Pair
    ) %>%
    summarise(
      count_raw = n()
    )
  
  # Final df:
  final_obs_df <-
    reduce(
      list(baseline_pairs, boot_pairs, raw_pairs), 
      left_join, 
      by = c("Pair", "cohort")
    ) %>%
    mutate(
      boot_mean_round = round(boot_mean),
      mean_bc = 2 * baseline_mean - boot_mean_round,
      CI_bc_l = 2 * baseline_mean - boot_CI_np_u,
      CI_bc_u = 2 * baseline_mean - boot_CI_np_l
    ) 
  return(final_obs_df)
}


## ==================== Function: # of obs visualization =================##

visualize_tm_num <- function(obs_data, group_num){
  group_vec <-paste0(group_num,"-",seq(1,5)) 
  pair_colors <- setNames(RColorBrewer::brewer.pal(length(group_vec), "Set2"), group_vec)
  df_filtered <- 
    obs_data %>%
    filter(
      Pair %in% group_vec)
  
  p <- plot_ly()
  
  for (pair in unique(df_filtered$Pair)) {
    
    df_sub <- df_filtered %>% filter(Pair == pair)  
    pair_color <- pair_colors[pair] 
    is_visible <- ifelse(pair == group_vec[1], TRUE, "legendonly")
    
    # Add Mean_bc line with circular markers
    p <- p %>%
      add_trace(
        data = df_sub, 
        x = ~cohort, 
        y = ~mean_bc, 
        type = 'scatter', 
        mode = 'lines+markers',
        marker = list(symbol = "circle", color = pair_color, size = 8),  
        line = list(color = pair_color),  
        name = paste0(pair, " - Mean"),
        legendgroup = pair,  
        visible = is_visible, 
        showlegend = TRUE) %>%
      
      add_ribbons(
        data = df_sub, 
        x = ~cohort, 
        ymin = ~CI_bc_l, 
        ymax = ~CI_bc_u, 
        fillcolor = toRGB(pair_color, alpha = 0.2),  
        line = list(color = "transparent"),  
        name = paste0(pair, " - CI"), 
        legendgroup = pair,  
        visible = is_visible,  
        showlegend = FALSE) %>%
      
      add_trace(
        data = df_sub, 
        x = ~cohort, 
        y = ~count_raw, 
        type = 'scatter', 
        mode = 'lines+markers',
        marker = list(symbol = "triangle-up", color = pair_color, size = 4),  
        line = list(color = pair_color, dash = "solid", width = 0.7),  
        name = paste0(pair, " - Raw"),
        legendgroup = pair,  
        visible = is_visible,  
        showlegend = TRUE)
  }
  
  p <- p %>%
    layout(title = paste("# of Obs of Class", group_num),
           xaxis = list(title = "Cohort"),
           yaxis = list(title = "# of Obs"),
           legend = list(title = list(text = "Pair")))
  return(p)
}


## =========================== Visualize D prime ============================##

visualize_d_prime_comb <- function(boot_rst){
  
  D_prime_percent <-
    bind_rows(
      lapply(
        boot_rst,
        function(sub_lst){
          sub_lst$d_prime_comb %>%
            filter(
              t == 1
            ) %>%
            arrange(
              percent
            ) %>%
            summarise_all(
              last
            ) %>%
            dplyr::select(
              combination,
              percent
            )
        }
      )) %>%
    mutate(
      cohort = seq(1945, 1992, 1)
    )
  
  p <-
    plot_ly(
      data = D_prime_percent, 
      x = ~cohort, 
      y = ~combination, 
      type = 'scatter', 
      mode = 'markers') %>%
    layout(title = "D prime Mix Level",
           xaxis = list(title = "cohort"),
           yaxis = list(title = "D Prime %"))
  return(p)
}

## =========================== Check ergodic  ============================##

check_ergodic <- function(boot_rst){
  MFP_tab <-
    bind_rows(
      lapply(
        boot_rst,
        function(cohort_lst){
          MFP <-
            bind_rows(
              lapply(
                cohort_lst,
                function(cohort){
                  if(all(is.na(cohort$MFP))){
                    return(NA)
                  }else{
                    return(nrow(cohort$MFP))
                  }
                }
              )
            )
        }
      )
    )
  # Generate a summary table:
  return(as.data.frame(which(is.na(MFP_tab), arr.ind = TRUE)))
}

## =========================== Bootstrap Validation  ============================##
valid_boot <- function(boot_rst, valid_rate, boot_times, valid_times){
  non_primitive_tm <-
    boot_rst %>%
    lapply(
      function(boot) 
        which(is.na(boot))
    ) %>%
    bind_rows()
  # Step 1: Generate a summary data for non-primitive rate:
  if(ncol(non_primitive_tm != 0)){
    non_primitive_tm <-
      non_primitive_tm %>%
      summarise(
        across(everything(), ~ sum(!is.na(.)))) %>%
      pivot_longer(
        everything(), 
        names_to = "cohort",
        values_to = "count") %>%
      mutate(
        cohort = as.integer(cohort),
        rate = count / boot_times,
        percent = round(100 * rate, 2),
        no_valid = ifelse(rate >= valid_rate, 1, 0)
      )}else{
        assert_that(
          nrow(check_ergodic(boot_rst)) == 0
        ) 
        return(list(
          boot_rst = boot_rst[c(1:valid_times)])
        )
      }
  
  # Step 2: Generate the current valid cohort:
  cohort_vec <- seq(1945, 1990)
  valid_cohort <- setdiff(
    cohort_vec,
    non_primitive_tm %>%
      filter(no_valid == 1) %>%
      pull(cohort)
  )
  
  assert_that(
    all(diff(valid_cohort) == 1)
  )
  
  # Step 3: Focus the bootstrap with this range:
  boot_rst_f <-
    lapply(
      boot_rst,
      function(cohort_lst){
        cohort_lst[
          names(cohort_lst) %in% valid_cohort
        ]
      }
    )
  
  # Step 4: Remove the whole bootstrap sample if there is any non-primitive steady state:
  boot_rst_f <- 
    boot_rst_f[
      vapply(boot_rst_f, function(boot) {
        vals <- unlist(boot, recursive = TRUE, use.names = FALSE)
        !any(is.na(vals))}, logical(1))
    ]
  # Step 4: Double check the irreducible(Primitive actually include ergodic):
  assert_that(
    nrow(check_ergodic(boot_rst_f)) == 0
  ) 
  # Step 5: Return the valid cohort vector along with the valid cohort vector:
  return(
    list(
      boot_rst = boot_rst_f[1:valid_times],
      non_validate_tm = non_primitive_tm,
      valid_cohort = valid_cohort
    )
  )
  
}


