####### Functions about Basic Operation:
####### @copyright: Weiqi Wang


#=================== Function 0: Cut the Cohort Cross-Sectionally =============#
cohort_cut <- function(data, interval_length){
  # Grab the cohort:
  cohort <- 
    data %>%
    dplyr::select(
      cohort
    ) %>%
    drop_na() %>%
    mutate(
      cohort = as.numeric(cohort)
    ) %>%
    pull()
  
  breaks <- seq(min(cohort), max(cohort), by = interval_length)
  if (max(cohort) > max(breaks)) {
    breaks <- c(breaks, max(cohort) + 1)
  }
  
  bins <-
    cut(
      cohort, 
      breaks = breaks, 
      include.lowest = TRUE, 
      right = FALSE)
  
  levels(bins) <- sapply(1:(length(breaks) - 1), function(i) {
    if(breaks[i + 1] < max(cohort)){
      paste0("[", breaks[i], ", ", breaks[i + 1], ")")}else{
        paste0("[", breaks[i], ", ", breaks[i + 1], "]")
      }
  })
  
  
  final_data <-
    cbind(
      bins,
      data
    ) %>%
    arrange(
      -desc(cohort)
    )
  return(final_data)
}
#=================== Function 1: Generate Cross-Sectional DFs==================#

cross_tm_df <- function(data, age_l=0, age_u=100, bin_width){
  cut_data <- cohort_cut(data, bin_width)
  cut_bins_lst <- 
    cut_data %>%
    dplyr::select(
      bins
    ) %>%
    unique() %>%
    pull()
  
  cross_tm_lst <- list()
  for(i in seq(length(cut_bins_lst))){
    cohort_bins <- cut_bins_lst[i]
    temp <-
      cut_data %>%
      filter(
        bins == cohort_bins
      ) %>%
      filter(
        age >= age_l & age <= age_u
      )
    cross_tm_lst[[i]] <- temp
  }
  return(cross_tm_lst)
}


#=================== Function 2: Fetch Non-Empty Cohorts ======================#
fetch_non_empty_cohort <- function(cs_tm_df_lst){
  non_emplty_lst <- list()
  index <- 1
  for(i in seq(length(cs_tm_df_lst))){
    temp <- cs_tm_df_lst[[i]]
    if(nrow(temp) == 0){
      next
    }else{
      non_emplty_lst [[index]] <- temp
      index <- index + 1
    }
  }
  return(non_emplty_lst)
}


#=================== Function 3: From tm to df ======================#
tm_to_df <- function(tm){
  tm <-
    tm %>%
    dplyr::select(
      -rowsum,
      -validation
    )
  data <- 
    data.frame(
    )
  for(i in seq(nrow(tm))){
    for(j in seq(ncol(tm[i,]))){
      temp <- as.data.frame(as.list(c( tm[i,j], i, j)))
      colnames(temp) <- c("cell_count", "parent_status", "children_status")
      data <- rbind(data, temp)
    }
  }
  return(data)
}

#=================== Function 4: 7 Class Transition Matrices ===================#
transition_matrix_7class_fun <- function(data){
  pc_status <- 
    data %>% 
    mutate(
      # Parent is in class 1:
      pc_11 = ifelse(status_p == 1 & status_c == 1, 1, 0),
      pc_12 = ifelse(status_p == 1 & status_c == 2, 1, 0),
      pc_13 = ifelse(status_p == 1 & status_c == 3, 1, 0),
      pc_14 = ifelse(status_p == 1 & status_c == 4, 1, 0),
      pc_15 = ifelse(status_p == 1 & status_c == 5, 1, 0),
      pc_16 = ifelse(status_p == 1 & status_c == 6, 1, 0),
      pc_17 = ifelse(status_p == 1 & status_c == 7, 1, 0),
      
      # Parent is in class 2:
      pc_21 = ifelse(status_p == 2 & status_c == 1, 1, 0),
      pc_22 = ifelse(status_p == 2 & status_c == 2, 1, 0),
      pc_23 = ifelse(status_p == 2 & status_c == 3, 1, 0),
      pc_24 = ifelse(status_p == 2 & status_c == 4, 1, 0),
      pc_25 = ifelse(status_p == 2 & status_c == 5, 1, 0),
      pc_26 = ifelse(status_p == 2 & status_c == 6, 1, 0),
      pc_27 = ifelse(status_p == 2 & status_c == 7, 1, 0),
      
      # Parent is in class 3:
      pc_31 = ifelse(status_p == 3 & status_c == 1, 1, 0),
      pc_32 = ifelse(status_p == 3 & status_c == 2, 1, 0),
      pc_33 = ifelse(status_p == 3 & status_c == 3, 1, 0),
      pc_34 = ifelse(status_p == 3 & status_c == 4, 1, 0),
      pc_35 = ifelse(status_p == 3 & status_c == 5, 1, 0),
      pc_36 = ifelse(status_p == 3 & status_c == 6, 1, 0),
      pc_37 = ifelse(status_p == 3 & status_c == 7, 1, 0),
      
      # Parent is in class 4:
      pc_41 = ifelse(status_p == 4 & status_c == 1, 1, 0),
      pc_42 = ifelse(status_p == 4 & status_c == 2, 1, 0),
      pc_43 = ifelse(status_p == 4 & status_c == 3, 1, 0),
      pc_44 = ifelse(status_p == 4 & status_c == 4, 1, 0),
      pc_45 = ifelse(status_p == 4 & status_c == 5, 1, 0),
      pc_46 = ifelse(status_p == 4 & status_c == 6, 1, 0),
      pc_47 = ifelse(status_p == 4 & status_c == 7, 1, 0),
      
      # Parent is in class 5:
      pc_51 = ifelse(status_p == 5 & status_c == 1, 1, 0),
      pc_52 = ifelse(status_p == 5 & status_c == 2, 1, 0),
      pc_53 = ifelse(status_p == 5 & status_c == 3, 1, 0),
      pc_54 = ifelse(status_p == 5 & status_c == 4, 1, 0),
      pc_55 = ifelse(status_p == 5 & status_c == 5, 1, 0),
      pc_56 = ifelse(status_p == 5 & status_c == 6, 1, 0),
      pc_57 = ifelse(status_p == 5 & status_c == 7, 1, 0),
      
      # Parent is in class 6:
      pc_61 = ifelse(status_p == 6 & status_c == 1, 1, 0),
      pc_62 = ifelse(status_p == 6 & status_c == 2, 1, 0),
      pc_63 = ifelse(status_p == 6 & status_c == 3, 1, 0),
      pc_64 = ifelse(status_p == 6 & status_c == 4, 1, 0),
      pc_65 = ifelse(status_p == 6 & status_c == 5, 1, 0),
      pc_66 = ifelse(status_p == 6 & status_c == 6, 1, 0),
      pc_67 = ifelse(status_p == 6 & status_c == 7, 1, 0),
      
      # Parent is in class 7:
      pc_71 = ifelse(status_p == 7 & status_c == 1, 1, 0),
      pc_72 = ifelse(status_p == 7 & status_c == 2, 1, 0),
      pc_73 = ifelse(status_p == 7 & status_c == 3, 1, 0),
      pc_74 = ifelse(status_p == 7 & status_c == 4, 1, 0),
      pc_75 = ifelse(status_p == 7 & status_c == 5, 1, 0),
      pc_76 = ifelse(status_p == 7 & status_c == 6, 1, 0),
      pc_77 = ifelse(status_p == 7 & status_c == 7, 1, 0)
    )
  
  
  trans_matrix_01_num <- 
    tibble(
      Class1 = c(
        sum(pc_status$pc_11, na.rm = TRUE),
        sum(pc_status$pc_21, na.rm = TRUE),
        sum(pc_status$pc_31, na.rm = TRUE),
        sum(pc_status$pc_41, na.rm = TRUE),
        sum(pc_status$pc_51, na.rm = TRUE),
        sum(pc_status$pc_61, na.rm = TRUE),
        sum(pc_status$pc_71, na.rm = TRUE)
      ),
      Class2 = c(
        sum(pc_status$pc_12, na.rm = TRUE),
        sum(pc_status$pc_22, na.rm = TRUE),
        sum(pc_status$pc_32, na.rm = TRUE),
        sum(pc_status$pc_42, na.rm = TRUE),
        sum(pc_status$pc_52, na.rm = TRUE),
        sum(pc_status$pc_62, na.rm = TRUE),
        sum(pc_status$pc_72, na.rm = TRUE)
      ),
      Class3 = c(
        sum(pc_status$pc_13, na.rm = TRUE),
        sum(pc_status$pc_23, na.rm = TRUE),
        sum(pc_status$pc_33, na.rm = TRUE),
        sum(pc_status$pc_43, na.rm = TRUE),
        sum(pc_status$pc_53, na.rm = TRUE),
        sum(pc_status$pc_63, na.rm = TRUE),
        sum(pc_status$pc_73, na.rm = TRUE)
      ),
      Class4 = c(
        sum(pc_status$pc_14, na.rm = TRUE),
        sum(pc_status$pc_24, na.rm = TRUE),
        sum(pc_status$pc_34, na.rm = TRUE),
        sum(pc_status$pc_44, na.rm = TRUE),
        sum(pc_status$pc_54, na.rm = TRUE),
        sum(pc_status$pc_64, na.rm = TRUE),
        sum(pc_status$pc_74, na.rm = TRUE)
      ),
      Class5 = c(
        sum(pc_status$pc_15, na.rm = TRUE),
        sum(pc_status$pc_25, na.rm = TRUE),
        sum(pc_status$pc_35, na.rm = TRUE),
        sum(pc_status$pc_45, na.rm = TRUE),
        sum(pc_status$pc_55, na.rm = TRUE),
        sum(pc_status$pc_65, na.rm = TRUE),
        sum(pc_status$pc_75, na.rm = TRUE)
      ),
      Class6 = c(
        sum(pc_status$pc_16, na.rm = TRUE),
        sum(pc_status$pc_26, na.rm = TRUE),
        sum(pc_status$pc_36, na.rm = TRUE),
        sum(pc_status$pc_46, na.rm = TRUE),
        sum(pc_status$pc_56, na.rm = TRUE),
        sum(pc_status$pc_66, na.rm = TRUE),
        sum(pc_status$pc_76, na.rm = TRUE)
      ),
      Class7 = c(
        sum(pc_status$pc_17, na.rm = TRUE),
        sum(pc_status$pc_27, na.rm = TRUE),
        sum(pc_status$pc_37, na.rm = TRUE),
        sum(pc_status$pc_47, na.rm = TRUE),
        sum(pc_status$pc_57, na.rm = TRUE),
        sum(pc_status$pc_67, na.rm = TRUE),
        sum(pc_status$pc_77, na.rm = TRUE)
      )
    ) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5 + Class6 + Class7,
           validation = c(data %>% filter(status_p == 1) %>% nrow(),
                          data %>% filter(status_p == 2) %>% nrow(), 
                          data %>% filter(status_p == 3) %>% nrow(), 
                          data %>% filter(status_p == 4) %>% nrow(),
                          data %>% filter(status_p == 5) %>% nrow(),
                          data %>% filter(status_p == 6) %>% nrow(),
                          data %>% filter(status_p == 7) %>% nrow())
    )
  
  rownames(trans_matrix_01_num) <- c("Class1",
                                     "Class2",
                                     "Class3",
                                     "Class4",
                                     "Class5",
                                     "Class6",
                                     "Class7")
  
  
  trans_matrix_01_perc <- 
    as.data.frame(trans_calc(
      trans_matrix_01_num %>%
        dplyr::select(
          Class1,
          Class2,
          Class3,
          Class4,
          Class5,
          Class6,
          Class7))) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5 + Class6 + Class7) 
  
  rownames(trans_matrix_01_perc) <- 
    c("Class1",
      "Class2",
      "Class3",
      "Class4",
      "Class5",
      "Class6",
      "Class7")
  
  
  return(list(trans_matrix_num = trans_matrix_01_num,
              trans_matrix_perc = trans_matrix_01_perc))
}


#=================== Function 4: 6 Class Transition Matrices ===================#

transition_matrix_6class_fun <- function(data){
  pc_status <- 
    data %>% 
    mutate(
      # Parent is in class 1:
      pc_11 = ifelse(status_p == 1 & status_c == 1, 1, 0),
      pc_12 = ifelse(status_p == 1 & status_c == 2, 1, 0),
      pc_13 = ifelse(status_p == 1 & status_c == 3, 1, 0),
      pc_14 = ifelse(status_p == 1 & status_c == 4, 1, 0),
      pc_15 = ifelse(status_p == 1 & status_c == 5, 1, 0),
      pc_16 = ifelse(status_p == 1 & status_c == 6, 1, 0),
      
      # Parent is in class 2:
      pc_21 = ifelse(status_p == 2 & status_c == 1, 1, 0),
      pc_22 = ifelse(status_p == 2 & status_c == 2, 1, 0),
      pc_23 = ifelse(status_p == 2 & status_c == 3, 1, 0),
      pc_24 = ifelse(status_p == 2 & status_c == 4, 1, 0),
      pc_25 = ifelse(status_p == 2 & status_c == 5, 1, 0),
      pc_26 = ifelse(status_p == 2 & status_c == 6, 1, 0),
      
      # Parent is in class 3:
      pc_31 = ifelse(status_p == 3 & status_c == 1, 1, 0),
      pc_32 = ifelse(status_p == 3 & status_c == 2, 1, 0),
      pc_33 = ifelse(status_p == 3 & status_c == 3, 1, 0),
      pc_34 = ifelse(status_p == 3 & status_c == 4, 1, 0),
      pc_35 = ifelse(status_p == 3 & status_c == 5, 1, 0),
      pc_36 = ifelse(status_p == 3 & status_c == 6, 1, 0),
      
      # Parent is in class 4:
      pc_41 = ifelse(status_p == 4 & status_c == 1, 1, 0),
      pc_42 = ifelse(status_p == 4 & status_c == 2, 1, 0),
      pc_43 = ifelse(status_p == 4 & status_c == 3, 1, 0),
      pc_44 = ifelse(status_p == 4 & status_c == 4, 1, 0),
      pc_45 = ifelse(status_p == 4 & status_c == 5, 1, 0),
      pc_46 = ifelse(status_p == 4 & status_c == 6, 1, 0),
      
      # Parent is in class 5:
      pc_51 = ifelse(status_p == 5 & status_c == 1, 1, 0),
      pc_52 = ifelse(status_p == 5 & status_c == 2, 1, 0),
      pc_53 = ifelse(status_p == 5 & status_c == 3, 1, 0),
      pc_54 = ifelse(status_p == 5 & status_c == 4, 1, 0),
      pc_55 = ifelse(status_p == 5 & status_c == 5, 1, 0),
      pc_56 = ifelse(status_p == 5 & status_c == 6, 1, 0),
      
      # Parent is in class 6:
      pc_61 = ifelse(status_p == 6 & status_c == 1, 1, 0),
      pc_62 = ifelse(status_p == 6 & status_c == 2, 1, 0),
      pc_63 = ifelse(status_p == 6 & status_c == 3, 1, 0),
      pc_64 = ifelse(status_p == 6 & status_c == 4, 1, 0),
      pc_65 = ifelse(status_p == 6 & status_c == 5, 1, 0),
      pc_66 = ifelse(status_p == 6 & status_c == 6, 1, 0)
    )
  
  
  trans_matrix_01_num <- 
    tibble(
      Class1 = c(
        sum(pc_status$pc_11, na.rm = TRUE),
        sum(pc_status$pc_21, na.rm = TRUE),
        sum(pc_status$pc_31, na.rm = TRUE),
        sum(pc_status$pc_41, na.rm = TRUE),
        sum(pc_status$pc_51, na.rm = TRUE),
        sum(pc_status$pc_61, na.rm = TRUE)
      ),
      Class2 = c(
        sum(pc_status$pc_12, na.rm = TRUE),
        sum(pc_status$pc_22, na.rm = TRUE),
        sum(pc_status$pc_32, na.rm = TRUE),
        sum(pc_status$pc_42, na.rm = TRUE),
        sum(pc_status$pc_52, na.rm = TRUE),
        sum(pc_status$pc_62, na.rm = TRUE)
      ),
      Class3 = c(
        sum(pc_status$pc_13, na.rm = TRUE),
        sum(pc_status$pc_23, na.rm = TRUE),
        sum(pc_status$pc_33, na.rm = TRUE),
        sum(pc_status$pc_43, na.rm = TRUE),
        sum(pc_status$pc_53, na.rm = TRUE),
        sum(pc_status$pc_63, na.rm = TRUE)
      ),
      Class4 = c(
        sum(pc_status$pc_14, na.rm = TRUE),
        sum(pc_status$pc_24, na.rm = TRUE),
        sum(pc_status$pc_34, na.rm = TRUE),
        sum(pc_status$pc_44, na.rm = TRUE),
        sum(pc_status$pc_54, na.rm = TRUE),
        sum(pc_status$pc_64, na.rm = TRUE)
      ),
      Class5 = c(
        sum(pc_status$pc_15, na.rm = TRUE),
        sum(pc_status$pc_25, na.rm = TRUE),
        sum(pc_status$pc_35, na.rm = TRUE),
        sum(pc_status$pc_45, na.rm = TRUE),
        sum(pc_status$pc_55, na.rm = TRUE),
        sum(pc_status$pc_65, na.rm = TRUE)
      ),
      Class6 = c(
        sum(pc_status$pc_16, na.rm = TRUE),
        sum(pc_status$pc_26, na.rm = TRUE),
        sum(pc_status$pc_36, na.rm = TRUE),
        sum(pc_status$pc_46, na.rm = TRUE),
        sum(pc_status$pc_56, na.rm = TRUE),
        sum(pc_status$pc_66, na.rm = TRUE)
      )
    ) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5 + Class6,
           validation = c(data %>% filter(status_p == 1) %>% nrow(),
                          data %>% filter(status_p == 2) %>% nrow(), 
                          data %>% filter(status_p == 3) %>% nrow(), 
                          data %>% filter(status_p == 4) %>% nrow(),
                          data %>% filter(status_p == 5) %>% nrow(),
                          data %>% filter(status_p == 6) %>% nrow())
    )
  
  rownames(trans_matrix_01_num) <- c("Class1",
                                     "Class2",
                                     "Class3",
                                     "Class4",
                                     "Class5",
                                     "Class6")
  
  
  trans_matrix_01_perc <- 
    as.data.frame(trans_calc(
      trans_matrix_01_num %>%
        dplyr::select(
          Class1,
          Class2,
          Class3,
          Class4,
          Class5,
          Class6))) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4 + Class5 + Class6) 
  
  rownames(trans_matrix_01_perc) <- 
    c("Class1",
      "Class2",
      "Class3",
      "Class4",
      "Class5",
      "Class6")
  
  
  return(list(trans_matrix_num = trans_matrix_01_num,
              trans_matrix_perc = trans_matrix_01_perc))
}


#=================== Function 4: 5 Class Transition Matrices ===================#

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


#=================== Function 4: 4 Class Transition Matrices ==================#
transition_matrix_4class_fun <- function(data){
  pc_status <- 
    data %>% 
    mutate(
      # Parent is in class 1:
      pc_11 = ifelse(status_p == 1 & status_c == 1, 1, 0),
      pc_12 = ifelse(status_p == 1 & status_c == 2, 1, 0),
      pc_13 = ifelse(status_p == 1 & status_c == 3, 1, 0),
      pc_14 = ifelse(status_p == 1 & status_c == 4, 1, 0),
      
      # Parent is in class 2:
      pc_21 = ifelse(status_p == 2 & status_c == 1, 1, 0),
      pc_22 = ifelse(status_p == 2 & status_c == 2, 1, 0),
      pc_23 = ifelse(status_p == 2 & status_c == 3, 1, 0),
      pc_24 = ifelse(status_p == 2 & status_c == 4, 1, 0),
      
      # Parent is in class 3:
      pc_31 = ifelse(status_p == 3 & status_c == 1, 1, 0),
      pc_32 = ifelse(status_p == 3 & status_c == 2, 1, 0),
      pc_33 = ifelse(status_p == 3 & status_c == 3, 1, 0),
      pc_34 = ifelse(status_p == 3 & status_c == 4, 1, 0),
      
      # Parent is in class 4:
      pc_41 = ifelse(status_p == 4 & status_c == 1, 1, 0),
      pc_42 = ifelse(status_p == 4 & status_c == 2, 1, 0),
      pc_43 = ifelse(status_p == 4 & status_c == 3, 1, 0),
      pc_44 = ifelse(status_p == 4 & status_c == 4, 1, 0)
    )
  
  
  trans_matrix_01_num <- 
    tibble(
      Class1 = c(
        sum(pc_status$pc_11, na.rm = TRUE),
        sum(pc_status$pc_21, na.rm = TRUE),
        sum(pc_status$pc_31, na.rm = TRUE),
        sum(pc_status$pc_41, na.rm = TRUE)
      ),
      Class2 = c(
        sum(pc_status$pc_12, na.rm = TRUE),
        sum(pc_status$pc_22, na.rm = TRUE),
        sum(pc_status$pc_32, na.rm = TRUE),
        sum(pc_status$pc_42, na.rm = TRUE)
      ),
      Class3 = c(
        sum(pc_status$pc_13, na.rm = TRUE),
        sum(pc_status$pc_23, na.rm = TRUE),
        sum(pc_status$pc_33, na.rm = TRUE),
        sum(pc_status$pc_43, na.rm = TRUE)
      ),
      Class4 = c(
        sum(pc_status$pc_14, na.rm = TRUE),
        sum(pc_status$pc_24, na.rm = TRUE),
        sum(pc_status$pc_34, na.rm = TRUE),
        sum(pc_status$pc_44, na.rm = TRUE)
      )) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4,
           validation = c(
             data %>% filter(status_p == 1) %>% nrow(),
             data %>% filter(status_p == 2) %>% nrow(), 
             data %>% filter(status_p == 3) %>% nrow(), 
             data %>% filter(status_p == 4) %>% nrow()
           ))
  
  rownames(trans_matrix_01_num) <- c("Class1",
                                     "Class2",
                                     "Class3",
                                     "Class4")
  
  
  trans_matrix_01_perc <- 
    as.data.frame(trans_calc(
      trans_matrix_01_num %>%
        dplyr::select(
          Class1,
          Class2,
          Class3,
          Class4))) %>%
    mutate(rowsum = Class1 + Class2 + Class3 + Class4) 
  
  rownames(trans_matrix_01_perc) <- 
    c("Class1",
      "Class2",
      "Class3",
      "Class4")
  
  
  return(list(trans_matrix_num = trans_matrix_01_num,
              trans_matrix_perc = trans_matrix_01_perc))
}


#=================== Function 4: 3 Class Transition Matrices ==================#
transition_matrix_3class_fun <- function(data){
  pc_status <- 
    data %>% 
    mutate(
      # Parent is in class 1:
      pc_11 = ifelse(status_p == 1 & status_c == 1, 1, 0),
      pc_12 = ifelse(status_p == 1 & status_c == 2, 1, 0),
      pc_13 = ifelse(status_p == 1 & status_c == 3, 1, 0),
      
      # Parent is in class 2:
      pc_21 = ifelse(status_p == 2 & status_c == 1, 1, 0),
      pc_22 = ifelse(status_p == 2 & status_c == 2, 1, 0),
      pc_23 = ifelse(status_p == 2 & status_c == 3, 1, 0),
      
      # Parent is in class 3:
      pc_31 = ifelse(status_p == 3 & status_c == 1, 1, 0),
      pc_32 = ifelse(status_p == 3 & status_c == 2, 1, 0),
      pc_33 = ifelse(status_p == 3 & status_c == 3, 1, 0)
    )
  
  
  trans_matrix_01_num <- 
    tibble(
      Class1 = c(
        sum(pc_status$pc_11, na.rm = TRUE),
        sum(pc_status$pc_21, na.rm = TRUE),
        sum(pc_status$pc_31, na.rm = TRUE)
      ),
      Class2 = c(
        sum(pc_status$pc_12, na.rm = TRUE),
        sum(pc_status$pc_22, na.rm = TRUE),
        sum(pc_status$pc_32, na.rm = TRUE)
      ),
      Class3 = c(
        sum(pc_status$pc_13, na.rm = TRUE),
        sum(pc_status$pc_23, na.rm = TRUE),
        sum(pc_status$pc_33, na.rm = TRUE)
      )
    ) %>%
    mutate(
      rowsum = Class1 + Class2 + Class3,
      validation = c(data %>% filter(status_p == 1) %>% nrow(),
                     data %>% filter(status_p == 2) %>% nrow(), 
                     data %>% filter(status_p == 3) %>% nrow())
    )
  
  rownames(trans_matrix_01_num) <- c("Class1",
                                     "Class2",
                                     "Class3")
  
  
  trans_matrix_01_perc <- 
    as.data.frame(trans_calc(
      trans_matrix_01_num %>%
        dplyr::select(
          Class1,
          Class2,
          Class3))) %>%
    mutate(rowsum = Class1 + Class2 + Class3) 
  
  rownames(trans_matrix_01_perc) <- 
    c("Class1",
      "Class2",
      "Class3")
  
  
  return(list(trans_matrix_num = trans_matrix_01_num,
              trans_matrix_perc = trans_matrix_01_perc))
}


###### Basic Operation:
calculate_percentage <- function(row) {
  row_percentage <- row / sum(row) 
  return(row_percentage)
}


trans_calc <- function(data){
  trans_matrix <- t(apply(data[, 1:ncol(data)], 1, calculate_percentage))
  return(trans_matrix)
}



#=================== Function 5: Generate Fitted TMs ==================#
fitted_TMs <- function(model, new_data, training_data){
  predictions <- predict(model, new_data, type = "response")
  round_predictions<- round(predictions)
  prediction_data <- 
    cbind(
      round_predictions,
      predictions, 
      training_data)
  return(prediction_data)
}


#=================== Function 6: Generate Regression Data ==================#
generate_model_data <- function(fitted, model, tm_df_org_lst, size, start_year, end_year){
  # 1. Get the regression data for comparison:
  regression_data <-bind_rows(lapply(tm_df_org_lst, tm_to_df))
  regression_data <-
    regression_data %>%
    mutate(
      cohort = rep(seq(start_year,end_year,1),each = size*size)
    ) %>%
    mutate(
      parent_children = as.factor(paste0(parent_status,"x", children_status)),
      parent_cohort = as.factor(paste0(parent_status,"x", cohort)),
      children_cohort = as.factor(paste0(children_status,"x", cohort))
    ) %>%
    mutate(
      parent_status = as.factor(parent_status),
      children_status = as.factor(children_status)
    )
  # 2. Generate the prediction data:
  year_length <- end_year - start_year + 1
  new_data <-
    data.frame(
      parent_status = rep(rep(seq(1,size,1),each=size),year_length),
      children_status = rep(rep(seq(1,size,1),size),year_length),
      cohort = rep(seq(start_year,end_year,1),each = size * size)
    ) %>%
    mutate(
      parent_children = as.factor(paste0(parent_status,"x", children_status)),
      parent_cohort = as.factor(paste0(parent_status,"x", cohort)),
      children_cohort = as.factor(paste0(children_status,"x", cohort))
    ) %>%
    mutate(
      parent_status = as.factor(parent_status),
      children_status = as.factor(children_status),
      cohort_f = as.factor(cohort)
    )
  
  if(fitted == TRUE){
    fitted_data <-
      fitted_TMs(
        model,
        new_data,
        regression_data
      )
    return(fitted_data)
  }else{
    return(regression_data)
  }
}

#=================== Function 7: Summary the model performance ==================#
summary_model_performance <- function(sp_param, basis, k_num, regression_data){
  # Section 1: Fit the model:
  control_params <- gam.control(trace = FALSE)
  model <-
    gam(
      cell_count ~
        parent_children +
        s(cohort,
          by = interaction(parent_status, children_status),
          bs = basis,
          k = k_num,
          sp = sp_param),
      family = poisson(link="log"),
      select=FALSE, # Keep this as FALSE, avoiding over-smoothing
      method = "REML",
      data = regression_data,
      control = control_params,
      gc.level = 2,
      nthreads = 3)
  # Section 2: Summary the model statistics
  summary_model <- summary(model)
  # Step 1: Grab the R Squared：
  r_sq <- summary_model$r.sq
  # Step 2: Grab the AIC:
  aic <- model$aic
  # Step 3: Grab the edf table:
  edf_table <-
    as.data.frame(summary_model$s.table) %>%
    dplyr::select(
      edf,
      Ref.df
    )
  rownames(edf_table) <- gsub(
    "^s\\(cohort\\):interaction\\(parent_status, children_status\\)(\\d)\\.(\\d)",
    "Parent Class: \\1 x Children Class: \\2",
    rownames(edf_table)
  )
  model_summary_lst <-
    list(
      AIC = aic,
      `R Sqaured` = r_sq,
      `EDF Table` = edf_table
    )
  return(list(model_summary = model_summary_lst, model = model))
}

#=================== Function 8: Generate Param Table for Param Sweep ==================#
generate_param_table <- function(k_vec,s_vec){
  # Generate the param_table:
  param_table <-
    tibble(
      sp_param = sp_vec,
      k = k_vec
    )
  # Cut it to row:
  param_lists <- split(param_table, seq(nrow(param_table)))
  return(param_lists)
}

#=================== Function 9: Model 5 Param Sweep ==================#
model5_param_sweep <- function(param_lists,regression_data){
  model_lst <-
    lapply(
      param_lists,
      function(param_lst){
        model_info <-
          summary_model_performance(
            param_lst$sp_param, 
            "tp",
            param_lst$k, 
            regression_data) 
      }
    )
  return(model_lst)
}

#=================== Function 10: Fetch model Results ==================#
model_result <- function(model_lst,tm_df_org_lst, size, start_year, end_year){
  result_list <- list()
  # Grab all the models:
  result_list[["model"]] <- 
    lapply(model_lst, function(x) x[["model"]])
  # Grab the model performance:
  result_list[["R2"]] <-
    lapply(model_lst, function(x) x[["model_summary"]][["R Sqaured"]])
  # Generate the fitted data:
  result_list[["fitted_data"]] <-
    lapply(
      result_list[["model"]],
      function(model){
        fitted_data <-
          generate_model_data(
            TRUE, 
            model, 
            tm_df_org_lst, 
            size, 
            start_year,
            end_year
          )
      }
    )
  return(result_list)
}


