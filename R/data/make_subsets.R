find_missing_bl_cpd <- function(df) {
  df %>%
    filter(visit == "baseline" & is.na(cpd))
}

count_missing_baseline <- function(df) {
  df %>%
    select(-c(quit, cpm, cpw)) %>%
    filter(visit == "baseline") %>%
    mutate(across(-subject_id, ~ifelse(is.na(.), 1, 0))) %>%
    tidyr::pivot_longer(-subject_id) %>%
    group_by(subject_id) %>%
    mutate(n_miss = sum(value))
}

filter_overall <- function(df, n_miss_thresh = 5) {
  
  # make subset with only subjects missing baseline CPD values removed
  # for LCA
  df_base <- df %>%
    filter(!subject_id %in% find_missing_bl_cpd(df)$subject_id)
  
  # remove largely-missing variables and subjects to make
  # "low-missing" data set -- will be used for imputation
  
  # variable level ----- 
  df_sub <- df_base %>% select(-c(sf_general, sf_physfunc))
  
  # subject level ------
  rm_sub <- df_sub %>% 
    count_missing_baseline() %>%
    filter(n_miss >= n_miss_thresh)
  
  df_lowmiss <- df_sub %>% filter(!subject_id %in% rm_sub$subject_id)
  
  # remove subjects missing any baseline variable
  rm_sub_all <- df_sub %>%
    count_missing_baseline() %>% 
    filter(n_miss > 0)
  
  df_nomiss <- df_sub %>% filter(!subject_id %in% rm_sub_all$subject_id)
  
  list(
    "df_base" = df_base,
    "df_lowmiss" = df_lowmiss,
    "df_nomiss" = df_nomiss
  )
}

filter_quit <- function(df) {
  df %>%
    filter(visit == "month 12") %>%
    select(subject_id, co, cpd, cpw, cpm)
}


write_subsets <- function(df) {
  data <- filter_overall(df)
  df_base <- data$df_base
  df_lowmiss <- data$df_lowmiss
  df_nomiss <- data$df_nomiss
  

  # baseline data sets -----
  # baseline variables for all subjects in LCA subset
  df_bl <- df_base %>% filter(visit == "baseline")
  save(df_bl, file = "../data/interim/df_bl.rda")
  
  # subjects with 5 or more missing baseline variables removed
  df_lowmiss_bl <- df_lowmiss %>% filter(visit == "baseline")
  save(df_lowmiss_bl, file = "../data/interim/df_lowmiss_bl.rda")
  
  # subjects with any missing baseline data removed
  df_nomiss_bl <- df_nomiss %>% filter(visit == "baseline")
  save(df_nomiss_bl, file = "../data/interim/df_nomiss_bl.rda")
  
  
  # LCA data set ----
  
  # The max number of subjects we'll use is the number in df_lowmiss_bl
  
  visit_levels <- c("baseline", "week 2", "week 10", "month 4", "month 6")
  visit_recode <-
    c(
      "baseline" = "baseline",
      "week 2" = 2,
      "week 10" = 10,
      "month 4" = 18,
      "month 6" = 26
    )
  
  df_lca <- df_base %>% 
    filter(subject_id %in% df_lowmiss_bl$subject_id) %>%
    filter(visit %in% visit_levels) %>%
    mutate(week = recode(as.character(visit), !!!visit_recode)) %>%
    select(subject_id, quit, week, cpd, cpw, cpm)
  
  save(df_lca, file = "../data/interim/df_lca.rda")
  
  
  # 1-year follow-up
  df_quit <- filter_quit(df)
  save(df_quit, file = "../data/interim/df_quit.rda")
 
}


  
load("../data/interim/compiled_data.rda")
write_subsets(df)

