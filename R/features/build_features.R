library(dplyr)
library(ggplot2)

make_avg_cpd <- function(df) {
  # looking for cigarettes per day as an average.
  # if quit --> 0 cpd; if missing cpd, use weekly average
  # if missing weekly average, use monthly average
  
  dpw <- 7
  dpm <- 30
  
  
  df %>%
    mutate(avg_cpd = case_when(
      quit == "YES" ~ as.numeric(0),
      !is.na(cpd) ~ as.numeric(cpd),
      is.na(cpd) & !is.na(cpw) ~ cpw/dpw,
      is.na(cpd) & is.na(cpw) & !is.na(cpm) ~ cpm/dpm
    ))
}

smoking_trajectories <- function(df, visit_levels, visit_recode) {
  
  df <- df %>%
    filter(week %in% visit_levels) %>%
    make_avg_cpd()
  
  # baseline cpd
  bl <- df %>%
    filter(week == "baseline") %>%
    mutate(bl_cpd = avg_cpd) %>%
    select(subject_id, bl_cpd)
  
  # construct proportion change
  traj <- df %>%
    filter(week != "baseline", !is.na(avg_cpd)) %>%
    left_join(bl, by = "subject_id") %>%
    mutate(prop_change = case_when(
      avg_cpd == 0 ~ -1,
      TRUE ~ avg_cpd / bl_cpd - 1
    )) %>%
    mutate(week = recode(as.character(week), !!!visit_recode)) %>%
    mutate(week = as.numeric(week)) %>%
    select(subject_id, week, avg_cpd, bl_cpd, prop_change)
  
  traj
}

rm_no_cpd <- function(df) {
  
  df_sub <- df %>%
    filter(visit %in% visit_levels)
  
  # only keep subjects with baseline cpd and at least one AVG cpd value
  has_bl <- df_sub %>%
    filter(visit == "baseline" & !is.na(cpd)) %>% 
    select(subject_id) %>%
    left_join(df_sub, by = "subject_id")
  
  has_bl %>%
    make_avg_cpd() %>%
    filter(visit != "baseline" & !is.na(avg_cpd)) %>%
    select(subject_id) %>%
    distinct() %>%
    left_join(has_bl, by = "subject_id") %>%
    rename("week" = visit)
}


rm_subj_var <- function(df) {
  miss_count <- df %>%
    filter(week == "baseline") %>%
    select(-c(cpw, cpm)) %>%
    mutate(across(-subject_id, ~ifelse(is.na(.), 1, 0))) %>%
    tidyr::pivot_longer(-subject_id) %>%
    group_by(subject_id) %>%
    mutate(subj_int = cur_group_id()) %>% 
    summarize(n_miss = sum(value)) %>%
    ungroup()

  rm_subj <- miss_count %>% filter(n_miss >= 7) %>% distinct(subject_id)

  df %>% filter(!subject_id %in% rm_subj$subject_id)
}


rm_var <- function(df) {
  df %>%
    select(-c(sf_general, sf_physfunc))
}

impute <- function(df_lowmiss_bl) {
  find_mode <- function(x) {
    val <- unique(x[!is.na(x)])
    mode_val <- val[which.max(tabulate(match(x, val)))]
    mode_val
  }
  

  # impute by study site, sex, and age blocks.
  # bin age into 2 groups --+ 17.9 to 51, 51 to 84.1
  # any fewer will cause 1 or 2 people per imputation box
  
  df_lowmiss_bl %>%
    mutate(age_bin = cut(age, breaks = 2)) %>%
    group_by(site, sex, age_bin) %>%
    mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    mutate(across(where(is.factor), ~dplyr::if_else(is.na(.), find_mode(.), .))) %>%
    ungroup() %>%
    select(-age_bin)

}


subset_bl <- function(df, traj) {

  # add end CPD for each subject
  w26 <- traj %>% 
    filter(week == 26) %>%
    select(subject_id, "week26_cpd" = avg_cpd)
  
  df %>%
    filter(visit == "baseline") %>%
    select(-c(visit_date, visit, cpw, cpm)) %>%
    left_join(w26, by = "subject_id")
}

subset_quit <- function(df, df_bl) {
  df_quitvars <- df %>%
    filter(visit == "month 12") %>%
    # mutate(quit = NA) %>%
    make_avg_cpd() %>%
    select(
      subject_id,
      "co_1year" = co,
      "avg_cpd_1year" = avg_cpd,
      quit
    )
  
  df_bl %>%
    select(-quit) %>%
    left_join(df_quitvars, by = "subject_id")
}



visit_levels <- c("baseline", "week 2", "week 10", "month 4", "month 6")
visit_recode <-
  c(
    "baseline" = "baseline",
    "week 2" = 2,
    "week 10" = 10,
    "month 4" = 18,
    "month 6" = 26
  )

# collapse categories with fewer than ~ 1% of the dataset

collapse_counts <- function(df) {
  
  # handle ts last quit attempt for never-quittters
  
  df <- df %>%
    mutate(
      # n_quit_attempts = recode(n_quit_attempts, "NEVER" = "ONCE"),
      anxiety = recode(anxiety, "EXTREMELY SO" = "VERY MUCH SO"),
      depression = recode(depression, "EXTREMELY SO" = "VERY MUCH SO")
    )
  
  return(df)
}


construct_ordinal <- function(df) {
  non_vars <- c("site", "visit", "trt_recode", 'trt_grp', "sex", "quit")
  zero_vars <- c("n_quit_attempts", "ts_last_quit_attempt", "anxiety", "depression", "int_to_quit")
  
  df %>%
    mutate(across((where(is.factor) & !non_vars), as.numeric)) %>%
    mutate(across(all_of(zero_vars), ~. - 1)) # certain ordinal variables start at 0
}


features <- function(df, nme) {
  
  df_lowmiss <- 
    rm_no_cpd(df) %>% 
    rm_subj_var() %>%
    rm_var()
  
  if (nme == "cat") {
    df_lowmiss <- df_lowmiss %>% collapse_counts()
  }
  
  traj <-  smoking_trajectories(df_lowmiss, visit_levels, visit_recode)
  
  df_lowmiss_bl <- df_lowmiss %>% filter(week == "baseline")
  
  df_imputed <- df_lowmiss_bl %>% 
    impute() %>% 
    make_avg_cpd() %>%
    select(-c(cpd, cpw, cpm))

  df_nomiss <- df_imputed %>% tidyr::drop_na()
  
  pred_imputed <- subset_quit(df, df_imputed)
  pred_nomiss <- subset_quit(df, df_nomiss)
  
  
  save(traj, file = "../data/processed/smoking_traj.rda")
  save(pred_imputed, file = paste0("../data/processed/pred_imputed", nme, ".rda"))
  save(pred_nomiss, file = paste0("../data/processed/pred_nomiss", nme, ".rda"))
}


# main --- 

load("../data/interim/compiled_data.rda")


# categorical
features(df, nme = "cat")

# ordinal
df_ord <- construct_ordinal(df)
features(df_ord, nme = "ord")
  