# McNeil LCA Project Data Management
# Tony Barrows
# Updated 2022-09-30

# all paths relative to ./R

library(dplyr)

# functions ----

load_raw_data <- function(path = "../data/raw/") {
  # Create one large vector for all trial data filled with vectors
  # for each trial data set. Name each with file names.
  
  dirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  df_vector <- list()
  
  for (i in 1:length(dirs)) {
    files <- list.files(dirs[[i]], pattern = ".txt", full.names = TRUE)
    tmp_lst <- list()
    for (j in 1:length(files)) {
      # read each data file, add blank fields if misaligned
      tmp_lst[[j]] <-
        read.table(
          files[[j]],
          sep = "\t",
          header = TRUE,
          quote = "",
          fill = TRUE
        )
      names(tmp_lst[[j]]) <- toupper(names(tmp_lst[[j]]))
    }
    names(tmp_lst) <- sub(".*/", "", files)
    df_vector[[i]] <- tmp_lst
    print(paste(i, "successfully loaded"))
  }
  names(df_vector) <- sub(".*/", "", dirs)
  df_vector
}

assign_id <- function(df_vector_sub) {
  # subject_id = the first and last 2 characters of STUDYID plus the subject_id
  
  for (i in 1:length(df_vector_sub)) {
    for (j in 1:length(df_vector_sub[[i]])) {
      study_id <- df_vector_sub[[i]][[j]]$STUDYID[[1]]
      c1 <- substr(study_id, 1, 2)
      c2 <- substr(study_id, nchar(study_id) - 1, nchar(study_id))
      
      tryCatch(
        df_vector_sub[[i]][[j]] <-
          df_vector_sub[[i]][[j]] %>%
          mutate(
            subject_id = paste(c1, c2, sep = ""),
            subject_id = as.numeric(paste(subject_id, SUBJID, sep = ""))
          ) %>%
          select(subject_id, everything()),
        error = function(e) {
          print(paste(
            "no_study_id:",
            names(df_vector_sub[i]),
            ":",
            names(df_vector_sub[[i]][j])
          ))
        }
      )
    }
  }
  df_vector_sub
}

subject_master_list <- function(df_vector_sub) {
  # "**Important** Only `SUBJID`s contained in the `DEM.txt`"
  # "files were included in the trial. Return master list."
  
  sid_df <- data.frame(subject_id = NA,
                       SUBJID = NA,
                       STUDYID = NA)
  
  for (i in 1:length(df_vector_sub)) {
    tmp <- df_vector_sub[[i]]$DEM.txt %>%
      select(subject_id, SUBJID, STUDYID)
    sid_df <- rbind(sid_df, tmp)
  }
  
  sid_df <-
    sid_df %>% distinct(.keep_all = TRUE) %>% tidyr::drop_na()
  
}

clean_trials <- function(df_vector) {
  # remove failed trial
  df_vector_sub <- df_vector[-2]
  
  # not data
  for (i in 1:length(df_vector_sub)) {
    df_vector_sub[[i]][which(startsWith(names(df_vector_sub[[i]]), "format"))] <-
      NULL
  }
  
  # parallelize
  for (i in 1:length(df_vector_sub)) {
    df_vector_sub[[i]]$QOLDATA.txt$STUDYID <-
      df_vector_sub[[i]]$AE.txt$STUDYID[1]
    df_vector_sub[[i]]$DEM.txt$STUDYID <-
      df_vector_sub[[i]]$AE.txt$STUDYID[1]
  }
  
  # assign subject IDs
  df_vector_sub <- assign_id(df_vector_sub)
  
  return(df_vector_sub)
}

rename_identical <- function(df) {
  # renaming function
  
  df %>%
    plyr::rename(c("ANXA1" = "ANXA",
                   "DEPR2" = "DEPR2X"),
                 warn_missing = FALSE)
}

recode_variables <- function(df) {
  # recode categorical
  
  df <-  df %>%
    mutate(
      ANXA1N = recode(
        ANXA1N,
        "0" = "NOT AT ALL",
        "1" = "SOMEWHAT",
        "2" = "MODERATELY SO",
        "3" = "VERY MUCH SO",
        "4" = "EXTREMELY SO",
        .missing = NA_character_
      ),
      PEPUP = cut(
        PEPUP,
        breaks = c(0, 2, 4, 6, 8, 10),
        labels = c("NOT AT ALL",
                   "MILD",
                   "MODERATE",
                   "STRONG",
                   "VERY STRONG")
      ),
      CALMEFFT = cut(
        CALMEFFT,
        breaks = c(0, 2, 4, 6, 8, 10),
        labels = c("NOT AT ALL",
                   "MILD",
                   "MODERATE",
                   "STRONG",
                   "VERY STRONG")
      )
    )
  
  # recode intention to quit
  
  quit_labels <- c("NOT AT ALL",
                   "A LITTLE",
                   "SOMEWHAT ",
                   "A LOT")
  
  df <- df %>%
    mutate(STOPSMK = factor(STOPSMK)) %>%
    mutate(
      QTNX6MO = cut(
        QTNX6MO,
        breaks = c(0, 2.5, 5, 7.5, 10),
        labels = quit_labels
      ),
      INTMOT = cut(
        INTMOT,
        breaks = c(0, 2.5, 5, 7.5, 10),
        labels = quit_labels
      )
    )
  
  # treatment group
  df <- df %>%
    mutate(
      trt_recode =
        ifelse(
          stringr::str_detect(tolower(TRTGRP), "active"),
          "active",
          "placebo"
        ),
      trt_recode = as.factor(trt_recode)
    )
  
  
  # unite across
  df %>%
    mutate(
      ANXA = ifelse(is.na(ANXA), ANXA1N, ANXA),
      CLMEFTA = ifelse(is.na(CLMEFTA), CALMEFFT, CLMEFTA),
      PEPUP = ifelse(is.na(PEPUP), PEPUPA, PEPUP),
      QTNX6MO = ifelse(is.na(QTNX6MO), STOPSMK, QTNX6MO),
      QTNX6MO = ifelse(is.na(QTNX6MO), INTMOT, QTNX6MO)
    ) %>%
    mutate(QTNX6MO = factor(QTNX6MO, labels = quit_labels))
}

rename_timepoints <- function(df) {
  df %>%
    mutate(VISIT = tolower(VISIT)) %>%
    mutate(VISIT = ifelse(VISIT == "week 0", "baseline", VISIT)) %>%
    mutate(VISIT = ifelse(VISIT == "visit 1", "baseline", VISIT)) %>%
    mutate(VISIT = as.factor(VISIT))
}


select_vars <- function(df_vector_sub, n, sid_df, site_name) {
  # need smkqty, smkqtyw, nucores, smkage, smknon, smkqt,
  # ftnscore, intmot, qtnx6mo, stopsmk, wantsmk, smkqtla,
  # anxa, depr, emowell, pepup, calmefft, lstcigex, treatment group
  
  vec <- df_vector_sub[[n]]
  names(vec) <- tolower(names(vec))
  names(vec) <- stringr::str_remove_all(names(vec), ".txt")
  
  # variables not in ger, den, usa
  if (n %in%  c(2, 3, 5)) {
    vec$craveq$CLMEFTA <- NA
    vec$wsq$ANXA1N <- NA
    vec$craveq$PEPUPA <- NA
    vec$ftn$FTNQTY <- NA
  }
  
  # variable only in swi
  if (n != 1) {
    vec$sminchk$SMKQTYM <- NA
    vec$intquit$STOPSMK <- NA
  }
  
  # align variables
  if (n == 1) {
    vec$craveq$CALMEFFT <- NA
    vec$craveq$PEPUP <- NA
  }
  
  if (n == 4) {
    vec$wsq$ANXA1N <- NA
  }
  
  if (n %in% c(2, 3, 4)) {
    vec$smkhist$CIGCQTY <- NA
  }
  
  if (n %in% c(1, 3)) {
    vec$intquit$QTNX6MO <- NA
  }
  
  if (n != 3) {
    vec$intquit$INTMOT <- NA
  }
  
  visit <- vec$visit %>% select(subject_id, VISIT, VISDTF)
  ftn <- vec$ftn %>% select(subject_id, VISIT, FTNSCORE)
  nu <- vec$nu %>% select(subject_id, VISIT, NUCORES)
  scq <-
    vec$scq %>% select(subject_id,
                       VISIT,
                       PHYSFUNC,
                       PHYSHEAL,
                       EMOPROB,
                       SOCFUNC,
                       PAIN,
                       GENERAL,
                       EMOWELL)
  sminchk <-
    vec$sminchk %>% select(subject_id, VISIT, SMKQTY, SMKQTYW, SMKQTYM, SMKSTP)
  smkhist <-
    vec$smkhist %>% select(subject_id, VISIT, SMKNON, SMKAGE, SMKQT, SMKQTLA, CIGCQTY)
  intquit <-
    vec$intquit %>% select(subject_id, VISIT, STOPSMK, QTNX6MO, INTMOT)
  
  if (n == 4) {
    wsq <-
      vec$wsq %>% rename_identical() %>% select(subject_id,
                                                VISIT,
                                                ANXA,
                                                ANXA1N,
                                                DEPR2X,
                                                PEPUP,
                                                CALMEFFT,
                                                LSTCGEX)
  } else {
    wsq <-
      vec$wsq %>% rename_identical() %>% select(subject_id, VISIT, ANXA, ANXA1N, DEPR2X)
  }
  
  dem <- vec$dem %>% select(subject_id, SEX, AGE)
  trt <- vec$trt %>% select(subject_id, TRTGRP)
  
  # join
  join_vars <- c("subject_id", "VISIT")
  df <- visit %>% rename_timepoints() %>%
    left_join(ftn %>% rename_timepoints(), by = join_vars) %>%
    left_join(nu %>% rename_timepoints(), by = join_vars) %>%
    left_join(scq %>% rename_timepoints(), by = join_vars) %>%
    left_join(sminchk %>% rename_timepoints(), by = join_vars) %>%
    left_join(smkhist %>% rename_timepoints(), by = join_vars) %>%
    left_join(wsq %>% rename_timepoints(), by = join_vars) %>%
    left_join(intquit %>% rename_timepoints(), by = join_vars) %>%
    left_join(dem, by = "subject_id") %>%
    left_join(trt, by = "subject_id")
  
  # file not in aus
  if (n != 4) {
    craveq <- vec$craveq %>%
      select(subject_id, VISIT, CLMEFTA, CALMEFFT, LSTCGEX, PEPUP, PEPUPA)
    df <-
      df %>% left_join(craveq %>% rename_timepoints(), by = join_vars)
  }
  
  fill_vars <- c(names(dem), names(smkhist), names(trt))
  fill_vars <- fill_vars[fill_vars != "subject_id"]
  
  df <- df %>%
    mutate(site = site_name) %>%
    select(subject_id, site, VISIT, TRTGRP, SEX, AGE, everything()) %>%
    filter(subject_id %in% sid_df$subject_id) %>%
    group_by(subject_id) %>%
    tidyr::fill(all_of(fill_vars)) %>%
    ungroup()
  
  
  df[df == ""] <- NA
  
  df
}

rename_recode <- function(df) {
  visit_levels <-
    c(
      "baseline",
      "week 1",
      "week 2",
      "week 6",
      "week 10",
      "month 3",
      "month 4",
      "month 6",
      "month 9",
      "month 12",
      "month 15",
      "month 18",
      "month 24"
    )
  nme_vec <- c(
    "TRTGRP" = "trt_grp",
    "VISIT" = "visit",
    "SEX" = "sex",
    "AGE" = "age",
    "VISDTF" = "visit_date",
    "FTNSCORE" = "ftnd",
    "NUCORES" = "co",
    "PHYSFUNC" = "sf_physfunc",
    "PHYSHEAL" = "sf_physheal",
    "EMOPROB" = "sf_emoprob",
    "SOCFUNC" = "sf_socfunc",
    "PAIN" = "sf_pain",
    "GENERAL" = "sf_general",
    "EMOWELL" = "sf_emowell",
    "SMKQTY" = "cpd",
    "SMKQTYW" = "cpw",
    "SMKQTYM" = "cpm",
    "SMKSTP" = "quit",
    "SMKQT" = "n_quit_attempts",
    "SMKNON" = "longest_period_wo_smoking",
    "SMKAGE" = "age_started_smoking",
    "SMKQTLA" = "ts_last_quit_attempt",
    "ANXA" = "anxiety",
    "DEPR2X" = "depression",
    "QTNX6MO" = "int_to_quit",
    "CLMEFTA" = "rsq_calming",
    "LSTCGEX" = "rsq_last_cig_exp",
    "PEPUP" = "rsq_pepping_up_eff"
  )
  
  df %>%
    plyr::rename(nme_vec) %>%
    mutate(visit = factor(visit, levels = visit_levels))
}

reorder_factors <- function(df) {
  trt_levels <- c("placebo", "active")
  quit_levels <-
    c("NEVER",
      "ONCE",
      "2 TO 5 TIMES",
      "6 TO 10 TIMES",
      "MORE THAN 10 TIMES")
  ts_levels <-
    c("never_quit",
      ">12 MONTHS",
      ">12-24 MONTHS",
      ">6-12 MONTHS",
      "0-6 MONTHS")
  anx_dep_levels <-
    c("NOT AT ALL",
      "SOMEWHAT",
      "MODERATELY SO",
      "VERY MUCH SO",
      "EXTREMELY SO")
  int_quit_levels <-
    c("NOT AT ALL", "A LITTLE", "SOMEWHAT ", "A LOT")
  lst_cig_levels <-
    c(
      "VERY UNPLEASANT",
      "SOMEWHAT UNPLEASANT",
      "NEUTRAL",
      "SOMEWHAT PLEASANT",
      "VERY PLEASANT"
    )
  
  longest_period <-
    c(
      "1" = "<1 week",
      "2" = "1week-1month",
      "3" = ">1month-3months",
      "4" = ">3months"
    )
  
  df %>%
    mutate(
      longest_period_wo_smoking = recode_factor(longest_period_wo_smoking,!!!longest_period),
      ts_last_quit_attempt = as.character(ts_last_quit_attempt),
      ts_last_quit_attempt = ifelse(
        n_quit_attempts == "NEVER",
        "never_quit",
        ts_last_quit_attempt
      ),
      trt_recode = forcats::fct_relevel(trt_recode, trt_levels),
      n_quit_attempts = forcats::fct_relevel(n_quit_attempts, quit_levels),
      ts_last_quit_attempt = forcats::fct_relevel(factor(ts_last_quit_attempt), ts_levels),
      anxiety = forcats::fct_relevel(anxiety, anx_dep_levels),
      depression = forcats::fct_relevel(depression, anx_dep_levels),
      int_to_quit = forcats::fct_relevel(int_to_quit, int_quit_levels),
      rsq_last_cig_exp = forcats::fct_relevel(rsq_last_cig_exp, lst_cig_levels)
    )
  
}


management_procedure <- function() {
  # main procedure
  
  df_vector <- load_raw_data()
  
  # clean
  df_vector_sub <- clean_trials(df_vector)
  sid_df <- subject_master_list(df_vector_sub)
  
  swi <- select_vars(df_vector_sub, 1, sid_df, site_name = "swi")
  ger <- select_vars(df_vector_sub, 2, sid_df, site_name = "ger")
  den <- select_vars(df_vector_sub, 3, sid_df, site_name = "den")
  aus <- select_vars(df_vector_sub, 4, sid_df, site_name = "aus")
  usa <- select_vars(df_vector_sub, 5, sid_df, site_name = "usa")
  #
  # add missing columns to aus
  aus[setdiff(names(swi), names(aus))] <- NA
  #
  # # join
  df <- rbind(swi, ger, den, aus, usa)
  
  df <- df %>%
    recode_variables() %>%
    mutate(across(
      c(
        site,
        TRTGRP,
        SEX,
        SMKQT,
        SMKQTLA,
        ANXA,
        DEPR2X,
        LSTCGEX,
        CLMEFTA,
        PEPUP,
        SMKSTP,
        SMKNON
      ),
      as.factor
    ),
    across(c(CLMEFTA, PEPUP), as.numeric)) %>%
    mutate(VISDTF = as.Date(VISDTF)) %>%
    select(-c(ANXA1N, CALMEFFT, PEPUPA, STOPSMK, INTMOT))
  
  # baseline smoking
  df <-
    df %>% mutate(SMKQTY = ifelse(is.na(SMKQTY), CIGCQTY, SMKQTY)) %>%
    select(-CIGCQTY)
  
  # rename time points
  df <- rename_timepoints(df)
  
  df <- rename_recode(df)
  
  df <- reorder_factors(df)
  
  df <- df %>%
    # select(-trt_grp) %>%
    select(subject_id, site, visit, trt_recode, everything())
  
  # make sure nothing is duplicated
  df <- df %>% distinct(subject_id, visit, .keep_all = TRUE)
  
  save(df, file = "../data/interim/compiled_data.rda")
  df
}

# main
df <- management_procedure()

