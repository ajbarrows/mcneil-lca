library(dplyr)


# load("../data/interim/ordinal_data.rda")

load("../data/interim/compiled_data.rda")
load("../data/predicted/lca3_predict.rda")
load("../data/processed/pred_imputedord.rda")


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


model_traj_df %>%
  distinct(subject_id) %>%
  nrow()



tab1_df <- model_traj_df %>%
  distinct(subject_id, class, .keep_all = TRUE) %>%
  select(subject_id, class, bl_cpd) %>%
  right_join(pred_imputed, by = "subject_id") %>%
  filter(!is.na(class) & week == "baseline") %>%
  select(-c(week, visit_date, co_1year, avg_cpd_1year, quit))


t1 <- tableone::CreateTableOne(
  strata = "class",
  data = tab1_df %>% select(-subject_id),
  test = FALSE,
  includeNA = TRUE
)

t1_overall <- tableone::CreateTableOne(
  data = tab1_df %>% select(-subject_id),
  test = FALSE,
  includeNA = TRUE
)

df_fullsample <- df %>%
  filter(visit == "baseline") %>%
  make_avg_cpd() %>%
  select(any_of(names(tab1_df)))

t1_fullsample <- tableone::CreateTableOne(
  data = df_fullsample %>% select(-subject_id),
  test = FALSE,
  includeNA = FALSE
)



t1_print <- print(t1, quote = FALSE, noSpaces = TRUE)
t1_overall_print <- print(t1_overall, quote = FALSE, noSpaces = TRUE)
t1_fullsample_print <- print(t1_fullsample, quote = FALSE, noSpaces = TRUE)

write.csv(t1_print, file = "../reports/tables/tableone.csv")
write.csv(t1_overall_print, file = "../reports/tables/tableone_overall.csv")
write.csv(t1_fullsample_print, file = "../reports/tables/tableone_fullsample.csv")