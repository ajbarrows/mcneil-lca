library(dplyr)


load("../data/interim/ordinal_data.rda")

load("../data/predicted/lca3_predict.rda")
load("../data/processed/pred_imputed.rda")

m3_class_df %>%
  distinct(subject_id) %>%
  nrow()



tab1_df <- m3_class_df %>%
  distinct(subject_id, class, .keep_all = TRUE) %>%
  select(subject_id, class, bl_cpd) %>%
  right_join(pred_imputed, by = "subject_id") %>%
  filter(!is.na(class) & week == "baseline") %>%
  select(subject_id, class, site, trt_recode, sex, age, ftnd, bl_cpd)


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
  select(subject_id, site, trt_recode, sex, age, ftnd, cpd)

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