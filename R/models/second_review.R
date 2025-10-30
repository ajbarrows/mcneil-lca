source("./models/predict_class_ordinal.R")

make_test_train_split <- function(df, prop = 0.8) {
  df_split <- initial_split(df, prop = prop)
  train <- training(df_split)
  test <- testing(df_split)

  list(
    'train_split' = split_df(train),
    'test_split' = split_df(test)
  )
}

source("visualization/visualize.R")

load("../data/processed/pred_imputedord.rda")
load("../data/processed/pred_nomissord.rda")
load("../data/predicted/lca3_predict.rda")

set.seed(42)

# parallelization

all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)


# main --- 
m3_class_df <- model_traj_df

df <- join_lca(m3_class_df, pred_imputed) %>%
  select(subject_id, class, all_of(predictors))



train_proportions <- c(0.9, 0.8, 0.7, 0.6)
fits <- list()


for (prop in train_proportions) {
  print(paste("Training on ", prop * 100, "% of the data."))
  split <- make_test_train_split(df, prop)
  fits[[as.character(prop)]] <- fit_procedure(
    split$train_split,
    split$test_split,
    include_nonrtfit = FALSE,
    include_comparison_fits = FALSE
  )
}

saveRDS(fits, "../models/data_split_sensitivity.rds")