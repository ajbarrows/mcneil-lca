library(dplyr)
library(tidymodels)
library(tictoc)

load("../data/processed/pred_imputedord.rda")
load("../data/predicted/lca3_predict.rda")

set.seed(123)

doParallel::registerDoParallel()


build_predictor_set <- function(pred_df, lca_df) {
  
  # merge with predicted Latent class
  df <- lca_df %>%
    select(subject_id, class) %>%
    distinct() %>%
    right_join(pred_df, by = "subject_id") %>%
    mutate(class = factor(class)) %>%
    select(-c(subject_id, avg_cpd_1year)) %>%
    as_tibble()
  
  # limit to those assigned a class
  has_class <- df %>% filter(!is.na(class))
  
  # limit to those with CO values at 1-year follow-up
  has_co <- has_class %>% filter(!is.na(co_1year))
  
  list(
    "has_class" = has_class,
    "has_co" = has_co
  )
}

summarize_followup <- function(has_class) {
  
  quit_threshold <- 11
  
  has_class <- has_class %>%
    mutate(quit = ifelse(co_1year <= 11, 1, 0))
  
  has_value_overall <- has_class %>%
    count(!is.na(co_1year)) %>%
    mutate(prp = n / sum(n))
  
  has_value_byclass <- has_class %>%
    group_by(class) %>%
    count(!is.na(co_1year)) %>%
    mutate(prp = n / sum(n))
  
  quit_overall <- has_class %>%
    filter(!is.na(co_1year)) %>%
    group_by(quit) %>%
    count() %>%
    ungroup() %>%
    mutate(prp = n / sum(n))
  
  quit_byclass <- has_class %>%
    filter(!is.na(co_1year)) %>%
    group_by(class, quit) %>%
    count() %>%
    group_by(class) %>%
    mutate(prp = n / sum(n))
  
  tab_list <- list(
    "has_co" = has_value_overall,
    "has_co_byclass" = has_value_byclass,
    "quit" = quit_overall,
    "quit_byclass" = quit_byclass
  )
  
  fpath <- "../reports/tables/quit_counts/"
  for (t in 1:length(tab_list)) {
    write.csv(tab_list[t], paste0(fpath, names(tab_list)[t], ".csv"), row.names = FALSE)
  }
}

prepare_data <- function(train, test) {
  # create recipe with training data
  
  
  enet_recipe <-
    recipe(co_1year ~ ., data = train) %>%
    step_center(all_numeric_predictors()) %>%
    step_scale(all_numeric_predictors()) %>%
    step_novel(all_nominal_predictors())
  
  if ("site" %in% names(train)) {
    enet_recpie <-
      enet_recipe %>%
      step_relevel(site, ref_level = "usa") %>%
      step_relevel(trt_recode, ref_level = "placebo") %>%
      step_relevel(sex, ref_level = "FEMALE")
  }

  enet_recipe <-
    enet_recipe %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    prep()
  # 
  train <- bake(enet_recipe, new_data = NULL)
  test <- bake(enet_recipe, new_data = test)
  
  list(
    "train" = train,
    "test" = test
  )
}

nested_cv_enet <- function(train, n_outer_folds = 5, n_inner_folds = 5) {
  
  outer_split <- caret::createFolds(
    train$co_1year, 
    k = n_outer_folds, 
    list = TRUE, 
    returnTrain = TRUE
    )
  
  metrics <- data.frame()
  models <- list()
  
  for (s in 1:length(outer_split)) {
    fit <- train[outer_split[[s]], ]
    eval <- train[-outer_split[[s]], ]
    
    params <- param_tune(fit, n_inner_folds = n_inner_folds)
    
    enet_model <- 
      linear_reg(penalty = params$tune_best$penalty,
                   mixture = params$tune_best$mixture) %>%
      set_engine("glmnet")
    
    # train model using fit data
    enet_fit <-
      enet_model %>%
      fit(co_1year ~ ., data = fit)
    
    # evaluate model using eval data
    result <- predict(enet_fit, eval)
    result <- augment(enet_fit, eval)
    r_squared <- rsq(result, truth = co_1year, estimate = .pred)
    
    # collect
    test_metrics <- r_squared
    test_metrics$fold <- names(outer_split)[s]
    test_metrics <- cbind(test_metrics, params$tune_best)
    metrics <- rbind(test_metrics, metrics)
    models[[s]] <- enet_fit
  }  
  
  list(
    "cv_metrics" = metrics,
    "cv_models" = models
  )
  
}

fit_ols <- function(train, coefs) {
  
  # get nonzero coefficients 
  features <- coefs %>% 
    filter(term != "(Intercept)" & estimate != 0)
  
  data <- train %>% 
    select(features$term, co_1year)
  
  linear_model <- linear_reg()
  
  
  reg_fit <- 
    linear_model %>%
    fit(co_1year ~ ., data = data)
  # 
  # tidy(log_fit, conf.int = TRUE, exponentiate = TRUE)
  
  return(reg_fit)
}

param_tune <- function(df_prepared, n_inner_folds) {
  tune_split <- vfold_cv(df_prepared, v = n_inner_folds)
  
  tune_spec <-
    linear_reg(penalty = tune(), mixture = tune()) %>%
    set_engine("glmnet")
  
  param_grid <- grid_regular(penalty(),
                             mixture(),
                             levels = list(penalty = 10,
                                           mixture = 10))
  workflow <-
    workflow() %>%
    add_model(tune_spec) %>%
    add_formula(co_1year ~ .)
  
  tune_result <- workflow %>%
    tune_grid(tune_split,
              grid = param_grid,
              metrics = metric_set(rsq))
  
  tune_best <-
    tune_result %>% select_by_one_std_err(mixture, metric = "rsq")
  
  list(
    "tune_result" = tune_result,
    "tune_best" = tune_best
  )
}


regression_pipeline <- function(train, test) {
  tic()
  
  dta <- prepare_data(train, test)
  train <- dta$train
  test <- dta$test
  
  # cross-validated fit
  cv <- nested_cv_enet(
    train,
    n_outer_folds = 5,
    n_inner_folds = 10
  )
  
  
  cvm <- cv$cv_metrics
  max_idx <- which(cvm$.estimate == max(cvm$.estimate))[1]
  
  best_penalty <- cvm$penalty[max_idx]
  best_mixture <- cvm$mixture[max_idx]
  
  print(paste("Best penalty =", best_penalty))
  print(paste("Best mixture =", best_mixture))
  
  best_model <- cv$cv_models[[max_idx]]
  coefs <- tidy(best_model)
  
  # obtain score from lockbox set
  # preds <- predict(best_model, test)
  preds <- augment(best_model, test)
  preds_full <- augment(best_model, rbind(train, test))
  
  r_squared <- rsq(preds, truth = co_1year, estimate = .pred)
  r_squared_full <- rsq(preds_full, truth = co_1year, estimate = .pred)
  
  
  # retrain best model using OLS and non-zero predictors
  # currently using all data
  reg_fit <- fit_ols(rbind(train, test), coefs)
  
  toc()
  
  list(
    "cv" = cv,
    "best_model" = best_model,
    "coef" = coefs,
    "test_score" = r_squared,
    "test_preds" = preds,
    "full_score" = r_squared_full,
    "full_preds" = preds_full,
    "reg_fit" = reg_fit
  )
}

fit_procedure <- function(train, test) {
  
  # only class
  onlyclass_train <- train %>% select(class, co_1year)
  onlyclass_test <- test %>% select(class, co_1year)
  
  # without class
  noclass_train <- train %>% select(-class)
  noclass_test <- test %>% select(-class)
  
  # within classes
  c1_train <- train %>% filter(class == 1) %>% select(-class)
  c1_test <- test %>% filter(class == 1) %>% select(-class)
  
  c2_train <- train %>% filter(class == 2) %>% select(-class)
  c2_test <- test %>% filter(class == 2) %>% select(-class)
  
  c3_train <- train %>% filter(class == 3) %>% select(-class)
  c3_test <- test %>% filter(class == 3) %>% select(-class)
  
  full_fit <- regression_pipeline(train, test)
  noclass_fit <- regression_pipeline(noclass_train, noclass_test)
  onlyclass_fit <- regression_pipeline(onlyclass_train, onlyclass_test)
  c1_fit <- regression_pipeline(c1_train, c1_test)
  c2_fit <- regression_pipeline(c2_train, c2_test)
  c3_fit <- regression_pipeline(c3_train, c3_test)
  
  list(
    "full_fit" = full_fit,
    "noclass_fit" = noclass_fit,
    "onlyclass_fit" = onlyclass_fit,
    "c1_fit" = c1_fit,
    "c2_fit" = c2_fit,
    "c3_fit" = c3_fit
  )
}


# main

# load data
m3_class_df <- model_traj_df

pred_set <- build_predictor_set(pred_imputed, m3_class_df)
has_class <- pred_set$has_class
has_co <- pred_set$has_co

# summarize quit counts
summarize_followup(has_class)


# main

has_co_split <- initial_split(has_co, prop = 0.8, strata = "co_1year")
has_co_train <- training(has_co_split)
has_co_test <- testing(has_co_split)



quit_fit <- fit_procedure(has_co_train, has_co_test)
save(quit_fit, file = "../models/quit_predict_imputedord.rda")


# 
# train <- has_co_train
# test <- has_co_test