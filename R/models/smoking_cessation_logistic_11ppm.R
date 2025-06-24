library(dplyr)
library(tidymodels)
library(tictoc)
library(doParallel)

predictors <- c(
  'class',
  'site',
  'trt_recode',
  'sex',
  'age',
  'ftnd',
  'co',
  'sf_physheal',
  'sf_emoprob',
  'sf_socfunc',
  'sf_pain',
  'sf_emowell',
  'longest_period_wo_smoking',
  'age_started_smoking',
  'n_quit_attempts',
  'ts_last_quit_attempt',
  'anxiety',
  'depression',
  'int_to_quit',
  'rsq_calming',
  'rsq_last_cig_exp',
  'rsq_pepping_up_eff',
  'avg_cpd'
)


load("../data/processed/pred_imputedord.rda")
load("../data/predicted/lca3_predict.rda")

set.seed(123)

# parallelization

all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)

impose_quit <- function(df, quit_thresh = 12) {
  df <- df %>%
    mutate(quit_verified = case_when(
      quit == "YES" & co_1year < quit_thresh ~ 1,
      co_1year >= quit_thresh ~ 0,
      quit == "NO" ~ 0,
      is.na(co_1year) | is.na(quit) ~ 0
    )) %>%
    select(-c(quit))
  
  return(df)
}


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

summarize_followup <- function(quit_set) {
  
  # quit_threshold <- 11
  
  # has_class <- has_class %>%
  #   mutate(quit = ifelse(co_1year <= quit_threshold, 1, 0))
  
  has_value_overall <- quit_set %>%
    count(!is.na(co_1year)) %>%
    mutate(prp = n / sum(n))
  
  has_value_byclass <- quit_set %>%
    group_by(class) %>%
    count(!is.na(co_1year)) %>%
    mutate(prp = n / sum(n))
  
  quit_overall <- quit_set %>%
    group_by(quit_verified) %>%
    count() %>%
    ungroup() %>%
    mutate(prp = n / sum(n))
  
  quit_byclass <- quit_set %>%
    group_by(class, quit_verified) %>%
    count() %>%
    group_by(class) %>%
    mutate(prp = n / sum(n))
  
  tab_list <- list(
    "has_co" = has_value_overall,
    "has_co_byclass" = has_value_byclass,
    "quit" = quit_overall,
    "quit_byclass" = quit_byclass
  )
  
  fpath <- "../reports/tables/11ppm/quit_counts/"
  for (t in 1:length(tab_list)) {
    write.csv(tab_list[t], paste0(fpath, names(tab_list)[t], ".csv"), row.names = FALSE)
  }
}

prepare_data <- function(train, test) {
  # create recipe with training data

  enet_recipe <-
    recipe(quit_verified ~ ., data = train)

  if ("site" %in% names(train)) {
    enet_recipe <-
      enet_recipe %>%
      step_relevel(site, ref_level = "usa") %>%
      step_relevel(trt_recode, ref_level = "placebo") %>%
      step_relevel(sex, ref_level = "FEMALE")
  }

  enet_recipe <-
    enet_recipe %>%
    step_center(all_numeric_predictors()) %>%
    step_scale(all_numeric_predictors()) %>%
    step_novel(all_nominal_predictors()) %>%
    step_bin2factor(all_outcomes(), ref_first = FALSE) %>%
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


return_counts <- function(train_split, test_split) {
  
  n_train <- nrow(train_split)
  n_test <- nrow(test_split)
  
  df <- data.frame(n_train, n_test)
  
  # names(df) <- colnames
  write.csv(
    df,
    "../reports/tables/model_setup/cessation_predict_split.csv",
    row.names = FALSE
  )
  return (df)
}

nested_cv_enet <- function(train, n_outer_folds = 5, n_inner_folds = 5) {
  
  outer_split <- caret::createFolds(
    train$quit_verified, 
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
      logistic_reg(penalty = params$tune_best$penalty,
                   mixture = params$tune_best$mixture) %>%
      set_engine("glmnet")
    
    # train model using fit data
    enet_fit <-
      enet_model %>%
      fit(quit_verified ~ ., data = fit)
    
    # evaluate model using eval data

    result <- augment(enet_fit, eval)

    auc <- roc_auc(
      result, 
      truth = quit_verified, 
      .pred_yes, 
      estimator = "binary")
    # f1 <- f_meas(result, truth = class, estimate = .pred_class)
    
    # collect
    # test_metrics <- rbind(auc, f1)
    test_metrics <- auc
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

param_tune <- function(df_prepared, n_inner_folds) {
  tune_split <- vfold_cv(df_prepared, v = n_inner_folds)
  
  tune_spec <-
    logistic_reg(penalty = tune(), mixture = tune()) %>%
    set_engine("glmnet")
  
  param_grid <- grid_regular(penalty(),
                             mixture(),
                             levels = list(penalty = 10,
                                           mixture = 10))
  workflow <-
    workflow() %>%
    add_model(tune_spec) %>%
    add_formula(quit_verified ~ .)
  
  tune_result <- workflow %>%
    tune_grid(tune_split,
              grid = param_grid,
              metrics = metric_set(roc_auc))
  
  tune_best <-
    tune_result %>% select_by_one_std_err(mixture, metric = "roc_auc")
  
  list(
    "tune_result" = tune_result,
    "tune_best" = tune_best
  )
}


auc2p <- function(auc, n1, n2, auc2 = 0.5) {
  # From Nick Allgaier
  w_sig = sqrt(n1*n2*(n1+n2+1)/12)
  z=n1*n2*(auc-auc2)/w_sig
  pval=1-pnorm(abs(z))
  pval
}

sample_null_auc <- function(enet_model, true_auc, train, n_samples = 100) {
  train_null <- train
  null_aucs <- c()
  for (i in 1:n_samples) {
    train_null$quit_verified <- sample(train_null$quit_verified)
    null_fit <- 
      enet_model %>%
      fit(quit_verified ~ ., data = train_null)
    result <- augment(null_fit, train_null)
    
    auc <-
      roc_auc(result,
              truth = quit_verified,
              .pred_yes,
              estimator = "binary")
    null_aucs <- append(null_aucs, auc$.estimate)
  }
  
  null_aucs
}

enet_fit <- function(data, penalty, mixture, permute_null = FALSE) {
  
  enet_model <- 
    logistic_reg(penalty = penalty,
                 mixture = mixture) %>%
    set_engine("glmnet")
  
  # model-fit
  enet_workflow <-
    workflow() %>%
    add_model(enet_model) %>%
    add_formula(quit_verified ~ .)
  
  # fit model 
  enet_fit <-
    enet_model %>%
    fit(quit_verified ~ ., data = data)
  
  coef <- tidy(enet_fit)
  
  # evaluate model
  result <- predict(enet_fit, data)
  result <- augment(enet_fit, data)
  auc <- roc_auc(result, truth = class, estimate = .pred_yes, estimator = "binary")
  roc <- roc_curve(result, truth = class)
  test_metrics <- auc
  
  if (permute_null) {
    # permute outcomes, test AUC against null AUC
    null_aucs <- sample_null_auc(
      enet_model,
      true_auc = auc$.estimate,
      data,
      n_samples = 100
    )
    
    auc_pval <- auc2p(
      auc = auc$.estimate,
      n1 = data %>% filter(class == 1) %>% nrow(),
      n2 = data %>% filter(class == 0) %>% nrow(),
      auc2 = mean(null_aucs)
    )
    
  } else {
    null_aucs <- auc_pval <- NULL
  }
  
  
  list(
    "coef" = coef,
    "test_metrics" = test_metrics,
    "test_roc" = roc,
    "predicted" = result,
    "full_fit" = enet_fit,
    "null_aucs" = null_aucs,
    "auc_pval" = auc_pval,
    "workflow" = enet_workflow
  )
}
enet_pipeline <- function(train, test) {
  # pipeline
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
  
  # Obtain test prediction
  preds <- augment(best_model, test)
  test_score <- roc_auc(preds, truth = quit_verified, .pred_yes, event_level = "second")
  test_roc <- roc_curve(preds, truth = quit_verified, .pred_yes, event_level = "second")
  # Obtain p-value from test prediction
  mod <-
    logistic_reg(penalty = best_penalty,
                 mixture = best_mixture) %>%
    set_engine("glmnet")
  
  null_aucs <- sample_null_auc(
    mod,
    true_auc = test_score$.estimate,
    test,
    n_samples = 100
  )
  
  auc_pval <- auc2p(
    auc = test_score$.estimate,
    n1 = test %>% filter(quit_verified == "yes") %>% nrow(),
    n2 = test %>% filter(quit_verified == "no") %>% nrow(),
    auc2 = mean(null_aucs)
  )
  
  
  # retrain best model using OLS and non-zero predictors
  # currently using all data
  # log_fit <- fit_ols(rbind(train, test), coefs)
  # 
  
  toc()
  
  list(
    "cv" = cv,
    "best_model" = best_model,
    "coef" = coefs,
    "test_score" = test_score,
    "test_preds" = preds,
    "null_aucs" = null_aucs,
    "auc_pval" = auc_pval,
    # "log_fit" = log_fit,
    "test_roc" = test_roc
  )
}

fit_procedure <- function(train, test) {
  
  # only class
  onlyclass_train <- train %>% select(class, quit_verified)
  onlyclass_test <- test %>% select(class, quit_verified)
  
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
  
  full_fit <- enet_pipeline(train, test)
  noclass_fit <- enet_pipeline(noclass_train, noclass_test)
  onlyclass_fit <- enet_pipeline(onlyclass_train, onlyclass_test)
  c1_fit <- enet_pipeline(c1_train, c1_test)
  c2_fit <- enet_pipeline(c2_train, c2_test)
  c3_fit <- enet_pipeline(c3_train, c3_test)
  
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

quit_set <- impose_quit(has_class)

# summarize quit counts
summarize_followup(quit_set)

quit_set <- quit_set %>% 
  select(quit_verified, all_of(predictors))

# main

quit_split <- initial_split(quit_set, prop = 0.8, strata = "quit_verified")
quit_train <- training(quit_split)
quit_test <- testing(quit_split)


return_counts(quit_train, quit_test)

# quit_test %>%
#   group_by(class) %>%
#   count(quit_verified)




quit_fit <- fit_procedure(quit_train, quit_test)
save(quit_fit, file = "../models/quit_predict_11ppm.rda")


# 
# train <- has_co_train
# test <- has_co_test