library(tidymodels)
library(tictoc)
library(doParallel)


predictors <- c(
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

# functions --- 

join_lca <- function(df_lca, df_bl) {
  as_tibble(df_lca) %>%
    distinct(subject_id, class) %>%
    right_join(df_bl, by = "subject_id") %>%
    select(-c(visit_date, week))
}

split_df <- function(df_joined) {
  # want to fit binomial models separately.
  # also want contrasts. this function should 
  # be run with training/testing data separately
  # not the most elegant code, but it's explicit
  
  # select_vars
  # df_joined <- df_joined %>%
  #   select(-c(co_1year, avg_cpd_1year, quit))
  # 
  # make one-hot encodings for class
  df_joined <- df_joined %>% mutate(class = as.factor(class))
  class_map <- df_joined %>% select(subject_id, class)
  class_map <- as_tibble(model.matrix(~.-1, data = class_map))
  

  class1 <- class_map %>% 
    select(subject_id, class = class1) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
  class2 <- class_map %>% 
    select(subject_id, class = class2) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
  class3 <- class_map %>% 
    select(subject_id, class = class3) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
  
  # contrasts
  
  class12 <- class_map %>%
    filter(class1 == 1 | class2 == 1) %>%
    select(subject_id, class = class1) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
  class13 <- class_map %>%
    filter(class1 == 1 | class3 == 1) %>%
    select(subject_id, class = class1) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
  class23 <- class_map %>%
    filter(class2 == 1 | class3 == 1) %>%
    select(subject_id, class = class2) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
    
  
  # make one-hot encodings for class
  # df_joined <- df_joined %>% mutate(class = as.factor(class))
  class_map <- df_joined %>% select(subject_id, class, trt_recode)
  class_map <- as_tibble(model.matrix(~.-1, data = class_map))
  
  
  # examine people who reduce w/o NRT vs. those who do not reduce w/o NRT.
  # limit sample to placebo condition, use class1 as binary outcome
  
  nonrt <- class_map %>%
    filter(trt_recodeactive == 0) %>%
    select(subject_id, class = class1, trt_recodeactive) %>%
    left_join(df_joined %>% select(-class), by = "subject_id")
  
  
  list(
    "class1" = class1,
    "class2" = class2,
    "class3" = class3,
    "class12" = class12,
    "class13" = class13,
    "class23" = class23,
    "nonrt" = nonrt
  )
  
}


nested_cv_enet <- function(train, n_outer_folds = 5, n_inner_folds = 5) {
  
  outer_split <- caret::createFolds(
    train$class, 
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
      fit(class ~ ., data = fit)
    
    # evaluate model using eval data
    result <- augment(enet_fit, eval)
    auc <- roc_auc(result, truth = class, .pred_yes, estimator = "binary")
    
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


prepare_data <- function(train, test) {
  
  enet_recipe <-
    recipe(class ~ ., data = train) %>%
    step_rm(subject_id) %>%
    step_bin2factor(all_outcomes(), ref_first = FALSE) %>%
    step_normalize(all_numeric_predictors())
  
  if ('site' %in% names(train)) {
    enet_recipe <- enet_recipe %>%
      step_relevel(site, ref_level = "usa")
    
  }
  enet_recipe <- enet_recipe %>%
    step_relevel(sex, ref_level = "FEMALE")
  
  if ('trt_recode' %in% names(train)) {
    enet_recipe <- enet_recipe %>%
      step_relevel(trt_recode, ref_level = "placebo")
  }
  enet_recipe <- enet_recipe %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors()) %>%
    prep()
  
  train <- bake(enet_recipe, new_data = NULL)
  
  # use mean/variance from training data for centering
  test <- bake(enet_recipe, new_data = test)
  
  print(summary(enet_recipe))
  print(enet_recipe)
  
  list(
    "train" = train,
    "test" = test
  )
}

param_tune <- function(train, n_inner_folds = 5) {
  tune_split <- vfold_cv(train, v = n_inner_folds)
  
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
    add_formula(class ~ .)
  
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
  w_sig = sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z = n1 * n2 * (auc - auc2) / w_sig
  pval = 1 - pnorm(abs(z))
  
  return (list(
    "z" = z,
    "p" = pval
  ))
}

sample_null_auc <- function(enet_model, true_auc, train, n_samples = 100) {
  train_null <- train
  null_aucs <- c()
  for (i in 1:n_samples) {
    train_null$class <- sample(train_null$class)
    null_fit <- 
      enet_model %>%
      fit(class ~ ., data = train_null)
    result <- augment(null_fit, train_null)
    
    auc <-
      roc_auc(result,
              truth = class,
              .pred_yes,
              estimator = "binary")
    null_aucs <- append(null_aucs, auc$.estimate)
  }

  null_aucs
}


return_counts <- function(train_split, test_split) {
  n_train_all <- nrow(train_split$class1)
  n_test_all <- nrow(test_split$class1)
  
  n_train_nonrt <- nrow(train_split$nonrt)
  n_test_nonrt <- nrow(test_split$nonrt)
  
  colnames <- c("model", "train", "test")
  model <- c("all", "nonrt")
  all <- c(n_train_all, n_test_all)
  nonrt <- c(n_train_nonrt, n_test_nonrt)
  
  df <- data.frame(
    model, all, nonrt
  )
  names(df) <- colnames
  write.csv(
    df, 
    "../reports/tables/model_setup/class_predict_split.csv", 
    row.names = FALSE
  )
  return (df)
  
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
    add_formula(class ~ .)
  
  # fit model 
  enet_fit <-
    enet_model %>%
    fit(class ~ ., data = data)
  
  coef <- tidy(enet_fit)
  
  # evaluate model
  result <- predict(enet_fit, data)
  result <- augment(enet_fit, data)
  auc <- roc_auc(result, truth = class, .pred_yes, estimator = "binary")
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
    
    mw_u <- auc2p(
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
    "mw_u" = mw_u,
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
  test_score <- roc_auc(preds, truth = class, .pred_yes, event_level = "second")
  test_roc <- roc_curve(preds, truth = class, .pred_yes, event_level = "second")
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
  
  mw_u <- auc2p(
    auc = test_score$.estimate,
    n1 = test %>% filter(class == "yes") %>% nrow(),
    n2 = test %>% filter(class == "no") %>% nrow(),
    auc2 = mean(null_aucs)
  )
  
  
  toc()
  
  list(
    "cv" = cv,
    "best_model" = best_model,
    "coef" = coefs,
    "test_score" = test_score,
    "test_preds" = preds,
    "null_aucs" = null_aucs,
    "mw_u" = mw_u,
    "test_roc" = test_roc
  )
}

fit_procedure <- function(train_splits, test_splits) {
  c1_fit <- enet_pipeline(train_splits$class1, test_splits$class1)
  c2_fit <- enet_pipeline(train_splits$class2, test_splits$class2)
  c3_fit <- enet_pipeline(train_splits$class3, test_splits$class3)
  
  c12_fit <- enet_pipeline(train_splits$class12, test_splits$class12)
  c13_fit <- enet_pipeline(train_splits$class13, test_splits$class13)
  c23_fit <- enet_pipeline(train_splits$class23, test_splits$class23)

  if ('site' %in% names(train_splits$nonrt)) {
    nonrt_fit <- NULL
  } else {
    nonrt_fit <- enet_pipeline(train_splits$nonrt, test_splits$nonrt)
  }
  
  list(
    "c1_fit" = c1_fit,
    "c2_fit" = c2_fit,
    "c3_fit" = c3_fit,
    "c12_fit" = c12_fit,
    "c13_fit" = c13_fit,
    "c23_fit" = c23_fit,
    "nonrt_fit" = nonrt_fit
  )
  
}



if (sys.nframe() == 0) {
  
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
  
  df_imputed <- join_lca(m3_class_df, pred_imputed) %>%
    select(subject_id, class, all_of(predictors))

  df_imputed_split <- initial_split(df_imputed, prop = 0.8)
  imputed_train <- training(df_imputed_split)
  imputed_test <- testing(df_imputed_split)
  
  # df_nomiss <- join_lca(m3_class_df, pred_nomiss)
  # df_nomiss_split <- initial_split(df_nomiss, prop = 0.8)
  # nomiss_train <- training(df_nomiss_split)
  # nomiss_test <- testing(df_nomiss_split)
  
  
  # train
  
  imputed_train_split <- split_df(imputed_train)
  imputed_test_split <- split_df(imputed_test)
  
  return_counts(imputed_train_split, imputed_test_split)
  
  # usa
  imputed_train_usa <- split_df(imputed_train %>% filter(site == "usa") %>% select(-site))
  imputed_test_usa <- split_df(imputed_test %>% filter(site == "usa") %>% select(-site))
  
  # aus
  imputed_train_aus <- split_df(imputed_train %>% filter(site == "aus") %>% select(-site))
  imputed_test_aus <- split_df(imputed_test %>% filter(site == "aus") %>% select(-site))
  
  # swi
  imputed_train_swi <- split_df(imputed_train %>% filter(site == "swi") %>% select(-site))
  imputed_test_swi <- split_df(imputed_test %>% filter(site == "swi") %>% select(-site))
  
  # den
  imputed_train_den <- split_df(imputed_train %>% filter(site == "den") %>% select(-site))
  imputed_test_den <- split_df(imputed_test %>% filter(site == "den") %>% select(-site))
  
  # ger
  imputed_train_ger <- split_df(imputed_train %>% filter(site == "ger") %>% select(-site))
  imputed_test_ger <- split_df(imputed_test %>% filter(site == "ger") %>% select(-site))
  
  
  
  imputed_fits <- fit_procedure(imputed_train_split, imputed_test_split)
  # usa <- fit_procedure(imputed_train_usa, imputed_test_usa)
  # aus <- fit_procedure(imputed_train_aus, imputed_test_aus)
  # swi <- fit_procedure(imputed_train_swi, imputed_test_swi)
  # den <- fit_procedure(imputed_train_den, imputed_test_den)
  # ger <- fit_procedure(imputed_train_ger, imputed_test_ger)
  
  # bysite
  
  
  
  # nomiss_fits <- fit_procedure(nomiss_train_split, nomiss_test_split)
  # # 
  save(imputed_fits, file = "../models/enet_imputed_ord.rda")
  # save(usa, file = "../models/predict_class_usa.rda")
  # save(aus, file = "../models/predict_class_aus.rda")
  # save(swi, file = "../models/predict_class_swi.rda")
  # save(den, file = "../models/predict_class_den.rda")
  # save(ger, file = "../models/predict_class_ger.rda")

}











# save(nomiss_fits, file = "../models/enet_nomiss.rda")
# 
# train <- imputed_train_split$class1
# test <- imputed_test_split$class