library(tidymodels)
library(ggplot2)

source("visualization/visualize.R")

load("../models/quit_predict_imputedord.rda")

theme_set(theme_classic(base_size = 15))


metrics <- data.frame()
rocs <- data.frame()

for (f in 1:length(quit_fit)) {
  
  # elastic net results
  nme <- names(quit_fit[f])
  cv_rsq_mean <- mean(quit_fit[[f]]$cv$cv_metrics$mean)
  cv_rsq_se <- sd(quit_fit[[f]]$cv$cv_metrics$std_err)
  test_rsq <- quit_fit[[f]]$test_score
  
  tmp <- data.frame(nme, cv_rsq_mean, cv_rsq_se, test_rsq)
  names(tmp)[names(tmp) == ".estimate"] <- "test_rsq"
  metrics <- rbind(metrics, tmp)
  
  
  
  t <- quit_fit[[f]]$test_preds %>%
    select(truth = co_1year, estimate = .pred)
  
  for (thresh in seq(0, 100, by = 1)) {
    
    t <- arrange(t, truth)
    
    # actually quit based on threshold
    t$quit_actual <- t$truth <= thresh
    
    # predicted quit based on threshold
    t$quit_test <- round(t$estimate) <= thresh
    
    # labeling
    t$tp <- t$quit_actual & t$quit_test
    t$fp <- !t$quit_actual & t$quit_test
    t$fn <- t$quit_actual & !t$quit_test
    t$tn <- !t$quit_actual & !t$quit_test
    
    
    tpr <- sum(t$tp) / sum(t$tp, t$fn)
    tnr <- sum(t$tn) / sum(t$fp, t$tn)
    
    fpr <- 1 - tnr
    
    tmp_df <- data.frame(thresh, tpr, fpr, nme, sum(t$tp), sum(t$fp), sum(t$fn), sum(t$tn))
    rocs <- rbind(rocs, tmp_df)
  }

  
}


gather_cv_coef <- function(imputed_fits, use_ref = TRUE) {
  term_feature <- read.csv(
    "../data/term_feature_map.csv",
    stringsAsFactors = FALSE
  )
  
  level_order <- read.csv(
    "../data/level_ordering.csv",
    stringsAsFactors = FALSE
  )
  
  reference <- read.csv(
    "../data/reference_term_map_cessation.csv",
    stringsAsFactors = FALSE
  )
  
  cv_coefs <- data.frame()
  for (m in 1:length(imputed_fits)) {
    cv_models <- imputed_fits[[m]]$cv$cv_models
    for (f in 1:length(cv_models)) {
      tmp <- tidy(cv_models[[f]])
      tmp$fold <- f
      tmp$model <- names(imputed_fits[m])
      cv_coefs <- rbind(cv_coefs, tmp)
    }
  }
  
  
  cv_coefs_long <- cv_coefs %>%
    filter(!model %in% c("c12_fit", "c13_fit", "c23_fit")) %>%
    group_by(model, term) %>%
    summarize(mean_est = mean(estimate), sd_est = sd(estimate)) %>%
    left_join(term_feature, by = "term")
  
  if (use_ref) {
    cv_coefs_long <- cv_coefs_long %>% rbind(reference)
  }
  cv_coefs_long <- cv_coefs_long %>% 
    mutate(
      model = case_when(
        model == "c1_fit" ~ "Within Class 1",
        model == "c2_fit" ~ "Within Class 2",
        model == "c3_fit" ~ "Within Class 3",
        model == "full_fit" ~ "Baseline Char. &\nLatent Class",
        model == "noclass_fit" ~ "Baseline Characteristics\n Alone",
        model == "onlyclass_fit" ~ "Latent Class Alone"
      )
    ) %>%
    mutate(
      feature = as.factor(feature),
      level = factor(level, levels = c(level_order$levels, ""))
    ) %>%
    arrange(model, feature, level)
  
  
  cv_coefs_long$level[is.na(cv_coefs_long$level)] <- ""
  
  cv_coefs <- cv_coefs_long %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(mean_est, sd_est),
      names_vary = "slowest"
    )
  
  list(
    "cv_coefs" = cv_coefs,
    "cv_coefs_long" = cv_coefs_long
  )
}



metrics
write.csv(metrics, "../reports/tables/ord/co_predict_perf.csv", row.names = FALSE)


cv_coef <- gather_cv_coef(quit_fit, use_ref = FALSE)

t <- cv_coef$cv_coefs_long

feature_imp(cv_coef$cv_coefs_long %>% filter(!stringr::str_detect(model, "Within"))) + 
  scale_color_brewer(palette = "Set2")
ggsave("../reports/figures/ord/predict_co_features.png", dpi = 600, width = 11, height = 10)

# write.csv(cv_coefs, "../reports/tables/co_predict_coef.csv", row.names = FALSE)
# select best model
full_fit <- quit_fit$full_fit
# full_fit$coef <- full_fit$coef %>%
#   filter(term != "(Intercept)") %>%
#   mutate(grp = stringr::str_extract(term, "^(.*?\\_)")) %>%
#   mutate(grp = ifelse(is.na(grp), term, grp)) %>%
#   group_by(grp) %>%
#   mutate(max_est = max(abs(estimate))) %>%
#   arrange(-max_est)
# 
# write.csv(full_fit$coef, "../reports/tables/enet_co_coef.csv")
# 
# write.csv(coefs, "../reports/tables/co_predict_ols_coefs.csv", row.names = FALSE)



# t <- full_fit$coef
# res <- quit_fit$full_fit$test_preds
# 
# res_count <- res %>% 
#   select(class_X1, class_X2, class_X3, co_1year, .pred) %>% 
#   mutate(
#     quit_actual = factor(ifelse(co_1year <= 11, "quit", "not_quit")),
#     quit_pred = factor(ifelse(.pred <= 11, "quit", "not_quit")),
#     success = quit_actual == quit_pred
#   ) %>%
#   mutate(latent_class = case_when(
#     class_X1 == 1 ~ 1,
#     class_X2 == 1 ~ 2,
#     class_X3 == 1 ~ 3)
#   )
# 
# res_count %>% 
#   group_by(latent_class, quit_pred) %>% 
#   count(quit_actual) %>%
#   tidyr::pivot_wider(
#     names_from = c(latent_class, quit_actual),
#     values_from = n,
#     id_expand = TRUE
#   )
# 
# 
# res_count %>%
#   group_by(latent_class) %>%
#   count(quit_actual)

res <- quit_fit$full_fit$full_preds
# 
# res_count <- res %>% 
#   select(class_X1, class_X2, class_X3, co_1year, .pred) %>% 
#   mutate(
#     quit_actual = factor(ifelse(co_1year <= 11, "quit", "not_quit")),
#     quit_pred = factor(ifelse(.pred <= 11, "quit", "not_quit")),
#     success = quit_actual == quit_pred
#   ) %>%
#   mutate(latent_class = case_when(
#     class_X1 == 1 ~ 1,
#     class_X2 == 1 ~ 2,
#     class_X3 == 1 ~ 3)
#   )

# res_count %>% 
#   group_by(latent_class, quit_pred) %>% 
#   count(quit_actual) %>%
#   tidyr::pivot_wider(
#     names_from = c(latent_class, quit_actual),
#     values_from = n,
#     id_expand = TRUE
#   )


# res_count %>%
#   group_by(latent_class) %>%
#   count(quit_actual)



## main

compare_metrics(metrics %>% filter(!nme %in% c("c1_fit", "c2_fit", "c3_fit")))
ggsave("../reports/figures/ord/predict_quit_compare.png")

quit_coef(full_fit)
ggsave("../reports/figures/ord/predict_quit_coefs.png")

quit <- count_quit(full_fit)
ggsave(plot = quit$plot, "../reports/figures/ord/count_quit.png")

png("../reports/figures/ord/predict_quit_rocs.png", height = 2000, width = 2000, res = 300)
roc_plot_quit(rocs)
dev.off()


# auc <- function(roc){
#   height = (roc$tpr[-1] + roc$tpr[-length(roc$tpr)])/2
#   width = -diff(roc$fpr)
#   sum(height * width)
# }
# 
# # auc
# 
# 
# full_roc <- rocs %>% filter(nme == "full_fit") %>% tidyr::drop_na()
# 

