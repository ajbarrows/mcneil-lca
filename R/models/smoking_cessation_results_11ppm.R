library(tidymodels)
library(ggplot2)

source("visualization/visualize.R")

# load("../models/quit_predict_imputedord.rda")
load("../models/quit_predict_11ppm.rda")
theme_set(theme_classic())

fig_path <- '../reports/figures/11ppm/'
table_path <- '../reports/tables/11ppm/'


metrics <- data.frame()
rocs <- data.frame()

for (f in 1:length(quit_fit)) {
  
  # elastic net results
  nme <- names(quit_fit[f])
  cv_roc_mean <- mean(quit_fit[[f]]$cv$cv_metrics$mean)
  cv_roc_se <- sd(quit_fit[[f]]$cv$cv_metrics$std_err)
  test_auc <- quit_fit[[f]]$test_score$.estimate
  auc_pval <- quit_fit[[f]]$auc_pval
  
  tmp <- data.frame(nme, cv_roc_mean, cv_roc_se, test_auc, auc_pval)
  names(tmp)[names(tmp) == ".estimate"] <- "test_auc"
  metrics <- rbind(metrics, tmp)
  
  
  # 
  # t <- quit_fit[[f]]$test_preds %>%
  #   select(truth = quit_verified, estimate = .pred)
  # 
  # for (thresh in seq(0, 100, by = 1)) {
  #   
  #   t <- arrange(t, truth)
  #   
  #   # actually quit based on threshold
  #   t$quit_actual <- t$truth <= thresh
  #   
  #   # predicted quit based on threshold
  #   t$quit_test <- round(t$estimate) <= thresh
  #   
  #   # labeling
  #   t$tp <- t$quit_actual & t$quit_test
  #   t$fp <- !t$quit_actual & t$quit_test
  #   t$fn <- t$quit_actual & !t$quit_test
  #   t$tn <- !t$quit_actual & !t$quit_test
  #   
  #   
  #   tpr <- sum(t$tp) / sum(t$tp, t$fn)
  #   tnr <- sum(t$tn) / sum(t$fp, t$tn)
  #   
  #   fpr <- 1 - tnr
  #   
  #   tmp_df <- data.frame(thresh, tpr, fpr, nme, sum(t$tp), sum(t$fp), sum(t$fn), sum(t$tn))
  #   rocs <- rbind(rocs, tmp_df)
  # }

  
}


gather_cv_coef <- function(imputed_fits, use_ref = TRUE, exponentiate = TRUE) {
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
  
  if (exponentiate) {
    cv_coefs <- cv_coefs %>%
      mutate(estimate = exp(estimate))
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

auc2p <- function(auc, n1, n2, auc2 = 0.5) {
  # From Nick Allgaier
  w_sig = sqrt(n1*n2*(n1+n2+1)/12)
  z=n1*n2*(auc-auc2)/w_sig
  pval=1-pnorm(abs(z))
  pval
}


compare_models <- function(metrics, quit_fit) {
  quitcount_test <- table(quit_fit$full_fit$test_preds$quit_verified)
  n_yes = as.integer(quitcount_test[2])
  n_no = as.integer(quitcount_test[1])
  
  auc_full <- metrics$test_auc[metrics$nme == "full_fit"]
  auc_noclass <- metrics$test_auc[metrics$nme == "noclass_fit"]
  
  auc_pval <- auc2p(auc_full, n_yes, n_no, auc_noclass)
  
  return(paste0("P(Full-Fit AUC > No-Class AUC) = ", auc_pval))
}


metrics
write.csv(metrics, file = paste0(table_path, 'quit_predict_perf.csv'), row.names = FALSE)

# compare full fit and no-class fit
compare_models(metrics, quit_fit)


cv_coef <- gather_cv_coef(quit_fit, use_ref = FALSE)
# 
# t <- cv_coef$cv_coefs_long %>%
#   mutate(
#     mean_odds = exp(mean_est),
#     sd_odds = exp(sd_est))


feature_imp(cv_coef$cv_coefs_long %>% filter(!stringr::str_detect(model, "Within|Latent Class Alone"))) + 
  scale_color_brewer(palette = "Set2")
ggsave(paste0(fig_path, "predict_quit_features.png"), dpi = 600)

# write.csv(cv_coefs, "../reports/tables/co_predict_coef.csv", row.names = FALSE)
# select best model
full_fit <- quit_fit$full_fit
noclass_fit <- quit_fit$noclass_fit
## main

compare_metrics(metrics %>% filter(!nme %in% c("c1_fit", "c2_fit", "c3_fit")))
ggsave(paste0(fig_path, "predict_quit_compare.png"))

quit_coef(full_fit)
ggsave(paste0(fig_path, "predict_quit_coefs.png"))

quit <- count_quit(full_fit)
write.csv(quit, paste0(table_path, "predict_quit_count.csv"), row.names = FALSE)




# ggsave(plot = quit$plot, "../reports/figures/ord/count_quit.png")

# compare_roc(full_fit$test_roc, noclass_fit$test_roc)
# ggsave("../reports/figures/ord/predict_quit_rocs.svg")

svg(paste0(fig_path, "predict_quit_rocs.svg"))
roc_plot_quit(full_fit$test_roc, noclass_fit$test_roc)
dev.off()



full_fit$test_roc

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

