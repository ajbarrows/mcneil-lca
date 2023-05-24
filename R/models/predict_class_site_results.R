library(tidymodels)
library(ggplot2)
library(stringr)

source("visualization/visualize.R")

theme_set(theme_classic(base_size = 15))
# 
# load("../models/enet_imputed_ord.rda")
# 
# 
# enet_results <- function(fits, plot_loc, prefix) {
#   
#   # main effect 
#   c1 <- fits$c1_fit$coef
#   c2 <- fits$c2_fit$coef
#   c3 <- fits$c3_fit$coef
#   
#   c1$class <- "class_1"
#   c2$class <- "class_2"
#   c3$class <- "class_3"
#   
#   main_plt <- rbind(c1, c2, c3)
#   
#   plot_coef_byclass(main_plt, threshold = 0, scales = "free_x")
#   ggsave(paste0(plot_loc, prefix, "full.png"), width = 10, height = 8, units = "in")
#   
#   # contrasts
#   c12 <- fits$c12_fit$coef
#   c12$class <- NA
#   c12$class[c12$estimate < 0] <- "class_2"
#   c12$class[c12$estimate >= 0] <- "class_1"
#   
#   plot_coef_byclass(c12, threshold = 0, facet = FALSE)
#   ggsave(paste0(plot_loc, prefix, "c12.png"))
#   
#   c13 <- fits$c13_fit$coef
#   c13$class <- NA
#   c13$class[c13$estimate < 0] <- "class_3"
#   c13$class[c13$estimate >= 0] <- "class_1"
#   
#   plot_coef_byclass(c13, threshold = 0, facet = FALSE)
#   ggsave(paste0(plot_loc, prefix, "c13.png"))
#   
#   c23 <- fits$c23_fit$coef
#   c23$class <- NA
#   c23$class[c23$estimate < 0] <- "class_2"
#   c23$class[c23$estimate >= 0] <- "class_3"
#   
#   plot_coef_byclass(c23, threshold = 0, facet = FALSE)
#   ggsave(paste0(plot_loc,  prefix, "c23.png"))
#   
#   contrap <- fits$contrap_fit$coef
#   contrap$class <- NA
#   contrap$class[contrap$estimate < 0] <- "class_1_placebo"
#   contrap$class[contrap$estimate >=0] <- "class_23_active"
#   
#   main_table <- main_plt %>% select(-penalty)
#   
#   main_table <- main_table %>%
#     mutate(grp = stringr::str_extract(term, "^(.*?\\_)")) %>%
#     mutate(grp = ifelse(is.na(grp), term, grp)) %>%
#     # mutate(term = stringr::str_remove(term, grp)) %>%
#     group_by(grp) %>%
#     mutate(max_est = max(abs(estimate))) %>%
#     ungroup()
#   
#   main_table <- 
#     main_table %>%
#     tidyr::pivot_wider(
#       names_from = "class",
#       values_from = "estimate"
#     ) %>%
#     select(grp, term, everything()) %>%
#     arrange(-max_est)
#   
#   c12$contr <- "c12"
#   c13$contr <- "c13"
#   c23$contr <- "c23"
#   contrap$contr <- "c1placebo_c23active"
# 
#   
#   contr_table <- rbind(c12, c13, c23, contrap) %>%
#     select(-c(penalty, class)) %>%
#     tidyr::pivot_wider(
#       names_from = "contr",
#       values_from = "estimate"
#     )
# 
#   
#   main_table <- main_table %>%
#     left_join(contr_table, by = "term")
#   
#   write.csv(
#     main_table, 
#     paste0("../reports/tables/ord/", prefix, "_coef.csv"), row.names = FALSE
#     )
# }
# 
# 
# # imputed ----- 
# 
# enet_results(imputed_fits, "../reports/figures/ord/", "enet_imputed")
# 
# metrics <- data.frame()
# coefs <- data.frame()
# rocs <- data.frame()
# options(scipen = 999)
# for (f in 1:length(imputed_fits)) {
#   # enet metrics
#   nme <- names(imputed_fits[f])
#   cv_auc_mean <- imputed_fits[[f]]$cv$cv_metrics$mean[2]
#   cv_auc_se <- imputed_fits[[f]]$cv$cv_metrics$std_err[2]
#   test_auc <- imputed_fits[[f]]$test_score$.estimate
#   null_auc_mean <- mean(imputed_fits[[f]]$null_aucs)
#   auc_pval <- imputed_fits[[f]]$auc_pval
#   
#   tmp <- data.frame(nme, cv_auc_mean, cv_auc_se, test_auc, null_auc_mean, auc_pval)
#   metrics <- rbind(metrics, tmp)
#   
#   # ols coefficients
#   # 
#   # coefs_df <- tidy(imputed_fits[[f]]$log_fit, exponentiate = TRUE, conf.int = TRUE)
#   # coefs_df$model <- nme
#   # 
#   # coefs <- rbind(coefs, coefs_df)
#   roc_tmp <- imputed_fits[[f]]$test_roc
#   roc_tmp$model <- nme
#   rocs <- rbind(rocs, roc_tmp)
#   
# }
# 
# 
gather_cv_coef <- function(imputed_fits, site = NULL, use_ref = TRUE,
                           use_site = FALSE, exponentiate = TRUE) {


  term_feature <- read.csv(
    "../data/term_feature_map.csv",
    stringsAsFactors = FALSE
  )

  level_order <- read.csv(
    "../data/level_ordering.csv",
    stringsAsFactors = FALSE
  )

  reference <- read.csv(
    "../data/reference_term_map.csv",
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

  levels <- c("Class 1 vs. All", "Class 2 vs. All", "Class 3 vs. All", "C1 Placebo vs.\nC2 & C3 Active")

  cv_coefs_long <- cv_coefs %>%
    filter(!model %in% c("c12_fit", "c13_fit", "c23_fit")) %>%
    group_by(model, term) %>%
    summarize(mean_est = mean(estimate), sd_est = sd(estimate)) %>%
    mutate(
      model = case_when(
        model == "c1_fit" ~ "Class 1 vs. All",
        model == "c2_fit" ~ "Class 2 vs. All",
        model == "c3_fit" ~ "Class 3 vs. All",
        model == "contrap_fit" ~ "C1 Placebo vs.\nC2 & C3 Active"
        # model == "c12_fit" ~ "Class 1 vs. Class 2"
      ),
      model = factor(model, levels = levels)
    ) %>%
    left_join(term_feature, by = "term")

  if (use_ref) {
    cv_coefs_long <- cv_coefs_long %>% rbind(reference)
  }
  cv_ceofs_long <- cv_coefs_long %>%
    mutate(
      feature = as.factor(feature),
      level = factor(level, levels = c(level_order$levels, ""))
    ) %>%
    arrange(model, feature, level)


  cv_coefs_long$level[is.na(cv_coefs_long$level)] <- ""

  if (use_site) {
    cv_coefs_long$site <- site
  }

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

# 
# 
# compare_metrics_class(metrics)
# ggsave("../reports/figures/ord/predict_class_compare.png")
# 
# cv_coef <- gather_cv_coef(imputed_fits, use_ref = FALSE)
# 
# 
# feature_imp(cv_coef$cv_coefs_long, use_ref = FALSE)
# # ggsave("../reports/figures/ord/explore/predict_class_features_contrast.png", dpi = 600)
# ggsave("../reports/figures/ord/predict_class_features.png", dpi = 600, height = 10)
# 
# 
# svg("../reports/figures/ord/predict_class_rocs.svg")
# roc_plot_class(rocs)
# dev.off()
# 
# 
# 
# write.csv(cv_coef$cv_coefs, "../reports/tables/ord/class_predic_coef.csv")
# 
# 
# 
# 
# write.csv(metrics, "../reports/tables/ord/class_predict_perf.csv", row.names = FALSE)

# write.csv(coefs, "../reports/tables/class_predict_ols_coefs.csv", row.names = FALSE)


# bysite

sites <- c("aus", "den", "ger", "swi", "usa")


coefs <- data.frame()


for (s in sites) {
  load(paste0("../models/predict_class_", s, ".rda"))
}


collect_metrics <- function(m, site) {
  metrics <- data.frame()
  rocs <- data.frame()
  for (f in 1:length(m)) {
    # enet metrics
    nme <- names(m[f])
    cv_auc_mean <- m[[f]]$cv$cv_metrics$mean[2]
    cv_auc_se <- m[[f]]$cv$cv_metrics$std_err[2]
    test_auc <- m[[f]]$test_score$.estimate
    null_auc_mean <- mean(m[[f]]$null_aucs)
    auc_pval <- m[[f]]$auc_pval

    tmp <- data.frame(nme, cv_auc_mean, cv_auc_se, test_auc, null_auc_mean, auc_pval)
    tmp$site <- site
    metrics <- rbind(metrics, tmp)

    roc_tmp <- m[[f]]$test_roc
    roc_tmp$model <- nme
    roc_tmp$site <- site
    rocs <- rbind(rocs, roc_tmp)
  }

  list(metrics, rocs)
}

aus_metrics <- collect_metrics(aus, "aus")
den_metrics <- collect_metrics(den, "den")
ger_metrics <- collect_metrics(ger, "ger")
swi_metrics <- collect_metrics(swi, "swi")
usa_metrics <- collect_metrics(usa, "usa")

site_metrics <- rbind(aus_metrics[[1]], den_metrics[[1]], ger_metrics[[1]], swi_metrics[[1]], usa_metrics[[1]])
site_rocs <- rbind(aus_metrics[[2]], den_metrics[[2]], ger_metrics[[2]], swi_metrics[[2]], usa_metrics[[2]])

site_metrics %>%
  filter(nme %in% c("c1_fit", "c2_fit", "c3_fit")) %>%
  ggplot(aes(x = site, y = cv_auc_mean, fill = nme)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(
    ymin = cv_auc_mean - cv_auc_se,
    ymax = cv_auc_mean + cv_auc_se
  ),
  position = position_dodge()) +
  labs(
    fill = "model",
    y = "Average CV AUC"
  )
ggsave("../reports/figures/ord/explore/site_metrics.png")

site_rocs %>%
  filter(model %in% c("c1_fit", "c2_fit", "c3_fit")) %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, color = model)) +
  geom_path() +
  geom_abline(linetype = "dotted") +
  coord_equal() +
  theme_bw(base_size = 15) +
  facet_wrap(~site)
ggsave("../reports/figures/ord/explore/site_roc.png")



aus_coef <- gather_cv_coef(aus, site = "aus", use_ref = FALSE, use_site = TRUE)
den_coef <- gather_cv_coef(den, site = "den", use_ref = FALSE, use_site = TRUE)
ger_coef <- gather_cv_coef(ger, site = "ger", use_ref = FALSE, use_site = TRUE)
swi_coef <- gather_cv_coef(swi, site = "swi", use_ref = FALSE, use_site = TRUE)
usa_coef <- gather_cv_coef(usa, site = "usa", use_ref = FALSE, use_site = TRUE)

site_coefs <- rbind(
  aus_coef$cv_coefs_long,
  den_coef$cv_coefs_long,
  ger_coef$cv_coefs_long,
  swi_coef$cv_coefs_long,
  usa_coef$cv_coefs_long
)


feature_imp_bysite(site_coefs)
ggsave("../reports/figures/ord/explore/site_features.png")





#





