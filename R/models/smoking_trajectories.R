library(lcmm)
library(dplyr)
library(ggplot2)

source("visualization/visualize.R")

theme_set(theme_classic(base_size = 15))
set.seed(42)

load("../data/processed/smoking_traj.rda")

cpd_hlme <- function(df, n_classes, B = NULL) {
  # specify "warm start" parameter estimates from the first model
  if (n_classes == 1) {
    hlme(
      prop_change ~ week,
      subject = "subject_id",
      data = df,
      var.time = "week",
      verbose = FALSE
    )
  } else if (n_classes > 1) {
    hlme(
      prop_change ~ week,
      mixture = ~ week,
      subject = "subject_id",
      data = df,
      ng = n_classes,
      B = B,
      var.time = "week",
      verbose = FALSE
    )
  }
}

fit_model <- function(df_lca) {
  m1 <- cpd_hlme(df_lca, n_classes = 1)
  m2 <- cpd_hlme(df_lca, n_classes = 2, B = m1)
  m3 <- cpd_hlme(df_lca, n_classes = 3, B = m1)
  m4 <- cpd_hlme(df_lca, n_classes = 4, B = m1)
  m5 <- cpd_hlme(df_lca, n_classes = 5, B = m1)
  m6 <- cpd_hlme(df_lca, n_classes = 6, B = m1)
  # m7 <- cpd_hlme(df_lca, n_classes = 7, B = m1)
  # m8 <- cpd_hlme(df_lca, n_classes = 8, B = m1)
  # m9 <- cpd_hlme(df_lca, n_classes = 9, B = m1)
  # # grid search
  
  # m4g <- gridsearch(
  #   hlme(
  #     prop_change ~ visit,
  #     mixture = ~ visit,
  #     subject = "subject_id",
  #     data = df_lca,
  #     ng = 4,
  #     verbose = FALSE
  #   ),
  #   rep = 100,
  #   maxiter = 30,
  #   minit = m1
  # )
  # 
  # m5g <- gridsearch(
  #   hlme(
  #     prop_change ~ visit,
  #     mixture = ~ visit,
  #     subject = "subject_id",
  #     data = df_lca,
  #     ng = 5,
  #     verbose = FALSE
  #   ),
  #   rep = 100,
  #   maxiter = 30,
  #   minit = m1
  # )
  # 

 list(m1, m2, m3, m4, m5, m6)
}


goodness_of_fit_df <- function(mod_list) {
  df <- data.frame("BIC" = NA, "AIC" = NA, "class" = NA)
  
  for (i in 1:length(mod_list)) {
    tmp <- data.frame(
      "BIC" = mod_list[[i]]$BIC,
      "AIC" = mod_list[[i]]$AIC,
      "class" = i)
    
    df <- rbind(df, tmp)
    
  }
  return(df %>% tidyr::drop_na())
}


merge_class <- function(model, df_lca) {
  class_probs <- as.data.frame(summarytable(model))
  df <- model$pprob %>% left_join(df_lca, by = "subject_id")
  
  list(
    "class_probs" = class_probs,
    "df" = df
  )
}


# main ----

# CPD distributions at each time point
tpt <- time_point_normality(traj)
ggsave(plot = tpt, file = "../reports/figures/lca/cpd_distributions.png", width = 7, height = 5, dpi=600)

# find smoking trajectories
cpd_lca <- fit_model(traj)
save(cpd_lca, file = "../models/lca_models.RData")

# goodness of fit
bic_df <- goodness_of_fit_df(cpd_lca)
bic <- plot_bic(bic_df)
ggsave(plot = bic, file = "../reports/figures/lca/lmm_bic.png", width = 7, height = 5, dpi=600)


fit_table <- data.frame()  
for (m in cpd_lca) {
    nme <- paste0("lca", m$ng)
    
    # print results
    sink(file = paste0("../reports/tables/lca_results/", nme, ".txt"))
    print(summary(m))
    sink()
    
    # assemble table
    tab <- data.frame(m$ng, m$BIC, m$AIC, m$loglik)
    fit_table <- rbind(fit_table, tab)
    # write class predictions
    model_traj <- merge_class(m, traj)
    model_traj_df <- model_traj$df
    save(model_traj_df, file = paste0("../data/predicted/", nme, "_predict.rda"))
    
    # generate trajectory plots
    plot_traj <- plot_trajectories(df_lca, model_traj, n_classes = m$ng)
    # ggsave(plot = plot_traj$plot, file = paste0("../reports/figures/lca/", nme, ".png"), width = 7, height = 5, dpi='retina')
    ggsave(plot = plot_traj$plot, file = paste0("../reports/figures/lca/", nme, ".pdf"), width = 7, height = 6)
    
    # generate cpd plot
    plot_cpd_traj <- plot_cpd_trajectories(df_lca, model_traj, n_classes = m$ng)
    ggsave(plot = plot_cpd_traj, file = paste0("../reports/figures/lca/cpd/", nme, ".png"), width = 7, height = 5, dpi=600)
    
}

write.csv(fit_table, "../reports/tables/lca_results/lca_fit.csv", row.names = FALSE)











# 
# traj$latent_sum %>%
#   distinct(class, week, avg_change, sd_change)






