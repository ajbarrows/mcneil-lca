time_point_normality <- function(df_lca, n_values) {
  
  n_values <- df_lca %>%
    group_by(week) %>%
    filter(!is.na(prop_change)) %>%
    count()
  
  n_vector <- paste0("n=", n_values$n)
  
  tmpts <- c("Week 2", "Week 10", "Week 18", "Week 26")

  dat_text <- data.frame(label = n_vector, week = factor(tmpts, levels = tmpts))
  
  p <- df_lca %>%
    mutate(week = factor(
      paste("Week", week),
      levels = tmpts),
      pct_change = (prop_change) * 100
    ) %>%
    
    ggplot(aes(x = pct_change, color = week, fill = week)) +
    # geom_density(alpha = 0.1) +
    geom_histogram(aes(y = stat(count / sum(count) * 100 )), show.legend = FALSE) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    facet_wrap(~week) +
    # annotate(geom = "text", x = 50, y = 400, label = "No Change") +
    geom_text(
      data = dat_text,
      aes(x = -Inf, y = -Inf, label = label),
      hjust = -0.5,
      vjust = -10,
      show.legend = FALSE) +
    labs(
      x = "Percent Change in Baseline CPD",
      y = "Percent of Subjects",
      color = "Trial\nFollow-Up Point",
      fill = "Trial\nFollow-Up Point",
      title = "Distribution of Changes in CPD at Each Follow-Up Point"
    )
  
  return(p)
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


plot_bic <- function(bic_df) {
  bic_df %>%
    mutate(class = factor(class)) %>%
    ggplot(aes(x = class, y = BIC)) +
    geom_point() +
    geom_line(group = 1) +
    labs(
      x = "Class",
      title = "Linear Mixed Model Goodness of Fit"
    )
}

plot_trajectories <- function(df_lca, class_probs, n_classes) {
  no_change <- 0
  scale_factor <- 100
  
  latent_df <- class_probs$class_probs %>% 
    select(starts_with("%")) %>%
    tidyr::pivot_longer(everything()) %>%
    mutate(
      class = as.numeric(substr(name, 7, 7)),
      value = paste(round(value, 2), "%", sep = "")
    ) %>%
    select(class, value) %>%
    left_join(class_probs$df, by = "class")
  
  latent_df <- latent_df %>%
    group_by(class) %>%
    distinct(subject_id) %>%
    count() %>%
    left_join(latent_df, by = "class") %>%
    mutate(
      week = factor(week),
      class = paste(class, ": ", value, " (n = ", n, ")", sep = ""),
      class = factor(class)
    )
  
  latent_sum <- latent_df %>%
    group_by(class, week) %>%
    summarize(
      pct_change = prop_change * scale_factor,
      avg_change = mean(pct_change, na.rm = TRUE),
      sd_change = sd(pct_change, na.rm = TRUE)
    )
  
  classes <- latent_sum %>%
    ungroup() %>%
    distinct(class)
  
  
  fill_values <- c(
   "#F8766D", 
   "#00BA38", 
   "#00A9FF",
   "#00BFC4",
   "#C77CFF",
   "#FF61CC"
  )
  
  rand_subjects <- latent_df %>% distinct(subject_id) %>% sample_frac(0.1)
  variance_plot <- latent_df %>% filter(subject_id %in% rand_subjects$subject_id)
  
  fill_values <- setNames(fill_values, classes$class)
  
  
  p <- latent_sum %>%
    ggplot(
      aes(x = week, y = avg_change, color = class, group = class
      )
    ) +
    geom_point() +
    geom_line(linewidth = 1) +
    theme(legend.position = "top", plot.title = element_text(size=15))+
    geom_hline(aes(yintercept = no_change), linetype = "dashed") +
    annotate("text", x = 1, y = no_change + 5, label = "No Change") +
    labs(
      x = "Follow-Up Week",
      y = "Average Percent Change from Baseline CPD",
      title = "Avg. Smoking Trajectories by Predicted Latent Class",
      color = "Class"
    )  +
    scale_color_manual(values = fill_values)
  p <- p + 
    geom_errorbar(
      aes(ymin = avg_change - sd_change, ymax = avg_change + sd_change),
      width = 0.1, size = 1, alpha = .05
    ) 
    # geom_line(
    #   data = latent_df,
    #   aes(x = week, y = prop_change * scale_factor, group = subject_id),
    #   alpha = 0.09
    # )
  
  list(
    "plot" = p,
    "latent_sum" = latent_sum
  )
}


plot_cpd_trajectories <- function(df_lca, class_probs, n_classes) {
  no_change <- 0
  scale_factor <- 100
  
  latent_df <- class_probs$class_probs %>% 
    select(starts_with("%")) %>%
    tidyr::pivot_longer(everything()) %>%
    mutate(
      class = as.numeric(substr(name, 7, 7)),
      value = paste(round(value, 2), "%", sep = "")
    ) %>%
    select(class, value) %>%
    left_join(class_probs$df, by = "class")
  
  latent_df <- latent_df %>%
    group_by(class) %>%
    distinct(subject_id) %>%
    count() %>%
    left_join(latent_df, by = "class") %>%
    mutate(
      week = factor(week),
      class = paste(class, ": ", value, " (n = ", n, ")", sep = ""),
      class = factor(class)
    )
  
  latent_sum <- latent_df %>%
    group_by(class, week) %>%
    summarize(
      pct_change = prop_change * scale_factor,
      avg_change = mean(pct_change, na.rm = TRUE),
      sd_change = sd(pct_change, na.rm = TRUE)
    )
  
  classes <- latent_sum %>%
    ungroup() %>%
    distinct(class)
  
  
  fill_values <- c(
    "#F8766D", 
    "#00BA38", 
    "#00A9FF",
    "#00BFC4",
    "#C77CFF",
    "#FF61CC"
  )
  
  rand_subjects <- latent_df %>% distinct(subject_id) %>% sample_frac(0.1)
  variance_plot <- latent_df %>% filter(subject_id %in% rand_subjects$subject_id)
  
  fill_values <- setNames(fill_values, classes$class)
  
  
  p <- latent_df %>%
    ggplot(aes(x = week, y = avg_cpd, color = class, group = subject_id)) +
    geom_line()
  
  # p <- latent_sum %>%
  #   ggplot(
  #     aes(x = week, y = avg_change, color = class, group = class
  #     )
  #   ) +
  #   geom_point() +
  #   geom_line(linewidth = 1) +
  #   geom_hline(aes(yintercept = no_change), linetype = "dashed") +
  #   annotate("text", x = 1, y = no_change + 5, label = "No Change") +
  #   labs(
  #     x = "Follow-Up Week",
  #     y = "Average Percent Change from Baseline CPD",
  #     title = "Average Smoking Trajectories by Predicted Latent Class",
  #     color = "Class"
  #   )  +
  #   scale_color_manual(values = fill_values)
  # p <- p + 
  #   # geom_errorbar(
  #   #   aes(ymin = avg_change - sd_change, ymax = avg_change + sd_change),
  #   #   width = 0.1, size = 1, alpha = .05
  #   # ) +
  #   # geom_line(
  #   #   data = latent_df,
  #   #   aes(x = week, y = prop_change * scale_factor, group = subject_id),
  #   #   alpha = 0.09
  #   # )
  # 
  # list(
  #   "plot" = p,
  #   "latent_sum" = latent_sum
  # )
  return(p)
}

plot_coef_byclass <- function(plt_df, threshold = 0, facet = TRUE, scales = "fixed") {
  
  fill_values <- c(
    "class_1" = "#F8766D", 
    "class_2" = "#00BA38", 
    "class_3" = "#619CFF"
  )
  
  
  plt_df <- plt_df %>%
    filter(abs(estimate) >= threshold) %>%
    filter(term != "(Intercept)") %>%
    filter(term != 'week26_cpd')
  
  # bygroup
  
  plt_df <- plt_df %>%
    mutate(grp = stringr::str_extract(term, "^(.*?\\_)")) %>%
    mutate(grp = ifelse(is.na(grp), term, grp)) %>%
    group_by(grp) %>%
    mutate(avg_est = mean(abs(estimate)))
  
  p <- ggplot(
    plt_df,
    aes(
      y = reorder(term, -avg_est),
      # y = term,
      x = estimate,
      fill = class)) +
    geom_col() +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    labs(
      x = "Beta",
      y = ""
    )
  
  if(facet) {
    p <- p + facet_wrap(~class, scales = scales)
  }
  
  p
}


compare_metrics <- function(metrics) {
  
  # levels <- c("full_fit", "noclass_fit", "c1_fit", "c2_fit", "c3_fit")
  
  rename_list <- c(
    "full_fit" = "Predictors + Class",
    "noclass_fit" = "Predictors Alone",
    "onlyclass_fit" = "Latent Class Alone",
    "c1_fit" = "Within Class 1",
    "c2_fit" = "Within Class 2",
    "c3_fit" = "Within Class 3"
  )
  
  df_plt <- metrics %>%
    mutate(nme = stringr::str_replace_all(nme, rename_list)) %>%
    tidyr::pivot_longer(
      c(cv_roc_mean, test_auc)
    ) %>%
    mutate(
      cv_roc_se = ifelse(name == "test_auc", NA, cv_roc_se),
      nme = factor(nme, levels = rename_list)
    )
  
  
  ggplot(df_plt, aes(x = nme, y = value, fill = name)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(
      ymin = value - cv_roc_se,
      ymax = value + cv_roc_se
    ),
    position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "model",
      y = "ROC AUC",
      fill = "",
      title = "Predicting Smoking Cessation at 1-year follow up"
    )
  
}

compare_metrics_class <- function(metrics) {
  
  rename_list <- c(
    "c1_fit" = "Class 1",
    "c2_fit" = "Class 2",
    "c3_fit" = "Class 3",
    "c12_fit" = "Class 1:2 Comp.",
    "c13_fit" = "Class 1:3 Comp.",
    "c23_fit" = "Class 2:3 Comp."
  )
  
  df_plt <- metrics %>%
    mutate(nme = stringr::str_replace_all(nme, rename_list)) %>%
    tidyr::pivot_longer(
      c(cv_auc_mean, test_auc)
    ) %>%
    mutate(
      cv_auc_se = ifelse(name == "test_auc", NA, cv_auc_se),
      nme = factor(nme, levels = rename_list)
      )
  
  ggplot(df_plt, aes(x = nme, y = value, fill = name)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(
      ymin = value - cv_auc_se,
      ymax = value + cv_auc_se
    ),
    position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "model",
      y = "ROC AUC",
      fill = "",
      title = "Predicting Class Membership"
    )
  
}



quit_coef <- function(full_fit) {
  
  plt_df <- full_fit$coef %>%
    filter(term != "(Intercept)") %>%
    mutate(grp = stringr::str_extract(term, "^(.*?\\_)")) %>%
    mutate(grp = ifelse(is.na(grp), term, grp)) %>%
    group_by(grp) %>%
    mutate(avg_est = mean(abs(estimate)))
  
  ggplot(plt_df, aes(y = reorder(term, -avg_est), x = estimate)) +
    geom_col() +
    labs(
      x = "Beta",
      y = "",
      title = "Predicting CO Values at 1-year follow up"
    )
}

count_quit <- function(full_fit, quit_threshold = 6) {
  
  # quit_thresh <- 11
  
  quit_df <- full_fit$test_preds %>%
    select(quit_verified, .pred_class) %>%
    mutate(
      quit_true = ifelse(quit_verified == 1, "yes", "no"),
      # quit_pred = ifelse(.pred_class == "", "yes", "no"),
      success = quit_true == .pred_class
    )
  # 
  quit_table <- quit_df %>%
    group_by(quit_verified) %>%
    count(success)
  # 
  # p <- ggplot(quit_df, aes(x = co_1year, y = .pred, color = success)) +
  #   geom_point() +
  #   geom_vline(aes(xintercept = quit_thresh), linetype = "dashed") +
  #   geom_hline(aes(yintercept = quit_thresh), linetype = "dashed") +
  #   labs(
  #     x = "Actual 1-year CO (ppm)",
  #     y = "Predicted 1-year CO (ppm)",
  #     title  = "Predicting 1-Year CO Using Latent Class + Baseline Characteristics"
  #   )
  # 
  # list(
  #   "plot" = p,
  #   "df" = quit_df
  # )
  
  quit_table
}

feature_imp <- function(cv_coefs_long, use_ref = TRUE, dashline = 1) {
  fill_values <- c(
    "Class 1 vs. All" = "#F8766D", 
    "Class 1 vs. All (Placebo NRT Only)" = "grey",
    "Class 2 vs. All" = "#00BA38", 
    "Class 3 vs. All" = "#619CFF",
    "Class 1 vs. Class 2" = "#C77CFF"
  )
  
  plt_df <- cv_coefs_long %>% 
    filter(term != "(Intercept)" | is.na(term)) %>%
    filter(!is.na(feature)) %>%
    mutate(
      model=forcats::fct_recode(
        model, 
        "Class 1 vs. All\n(Placebo NRT Only)"="Class 1 vs. All (Placebo NRT Only)")
    )

  plt <- plt_df %>%
    ggplot(aes(y = level, x = mean_est, color = model)) +
    geom_pointrange(
      aes(xmin = mean_est - sd_est, xmax = mean_est + sd_est),
      show.legend = FALSE
    ) 
  
  if (use_ref) {
    plt <- plt + 
      geom_point(
        data = plt_df %>% filter(is.na(term)),
        color = "black",
        size = 2.75,
        shape = 15)
  }
  
  plt <- plt + 
    geom_vline(aes(xintercept = dashline), linetype = "dashed") +
    ggh4x::facet_nested(
      feature ~ model, 
      scales = "free_y", 
      space = "free_y",
      switch = "y",
      drop = TRUE) +
    theme_bw() +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0)) +
    scale_color_manual(values = fill_values) +
    labs(
      y = "", 
      x = "Average Feature Importance",
      color = ""
    ) +
    xlim(c(-1, 3))
    # scale_x_continuous(breaks=scales::pretty_breaks())
  
  plt
  
}

feature_imp_bysite <- function(cv_coefs_long, use_ref = TRUE) {
  plt_df <- cv_coefs_long %>% 
    filter(term != "(Intercept)" | is.na(term)) %>%
    filter(!is.na(feature))
  
  plt <- plt_df %>%
    ggplot(aes(y = level, x = mean_est, color = model)) +
    geom_pointrange(
      aes(xmin = mean_est - sd_est, xmax = mean_est + sd_est),
      show.legend = FALSE
    ) 
  
  if (use_ref) {
    plt <- plt + 
      geom_point(
        data = plt_df %>% filter(is.na(term)),
        color = "black",
        size = 2.75,
        shape = 15)
  }
  
  plt <- plt + 
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    ggh4x::facet_nested(
      feature ~ site, 
      scales = "free_y", 
      space = "free_y",
      switch = "y",
      drop = TRUE) +
    theme_bw(base_size = 15) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0)) +
    labs(
      y = "", 
      x = "Average Feature Importance",
      color = ""
    ) 

  
}


roc_plot_class <- function(rocs) {

  c1 <- rocs %>% filter(model == "c1_fit")
  c2 <- rocs %>% filter(model == "c2_fit")
  c3 <- rocs %>% filter(model == "c3_fit")

  
  plot(
    1 - c1$specificity, c1$sensitivity, 
    type = "l",
    col = "#F8766D",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    asp = 1
  )
  lines(1 - c2$specificity, c2$sensitivity, type = "l", col = "#00BA38")
  lines(1 - c3$specificity, c3$sensitivity, type = "l", col = "#00A9FF")
  abline(coef = c(0, 1), lty = 2)
  legend(
    "bottomright",
    legend = c(
      "Class 1",
      "Class 2",
      "Class 3"
    ),
    col = c("#F8766D", "#00BA38", "#00A9FF"),
    lty = 1
  )
  
  title("Predicting Latent Class Membership")
  
}
# 
roc_plot_quit <- function(roc1, roc2) {


  plot(
    1 - roc1$specificity, roc1$sensitivity,
    type = "S",
    col = "red",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    asp = 1
  )
  lines(1 - roc2$specificity, roc2$sensitivity, type = "S", col = "blue")
  abline(coef = c(0, 1), lty = 2)

  legend(
    "bottomright",
    legend = c(
      "Baseline Characteristics & Latent Class",
      "Baseline Characteristics Alone"
    ),
    col = c("red", "blue"),
    lty = 1
  )
  title("Predicting Smoking Cessation")
}




