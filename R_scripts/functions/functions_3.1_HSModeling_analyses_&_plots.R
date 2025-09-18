# Functions of model predictions and evaluations analyses and plots 

# Function to apply same scaling
scale_with_params <- function(df, params) {
  df %>%
    mutate(across(where(is.numeric),
                  ~ (.-params[[paste0(cur_column(), "_mean")]]) /
                    params[[paste0(cur_column(), "_sd")]]))
}


### Function to get variable importance and plot the results ===================
get_var_importance <- function(model, type) {
  if(type == "GLM") type <- "glm"
  if(type == "GAM") type <- "gam"
  if(type == "RF")  type <- "randomForest"
  if(type == "BRT") type <- "gbm" 
  
  if (type == "glm") {
    coefs <- coef(model)
    if ("(Intercept)" %in% names(coefs)) coefs <- coefs[-which(names(coefs) == "(Intercept)")]
    imp <- abs(coefs)
    df <- data.frame(
      Variable = names(imp),
      Importance = as.numeric(imp)
    )
  } else if (type == "gam") {
    # Extract summary of smooth terms
    s <- summary(model)
    # The "s.table" element has the smooth terms and their significance
    # We'll get the edf (effective degrees of freedom) or chi.sq as proxy of importance
    # Here, I'll use the absolute t-values or edf as importance (you can decide)
    sm_terms <- s$s.table
    if (is.null(sm_terms)) {
      stop("No smooth terms found in GAM summary")
    }
    # Clean variable names by removing s(...) wrapper
    vars <- rownames(sm_terms)
    vars_clean <- gsub("s\\((.*)\\)", "\\1", vars)
    
    # Use edf as importance (or use abs(t value) - pick what fits your goal)
    importance <- sm_terms[, "edf"]
    
    df <- data.frame(
      Variable = vars_clean,
      Importance = as.numeric(importance)
    )
  } else if (type == "rf") {
    imp <- randomForest::importance(model, type = 2) # MeanDecreaseGini
    df <- data.frame(
      Variable = rownames(imp),
      Importance = imp[, 1]
    )
  } else if (type == "brt") {
    imp <- summary(model, plotit = FALSE)
    df <- data.frame(
      Variable = imp$var,
      Importance = imp$rel.inf
    )
  } else {
    stop("Unknown model type")
  }
  
  df$Variable <- as.character(df$Variable)
  df$Importance <- as.numeric(df$Importance)
  
  df %>% arrange(desc(Importance))
  
  
}


plot_var_importance <- function(models, types, top_n = 10) {
  # models: list of model objects
  # types: vector of model types ("glm", "gam", "rf", "brt")
  # top_n: how many top variables to show per model
  scale_importance <- function(imp) {
    imp / sum(imp) * 100
  }
  
  all_imp <- lapply(seq_along(models), function(i) {
    imp <- get_var_importance(models[[i]], types[i])
    imp$Model <- toupper(types[i])
    imp
    imp$Importance_std<- scale_importance(imp$Importance)
    return(imp) 
  }) %>% bind_rows()
  
  # Keep top_n variables per model
  top_vars <- all_imp %>%
    group_by(Model) %>%
    top_n(top_n, wt = Importance_std) %>%
    ungroup() %>%
    distinct(Variable)
  
  # Filter all_imp to keep only these variables for consistent comparison
  plot_df <- all_imp %>%
    filter(Variable %in% top_vars$Variable) %>%
    mutate(Variable = factor(Variable, levels = rev(unique(Variable))))
  
  ggplot(plot_df, aes(x = Variable, y = Importance_std, fill = Model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() +
    labs(title = "Variable Importance Comparison",
         y = "Importance",
         x = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


plot_var_importance_v2 <- function(models, types, top_n = 10, stacked = FALSE) {
  # models: list of model objects
  # types: vector of model types ("glm", "gam", "rf", "brt")
  # top_n: number of top variables to show
  # stacked: if TRUE, bars are stacked; if FALSE, side-by-side
  
  scale_importance <- function(imp) {
    imp / sum(imp) * 100
  }
  
  # 1. Compute importance for each model
  all_imp <- lapply(seq_along(models), function(i) {
    imp <- get_var_importance(models[[i]], types[i])
    
    # Aggregate GLM/GAM quadratic terms
    if (types[i] %in% c("glm", "gam")) {
      imp$Variable <- gsub("^I\\((.*)\\^2\\)$", "\\1", imp$Variable)  # X^2 → X
      imp <- imp %>%
        group_by(Variable) %>%
        summarise(Importance = sum(Importance)) %>%
        ungroup()
    }
    
    imp$Importance_std <- scale_importance(imp$Importance)
    imp$Model <- toupper(types[i])
    return(imp)
  }) %>% bind_rows()
  
  # 2. Determine top_n variables across all models
  top_vars <- all_imp %>%
    group_by(Variable) %>%
    summarise(MaxImp = max(Importance_std)) %>%
    top_n(top_n, wt = MaxImp) %>%
    pull(Variable)
  
  # 3. Filter and order for plotting
  plot_df <- all_imp %>%
    filter(Variable %in% top_vars) %>%
    mutate(Variable = factor(Variable, levels = rev(top_vars)))
  
  # 4. Plot
  pos <- if (stacked) "stack" else position_dodge(width = 0.8)
  
  p<- ggplot(plot_df, aes(x = Variable, y = Importance_std, fill = Model)) +
    geom_col(position = pos, width = 0.7) +
    coord_flip() +
    labs(title = "Variable Importance Comparison",
         y = "Importance (%)",
         x = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(p, plot_df))
}



#### Function to get ROC curve data -------------------------------------------

get_roc_data <- function(model= glm_mod_s, test= test_sel_scaled, type="GLM") { #, quad_vars= quad_vars_sel
  type <- tolower(type)  # make case-insensitive
  
  if (type %in% c("glm", "gam")) {
    # test_scaled <- test %>%
    #   mutate(across(where(is.numeric), scale))
    test_glm<- test
    for (v in quad_vars) {
      test_glm[[paste0("I(", v, "^2)")]] <- test_glm[[v]]^2
    }
    # test_glm<- test_glm[,-1]
    pred <- predict(model, newdata = test_glm, type = "response")
  } else if (type == "rf") {
    pred <- predict(model, newdata = test, type = "prob")[, 2]
  } else if (type == "brt") {
    pred <- predict(model, newdata = test, n.trees = model$n.trees, type = "response")
  } else {
    stop("Unknown model type: ", type)
  }
  
  roc_obj <- roc(test$kelp, pred)
  coords_df <- coords(roc_obj, ret = c("specificity", "sensitivity"), transpose = FALSE)
  
  data.frame(
    Specificity = coords_df$specificity,
    Sensitivity = coords_df$sensitivity,
    Model = toupper(type)
  )
}
### # Variables importance RANKING Select the most important variables =========
# # Step 1: Keep variables important in >= 3 models
# # Variables importance scaled based on each model values
# 
# all_imp <- lapply(seq_along(models), function(i) {
#   imp <- get_var_importance(models[[i]], types[i])
#   imp$Model <- toupper(types[i])
#   imp$Importance_std <- scale_importance(imp$Importance)
#   imp
# }) %>% bind_rows()
# 
# 
# # Summarize across models
# var_summary <- all_imp %>%
#   group_by(Variable) %>%
#   summarise(
#     n_models = n(),
#     min_importance = min(Importance_std, na.rm = TRUE),
#     avg_importance = mean(Importance_std, na.rm = TRUE),
#     max_importance = max(Importance_std, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   filter(n_models >= 3) %>%        # Keep vars in ≥ 3 models
#   arrange(desc(avg_importance)) %>%
#   mutate(
#     Rank = row_number(),
#     Category = cut(
#       avg_importance,
#       breaks = seq(0, max(avg_importance, na.rm = TRUE) + 2, by = 2),
#       include.lowest = TRUE,
#       right = FALSE
#     )
#   )
# 
# print(var_summary, n = Inf)
# 
# 
# 
# ## Another approach =====
# # So you want to rescale the importance of variables within each model so that 
# # the maximum importance per model is 100%, then combine the models to summarize
# # importance across models and rank variables into categories (like the 2% bins 
# # we were trying before). 
# 
# 
# # Variables importance standardized based on each model maximum value (0-100)
# # Step 1: Standardize importance per model relative to the maximum in that model
# all_imp_std <- all_imp %>%
#   group_by(Model) %>%
#   mutate(Importance_std = Importance / max(Importance) * 100) %>%
#   ungroup()
# 
# # Step 2: Summarize across models
# var_summary <- all_imp_std %>%
#   group_by(Variable) %>%
#   summarise(
#     n_models = n(),
#     min_importance = min(Importance_std),
#     avg_importance = mean(Importance_std),
#     max_importance = max(Importance_std)
#   ) %>%
#   arrange(desc(avg_importance))
# 
# # Step 3: Rank variables by avg_importance
# var_summary <- var_summary %>%
#   mutate(
#     Rank = row_number(),
#     Category = cut(
#       avg_importance,
#       breaks = seq(0, 100, by = 2),   # 2% intervals
#       include.lowest = TRUE,
#       labels = FALSE
#     )
#   )


## The current speed is still not in the selected variables so we want to explore if considering the variability among models we get a more realistic model structure
# Explanation:
# Standard deviation of the standardized importance (Importance_std) across models for each variable.
# Importance_std = Importance / max(Importance) * 100 rescales each model’s variable importance to a 0–100% scale.
# Then you summarize min, mean, max importance across models.
# Importance_SD → variability of importance across models.
# Weighted_Score → high if a variable is both important and consistent.
# Keep → safeguards against variables that are consistently low.
# Sorting by Weighted_Score gives the final ranking.


rank_variables <- function(models, types, avg_threshold = NULL) {
  
  # 1. Get raw importance per model
  all_imp <- lapply(seq_along(models), function(i) {
    imp <- get_var_importance(models[[i]], types[i])
    
    # Aggregate linear + quadratic terms for GLM/GAM
    imp$Variable <- gsub("^I\\((.*)\\^2\\)$", "\\1", imp$Variable)
    imp <- imp %>%
      group_by(Variable) %>%
      summarise(Importance = sum(Importance, na.rm = TRUE), .groups = "drop")
    
    # Standardize per model so top variable = 100
    imp$Importance_std <- (imp$Importance / max(imp$Importance, na.rm = TRUE)) * 100
    imp$Model <- toupper(types[i])
    imp
  }) %>% bind_rows()
  
  # 2. Summarize across models
  var_summary <- all_imp %>%
    group_by(Variable) %>%
    summarise(
      n_models       = n(),
      min_importance = min(Importance_std, na.rm = TRUE),
      avg_importance = mean(Importance_std, na.rm = TRUE),
      max_importance = max(Importance_std, na.rm = TRUE),
      sd_importance  = sd(Importance_std, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3. Weighted score with min-importance safeguard
  overall_sd <- mean(var_summary$sd_importance, na.rm = TRUE)
  var_summary <- var_summary %>%
    mutate(
      consistency_factor = 1 / (1 + (sd_importance / overall_sd)),
      weighted_score     = avg_importance * consistency_factor,
      weighted_score     = ifelse(min_importance < 5,
                                  weighted_score * 0.8,
                                  weighted_score)
    ) %>%
    arrange(desc(weighted_score)) %>%
    mutate(
      Rank = row_number(),
      Category = cut(
        avg_importance,
        breaks = seq(0, 100, by = 2),
        include.lowest = TRUE,
        labels = paste0(seq(0, 98, by = 2), "-", seq(2, 100, by = 2), "%")
      )
    )
  
  # 4. Apply threshold filter if specified
  if (!is.null(avg_threshold)) {
    var_summary <- var_summary %>%
      filter(weighted_score > avg_threshold)
  }
  
  return(list(
    summary_table = var_summary,
    raw_importance = all_imp
  ))
}


plot_variable_ranking <- function(summary_table, threshold) {
  if (!is.null(threshold)) {
    summary_table <- summary_table %>%
    mutate(
      Variable = forcats::fct_reorder(Variable, weighted_score),
      ColorFlag = ifelse(is.na(weighted_score), "Below",
                         ifelse(weighted_score >= threshold, "Above", "Below")))
  } else{
    summary_table <- summary_table %>%
      mutate(
        Variable = forcats::fct_reorder(Variable, weighted_score),
        ColorFlag = "Above")
  }
  
  
  p<- ggplot(summary_table, aes(x = weighted_score, y = Variable, fill = ColorFlag)) +
    geom_col() +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(
      values = c("Above" = "steelblue", "Below" = "grey70"),
      guide = "none" ) +
    labs(
      x = "Weighted Score",
      y = NULL,
      title = "Variable Importance Ranking with Threshold"
    ) +
    theme(#legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 7),   # facet label size
          axis.title = element_text(size = 9),                 # axis title size
          axis.text = element_text(size = 7),                  # axis tick labels
          plot.title = element_text(size = 9, face = "bold")   # plot title
    )+
    theme_bw()
  
  return(list(summary_table, p))
}



# Plot Variables Importance Profile (ranking of variables based on average importance and weighted score)
plot_importance_profiles <- function(summary_table, threshold = NULL) {
  
  summary_table <- summary_table %>%
    mutate(Variable = forcats::fct_reorder(Variable, weighted_score,.desc = TRUE),
           x_num = as.numeric(Variable))
  
  ggplot(summary_table, aes(x = x_num)) +
    # Error bars from min to max importance
    geom_errorbar(
      aes(ymin = min_importance, ymax = max_importance),
      width = 0.5, color = "grey70", alpha = 0.8
    ) +
    # Average importance line and points
    geom_line(aes(y = avg_importance, color = "Average Importance"), size = 1) +
    geom_point(aes(y = avg_importance, color = "Average Importance"), size = 2) +
    
    {if (!is.null(threshold)) geom_hline(aes(yintercept = threshold, linetype = "Threshold"),
               color = "darkgreen", linewidth = 1) } +
    scale_color_manual(values = c("Average Importance" = "blue", 
                                  "Weighted Score" = "darkgreen")) +
    scale_linetype_manual(values = c("Threshold" = "dashed")) +
    # Weighted score line and points
    # {if (!is.null(threshold)) 
    # geom_line(aes(y = weighted_score, color = "Weighted Score"), size = 1) +
    # geom_point(aes(y = weighted_score, color = "Weighted Score"), size = 2)  } else{
      geom_line(aes(y = weighted_score, color = "blue"), size = 1) +
        geom_point(aes(y = weighted_score, color = "blue"), size = 2) +
  
    # Optional threshold
    # {if (!is.null(threshold)) geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgreen")} +
    scale_x_continuous(
      breaks = summary_table$x_num,
      labels = summary_table$Variable
    ) +
    scale_color_manual(values = c("Average Importance" = "blue", "Weighted Score" = "red")) +
    labs(
      x = NULL,
      y = "Importance Value",
      title = "Variable Importance Profiles",
      subtitle = "Gray error bars = min to max importance per variable",
      color = "",
      linetype=""
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = c(0.97, 0.97),
          legend.justification = c("right", "top"),
          # axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 8),
          axis.text.x = element_text(size = 8,angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          axis.title = element_text(size = 9),
          legend.key.size = unit(0.6, "lines"))
}



#-----------------------------
# 1. Bar plot of weighted scores with threshold
#-----------------------------
plot_variable_ranking_v2 <- function(summary_table, threshold = NULL) {
  
  summary_table <- summary_table %>%
    mutate(
      Variable = forcats::fct_reorder(Variable, weighted_score),
      ColorFlag = case_when(
        is.na(weighted_score) ~ "Below",
        !is.null(threshold) & weighted_score >= threshold ~ "Above",
        TRUE ~ "Below"
      )
    )
  
  ggplot(summary_table, aes(x = Variable, y = weighted_score, fill = ColorFlag)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = c("Above" = "steelblue", "Below" = "grey70"), guide = "none") +
    coord_flip() +
    labs(title = "Variable Importance Ranking with Threshold",
         y = "Weighted Score",
         x = NULL) +
    theme_minimal(base_size = 14)
}

#-----------------------------
# 2. Profile plot of weighted vs avg importance
#-----------------------------
plot_importance_profiles_v2 <- function(summary_table, threshold = NULL) {
  
  summary_table <- summary_table %>%
    mutate(
      Variable = forcats::fct_reorder(Variable, weighted_score, .desc = TRUE),
      x_num = as.numeric(Variable),
      ColorFlag = case_when(
        is.na(weighted_score) ~ "Below",
        !is.null(threshold) & weighted_score >= threshold ~ "Above",
        TRUE ~ "Below"
      )
    )
  
  ggplot(summary_table, aes(x = x_num)) +
    geom_errorbar(aes(ymin = min_importance, ymax = max_importance), width = 0.5, color = "grey70", alpha = 0.8) +
    geom_line(aes(y = avg_importance, color = "Average Importance"), size = 1) +
    geom_point(aes(y = avg_importance, color = "Average Importance"), size = 2) +
    geom_line(aes(y = weighted_score, color = "Weighted Score"), size = 1) +
    geom_point(aes(y = weighted_score, color = "Weighted Score"), size = 2) +
    {if (!is.null(threshold)) geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgreen")} +
    scale_x_continuous(breaks = summary_table$x_num, labels = summary_table$Variable) +
    scale_color_manual(values = c("Average Importance" = "blue", "Weighted Score" = "steelblue")) +
    labs(x = NULL,
         y = "Importance Value",
         title = "Variable Importance Profiles",
         subtitle = "Gray error bars = min to max importance per variable",
         color = "") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = c(0.97, 0.97),
      legend.justification = c("right", "top"),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 9),
      legend.key.size = unit(0.6, "lines")
    )
}












### Get metrics of models performance for model evaluation =====================
# Function to compute metrics
# get_metrics <- function(model, test, type = c("glm", "gam", "rf", "brt")) {
#   if (type %in% c("glm", "gam")) {
#     pred <- predict(model, newdata = test, type = "response")
#   } else if (type == "rf") {
#     pred <- predict(model, newdata = test, type = "prob")[, 2]
#   } else if (type == "brt") {
#     pred <- predict(model, newdata = test, n.trees = model$n.trees, type = "response")
#   }
#   
#   roc_obj <- pROC::roc(test$kelp, pred)
#   auc_val <- as.numeric(pROC::auc(roc_obj))
#   ci_val <- pROC::ci.auc(roc_obj)
#   
#   pred_class <- ifelse(pred > 0.5, 1, 0)
#   cm <- caret::confusionMatrix(as.factor(pred_class), as.factor(test$kelp), positive = "1")
#   
#   data.frame(
#     Model = type,
#     AUC = auc_val,
#     AUC_lower = ci_val[1],
#     AUC_upper = ci_val[3],
#     Sensitivity = cm$byClass["Sensitivity"],
#     Specificity = cm$byClass["Specificity"]
#   )
# }


# pROC::ci.auc() uses the single ROC curve you give it (from your test set predictions) and estimates the confidence interval for AUC statistically using resampling or an asymptotic formula.
# By default, ci.auc() uses DeLong’s method (DeLong, DeLong & Clarke-Pearson, 1988), which:
#   Treats your test predictions and true labels as the full dataset.
# Estimates the variability of the AUC based on the distribution of ranks of predicted scores between positive and negative cases.
# Produces a 95% CI (by default) without retraining the model.

get_metrics_optimized2 <- function(model=gam_mod_s, test_data=test_sel, 
                                   model_name, quad_vars = NULL, scale_params = NULL,
                                   threshold_type = c("youden", "10pct_omission")) {
  cls <- tolower(model_name)
  threshold_type <- match.arg(threshold_type)  # ensure valid choice
  
  if (cls == "glm") {
    # Scale if scale_params provided
    if (!is.null(scale_params)) {
      
      # Extract mean and sd values from your data frame
      mean_vals <- scale_params[1, grep("_mean$", names(scale_params))]
      sd_vals   <- scale_params[1, grep("_sd$", names(scale_params))]
      
      # Remove the _mean/_sd suffix to match variable names in test_data
      names(mean_vals) <- sub("_mean$", "", names(mean_vals))
      names(sd_vals)   <- sub("_sd$", "", names(sd_vals))
      
      # Scale numeric variables in test_data that have a corresponding mean/sd
      test_data_scaled <- as.data.frame(lapply(names(test_data), function(nm) {
        if (nm %in% names(mean_vals)) {
          (test_data[[nm]] - mean_vals[[nm]]) / sd_vals[[nm]]
        } else {
          test_data[[nm]]  # leave untouched if no scaling info
        }
      }))
      }else{test_data_scaled = test_data}
      
      # Restore names
      names(test_data_scaled) <- names(test_data)
      
      # Make sure response is factor
      test_data_scaled$kelp <- as.factor(test_data$kelp)
    
    # Add quadratic terms only for GLM
    if (!is.null(quad_vars)) {
      for (v in quad_vars) {
        if (v %in% names(test_data_scaled)) {
          test_data_scaled[[paste0("I(", v, "^2)")]] <- test_data_scaled[[v]]^2
        }
      }
    }
    
    pred <- predict(model, newdata = test_data_scaled, type = "response")
    
  } else if (cls == "gam") {
    # GAM: just scale if needed, no quadratic terms added
    if (!is.null(scale_params)) {
      
      # Extract mean and sd values from your data frame
      mean_vals <- scale_params[1, grep("_mean$", names(scale_params))]
      sd_vals   <- scale_params[1, grep("_sd$", names(scale_params))]
      
      # Remove the _mean/_sd suffix to match variable names in test_data
      names(mean_vals) <- sub("_mean$", "", names(mean_vals))
      names(sd_vals)   <- sub("_sd$", "", names(sd_vals))
      
      # Scale numeric variables in test_data that have a corresponding mean/sd
      test_data_scaled <- as.data.frame(lapply(names(test_data), function(nm) {
        if (nm %in% names(mean_vals)) {
          (test_data[[nm]] - mean_vals[[nm]]) / sd_vals[[nm]]
        } else {
          test_data[[nm]]  # leave untouched if no scaling info
        }
      }))
    }else{test_data_scaled = test_data}
    
     # Restore names
      names(test_data_scaled) <- names(test_data)
      
      # Make sure response is factor
      test_data_scaled$kelp <- as.factor(test_data_scaled$kelp)
      
    pred <- predict(model, newdata = test_data_scaled, type = "response")
    
  } else if (cls == "rf") {
    pred <- predict(model, newdata = test_data, type = "prob")[, 2]
  } else if (cls == "brt") {
    pred <- predict(model, newdata = test_data, n.trees = model$n.trees, type = "response")
  } else {
    stop("Unknown model type")
  }
  
  # Keep complete pairs only
  truth <- test_data$kelp
  keep <- is.finite(pred) & !is.na(truth)
  pred  <- pred[keep]
  truth <- truth[keep]
  
  if (length(pred) != length(truth)) stop("Length mismatch: pred=", length(pred), " truth=", length(truth))
  if (length(pred) == 0) stop("No valid predictions after NA filtering.")
  
  # ROC
  roc_obj <- pROC::roc(truth, pred, quiet = TRUE)
  
  # Sensitivities, specificities, thresholds
  sens_vec <- roc_obj$sensitivities
  spec_vec <- roc_obj$specificities
  thr_vec  <- roc_obj$thresholds
  
  # Determine threshold based on chosen type
  if (threshold_type == "youden") {
    J <- sens_vec + spec_vec - 1
    i_best <- which.max(J)
    opt_thresh <- thr_vec[i_best]
  } else if (threshold_type == "10pct_omission") {
    # Find threshold where sensitivity >= 0.9 (10% omission)
    finite_idx <- which(is.finite(thr_vec))
    sens_vec_f <- sens_vec[finite_idx]
    thr_vec_f  <- thr_vec[finite_idx]
    
    # Threshold for 10% omission (sensitivity >= 0.9)
    i_best <- which(sens_vec_f >= 0.9)[1]  # first threshold meeting condition
    if (is.na(i_best)) {
      warning("No threshold achieves 90% sensitivity; using maximum sensitivity threshold instead.")
      i_best <- which.max(sens_vec_f)
    }
    opt_thresh <- thr_vec_f[i_best]
  }
  
  pred_class <- as.integer(pred >= opt_thresh)
  tp <- sum(pred_class == 1 & truth == 1)
  tn <- sum(pred_class == 0 & truth == 0)
  fp <- sum(pred_class == 1 & truth == 0)
  fn <- sum(pred_class == 0 & truth == 1)
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  auc_val     <- as.numeric(pROC::auc(roc_obj))
  tss         <- as.numeric(sensitivity + specificity - 1)
  
  tibble::tibble(
    Model = model_name,
    ThresholdType = threshold_type, 
    Threshold = opt_thresh,
    AUC = auc_val,
    Sensitivity = sensitivity,
    Specificity = specificity,
    TSS = tss
  )
}



get_metrics_optimized <- function(model, test_data, model_name ) {
  # Predict probabilities
  if (model_name %in% c("glm", "gam")) {
    # if (!is.null(quad_vars)) {
    #   for (v in quad_vars) {
    #     if (v %in% names(test_data)) {
    #       test_data[[paste0("I(", v, "^2)")]] <- test_data[[v]]^2
    #     }
    #   }
    # }
    pred <- predict(model, newdata = test_data, type = "response")
  } else if (model_name == "rf") {
    pred <- predict(model, newdata = test_data, type = "prob")[, 2]
  } else if (model_name == "brt") {
    pred <- predict(model, newdata = test_data, n.trees = model$n.trees, type = "response")
  } else {
    stop("Unknown model type")
  }
  
  # Keep complete pairs only
  # True labels from *this* dataset
  truth <- test_data$kelp
  keep <- is.finite(pred) & !is.na(truth)
  pred  <- pred[keep]
  truth <- truth[keep]
  
  if (length(pred) != length(truth)) {
    stop("Length mismatch: pred=", length(pred), " truth=", length(truth))
  }
  if (length(pred) == 0) stop("No valid predictions after NA filtering.")
  
  # ROC
  roc_obj <- pROC::roc(truth, pred, quiet = TRUE)
  
  # Youden's J across all thresholds in the roc object
  sens_vec <- roc_obj$sensitivities
  spec_vec <- roc_obj$specificities
  thr_vec  <- roc_obj$thresholds
  J <- sens_vec + spec_vec - 1
  i_best <- which.max(J)
  opt_thresh <- thr_vec[i_best]
  
  
  # Binary predictions at optimal threshold
  pred_class <- as.integer(pred >= opt_thresh)
  
  # Manual metrics (avoid caret factor issues)
  tp <- sum(pred_class == 1 & truth == 1)
  tn <- sum(pred_class == 0 & truth == 0)
  fp <- sum(pred_class == 1 & truth == 0)
  fn <- sum(pred_class == 0 & truth == 1)
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  auc_val     <- as.numeric(pROC::auc(roc_obj))
  tss  <- as.numeric(sensitivity + specificity - 1)
  
  tibble::tibble(
    Model = model_name,
    Threshold = opt_thresh,
    AUC = auc_val,
    Sensitivity = sensitivity,
    Specificity = specificity,
    TSS = tss
  )
}




#### Function to run all models and compute 10-fold CV AUC  ====================
## Run models
## Select variables based on threshold of weighted score
## Calculate performance with cross validation AUC 
## Store response curves of each model (n=10 CV)

run_models_with_threshold <- function(thresh, train_data, rank_fun, scaling_params = NULL, quad_vars = NULL) {

  # 1. Rank variables using your existing function
  res_rank <- rank_fun(models, types, avg_threshold = thresh)
  vars_selected <- res_rank$summary_table$Variable

  # 2. Prepare training datasets
  # GLM
  quad_vars_sel <- quad_vars[quad_vars %in% vars_selected]
  linear_terms_sel <- setdiff(vars_selected, quad_vars_sel)

  terms_quad_sel <- paste0(quad_vars_sel, " + I(", quad_vars_sel, "^2)")
  glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel, linear_terms_sel), collapse = " + ")))

  train_sel_scaled <- scale_with_params(train_data, scaling_params)

  # GAM
  gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))

  # 3. Fit models
  glm_mod <- glm(glm_formula, data = train_sel_scaled, family = binomial)
  gam_mod <- gam(gam_formula, data = train_sel_scaled, family = binomial)

  rf_mod <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                         data = train_data, ntree = 500, importance = TRUE)

  train_brt <- train_data
  train_brt$kelp <- as.numeric(as.character(train_brt$kelp))
  brt_mod <- dismo::gbm.step(data = train_brt,
                             gbm.x = which(names(train_brt) != "kelp"),
                             gbm.y = which(names(train_brt) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5)

  # 4. Cross-validation (10-fold)
  folds <- caret::createFolds(train_data$kelp, k = 10)

  auc_cv <- function(mod, mod_type, folds, data, formula = NULL) {
    aucs <- numeric(length(folds))
    for(i in seq_along(folds)) {
      train_idx <- folds[[i]]
      test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
      train_fold <- data[train_idx, ]
      test_fold  <- data[test_idx, ]

      # Scale
      train_fold <- scale_with_params(train_fold, scaling_params)
      test_fold  <- scale_with_params(test_fold, scaling_params)

      if(mod_type == "GLM") {
        mod_fold <- glm(formula, data = train_fold, family = binomial)
        pred <- predict(mod_fold, newdata = test_fold, type = "response")
      } else if(mod_type == "GAM") {
        mod_fold <- gam(formula, data = train_fold, family = binomial)
        pred <- predict(mod_fold, newdata = test_fold, type = "response")
      } else if(mod_type == "RF") {
        mod_fold <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                                 data = train_fold, ntree = 500)
        pred <- predict(mod_fold, newdata = test_fold, type = "prob")[,2]
      } else if(mod_type == "BRT") {
        train_brt_fold <- train_fold
        train_brt_fold$kelp <- as.numeric(as.character(train_brt_fold$kelp))
        mod_fold <- dismo::gbm.step(data = train_brt_fold,
                                    gbm.x = which(names(train_brt_fold) != "kelp"),
                                    gbm.y = which(names(train_brt_fold) == "kelp"),
                                    family = "bernoulli",
                                    tree.complexity = 3,
                                    learning.rate = 0.01,
                                    bag.fraction = 0.5,
                                    verbose = FALSE)
        pred <- predict(mod_fold, newdata = test_fold, type = "response", n.trees = mod_fold$n.trees)
      }

      aucs[i] <- as.numeric(pROC::roc(test_fold$kelp, pred)$auc)
    }
    return(list(mean_auc = mean(aucs), se_auc = sd(aucs)/sqrt(length(aucs)), aucs = aucs))
  }

  aucs <- list(
    GLM = auc_cv(glm_mod, "GLM", folds, train_sel_scaled, formula = glm_formula),
    GAM = auc_cv(gam_mod, "GAM", folds, train_sel_scaled, formula = gam_formula),
    RF  = auc_cv(rf_mod, "RF", folds, train_data),
    BRT = auc_cv(brt_mod, "BRT", folds, train_data)
  )

  return(list(threshold = thresh, selected_vars = vars_selected, aucs = aucs))
}



# This function stores the response curves from each model run for the cross validation (fold=10)
run_models_with_threshold_curves <-function(thresh, train_data, rank_fun, scaling_params = NULL, quad_vars = NULL) {
  
  # ---
  res_rank <- rank_fun(models, types, avg_threshold = thresh)
  vars_selected <- res_rank$summary_table$Variable
  
  quad_vars_sel <- quad_vars[quad_vars %in% vars_selected]
  linear_terms_sel <- setdiff(vars_selected, quad_vars_sel)
  
  terms_quad_sel <- paste0(quad_vars_sel, " + I(", quad_vars_sel, "^2)")
  glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel, linear_terms_sel), collapse = " + ")))
  
  train_data$kelp<- as.factor(train_data$kelp)
  train_sel_scaled <- scale_with_params(train_data, scaling_params)
  
  gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))
  
  glm_mod <- glm(glm_formula, data = train_sel_scaled, family = binomial)
  gam_mod <- gam(gam_formula, data = train_sel_scaled, family = binomial)
  
  rf_mod <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                         data = train_data, ntree = 500, importance = TRUE)
  
  train_brt <- train_data
  train_brt$kelp <- as.numeric(as.character(train_brt$kelp))
  brt_mod <- dismo::gbm.step(data = train_brt, 
                             gbm.x = which(names(train_brt) != "kelp"),
                             gbm.y = which(names(train_brt) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5)
  
  folds <- caret::createFolds(train_data$kelp, k = 10)
  
  auc_cv <- function(mod, mod_type, folds, data, formula = NULL) {
    aucs <- numeric(length(folds))
    response_curves <- list()   # <--- add storage for curves
    for(i in seq_along(folds)) {
      train_idx <- folds[[i]]
      test_idx  <- setdiff(seq_len(nrow(data)), train_idx)
      train_fold <- data[train_idx, ]
      test_fold  <- data[test_idx, ]
      
      # Scale
      train_fold <- scale_with_params(train_fold, scaling_params)
      test_fold  <- scale_with_params(test_fold, scaling_params)
      
      if(mod_type == "GLM") {
        mod_fold <- glm(formula, data = train_fold, family = binomial)
        pred <- predict(mod_fold, newdata = test_fold, type = "response")
      } else if(mod_type == "GAM") {
        mod_fold <- gam(formula, data = train_fold, family = binomial)
        pred <- predict(mod_fold, newdata = test_fold, type = "response")
      } else if(mod_type == "RF") {
        mod_fold <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                                 data = train_fold, ntree = 500)
        pred <- predict(mod_fold, newdata = test_fold, type = "prob")[,2]
      } else if(mod_type == "BRT") {
        train_brt_fold <- train_fold
        train_brt_fold$kelp <- as.numeric(as.character(train_brt_fold$kelp))
        mod_fold <- dismo::gbm.step(data = train_brt_fold, 
                                    gbm.x = which(names(train_brt_fold) != "kelp"),
                                    gbm.y = which(names(train_brt_fold) == "kelp"),
                                    family = "bernoulli",
                                    tree.complexity = 3,
                                    learning.rate = 0.005,
                                    bag.fraction = 0.5,
                                    verbose = FALSE)
        pred <- predict(mod_fold, newdata = test_fold, type = "response", n.trees = mod_fold$n.trees)
      }
      
      aucs[i] <- as.numeric(pROC::roc(test_fold$kelp, pred)$auc)
      
      # --- compute response curves for this fold ---
    #   for(var in vars_selected) {
    #     curve_df <- get_curve(mod_fold, var, create_var_grid(var, train_fold), var,
    #                           scale_params = if(mod_type %in% c("GLM","GAM")) scaling_params else NULL)
    #     curve_df$Fold <- i
    #     curve_df$Model <- mod_type
    #     response_curves[[paste(mod_type, var, i, sep = "_")]] <- curve_df
    #   }
    }
    return(list(mean_auc = mean(aucs), se_auc = sd(aucs)/sqrt(length(aucs)), aucs = aucs
                # response_curves = response_curves
                ))
  }
  
  aucs <- list(
    GLM = auc_cv(glm_mod, "GLM", folds, train_sel_scaled, formula = glm_formula),
    GAM = auc_cv(gam_mod, "GAM", folds, train_sel_scaled, formula = gam_formula),
    RF  = auc_cv(rf_mod, "RF", folds, train_data),
    BRT = auc_cv(brt_mod, "BRT", folds, train_data)
  )
  
  # Combine all response curves
  # response_curves_all <- do.call(rbind, unlist(lapply(aucs, `[[`, "response_curves"), recursive = FALSE))
  
  return(list(threshold = thresh, selected_vars = vars_selected, aucs = aucs))
              # response_curves = response_curves_all))
}

  
  

 
# analyze_threshold_performance <- function(results) {
#   
#   # 1. Extract per-fold AUCs into long format
#   all_aucs <- list()
#   
#   for(i in seq_along(results)) {
#     thr <- results[[i]]$threshold
#     for(mod in c("GLM","GAM","RF","BRT")) {
#       df <- data.frame(
#         Threshold = thr,
#         Model = mod,
#         Fold = seq_along(results[[i]]$aucs[[mod]]$aucs),
#         AUC = results[[i]]$aucs[[mod]]$aucs
#       )
#       all_aucs[[length(all_aucs)+1]] <- df
#     }
#   }
#   
#   cv_df <- bind_rows(all_aucs)
#   
#   # 2. Perform Friedman test for each model
#   friedman_results <- cv_df %>%
#     group_by(Model) %>%
#     group_modify(~ {
#       ft <- friedman.test(AUC ~ Threshold | Fold, data = .x)
#       data.frame(Friedman_Chi2 = ft$statistic,
#                  df = ft$parameter,
#                  p.value = ft$p.value)
#     }) %>%
#     ungroup()
#   
#   # 3. Boxplot of AUC by threshold for each model
#   p <- ggplot(cv_df, aes(x = factor(Threshold), y = AUC, fill = factor(Threshold))) +
#     geom_boxplot() +
#     facet_wrap(~Model) +
#     theme_bw(base_size = 14) +
#     labs(x = "Variable Selection Threshold", y = "Cross-validated AUC",
#          fill = "Threshold")
#   
#   return(list(
#     cv_df = cv_df,
#     friedman_results = friedman_results,
#     plot = p
#   ))
# }
# 
#




### Predict the models in the entire area ======================================

# Function to align raster stack with model predictors
align_stack_to_model <- function(rstack, model) {
  # Extract predictor names from the model
  # Works for glm/gam; drops the response
  predictors <- attr(terms(model), "term.labels")
  
  # Get stack names
  rnames <- names(rstack)
  
  # Check for missing variables
  missing_vars <- setdiff(predictors, rnames)
  if (length(missing_vars) > 0) {
    stop("These predictors are missing from the raster stack: ", 
         paste(missing_vars, collapse = ", "))
  }
  
  # Reorder the raster stack to match model predictors
  rstack <- rstack[[predictors]]
  
  return(rstack)
}

# Example usage:
# raster_stack_aligned <- align_stack_to_model(raster_stack, model)


predict_raster_prob <- function(model=gam_mod_s, raster_stack= raster_stack_predict_scaled,  model_name= "gam") {
  model_name <- tolower(model_name)
 
  
  if (model_name %in% c("glm", "gam")) {
    ### Conver raster stack into dataframe to predict the habitat suitability 
    # raster_vals <- as.data.frame(raster_stack, xy = F, na.rm = FALSE)
    # pred_df[] <- lapply(pred_df, function(x) as.numeric(x)) # Ensure all columns are numeric (GLM does not accept matrix/array)
    
    prob <- terra::predict(newdata = raster_stack, model,  type = "response", na.rm = TRUE)
    
    return(prob)
    
  } else if (model_name == "rf") {
    # RF can return a matrix of probabilities; take the 2nd column ("1")
    prob <- raster::predict(raster_stack, model, type = "prob", index = 2)
    return(prob)
    
  } else if (model_name == "brt") {
    prob <- terra::predict(
      raster_stack,
      model,
      fun = function(model, data) {
        # Identify rows with any NA values
        na_rows <- apply(data, 1, function(x) any(is.na(x)))
        
        # Prepare output vector
        pred <- rep(NA, nrow(data))
        
        # Predict only for rows without NA
        pred[!na_rows] <- predict(model,
                                  data[!na_rows, , drop = FALSE],
                                  n.trees = model$n.trees,
                                  type = "response")
        return(pred)
      }
    )
  } else {
    stop("Unknown model type: ", model_name)
  }
}


### Function to binarize and plot a prediction raster using threshold ===================
binarize_raster <- function(pred_raster= gam_pred_raster, model_name="gam", thresholds_table) {
  
  pred_raster<- raster(pred_raster)
  # Get model-specific threshold
  thr <- thresholds_table$Threshold[thresholds_table$Model == model_name]
  if (length(thr) == 0) {
    stop(paste("No threshold found for model", model_name))
  }
  
  # Apply threshold -> TRUE/FALSE
  bin_raster <- pred_raster >= thr
  
  # Convert to integer 0/1
  bin_raster <- as.integer(bin_raster)
  
  # Keep names clear
  names(bin_raster) <- paste0(model_name, "_binary")
  
  return(bin_raster)
}




plot_raster_gg <- function(rast= glm_bin, title_plot="GLM Binary Prediction (Study Area)") {
  
  # rast<- raster(rast)
  # convert raster to dataframe
  rast_df <- as.data.frame(rast, xy = TRUE, na.rm = TRUE)
  colnames(rast_df) <- c( "value", "x", "y")
  rast_df$value<- as.factor(rast_df$value)
  # summary(rast_df)
  
  p<- ggplot() +
    # geom_sf(data = world_proj, fill = "grey90", color = "grey50", size = 0.3) +
    # geom_tile(data = land_highres, aes(x = x, y = y, fill = value)) +
    # geom_sf(data = coastline_proj, color = "black", size = 0.4) +
    geom_sf(data = land_highres, fill = "grey85", color = NA) +
    geom_tile(data = rast_df, aes(x = x, y = y, fill = value)) +  # BINARY RASTER
    geom_sf(data = coastline, color = "grey40", size = 0.3) +
    coord_sf(xlim = c(study_extent["xmin"], study_extent["xmax"]),
             ylim = c(study_extent["ymin"], study_extent["ymax"]),
             expand = FALSE) +
    scale_x_continuous(breaks = seq(-126, -122, by = 1)) +
    scale_y_continuous(breaks = seq(47, 51, by = 1))+
    scale_fill_manual(values = c("0" = "blue", "1" = "red"),
                      labels = c("0" = "Unsuitable", "1" = "Suitable"),
                      na.value = "transparent") +
    # scale_color_viridis_c(option = "plasma", name = "Standardized\nValue") +
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(size = 10, face = "bold"),
          # plot.text = element_text(size = 9, face = "bold"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          panel.grid = element_blank(),
          legend.position = c(0.86, 0.90))+
    labs(title = title_plot,
         fill = "Bull Kelp Habitat", #Suitable/Unsuitable\n
         x="",
         y="")
  
  return(p)
}






ensemble_raster <- function(pred_rasters_list, weights_df) {
  
  common_models <- intersect(names(pred_rasters_list), weights_df$Model)
  
  rasters <- pred_rasters_list[common_models]
  weights <- weights_df %>% filter(Model %in% names(rasters))
  
  # 2. Normalize weights so they sum to 1
  weights$Weight <- weights$Weight / sum(weights$Weight)
  
  # 3. Start with an empty raster (all 0s, same extent/resolution as first raster)
  ens <- rasters[[1]] * 0
  
  # 4. Weighted average = sum(predict_i * weight_i)
  for (i in seq_along(rasters)) {
    ens <- ens + rasters[[i]] * weights$Weight[i]
  }
  
  return(list(raster = ens, weights = weights))
}

# How to use:
# results_selected <- bind_rows(
#   get_metrics_optimized(model= glm_mod_s, test_data=test_sel, model_name = "glm"),
#   get_metrics_optimized(gam_mod_s, test_sel, "gam"),
#   get_metrics_optimized(rf_mod_s, test_sel, "rf"),
#   get_metrics_optimized(brt_mod_s, test_sel, "brt")
# )
# 
# # Pivot the results for plotting
# results_long_selected <- results_selected %>%
#   pivot_longer(cols = c("AUC", "Sensitivity", "Specificity"),
#                names_to = "Metric", values_to = "Value")
# head(results_long_selected)
# 
# 
# pred_rasters_list <- list(
#   glm = glm_pred_raster,
#   gam = gam_pred_raster,
#   rf  = rf_pred_raster,
#   brt = brt_pred_raster
# )
# 
# # Calculate TSS + weights
# tss_weights <- calc_tss_points(
#   models = list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s),
#   types  = c("glm", "gam", "rf", "brt"),
#   test_df = test_sel,
#   thresholds_table = thresholds_table
# )
# 
# print(tss_weights)
# # A tibble with Model, TSS, Weight
# 
# # Build ensemble raster
# ens_out <- ensemble_raster(pred_rasters_list, tss_weights)
# 
# # Save raster
# terra::writeRaster(ens_out$raster, "ensemble_suitability.tif", overwrite=TRUE)
# 
# # Save weights + TSS table
# readr::write_csv(ens_out$weights, "model_tss_weights.csv")



uncertainty_raster <- function(pred_rasters_list, weights_df) {
  
  # Only keep models that are in both lists
  common_models <- intersect(names(pred_rasters_list), weights_df$Model)
  rasters <- pred_rasters_list[common_models]
  weights <- weights_df %>% filter(Model %in% common_models)
  
  # Convert to SpatRaster if not already
  # rasters <- lapply(rasters, rast)
  
  # Stack rasters
  raster_stack <- rast(rasters)
  
  # Weighted mean
  # w_mean <- app(raster_stack, fun = function(x) sum(x * weights$Weight, na.rm = TRUE) / sum(weights$Weight))
  weights$Weight <- weights$Weight / sum(weights$Weight)# Normalize weights so they sum to 1
  w_mean <- rasters[[1]] * 0 # 3. Start with an empty raster (all 0s, same extent/resolution as first raster)
  
  # Weighted average = sum(predict_i * weight_i)
  for (i in seq_along(rasters)) {
    w_mean <- w_mean + rasters[[i]] * weights$Weight[i]
  }
  
  # Weighted standard deviation
  w_sd <- app(raster_stack, fun = function(x) sqrt(sum(weights$Weight * (x - sum(x * weights$Weight / sum(weights$Weight)))^2) / sum(weights$Weight)), 
              cores = 1)
  
  # Min and Max across models
  w_min <- app(raster_stack, min, na.rm = TRUE)
  w_max <- app(raster_stack, max, na.rm = TRUE)
  
  return(list(mean = w_mean, sd = w_sd, min = w_min, max = w_max, weights = weights))
}


### Function to get the different threshold values for a specific model =======
get_thresh_table <- function(roc_obj, pred, truth, model_name = "model") {
  sens_vec <- roc_obj$sensitivities
  spec_vec <- roc_obj$specificities
  thr_vec  <- roc_obj$thresholds
  auc_val  <- as.numeric(pROC::auc(roc_obj))
  
  # --- Helper to compute confusion-matrix-based metrics ---
  compute_metrics <- function(opt_thresh) {
    pred_class <- as.integer(pred >= opt_thresh)
    tp <- sum(pred_class == 1 & truth == 1)
    tn <- sum(pred_class == 0 & truth == 0)
    fp <- sum(pred_class == 1 & truth == 0)
    fn <- sum(pred_class == 0 & truth == 1)
    
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    tss         <- sensitivity + specificity - 1
    c(sensitivity, specificity, tss)
  }
  
  # --- Max TSS (Youden’s J) ---
  J <- sens_vec + spec_vec - 1
  i_best <- which.max(J)
  thr_best <- thr_vec[i_best]
  m_best <- compute_metrics(thr_best)
  
  # --- No omission threshold (sensitivity = 1) ---
  i_no_omiss <- which(sens_vec == 1)[1]
  thr_no_omiss <- thr_vec[i_no_omiss]
  m_no_omiss <- compute_metrics(thr_no_omiss)
  
  # --- 10% omission threshold (sensitivity >= 0.9) ---
  finite_idx <- which(is.finite(thr_vec))
  sens_vec_f <- sens_vec[finite_idx]
  thr_vec_f  <- thr_vec[finite_idx]
  i_10 <- which(sens_vec_f >= 0.9)[1]
  if (is.na(i_10)) {
    warning("No threshold achieves 90% sensitivity; using maximum sensitivity threshold instead.")
    i_10 <- which.max(sens_vec_f)
  }
  thr_10 <- thr_vec_f[i_10]
  m_10 <- compute_metrics(thr_10)
  
  # --- Return tidy table ---
  tibble::tibble(
    Model = model_name,
    Criterion = c("Max TSS", "No omission", "10% omission"),
    Threshold = c(thr_best, thr_no_omiss, thr_10),
    Sensitivity = c(m_best[1], m_no_omiss[1], m_10[1]),
    Specificity = c(m_best[2], m_no_omiss[2], m_10[2]),
    TSS = c(m_best[3], m_no_omiss[3], m_10[3]),
    AUC = auc_val
  )
}





### Function to convert model predictions into unsuitable and re-scale suitable areas =====
threshold_rescale <- function(ens_raster= ensemble_raster_pred, opt_thresh = TSS_threshold_ens) {
  require(terra)
  
  # 1. Binary raster: below threshold = 0, above = 1
  binary_rast <- ifel(ens_raster >= opt_thresh, 1, 0)
  
  # 2. Extract only values above threshold
  vals_above <- values(ens_raster, na.rm = TRUE)
  vals_above <- vals_above[vals_above > opt_thresh]
  
  if (length(vals_above) == 0) {
    warning("No values above threshold — only binary raster returned.")
    return(list(binary = binary_rast, rescaled = binary_rast))
  }
  
  min_above <- min(vals_above, na.rm = TRUE)
  max_above <- max(vals_above, na.rm = TRUE)
  
  # 3. Rescale values above threshold to [0,1]
  rescaled <- (ens_raster - min_above) / (max_above - min_above)
  
  # 4. Mask out below threshold → set to 0
  rescaled <- ifel(ens_raster < opt_thresh, 0, rescaled)
  
  return(list(binary = binary_rast, rescaled = rescaled))
}




