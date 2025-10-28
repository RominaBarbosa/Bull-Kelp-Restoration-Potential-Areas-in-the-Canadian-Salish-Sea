
# ---- Load packages ----
library(blockCV)         
library(caret)           
library(randomForest)    
library(mgcv)            
library(gbm)             
library(pROC)            
library(PresenceAbsence) 
library(dplyr)
library(tidyr)
library(raster)   # or terra if you use SpatRaster
library(gstat)
library(sp)



source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/plot_response_curves.R")
source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")
variables_selection_path<- "/Volumes/Romina_PSF/PSF/SDM/Variables_selection"

model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres"
setwd(model_results_path)
train_sel<- read.csv("train_selected_table_FINALMODELS_M7.csv")
train_sel<- train_sel[,-1]

# train_sel_scaled<- read.csv("train_selected_scaled_table_FINALMODELS.csv")
# train_sel_scaled<- train_sel_scaled[,-1]

test_sel<- read.csv("test_selected_table_FINALMODELS_M7.csv")
test_sel<- test_sel[,-1]

# test_sel_scaled<- read.csv("test_selected_scaled_table_FINALMODELS.csv")
# test_sel_scaled<- test_sel_scaled[,-1]

### Selected predictors
glm_mod_se<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))
vars <- attr(terms(glm_mod_se), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)

# Save scaling parameters for using to scale variables in predictions 
train<- train_sel
test<- test_sel

# Function to apply same scaling
scale_with_params <- function(df, params) {
  df %>%
    mutate(across(where(is.numeric),
                  ~ (.-params[[paste0(cur_column(), "_mean")]]) /
                    params[[paste0(cur_column(), "_sd")]]))
}

scaling_params <- train[, colnames(train)%in%vars_selected] %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd), na.rm = TRUE))


train$kelp<- as.factor(train$kelp)
test$kelp <- as.factor(test$kelp)


# 
explore_k_clusters <- function(data, k_max = 15, nstart = 25, seed = 42) {
  set.seed(seed)
  require(cluster)
  
  # ----- 1. Elbow Method -----
  wss <- numeric(k_max)
  for (k in 1:k_max) {
    km <- kmeans(data, centers = k, nstart = nstart)
    wss[k] <- km$tot.withinss
  }
  
  plot(1:k_max, wss, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters K",
       ylab = "Total within-clusters sum of squares",
       main = "Elbow Method")
  
  # ----- 2. Silhouette Method -----
  sil_width <- numeric(k_max)
  for (k in 2:k_max) {
    km <- kmeans(data, centers = k, nstart = nstart)
    ss <- silhouette(km$cluster, dist(data))
    sil_width[k] <- mean(ss[, 3])
  }
  
  plot(2:k_max, sil_width[2:k_max], type = "b", pch = 19,
       xlab = "Number of clusters K",
       ylab = "Average silhouette width",
       main = "Silhouette Method")
  
  # ----- 3. Gap Statistic -----
  gap_stat <- clusGap(data, FUN = kmeans, nstart = nstart, K.max = k_max, B = 50)
  fviz_gap_stat(gap_stat)
  
  return(list(wss = wss, sil_width = sil_width, gap_stat = gap_stat))
}

summary(train[, colnames(train)%in%vars_selected])
results <- explore_k_clusters(train[, colnames(train)%in%vars_selected], k_max = 7, nstart = 25, seed = 42)


K <-  5  # <-- set number of clusters
env_pres <- train[, colnames(train)%in%vars_selected] 

set.seed(123)
cl <- kmeans(env_pres, centers = K, nstart = 25)

# assign cluster IDs back to full dataset (presences + absences)
unique(train$weight)
train_cluster<- train
train_cluster$cluster <- NA
train_cluster$cluster<- cl$cluster

train_scaled_cluster<- scale_with_params(train_sel[,colnames(train)%in%vars_selected], scaling_params)
train_scaled_cluster$kelp<- train_sel$kelp
train_scaled_cluster$cluster <- NA
train_scaled_cluster$cluster <- cl$cluster
train_scaled_cluster$weight<- train_sel$weight
train_scaled_cluster$Cluster<- train_sel$Cluster

test_sel_scaled<- scale_with_params(test_sel[,colnames(train)%in%vars_selected], scaling_params)
test_sel_scaled$kelp<- test_sel$kelp
# ==============================================================================
#  Models Evaluation
# ==============================================================================
# ---- Calculate metrics ----
#=====================================================================
glm_mod_se<- readRDS("glm_mod_s.rds")
gam_mod_se<- readRDS("gam_mod_s.rds")
rf_mod_se<- readRDS("rf_mod_s.rds")
brt_mod_se<- readRDS("brt_mod_s.rds")

# ---- “How does the model trained globally perform on each area?”----
models <- c("glm","gam","rf","brt")
all_results <- list()
all_results_long<- list()
folds<- as.integer(train_scaled_cluster$Cluster)
models_list<- list(glm_mod_se, gam_mod_se, rf_mod_se, brt_mod_se)


for (m in models){
  metrics_list <- list()  # initialize per model
  
  for (i in 2:5){
    if(m == "glm"){
      train_data <- train_scaled_cluster[folds != i, ]
      test_data  <- train_scaled_cluster[folds == i, ]
      mod <- models_list[[1]]
      # pred <- predict(mod, newdata=test_data, type="response")
      fold_metrics <- get_metrics_optimized2(mod, test_data, scale_params= NULL, m,
                                             threshold_type = "youden")
    }
    if(m == "gam"){
      train_data <- train_scaled_cluster[folds != i, ]
      test_data  <- train_scaled_cluster[folds == i, ]
      mod <- models_list[[2]]
      # pred <- predict(mod, newdata=test_data, type="response")
      fold_metrics <- get_metrics_optimized2(mod, test_data, scale_params= NULL, m,
                                             threshold_type = "youden")
    }
    if(m == "rf"){
      train_data <- train_cluster[folds != i, ]
      test_data  <- train_cluster[folds == i, ]
      mod <- models_list[[3]]
      # pred <- predict(mod, newdata=test_data, type="prob")[,2]
      fold_metrics <- get_metrics_optimized2(mod, test_data, m,threshold_type = "youden")
    }
    if(m == "brt"){
      train_data <- train_cluster[folds != i, ]
      test_data  <- train_cluster[folds == i, ]
      mod <- models_list[[4]]
      # pred <- predict(mod, newdata=test_data, type="response", n.trees=2000)
      fold_metrics <- get_metrics_optimized2(mod, test_data, m,threshold_type = "youden")
    }
    
    metrics_list[[i]] <- fold_metrics
  }
  
  # Combine all folds into one data frame
  metrics_df <- do.call(rbind, metrics_list)
  
  # Compute mean and SD across folds
  summary_df <- data.frame(
    Model = m,
    Metric = c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS"),
    Mean = colMeans(metrics_df[, c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS")], na.rm = TRUE),
    SD   = apply(metrics_df[, c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS")], 2, sd, na.rm = TRUE)
  )
  
  all_results[[m]]<-  summary_df
  all_results_long[[m]]<- metrics_df
}

  # Combine all models
final_long_summary <- do.call(rbind, all_results)
final_long <- do.call(rbind, all_results_long)

ggplot(final_long, aes(x= Model, y= AUC, color= Sensitivity))+
  geom_boxplot()+
  ylim(c(0,1))

# Pivot to wide format
final_wide <- final_long_summary %>%
  pivot_wider(names_from = Metric,
              values_from = c(Mean, SD),
              names_sep = "_")

final_wide
final_wide_formatted <- final_wide %>%
  mutate(
    Threshold = sprintf("%.3f (%.3f)", Mean_Threshold, SD_Threshold),
    AUC       = sprintf("%.3f (%.3f)", Mean_AUC, SD_AUC),
    Sensitivity = sprintf("%.3f (%.3f)", Mean_Sensitivity, SD_Sensitivity),
    Specificity = sprintf("%.3f (%.3f)", Mean_Specificity, SD_Specificity),
    TSS = sprintf("%.3f (%.3f)", Mean_TSS, SD_TSS)
  ) %>%
  dplyr::select(Model, Threshold, AUC, Sensitivity, Specificity, TSS)

final_wide_formatted


# Plot of metrics
plot_model_metrics <- function(data, metric_order = c("Threshold","AUC","TSS","Sensitivity","Specificity"), title = "") {
  
  # Pivot Mean columns
  mean_long <- data %>%
    dplyr::select(Model, starts_with("Mean_")) %>%
    pivot_longer(
      cols = -Model,
      names_to = "Metric",
      names_prefix = "Mean_",
      values_to = "Mean"
    )
  
  # Pivot SD columns
  sd_long <- data %>%
    dplyr::select(Model, starts_with("SD_")) %>%
    pivot_longer(
      cols = -Model,
      names_to = "Metric",
      names_prefix = "SD_",
      values_to = "SD"
    )
  
  # Join Mean and SD
  plot_data <- left_join(mean_long, sd_long, by = c("Model", "Metric"))
  
  # Convert Metric to factor for plotting order
  plot_data <- plot_data %>%
    mutate(Metric = factor(Metric, levels = metric_order)) %>%
    as.data.frame()
  
  # Create the plot
  p <- plot_data %>%
    filter(Metric != "Threshold") %>%
    ggplot(aes(x = Metric, y = Mean, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                  position = position_dodge(width = 0.8), width = 0.2) +
    theme_bw() +
    ylim(c(-0.25, 1.25))+
    labs(y = "Metric Value", x = "Metric", title = title) +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}

plot_model_metrics(final_wide)



#=====================================================================
# ---- “Can the model generalize to new, unseen areas?” ----
# Cross validation by spatial environmental clusters =================
models <- c("glm","gam","rf","brt")
all_results_b <- list()
all_results_long_b<- list()
folds<- as.integer(train_cluster$Cluster)

train_weight<- train_cluster
train_scaled_weight<- train_scaled_cluster

for (m in models){
  metrics_list <- list()  # initialize per model
  
  for (i in 2:5){
    if(m == "glm"){
      train_data <- train_scaled_weight[folds != i, ]
      test_data  <- train_scaled_weight[folds == i, ]
      
      terms_quad_sel <- paste0(vars_selected, " + I(", vars_selected, "^2)")
      glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel), collapse = " + ")))
      mod <- glm_mod_se <- glm(glm_formula, 
                               data = train_data, 
                               family = binomial,
                               weights = weight)
      fold_metrics <- get_metrics_optimized2(mod, test_data, scale_params= NULL, m,
                                             threshold_type = "youden")
    }
    if(m == "gam"){
      train_data <- train_scaled_weight[folds != i, ]
      test_data  <- train_scaled_weight[folds == i, ]
      gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))
      mod <- gam(gam_formula, 
                 data = train_data, 
                 family = binomial,
                 weights = weight)
      
      fold_metrics <- get_metrics_optimized2(mod, test_data, scale_params= NULL, m,
                                             threshold_type = "youden")
    }
    if(m == "rf"){
      train_data <- train_weight[folds != i, ]
      test_data  <- train_weight[folds == i, ]
      train_data$kelp <- as.factor(train_data$kelp)
      mod <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                          data = train_data,
                          ntree = 500,
                          importance = TRUE,
                          sampsize = nrow(train_data),
                          case.weights = train_data$weight)
      
      fold_metrics <- get_metrics_optimized2(mod, test_data, m,threshold_type = "youden")
    }
    if(m == "brt"){
      train_data <- train_weight[folds != i, ]
      test_data  <- train_weight[folds == i, ]
      
      train_data$kelp <- as.numeric(as.character(train_data$kelp))
      
      mod <- dismo::gbm.step(data = train_data, 
                             gbm.x = which(names(train_data) %in% vars_selected),
                             gbm.y = which(names(train_data) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5,
                             site.weights = train_data$weight)
      
      fold_metrics <- get_metrics_optimized2(mod, test_data, m,threshold_type = "youden")
    }
    
    metrics_list[[i]] <- fold_metrics
  }
  
  # Combine all folds into one data frame
  metrics_df <- do.call(rbind, metrics_list)
  
  # Compute mean and SD across folds
  summary_df <- data.frame(
    Model = m,
    Metric = c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS"),
    Mean = colMeans(metrics_df[, c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS")], na.rm = TRUE),
    SD   = apply(metrics_df[, c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS")], 2, sd, na.rm = TRUE)
  )
  
  all_results_b[[m]]<-  summary_df
  all_results_long_b[[m]]<- metrics_df
}

# Combine all models
final_long_summary_b <- do.call(rbind, all_results_b)
final_long_b <- do.call(rbind, all_results_long_b)

ggplot(final_long_b, aes(x= Model, y= AUC, color= Sensitivity))+
  geom_boxplot()+
  ylim(c(0,1))

# Pivot to wide format
final_wide_b <- final_long_summary_b %>%
  pivot_wider(names_from = Metric,
              values_from = c(Mean, SD),
              names_sep = "_")

final_wide_b
final_wide_formatted_b <- final_wide_b %>%
  mutate(
    Threshold = sprintf("%.3f (%.3f)", Mean_Threshold, SD_Threshold),
    AUC       = sprintf("%.3f (%.3f)", Mean_AUC, SD_AUC),
    Sensitivity = sprintf("%.3f (%.3f)", Mean_Sensitivity, SD_Sensitivity),
    Specificity = sprintf("%.3f (%.3f)", Mean_Specificity, SD_Specificity),
    TSS = sprintf("%.3f (%.3f)", Mean_TSS, SD_TSS)
  ) %>%
  dplyr::select(Model, Threshold, AUC, Sensitivity, Specificity, TSS)



library(RColorBrewer)






# # 1. Define spatial blocks (e.g., 50 km side length)
# sb <- spatialBlock(speciesData   = species,
#                    species       = "y",          # response column
#                    rasterLayer   = predictors,   # predictor stack
#                    theRange      = 50000,        # block size (in map units, e.g., meters)
#                    k             = 5,            # number of folds
#                    selection     = "random",     # or "systematic"
#                    iteration     = 100,          # try 100 random partitions
#                    biomod2Format = TRUE)         # output compatible with biomod2 if needed
# 
# # 2. Inspect blocks
# plot(sb$blocks)   # see how blocks look in space
# 
# # 3. Access folds for CV
# folds <- sb$folds



### Spatial random cross validation ============================================
# Create k folds (random)
k <- 6
folds_rdom <- createFolds(train_sel$kelp, k = k, list = TRUE, returnTrain = FALSE)
folds <- integer(length(train_sel$kelp))

# Loop through each fold and assign the fold number
for (i in seq_along(folds_rdom)) {
  folds[folds_rdom[[i]]] <- i
}


models <- c("glm","gam","rf","brt")
all_results_CVrdm <- list()
all_results_long_CVrdm<- list()
train_weight<- train_cluster
train_scaled_weight<- train_scaled_cluster

for (m in models){
  metrics_list <- list()  # initialize per model
  
  for (i in 2:5){
    if(m == "glm"){
      train_data <- train_scaled_weight[folds != i, ]
      test_data  <- train_scaled_weight[folds == i, ]
      
      terms_quad_sel <- paste0(vars_selected, " + I(", vars_selected, "^2)")
      glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel), collapse = " + ")))
      mod <- glm_mod_se <- glm(glm_formula, 
                               data = train_data, 
                               family = binomial,
                               weights = weight)
      fold_metrics <- get_metrics_optimized2(mod, test_data, scale_params= NULL, m,
                                             threshold_type = "youden")
    }
    if(m == "gam"){
      train_data <- train_scaled_weight[folds != i, ]
      test_data  <- train_scaled_weight[folds == i, ]
      gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))
      mod <- gam(gam_formula, 
                 data = train_data, 
                 family = binomial,
                 weights = weight)
      
      fold_metrics <- get_metrics_optimized2(mod, test_data, scale_params= NULL, m,
                                             threshold_type = "youden")
    }
    if(m == "rf"){
      train_data <- train_weight[folds != i, ]
      test_data  <- train_weight[folds == i, ]
      train_data$kelp <- as.factor(train_data$kelp)
      mod <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                          data = train_data,
                          ntree = 500,
                          importance = TRUE,
                          sampsize = nrow(train_data),
                          case.weights = train_data$weight)
      
      fold_metrics <- get_metrics_optimized2(mod, test_data, m,threshold_type = "youden")
    }
    if(m == "brt"){
      train_data <- train_weight[folds != i, ]
      test_data  <- train_weight[folds == i, ]
      
      train_data$kelp <- as.numeric(as.character(train_data$kelp))
      
      mod <- dismo::gbm.step(data = train_data, 
                             gbm.x = which(names(train_data) %in% vars_selected),
                             gbm.y = which(names(train_data) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5,
                             site.weights = train_data$weight)
      
      fold_metrics <- get_metrics_optimized2(mod, test_data, m,threshold_type = "youden")
    }
    
    metrics_list[[i]] <- fold_metrics
  }
  
  # Combine all folds into one data frame
  metrics_df <- do.call(rbind, metrics_list)
  
  # Compute mean and SD across folds
  summary_df <- data.frame(
    Model = m,
    Metric = c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS"),
    Mean = colMeans(metrics_df[, c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS")], na.rm = TRUE),
    SD   = apply(metrics_df[, c("Threshold", "AUC", "Sensitivity", "Specificity", "TSS")], 2, sd, na.rm = TRUE)
  )
  
  all_results_CVrdm[[m]]<-  summary_df
  all_results_long_CVrdm[[m]]<- metrics_df
}

# Combine all models
final_long_summary_CVrdm <- do.call(rbind, all_results_CVrdm)
final_long_CVrdm <- do.call(rbind, all_results_long_CVrdm)

ggplot(final_long_CVrdm, aes(x= Model, y= AUC, color= Sensitivity))+
  geom_boxplot()+
  ylim(c(0,1))

# Pivot to wide format
final_wide_CVrdm <- final_long_summary_CVrdm %>%
  pivot_wider(names_from = Metric,
              values_from = c(Mean, SD),
              names_sep = "_")

final_wide_CVrdm
final_wide_formatted_CVrdm <- final_wide_CVrdm %>%
  mutate(
    Threshold = sprintf("%.3f (%.3f)", Mean_Threshold, SD_Threshold),
    AUC       = sprintf("%.3f (%.3f)", Mean_AUC, SD_AUC),
    Sensitivity = sprintf("%.3f (%.3f)", Mean_Sensitivity, SD_Sensitivity),
    Specificity = sprintf("%.3f (%.3f)", Mean_Specificity, SD_Specificity),
    TSS = sprintf("%.3f (%.3f)", Mean_TSS, SD_TSS)
  ) %>%
  dplyr::select(Model, Threshold, AUC, Sensitivity, Specificity, TSS)




# Save results
plot_CVrdm<- plot_model_metrics(final_wide_CVrdm)
plot_CVclusters<- plot_model_metrics(final_wide_b)
cowplot::plot_grid(plot_CVrdm, plot_CVclusters, ncol=2, valign=T, labels = c("A)", "B)"))

ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/CrossValidation_performance_metrics_M7.pdf", width = 20, height = 15, dpi= 300, units="cm")
ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/CrossValidation_performance_metrics_M7.png", width = 20, height = 12, dpi= 300, units="cm")

cv_results<- rbind(final_wide_formatted, final_wide_formatted_b, final_wide_formatted_CVrdm)
cv_results$approach<- rep(c("cluster_evaluation", "model_evaluation_cluster", "model_evaluation_random"),1, each=4)


# write.csv(cv_results, "CValidation_results_M7.csv")




# Get predictions
pred_df_sel <- bind_rows(
  data.frame(Probability = predict(glm_mod_se, newdata = test_sel_scaled, type = "response"), Model = "GLM", Truth = test_sel_scaled$kelp),
  data.frame(Probability = predict(gam_mod_se, newdata = test_sel_scaled, type = "response"), Model = "GAM", Truth = test_sel_scaled$kelp),
  data.frame(Probability = predict(rf_mod_se, newdata = test_sel, type = "prob")[, 2], Model = "RF", Truth = test_sel$kelp),
  data.frame(Probability = predict(brt_mod_se, newdata = test_sel, n.trees = brt_mod_se$n.trees, type = "response"), Model = "BRT", Truth = test_sel$kelp)
)

# Plot distributions
ggplot(pred_df_sel, aes(x = Probability, fill = as.factor(Truth))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = c("0" = "red", "1" = "blue"), name = "Observed") +
  labs(title = "Predicted Probability Distributions by Model",
       x = "Predicted Probability", y = "Density") +
  theme_minimal()


models_selected <- list(glm_mod_se, gam_mod_se, rf_mod_se, brt_mod_se)
types <- c("glm", "gam", "rf", "brt")

plot_var_importance_v2(models_selected, types,  top_n = 10, stacked = F)

variablesImportance_models_selected <- rank_variables(models_selected, types, avg_threshold = NULL)
plot<- plot_variable_ranking(variablesImportance_models_selected$summary_table, threshold = NULL)
plot[[2]]


summary_table<- plot[[1]]
summary_table <- summary_table %>%
  arrange(desc(weighted_score)) %>%
  mutate(Variable = factor(Variable, levels = Variable))  # preserve this order

# 2️⃣ Define expression labels
# var_labels_expr <- c(
#   "temperature_summer_mean"   = expression("Mean Summer Temperature ("*degree*C*")"),
#   "slope_5x5"                 = expression("Slope ("*degree*"; 100 x 100 m)"),
#   "turbidity_summer_mean"     = expression("Mean Summer Turbidity"),
#   "ammonium_spring_SD"        = expression("SD of Spring Ammonium"),
#   "PAR_summer_mean"           = expression("Mean Summer PAR"),
#   "currentSpeed_summer_mean"  = expression("Mean Summer Current Speed (m s^-1)"),
#   "salinity_summer_SD"        = expression("SD of Summer Salinity")
# )

var_labels_expr <- c(
  # "nitrate_summer_minimum"    = expression("Min of Summer Nitrate (uM.L^-1)"),
  "temperature_summer_mean"   = expression("Mean Summer SST ("*degree*C*")"),
  "slope_5x5"                 = expression("Slope ("*degree*"; 100 x 100 m)"),
  "turbidity_summer_mean"     = expression("Mean Summer Turbidity"),
  "ammonium_spring_mean"      = expression("Mean Spring Ammonium"),
  "PAR_summer_mean"           = expression("Mean Summer PAR"),
  "currentSpeed_summer_mean"  = expression("Mean Summer Current Speed (m s^-1)"),
  "nitrate_winter_mean"       = expression("Mean Winter Nitrate"),
  "salinity_summer_mean"      = expression("Mean Winter SSS")
)

# Extract the labels in the same order as the factor levels
labels_to_use <- var_labels_expr[levels(summary_table$Variable)]

summary_table$Variable <- factor(summary_table$Variable, 
                                 levels = summary_table$Variable[order(summary_table$avg_importance)])

# Plot
ggplot(summary_table, aes(x = avg_importance, y = Variable, fill = ColorFlag)) +
  geom_col() +
  geom_errorbarh(aes(xmin = avg_importance - sd_importance,
                     xmax = avg_importance + sd_importance),
                 height = 0.2, color = "black") +
  scale_fill_manual(values = c("Above" = "grey70", "Below" = "grey80"), guide = "none") +
  scale_y_discrete(labels = labels_to_use) +
  labs(x = "Weighted Score", y = NULL, title = "") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 9, face = "bold")
  )

# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/VariablesImportance_SelectedModels_M7.pdf", width = 11, height = 6, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/VariablesImportance_SelectedModels_M7.png", width = 11, height = 6, dpi= 300, units="cm")


plot_importance_profiles(variablesImportance_models_selected$summary_table)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_profile_selectedModels_Sep2025_M7.png", width = 18, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_profile_selectedModels_M7.pdf", width = 18, height = 12, dpi= 300, units="cm")



### PLOT RESPONSE CURVES ======================================================
# Combine response curves across models
create_full_grid <- function(var="PAR_summer_mean", train_sel, model_type = "glm", 
                             quad_vars = NULL, scale_params = NULL, n = 50) {
  
  predictors <- setdiff(names(train_sel), "kelp")
  
  # Sequence for the variable of interest
  x_unscaled <- seq(min(train_sel[[var]], na.rm = TRUE),
                    max(train_sel[[var]], na.rm = TRUE),
                    length.out = n)
  
  # Create a typical row: mean for numeric, mode for factor
  typical_row <- data.frame(matrix(nrow = 1, ncol = 0))
  for (col in predictors) {
    if (col == var) next
    if (is.numeric(train_sel[[col]])) {
      typical_row[[col]] <- mean(train_sel[[col]], na.rm = TRUE)
    } else if (is.factor(train_sel[[col]])) {
      typical_row[[col]] <- names(sort(table(train_sel[[col]]), decreasing = TRUE))[1]
    }
  }
  
  # Replicate the typical row n times and vary the variable of interest
  grid <- typical_row[rep(1, times = n), , drop = FALSE]
  
  grid[[var]] <- x_unscaled
  
  # Handle scaling and quadratic terms
  if (model_type == "glm") {
    if (is.null(scale_params)) stop("GLM requires scale_params.")
    
    scale_grid <- function(grid, params) {
      grid %>%
        mutate(across(where(is.numeric),
                      ~ (.-params[[paste0(cur_column(), "_mean")]]) /
                        params[[paste0(cur_column(), "_sd")]]))
    }
    
    grid_scaled <- scale_grid(grid, scale_params)
    
    # Add quadratic terms for all variables in quad_vars
    if (!is.null(quad_vars)) {
      for (qv in quad_vars) {
        if (qv %in% names(grid_scaled)) {
          if (qv == var) {
            # variable being varied
            grid_scaled[[paste0("I(", qv, "^2)")]] <- grid_scaled[[qv]]^2
          } else {
            # constant variable → square of its mean
            grid_scaled[[paste0("I(", qv, "^2)")]] <- (grid_scaled[[qv]][1])^2
          }
        }
      }
    }
    
    grid_to_predict <- grid_scaled
    
  } else if (model_type == "gam") {
    if (is.null(scale_params)) stop("GAM requires scale_params.")
    grid_to_predict <- as.data.frame(scale_grid(grid, scale_params))
    
  } else if (model_type %in% c("randomForest", "gbm")) {
    # no scaling, no quadratic terms
    grid_to_predict <- grid
    
  } else {
    stop("Unsupported model type")
  }
  
  return(list(grid = grid_to_predict, x_unscaled = x_unscaled))
}



get_curve2 <- function(model, var, grid_info, model_name) {
  
  x_unscaled <- grid_info$x_unscaled
  grid <- grid_info$grid
  cls <- class(model)[1]
  
  if (cls %in% c("glm", "gam")) {
    pred <- predict(model, newdata = grid, type = "response")
    df <- data.frame(
      x = x_unscaled,
      fit = pred,
      model = model_name,
      var = var
    )
    
  } else if (cls %in% c("randomForest", "randomForest.formula")) {
    pd <- pdp::partial(object = model,
                       pred.var = var,
                       train = grid,
                       prob = TRUE)
    df <- data.frame(
      x = pd[[var]],
      fit = 1 - pd$yhat,  # adjust depending on your positive class
      model = model_name,
      var = var
    )
    
  } else if (cls == "gbm") {
    pd <- pdp::partial(object = model,
                       pred.var = var,
                       train = grid,
                       pred.fun = function(object, newdata) {
                         predict(object, newdata, n.trees = model$n.trees, type = "response")
                       },
                       recursive = FALSE)
    df <- data.frame(
      x = pd[[var]],
      fit = pd$yhat,
      model = model_name,
      var = var
    )
    
  } else {
    stop("Unsupported model class")
  }
  
  return(df)
}



# Variables to plot
vars_to_plot <- names(train)[names(train)%in%vars_selected]  # all variables except response

# Models list
models_list <- list(
  list(model = glm_mod_se, type = "glm", name = "GLM", quad_vars = vars_selected),
  list(model = gam_mod_se, type = "gam", name = "GAM", quad_vars = NULL),
  list(model = rf_mod_se, type = "randomForest", name = "RF", quad_vars = NULL),
  list(model = brt_mod_se, type = "gbm", name = "GBM", quad_vars = NULL)
)

# Loop over variables and models
all_curves <- lapply(vars_to_plot, function(v) {
  do.call(rbind, lapply(models_list, function(m) {
    # Create the grid
    grid_info <- create_full_grid(
      var = v,
      train_sel = train_sel[,colnames(train)%in%vars_selected],
      model_type = m$type,
      quad_vars = m$quad_vars,
      scale_params = scaling_params
    )
    # Predict
    get_curve2(
      model = m$model,
      var = v,
      grid_info = grid_info,
      model_name = m$name
    )
  }))
})

# Combine all variables into a single data frame
all_curves_df <- do.call(rbind, all_curves)

var_labels <- c(
  # "nitrate_summer_minimum"    = "Min of Summer Nitrate (uM.L-1)",
  "temperature_summer_mean"   = "Mean Summer Temperature ('°'C)",
  "slope_5x5"                 = "Slope~('°')~; 100 x 100 m)",
  "turbidity_summer_mean"     = "Mean Summer Turbidity",
  "ammonium_spring_mean"        = "Mean Spring Ammonium",
  "PAR_summer_mean"           = "Mean Summer PAR",
  "currentSpeed_summer_mean"  = "Mean Summer Current Speed (m s-1)",
  "nitrate_winter_mean"        = "Mean Winter Nitrate",
  "salinity_summer_mean"        = "Mean Summer SSS"
)


# Add Raw training data points
df_rug <- train %>%
  dplyr::select(all_of(vars_to_plot), kelp) %>%
  filter(kelp=="1")%>%
  tidyr::pivot_longer(
    cols = all_of(vars_to_plot),
    names_to = "var",
    values_to = "x"
  )

df_rug_abs <- train %>%
  dplyr::select(all_of(vars_to_plot), kelp) %>%
  # filter(kelp=="0")%>%
  tidyr::pivot_longer(
    cols = all_of(vars_to_plot),
    names_to = "var",
    values_to = "x"
  )


# Faceted plot for all variables
curves_plot<- ggplot(all_curves_df, aes(x = x, y = fit, color = model)) +
  geom_line(size = 1) +
  facet_wrap(~var, scales = "free_x", ncol=4,
             labeller = labeller(var = as_labeller(var_labels_expr, label_parsed))) +
  labs(
    x = "Predictor value",
    y = "Predicted response",
    color = "Model",
    title = "Response Curves for All Variables and Models"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 7),   # facet label size
    axis.title = element_text(size = 9),                 # axis title size
    axis.text = element_text(size = 7),                  # axis tick labels
    plot.title = element_text(size = 9, face = "bold")   # plot title
  )


curves_plot
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/ResponseCurves_SelectedVariables_M7.png", width = 19, height = 9, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/ResponseCurves_SelectedVariables_M7.pdf", width = 19, height = 9, dpi= 300, units="cm")


# Plot summary average curves among models: 
all_curves_df$var<- as.factor(all_curves_df$var)
summary_curves <- all_curves_df %>%
  group_by(var, x, model) %>%
  summarise(
    fit = mean(fit))%>%
  group_by(var, x) %>%
  summarise(
    mean_fit = mean(fit),
    sd_fit   = sd(fit),
    ymin = min(fit),   # minimum across models
    ymax = max(fit),   # maximum across models
    .groups = "drop"
  )


# Plots
curves_plot2<- ggplot(summary_curves, aes(x = x, y = mean_fit)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.3) +  # shaded area = full range
  geom_line(color = "black", size = 0.6) +       # average curve
  facet_wrap(~var, scales = "free_x", ncol=4) +
  # labeller = labeller(var = as_labeller(var_labels, label_parsed))) +
  labs(
    x = "",
    y = "Predicted Habitat Suitability"#,
    # title = "Average Response Curve with Full Range Across Models"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 7),   # facet label size
    axis.title = element_text(size = 9),                 # axis title size
    axis.text = element_text(size = 7),                  # axis tick labels
    plot.title = element_text(size = 9, face = "bold")   # plot title
  )

curves_plot2 +
  geom_rug(
    data = df_rug_abs,
    aes(x = x, color = factor(kelp)), 
    inherit.aes = FALSE,
    sides = "b", 
    alpha = 0.3 
    # color = "black"
  ) +
  # Add raw predictor values along x-axis
  geom_rug(
    data = df_rug,
    aes(x = x), 
    inherit.aes = FALSE, 
    sides = "b", 
    alpha = 0.3, 
    color = "green"
  ) +
  scale_color_manual(
    values = c("0" = "blue", "1" = "green"),
    labels = c("0" = "Absence", "1" = "Presence")
  ) 


# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/ResponseCurves_SelectedVariables_average_Sep2025_M7.png", width = 19, height = 10, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/ResponseCurves_SelectedVariables_average_Sep2025_M7.pdf", width = 19, height = 10, dpi= 300, units="cm")



















### Exploring variogram - spatial 
# ---- Example data (replace with your own) ----
# data: dataframe with predictors (X1, X2, ...) and response (binary 0/1 in "y")
# coords: matrix/data.frame of coordinates (longitude, latitude)
# coords <- data[,c("lon","lat")]

set.seed(123)

# ---- Example data ----
# predictors: RasterStack or SpatRaster of environmental variables
# coords: matrix/data.frame of coordinates (longitude, latitude)
# species: dataframe with presence/absence, column "y"

# ---- Option 1: Interactive exploration using blockCV ----
# Opens a Shiny app to explore spatial autocorrelation range
predictors_rast<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_FINAL2.tif")


# ---- Compute variogram for one predictor manually ----
# Example using the first raster layer
var_data <- train.xy

coordinates(var_data) <- ~x+y

# Compute empirical variogram
vgm_model <- variogram(temperature_summer_mean ~ 1, var_data)

# Fit a theoretical model (e.g., spherical)
fit <- fit.variogram(vgm_model, model = vgm(psill=0.5, model="Sph", nugget=1, range=500))

# Plot the variogram
plot(vgm_model, fit)

# ---- Identify block size ----
# The 'range' of the fitted model is the distance at which autocorrelation
# essentially disappears. Use this as 'theRange' in spatialBlock:
block_size <- fit$range[2]   # or read from the plot manually
print(paste("Suggested block size (m):", block_size))


