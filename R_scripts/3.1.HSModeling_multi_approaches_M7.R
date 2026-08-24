###==================================================================
### Species Distribution models        SDMs          ################
###                                                  ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Last version: 07-Oct-2025                        ################
###==================================================================
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf)

source("/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")
variables_selection_path<- "/SDM/Variables_selection/M6"

# Load kelp records ============================================================
# kelp_presabs_df<- read.csv("/SDM/Presence_absences_kelp/presence_2014_to_2019_filtered_plus_AbsencesFINAL.csv")
# kelp_presabs_df<- kelp_presabs_df[,-1]

# kelp_data_with_variables<- read.csv("/Presence_absences_kelp_variablesSelectedFINAL_M7.csv")
# kelp_data_with_variables<- kelp_data_with_variables[,-1]


#### PERFORM MODELS ############################################################
# Load necessary libraries
library(dplyr)
library(caret)
library(pROC)
library(mgcv)      # GAM
library(randomForest) # RF
library(gbm)       # BRT
library(terra)  
library(gstat)
library(sp)
library(raster)   # or terra if you use SpatRaster

# ---- Example data ----
# predictors: RasterStack or SpatRaster of environmental variables
# coords: matrix/data.frame of coordinates (longitude, latitude)
# species: dataframe with presence/absence, column "y"


# ================================
# STEP 1: Balance dataset
# ================================
# Keep all presences
df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M7.csv")
df<- df[,-1]
colnames(df)[18]<- "bathymetry_20"

# Remove absences from period 2020-2022
df[which(df$period =="2020_2022" & df$kelp == 0), ]<- NA

df<- na.exclude(df)
  
# There are much more absences than presence records
df%>%
  group_by(kelp, period)%>%
  summarize(n= length(kelp))
# kelp period        n
# 1     0 2014_2019  2509
# 2     1 2014_2019  1183
# 3     1 2020_2022   854


### Downsample absences to match presence count =====
presences <- df %>% filter(kelp == 1)
absences <- df %>% filter(kelp == 0) %>% sample_n(nrow(presences))

summary(presences$bathymetry_20)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -9.6917 -0.8142  1.3000  2.3752  4.4972 26.4123 

# Combine balanced data
df_balanced <- bind_rows(presences, absences)
df_balanced%>%
  group_by(kelp, period)%>%
  summarize(n= length(kelp))
# M7
# kelp period        n
# 1     0 2014_2019  2037
# 2     1 2014_2019  1183
# 3     1 2020_2022   854

ggplot(df_balanced, aes(x=as.factor(kelp), y=bathymetry_20, fill=as.factor(kelp)))+
  geom_violin()+
  geom_boxplot(width=0.1)+ 
  scale_fill_manual(values=c("blue", "green"))+
  labs(x="", y= "Depth (m)", fill="Kelp Record")+
  theme_bw()



### Explore the max depth inhabited by kelp based on Mora-Soto et al (2024) data ====
quantile(presences$bathymetry_20, 0.9)   # 8.395664
quantile(presences$bathymetry_20, 0.95)  # 11.87806
quantile(presences$bathymetry_20, 0.93)  # 10.28292
quantile(presences$bathymetry_20, 0.97)  # 14.591 m
quantile(presences$bathymetry_20, 0.99)  # 18.35999

df$kelp<- as.factor(df$kelp)
df.p<- df
levels(df.p$kelp)<- c("Absence", "Presence")

ggplot(df.p, aes(x=as.factor(kelp), y=bathymetry_20, fill=as.factor(kelp)))+
  geom_violin()+
  geom_boxplot(width=0.1)+ 
  scale_fill_manual(values=c("blue", "green"))+
  labs(x="", y= "Depth (m)", fill="Kelp Record")+
  theme_bw()

ggplot(df.p, aes(x= kelp, y=slope_5x5, color=as.factor(kelp)))+
  geom_violin()+
  # scale_color_manual(values=c("blue", "green"))+
  labs(x="", y= "Depth (m)", fill="Kelp Record")+
  theme_bw()



# ================================
# STEP 2: Train-test split (70/30)
# ================================
library(caret)
set.seed(1234) 
train_index <- createDataPartition(df_balanced$kelp, p = 0.7, list = FALSE)
train <- df_balanced[train_index, ]
test  <- df_balanced[-train_index, ]

train_df<- train
train_df$set<- "train"
train_df$cell<- train_index
test_df<- test
test_df$set<- "test"
test_df$cell<- NA
  
train_test_dataset<- rbind(train_df, test_df)

train_test_dataset %>%
  group_by(set, period, kelp)%>%
  summarize(n= length(set))
# M7
# set   period        n
# 1 test  2014_2019   942
# 2 test  2020_2022   280
# 3 train 2014_2019  2278
# 4 train 2020_2022   574


# write.csv(train_test_dataset, "/training_testing_datasets_blob_M7.csv")
# train_test_dataset<- read.csv("/training_testing_datasets_blob_M7.csv")
# train_test_dataset<- (train_test_dataset[,-1])

train_test_dataset%>%
  group_by(set, period, kelp)%>%
  summarize(n= length(kelp))
# set   period     kelp     n
# 1 test  2014_2019     0   611
# 2 test  2014_2019     1   331
# 3 test  2020_2022     1   280
# 4 train 2014_2019     0  1426
# 5 train 2014_2019     1   852
# 6 train 2020_2022     1   574


# Exclude columns of variables not used
colnames(train_test_dataset)
predictors <- setdiff(
  names(train_test_dataset),
  c("kelp", "depth","x", "y", "set", "cell", "substrate_aligned", "period", "substrate", "bathymetry_20", "PAR_summer_mean") 
  )

library(vegan)   # or stats::prcomp
train_test_dataset_clean <- train_test_dataset %>%
  filter(if_all(all_of(predictors), ~ is.finite(.))) %>%
  drop_na(all_of(predictors))
summary(train_test_dataset_clean)

# Extract the cleaned environmental matrix
env_clean <- train_test_dataset_clean %>%
  dplyr::select(all_of(predictors))



# Subset the data to only predictors
predictors_data <- train_test_dataset %>% dplyr::select(all_of(predictors))
colnames(predictors_data)

# Calculate correlation matrix (use pairwise complete obs to handle NAs)
cor_mat <- cor(predictors_data, use = "pairwise.complete.obs")

# Plot 
corrplot::corrplot(cor_mat, method = "color", type = "upper", 
                   tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, diag = FALSE)


# Select predictor variables again based on correlations (>0.8 pearson cor)
library(tidyverse)
predictors <- setdiff(
  predictors,
  c("ammonium_winter_mean",
    "ammonium_summer_minimum",
    "nitrate_summer_minimum"))

# Subset the data to only predictors
predictors_data <- train_test_dataset %>% dplyr::select(all_of(predictors))
colnames(predictors_data)
cor_mat <- cor(predictors_data, use = "pairwise.complete.obs")

# Original long names
original_names <- colnames(cor_mat)

# Define shorter names in the same order
short_names <- c(  "Mean spring ammonium", 
                   "Mean summer current speed",
                   # "Mean summer PAR", 
                   "Mean summer SST",    
                   "Mean summer turbidity",
                   "Mean summer SSS",
                    "Slope" ,  
                   "TPI"
)  # adjust as needed

# Assign shorter names to correlation matrix dimnames
colnames(cor_mat) <- short_names
rownames(cor_mat) <- short_names

# pdf("/correlation_plot.pdf", width = 10, height = 8)
# corrplot::corrplot(cor_mat, method = "color", type = "upper", 
#                    tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, diag = FALSE)
# dev.off()



### Scale predictors for GLM and GAM  ==========================================
# Save scaling parameters for using to scale variables in predictions 
train<- train_test_dataset %>% 
  filter(set== "train")

test<- train_test_dataset %>%
  filter(set== "test")


# Function to apply same scaling
# Function to scale data using precomputed params
scale_with_params <- function(df, params) {
  df %>%
    mutate(across(where(is.numeric),
                  ~ (.-params[[paste0(cur_column(), "_mean")]]) /
                    params[[paste0(cur_column(), "_sd")]]))
}

# Compute scaling parameters correctly (no deprecation warning)
scaling_params <- train %>%
  summarise(across(where(is.numeric),
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x, na.rm = TRUE))))


train <- train %>% dplyr::select(kelp, all_of(predictors))
test <- test %>% dplyr::select(kelp, all_of(predictors))
train$kelp<- as.factor(train$kelp)
test$kelp <- as.factor(test$kelp)


train_scaled <- scale_with_params(train, scaling_params)
test_scaled  <- scale_with_params(test,  scaling_params)

train %>%
  pivot_longer(
    cols = -kelp,               # all columns except kelp
    names_to = "variable", 
    values_to = "value"
  ) %>%
  ggplot(aes(x = factor(kelp), y = value)) +
  geom_violin()+
  geom_boxplot(width=0.4) +
  facet_wrap(~ variable, scales = "free_y", ncol=3) +
  labs(x = "Kelp presence (0 = absent, 1 = present)",
       y = "Value") +
  theme_bw()


# ==============================================================================
# STEP 3: Fit models
# ==============================================================================
library(dplyr)
library(cluster)
library(mgcv)
library(randomForest)
library(dismo)
library(ggplot2)
library(cluster)
library(factoextra)

explore_k_clusters <- function(data, k_max = 15, nstart = 25, seed = 42) {
  set.seed(seed)
  
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

results <- explore_k_clusters(train[,2:ncol(train)], k_max = 7, nstart = 25, seed = 42)


K <-  5  # <-- set number of clusters
env_pres <- train_scaled 

set.seed(123)
cl <- kmeans(env_pres, centers = K, nstart = 25)

# assign cluster IDs back to full dataset (presences + absences)
train$cluster <- NA
train$cluster<- cl$cluster


# Visualize in biplot:
# Select only numeric columns for PCA
num_vars <- names(train)[sapply(train, is.numeric)]
pca_data <- scale(train[, num_vars[1:8]])

# Run PCA
pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_res$x)  # PCA coordinates
pca_scores$cluster <- as.factor(train$cluster)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA Biplot of Clusters",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Adding variables as arrows
library(factoextra)
fviz_pca_biplot(pca_res,
                label = "var",       # show variable names
                habillage = train$cluster, # color by cluster
                addEllipses = TRUE,  # optional: confidence ellipse
                ellipse.level = 0.95,
                palette = "Set1")


train.xy<- train_test_dataset %>% 
  filter(set== "train")

train.xy$Cluster<- train$cluster
str(train.xy)
ggplot(train.xy, aes(y= y, x= x, color=as.factor(Cluster)))+
  geom_point()+ theme_bw()

# -----------------------------
# 4. Compute cluster weights  
# -----------------------------
train[which(train$cluster == 4), "cluster"]<- 1 # all cluster 4 is integrated to cluster 1 (this changes everytime you run the pca so make sure to check the map)

freq <- table(train[which(train$kelp==1), "cluster"])
inv_freq <- 1 / as.numeric(freq)
names(inv_freq) <- names(freq)

# Now lookup by cluster number as character
pres_weights <- inv_freq[as.character(train$cluster[train$kelp==1])]
# pres_weights <- inv_freq[train[which(train$kelp==1), "cluster"]]

# normalize to mean = 1
pres_weights <- pres_weights * (length(pres_weights) / sum(pres_weights))

# add weights to dataset
train_weight<- train
train_weight$weight <- 1
train_weight$weight[train_weight$kelp == 1] <- pres_weights

train_scaled_weight<- train_scaled
train_scaled_weight$weight <- 1
train_scaled_weight$weight[train_scaled_weight$kelp == 1] <- pres_weights

# -----------------------------
# 5. Fit models with weights
# -----------------------------
## GLM
# ==============================================================================
library(car)        # vif
library(broom)      # tidy model output
library(visreg)     # partial effect plots


# Specify the predictors for quadratic terms 
quad_vars <- colnames(train_scaled_weight)[2:8]

# Build the formula with quadratic terms using I(x^2)
#    Use I(x^2) rather than poly() so prediction on new data is straightforward.
terms_quad <- paste0(quad_vars, " + I(", quad_vars, "^2)")
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad), collapse = " + ")))

# Fit the GLM (binomial)
glm_mod_s <- glm(glm_formula,
                 data = train_scaled_weight,
                 family = binomial(link = "logit"),
                 weights = weight)
# Basic summary & multicollinearity check
summary(glm_mod_s)
vif_vals <- vif(glm_mod_s)     # check for high VIFs (>5 or >10)
print(vif_vals)

#M7 
#  ammonium_spring_mean     I(ammonium_spring_mean^2)      currentSpeed_summer_mean 
# 1.774754                      1.776926                      3.694301 
# I(currentSpeed_summer_mean^2)               PAR_summer_mean          I(PAR_summer_mean^2) 
# 2.362984                      2.548872                      1.329940 
# temperature_summer_mean  I(temperature_summer_mean^2)         turbidity_summer_mean 
# 10.160168                      2.507361                      7.299808 
# I(turbidity_summer_mean^2)          salinity_summer_mean     I(salinity_summer_mean^2) 
# 4.010847                     19.617388                      7.369586 
# slope_5x5                I(slope_5x5^2)                       TPI_3x3 
# 2.371280                      2.269577                      1.589386 
# I(TPI_3x3^2) 
# 1.754806 


## GAM
# ==============================================================================
gam_formula <- as.formula(paste("kelp ~", paste0("s(", predictors, ")", collapse = " + ")))
gam_mod_s <- gam(gam_formula,
                                  data = train_scaled_weight,
                                  family = binomial,
                                  weights = weight)
                 
                 
visreg(gam_mod_s, "salinity_summer_mean", scale = "response")

## Random Forest
# ==============================================================================
train_weight$kelp <- as.factor(train_weight$kelp)
rf_mod <-randomForest(x = train_weight[, predictors], y = train_weight$kelp, ntree = 500)
rf_mod_s <- randomForest(as.formula(paste("kelp ~", paste(predictors, collapse = " + "))),
                         data = train_weight,
                         ntree = 1500,
                         importance = TRUE,
                         sampsize = floor(0.7 * nrow(train_weight)),  # 70% subsample per tree
                         # sampsize = nrow(train_weight),
                         case.weights = train_weight$weight)

# Call:
#   randomForest(formula = as.formula(paste("kelp ~", paste(predictors,      collapse = " + "))), data = train, ntree = 500, importance = TRUE,      sampsize = nrow(train), case.weights = train$weight) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 13.87%
# Confusion matrix:
#   0   1 class.error
# 0 683 146   0.1761158
# 1  84 745   0.1013269

## BRT
# ==============================================================================
train_brt <- train_weight
colnames(train_brt)
train_brt$kelp <- as.numeric(as.character(train_brt$kelp))
brt_mod_s <- dismo::gbm.step(data = train_brt,
                             gbm.x = which(names(train_brt) %in% predictors),
                             gbm.y = which(names(train_brt) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 1,
                             learning.rate = 0.005,
                             bag.fraction = 0.5,
                             site.weights = train_brt$weight)


# M7 weighted records, including presences from both periods and absences only from period 1.
# fitting final gbm model with a fixed number of 3800 trees for kelp
# mean total deviance = 1.386 
# mean residual deviance = 0.381 
# 
# estimated cv deviance = 0.617 ; se = 0.024 
# 
# training data correlation = 0.882 
# cv correlation =  0.786 ; se = 0.006 
# 
# training data AUC score = 0.983 
# cv AUC score = 0.943 ; se = 0.003 
# 
# elapsed time -  0.97 minutes 


# ================================
# STEP 4: Predictions & Evaluation
# ================================
#-----------------------------
# Collect results
results <- bind_rows(
  get_metrics_optimized(glm_mod_s, train_scaled, "glm"),
  get_metrics_optimized(gam_mod_s, train_scaled, "gam"),
  get_metrics_optimized(rf_mod_s, train, "rf"),
  get_metrics_optimized(brt_mod_s, train, "brt")
)

# M7
# Model Threshold   AUC Sensitivity Specificity   TSS
# 1 glm       0.548 0.907       0.818       0.816 0.634
# 2 gam       0.491 0.928       0.867       0.827 0.694
# 3 rf        0.504 1           1           1     1    
# 4 brt       0.458 0.983       0.943       0.921 0.865


# Pivot the results for plotting
results_long <- results %>%
  pivot_longer(cols = c("AUC", "Sensitivity", "Specificity"),
               names_to = "Metric", values_to = "Value")
head(results_long)


# Plot
ggplot(results_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set2") +
  ylim(c(0,1))+
  labs(title = "Model Performance Comparison (testing data)",
       y = "Score",
       x = "Model") +
  theme_bw()


# Metrics with testing data
results_test <- bind_rows(
  get_metrics_optimized(glm_mod_s, test_scaled, "glm"),
  get_metrics_optimized(gam_mod_s, test_scaled, "gam"),
  get_metrics_optimized(rf_mod_s, test, "rf"),
  get_metrics_optimized(brt_mod_s, test, "brt")
)
# Model Threshold   AUC Sensitivity Specificity   TSS
# 1 glm       0.580 0.909       0.805       0.833 0.638
# 2 gam       0.537 0.917       0.835       0.831 0.666
# 3 rf        0.567 0.961       0.903       0.884 0.787
# 4 brt       0.461 0.949       0.900       0.845 0.745


# =============================================================================
# STEP 5: Variable Selection based on Variables Importance
# =============================================================================
# each model family measures "importance" differently, so we’ll need to standardize how we extract and plot them:
# GLM & GAM → We can use caret::varImp() which ranks predictors based on their absolute t-statistics (or z-values).
# RF → Use randomForest::importance() or ranger::importance().
# BRT → Use summary() from gbm which returns relative influence of predictors.
# Then we can combine them into one table for side-by-side comparison.


# Use functions to Get variable importance and plot the results
models <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")

p_vars<- plot_var_importance_v2(models, types,  top_n = 16, stacked = TRUE)
p_vars2<- plot_var_importance_v2(models, types,  top_n = 16)


# Explanation for the ranking of variables strategy:
# Standard deviation of the standardized importance (Importance_std) across models for each variable.
# Importance_std = Importance / max(Importance) * 100 rescales each model’s variable importance to a 0–100% scale.
# Then you summarize min, mean, max importance across models.
# Importance_SD → variability of importance across models.
# Weighted_Score → high if a variable is both important and consistent.
# Keep → safeguards against variables that are consistently low.
# Sorting by Weighted_Score gives the final ranking.

# Ranking variables and exclude less important ones using a threshold of 10 weighted score
res <- rank_variables(models, types, avg_threshold = NULL)
res$summary_table 
plot_variable_ranking_v2(res$summary_table, threshold = 10)
# ggsave("/VariablesImportance_for_SelecBasedOnImportance.pdf", width = 15, height = 13, dpi= 300, units="cm")


plot_importance_profiles_v2(res$summary_table, threshold=10)
# ggsave("/VariablesImportance_profile_VarsSelection.pdf", width = 18, height = 12, dpi= 300, units="cm")



# ==============================================================================
# STEP 7: SDMs with selected variables 
# ==============================================================================
setwd("/SDM/SDM_results/Sep2025_M7_weightedPres")
thresh= 10

vars_selected <- rank_variables(models, types, avg_threshold = 10)
vars_selected <- vars_selected$summary_table$Variable
# M7
# "turbidity_summer_mean"    "temperature_summer_mean"  "PAR_summer_mean"         
# [4] "salinity_summer_mean"     "slope_5x5"                "currentSpeed_summer_mean"
# [7] "ammonium_spring_mean"

# Refit models using only selected vars
train_sel <- train_weight %>% dplyr::select(all_of(c("kelp", vars_selected)))
test_sel <- test %>% dplyr::select(all_of(c("kelp", vars_selected)))

train_sel$kelp<- as.factor(train_sel$kelp)

scaling_params_2 <- train_sel %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd), na.rm = TRUE))

train_sel_scaled <- scale_with_params(train_sel, scaling_params_2)
test_sel_scaled  <- scale_with_params(test_sel,  scaling_params_2)


K <-  5  # <-- set number of clusters
env_pres <- train_sel_scaled 
set.seed(42)
cl <- kmeans(env_pres[,2:8], centers = K, nstart = 25)

# assign cluster IDs back to full dataset (presences + absences)
train_sel_scaled$cluster <- NA
train_sel_scaled$cluster<- cl$cluster

# assign cluster to testing points
predict_kmeans <- function(kmod, newdata, 
                           train_for_scaling = NULL, 
                           center = NULL, scale_values = NULL, 
                           do_scale = FALSE, 
                           na.action = c("error","omit")) {
  # kmod: object returned by stats::kmeans (must have $centers)
  # newdata: data.frame or matrix with same variables (same order) as used to fit kmod
  # train_for_scaling: optional training data (used to compute center/scale if do_scale=TRUE and center/scale_values not supplied)
  # center/scale_values: numeric vectors of length p (optional alternative to train_for_scaling)
  # do_scale: logical; if TRUE will scale newdata and centers using provided or derived center/scale
  # na.action: "error" (default) or "omit" (will drop rows with NA and return NA for them)
  
  # --- checks & coercions
  if (is.null(kmod$centers)) stop("kmod must be a kmeans object with $centers.")
  centers <- as.matrix(kmod$centers)                # k x p
  if (is.data.frame(newdata)) newmat <- as.matrix(newdata)
  else newmat <- as.matrix(newdata)
  if (!is.numeric(newmat)) stop("newdata must be numeric (no factors/characters).")
  
  p_centers <- ncol(centers)
  p_new <- ncol(newmat)
  if (p_centers != p_new) stop(sprintf("Dim mismatch: centers have %d cols but newdata has %d cols.", p_centers, p_new))
  
  na.action <- match.arg(na.action)
  if (na.action == "omit") {
    nas <- apply(newmat, 1, function(r) any(is.na(r)))
    newmat_na <- newmat
    newmat <- newmat[!nas, , drop = FALSE]
  } else {
    if (any(is.na(newmat))) stop("newdata contains NA values; set na.action='omit' or impute before calling.")
  }
  
  # --- scaling if requested
  if (do_scale) {
    if (is.null(center) || is.null(scale_values)) {
      if (is.null(train_for_scaling)) stop("When do_scale=TRUE provide either center/scale_values or train_for_scaling.")
      trainm <- as.matrix(train_for_scaling)
      if (ncol(trainm) != p_new) stop("train_for_scaling must have same columns/order as newdata/centers.")
      center <- colMeans(trainm, na.rm = TRUE)
      scale_values <- apply(trainm, 2, sd, na.rm = TRUE)
      scale_values[scale_values == 0] <- 1    # avoid division by zero
    }
    # apply scaling to both new data and centers (centers are in the same space as training)
    newmat <- sweep(newmat, 2, center, FUN = "-")
    newmat <- sweep(newmat, 2, scale_values, FUN = "/")
    centers <- sweep(centers, 2, center, FUN = "-")
    centers <- sweep(centers, 2, scale_values, FUN = "/")
  }
  
  # --- efficient squared Euclidean distance calc
  # D_ij = sum_k (x_i_k - c_j_k)^2 = sum_k x_i_k^2 + sum_k c_j_k^2 - 2 * x_i %*% t(c_j)
  X2 <- rowSums(newmat * newmat)
  C2 <- rowSums(centers * centers)
  # tcrossprod gives X %*% t(C)
  # D is n x k
  D <- outer(X2, C2, "+") - 2 * tcrossprod(newmat, centers)
  
  # numerical safety (tiny negatives due to floating point)
  D[D < 0 & D > -1e-12] <- 0
  
  # get index of min distance (max.col on -D is fast)
  clusters_assigned <- max.col(-D, ties.method = "first")
  
  # create output with NA rows if needed
  if (exists("nas") && any(nas)) {
    out <- rep(NA_integer_, length(nas))
    out[!nas] <- clusters_assigned
    clusters_assigned <- out
  }
  
  return(clusters_assigned)
}



test_sel_scaled$cluster<-  predict_kmeans(cl, test_sel_scaled[,2:8], train_for_scaling = train, do_scale = F)
summary(as.factor(test_sel_scaled$cluster))
# 1   2   3   4   5 
# 146 197 327 447 105 


# Visualize in biplot:
# Select only numeric columns for PCA
num_vars <- names(train_sel_scaled)[sapply(train_sel_scaled, is.numeric)]
pca_data <- scale(train_sel_scaled[, num_vars[1:7]])

# Run PCA
pca_res <- prcomp(pca_data, center = TRUE, scale. = F)
pca_scores <- as.data.frame(pca_res$x)  # PCA coordinates
pca_scores$cluster <- as.factor(train_sel_scaled$cluster)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA Biplot of Clusters",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Adding variables as arrows
factoextra::fviz_pca_biplot(pca_res,
                label = "var",       # show variable names
                habillage = train_sel_scaled$cluster, # color by cluster
                addEllipses = F,  # optional: confidence ellipse
                ellipse.level = 0.95,
                palette = "Set1")


# Plot points in map
train.xy<- train_test_dataset %>% 
  filter(set== "train")

train.xy$Cluster<- train_sel_scaled$cluster
str(train.xy)

train_map<- ggplot(train.xy, aes(y= y, x= x, color=as.factor(Cluster)))+
  geom_point()+ theme_bw()+ 
  labs(title = "Training records (70%)", fill= "Env. Cluster")

test.xy<- train_test_dataset %>% 
  filter(set== "test")

test.xy$Cluster<- test_sel_scaled$cluster
testing_map<- ggplot(test.xy, aes(y= y, x= x, color=as.factor(Cluster)))+
  geom_point()+ 
  labs(title = "Testing records (30%)", color= "Env. Cluster")+ 
  theme_bw()

cowplot::plot_grid(train_map, testing_map)

# merge cluster 4 with 5 because 4 has few records and is widely distributed
# Carfull that every run of the kmeans give different order of cluster numbers 
train_sel_scaled[which(train_sel_scaled$cluster==4), "cluster"]<- 1 
test_sel_scaled[which(test_sel_scaled$cluster==4), "cluster"]<- 1 #

# Compute cluster weights  
# raw weights for presences
# Frequency table already has cluster IDs as names
freq <- table(train_sel_scaled$cluster[train_sel_scaled$kelp == 1])

# Inverse frequency, preserving names
inv_freq <- 1 / freq

# Extract weights for presence rows
pres_weights <- inv_freq[as.character(train_sel_scaled$cluster[train_sel_scaled$kelp == 1])]

# normalize to mean = 1
pres_weights <- pres_weights * (length(pres_weights) / sum(pres_weights))

# add weights to dataset
train_sel_scaled$weight <- 1
train_sel_scaled$weight[train_sel_scaled$kelp == 1] <- pres_weights
head(train_sel_scaled)


train_sel$weight <- 1
train_sel$weight[train_sel$kelp == 1] <- pres_weights
head(train_sel)


train.xy$weight <- 1
train.xy$weight[train.xy$kelp == 1] <- pres_weights
head(train.xy)

# ==============================================================================
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres")

# write.csv(train.xy, "train_selected_table_FINALMODELS_M7.csv")
# write.csv(test.xy, "test_selected_table_FINALMODELS_M7.csv")
# write.csv(test_sel_scaled, "test_selected_scaled_table_FINALMODELS.csv")
# write.csv(train_sel_scaled, "train_selected_scaled_table_FINALMODELS.csv")

## GLM
quad_vars_sel <- quad_vars[quad_vars %in% vars_selected]
terms_quad_sel <- paste0(quad_vars_sel, " + I(", quad_vars_sel, "^2)")
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel), collapse = " + ")))
glm_mod_se <- glm(glm_formula, 
                 data = train_sel_scaled, 
                 family = binomial,
                 weights = weight)

# saveRDS(glm_mod_se, "glm_mod_s.rds")

## GAM
gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))
gam_mod_se <- gam(gam_formula, 
                 data = train_sel_scaled, 
                 family = binomial,
                 weights = weight)

# saveRDS(gam_mod_se, "gam_mod_s.rds")

## Random Forest
train_sel$kelp <- as.factor(train_sel$kelp)
rf_mod_se <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                         data = train_sel,
                         ntree = 500,
                         importance = TRUE,
                         sampsize = floor(0.7 * nrow(train_weight)),  # 70% subsample per tree NEW
                         # ntree = 500,       # M7 
                         # importance = TRUE, # M7 
                         # sampsize = nrow(train_sel), # M7
                         case.weights = train_sel$weight)

# saveRDS(rf_mod_se, "rf_mod_s.rds")

library(pdp)
varImpPlot(rf_mod_se) 


## BRT
train_brt_sel <- train_sel# %>% select(-env_cluster)
colnames(train_brt_sel)
train_brt_sel$kelp <- as.numeric(as.character(train_brt_sel$kelp))
brt_mod_se <- dismo::gbm.step(data = train_brt_sel, 
                             gbm.x = which(names(train_brt_sel) %in% vars_selected),
                             gbm.y = which(names(train_brt_sel) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 1,  # M7 = 3
                             learning.rate = 0.005, # M7 = 0.01
                             bag.fraction = 0.5,
                             site.weights = train_brt_sel$weight)

# saveRDS(brt_mod_se, "brt_mod_s.rds")


# M7 weighted presences with merged presences from both periods and only 2014-2019 absences
# fitting final gbm model with a fixed number of 4650 trees for kelp
# 
# mean total deviance = 1.386 
# mean residual deviance = 0.359 
# 
# estimated cv deviance = 0.601 ; se = 0.022 
# 
# training data correlation = 0.889 
# cv correlation =  0.793 ; se = 0.011 
# 
# training data AUC score = 0.984 
# cv AUC score = 0.946 ; se = 0.005 
# 
# elapsed time -  0.02 minutes 

# save.image("3.1.HSModeling_multi_approaches_M7.RData")

