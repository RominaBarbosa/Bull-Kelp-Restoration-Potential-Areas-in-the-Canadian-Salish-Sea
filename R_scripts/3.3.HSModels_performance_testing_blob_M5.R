###==================================================================
### Species Distribution models        SDMs          ################
###     Ensemble model  - Blob conditions            ################
### Testing performance (30% of dataset)             ################
###        Table and barplots                        ################
### Author: Romina Barbosa                           ################
### Last version: 03-Sep-2025                        ################
###==================================================================
library(tidyr)
library(stringr)



# Predicted rasters were masked by substrate availability to get the final prediction before evaluation of predictions 
# i.e., the testing points at non suitable substrate were not considered when evaluating the predictions accuracy to make it comparable with a
#the evaluation of the ensemble model

setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres")
# suit_t1 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres/tifs/ensemble_Average_suitability_blob_M5.tif")
# suit_t2 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M5/ensemble_ave_postBlob_M5.tif")

ens_average<- terra::rast("tifs/ensemble_Average_suitability_blob_M5.tif")
ens_pca<- terra::rast("tifs/ensemble_PCA_blob_M5.tif")
gam_pred_raster<- rast("tifs/gam_pred_raster.tif")
glm_pred_raster<- rast("tifs/glm_pred_raster.tif")
rf_pred_raster<- rast("tifs/rf_pred_raster.tif")
brt_pred_raster<- rast("tifs/brt_pred_raster.tif")

substrate_aligned<- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Aug2025/tifs/substrate_SOG_aligned.tif")

ens_average_masked <- terra::mask(ens_average, substrate_aligned)
ens_pca_masked <- terra::mask(ens_pca, substrate_aligned)

gam_pred_masked <- terra::mask(gam_pred_raster, substrate_aligned)
glm_pred_masked <- terra::mask(glm_pred_raster, substrate_aligned)
rf_pred_masked <- terra::mask(rf_pred_raster, substrate_aligned)
brt_pred_masked <- terra::mask(brt_pred_raster, substrate_aligned)


vars_selected<- c("ammonium_spring_mean",     "currentSpeed_summer_mean", "nitrate_winter_mean",      "PAR_summer_mean" ,
                  "temperature_summer_mean",  "turbidity_summer_mean",    "slope_5x5" )      


# train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_FINAL.csv")
# train_test_dataset<- (train_test_dataset[,-1])
# names(train_test_dataset)

train<- train_test_dataset %>%
  filter(set == "train")

test<- train_test_dataset %>%
  filter(set == "test")


train$kelp<- as.factor(train$kelp)
train_points <- vect(train, geom = c("x", "y"), crs = "EPSG:3005")

test$kelp<- as.factor(test$kelp)
test_points <- vect(test, geom = c("x", "y"), crs = "EPSG:3005")


plot(train_points, cex=0.4)
plot(test_points, cex=0.4, add=T, col="red")

train$ensemble_ave_pred <- terra::extract(ens_average_masked_all, train_points)[,2]
train$ensemble_ave_pred
# train$ensemble_wAve_pred <- terra::extract(ens_wAverage_masked_all, train_points)[,2]

### Load model performance table for training dataset =========
# models_performance<- read.csv("models_performance_calibration_Table.csv")

train_masked<- train[-which(is.na(train$ensemble_ave_pred)),]
roc_ave_ens <- pROC::roc(response = train_masked$kelp,
                         predictor = train_masked$ensemble_ave_pred)
plot(roc_ave_ens, col = "green", main = "ROC curve: Ensemble model")

# roc_wAve_ens <- pROC::roc(response = train_masked$kelp,
#                           predictor = train_masked$ensemble_wAve_pred)
# plot(roc_wAve_ens, col = "blue", main = "", add=T)

# tab_ens <- get_thresh_table(
#   roc_obj = roc_wAve_ens,
#   pred = train[-which(is.na(train$ensemble_wAve_pred)),"ensemble_wAve_pred"],
#   truth = train[-which(is.na(train$ensemble_wAve_pred)),"kelp"],
#   model_name = "Ensemble_wAve"
# )

ens_performance <- get_thresh_table(
  roc_obj = roc_ave_ens,
  pred = train[-which(is.na(train$ensemble_ave_pred)),"ensemble_ave_pred"],
  truth = train[-which(is.na(train$ensemble_ave_pred)),"kelp"],
  model_name = "Ensemble_Ave"
)

# With training data: 
# Model         Criterion     Threshold Sensitivity Specificity     TSS   AUC
# 1 Ensemble_wAve Max TSS         0.638         0.893     0.964   0.857   0.982
# 2 Ensemble_wAve No omission  -Inf             1         0       0       0.982
# 3 Ensemble_wAve 10% omission    0.00747       1         0.00331 0.00331 0.982
# 4 Ensemble_Ave  Max TSS         0.694         0.841     0.983   0.824   0.974
# 5 Ensemble_Ave  No omission  -Inf             1         0       0       0.974
# 6 Ensemble_Ave  10% omission    0.00707       1         0.00331 0.00331 0.974

# With testing dataset
# Model        Criterion     Threshold Sensitivity Specificity     TSS   AUC
# 1 Ensemble_Ave Max TSS         0.694         0.841     0.984   0.824   0.974
# 2 Ensemble_Ave No omission  -Inf             1         0       0       0.974
# 3 Ensemble_Ave 10% omission    0.00707       1         0.00327 0.00327 0.974

# ================================================================
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/tifs"
setwd(model_results_path)
dir()

# glm_mod_s<- readRDS("glm_mod_s.rds")
# gam_mod_s<- readRDS("gam_mod_s.rds")
# rf_mod_s<-  readRDS("rf_mod_s.rds")
# brt_mod_s<- readRDS("brt_mod_s.rds")


# Using thresholds derived from training
thresholds <- c("0.1%_Omission" = 0.2054824,  #  from training - # 0.1% omission 
                # "MaxF1" = 0.46898,
                "MaxTSS" = 0.6940952 )##  from training


# # Classify test points
# test_class <- cut(test_pred,
#                   breaks = c(-Inf, thr_0.1pct, thr_maxTSS, Inf),
#                   labels = c("Inadequate", "Adequate", "Optimal"))

# Function to evaluate predictions using specified thresholds
pred_vals <- terra::extract(pred_raster, test_pts)[,2]
truth_vals <- test_pts[[truth_col]]

# Remove NAs
valid_idx <- which(!is.na(pred_vals))
pred_vals <- pred_vals[valid_idx]
truth_vals <- truth_vals[valid_idx,]


evaluate_sdm_thresholds <- function(pred_raster, test_pts, truth_col = "kelp", thresholds, model_name) {
  # Extract predicted values at test points
  pred_vals <- terra::extract(pred_raster, test_pts)[,2]
  truth_vals <- test_pts[[truth_col]]
  
  # Remove NAs
  valid_idx <- which(!is.na(pred_vals))
  pred_vals <- pred_vals[valid_idx]
  truth_vals <- truth_vals[valid_idx,]
  
  # Compute AUC
  roc_obj <- pROC::roc(response = truth_vals, predictor = pred_vals)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Prepare results table
  results_list <- list()
  
  for (thr_name in names(thresholds)) {
    thr <- thresholds[[thr_name]]
    
    # Classify predictions
    pred_class <- as.integer(pred_vals >= thr)
    
    # Confusion matrix
    tp <- sum(pred_class == 1 & truth_vals == 1)
    tn <- sum(pred_class == 0 & truth_vals == 0)
    fp <- sum(pred_class == 1 & truth_vals == 0)
    fn <- sum(pred_class == 0 & truth_vals == 1)
    
    sensitivity <- if((tp + fn) > 0) tp / (tp + fn) else NA
    specificity <- if((tn + fp) > 0) tn / (tn + fp) else NA
    TSS_val     <- sensitivity + specificity - 1
    
    results_list[[thr_name]] <- data.frame(
      ThresholdName = thr_name,
      Threshold = thr,
      AUC = auc_val,
      Sensitivity = sensitivity,
      Specificity = specificity,
      TSS = TSS_val,
      Model= model_name
    )
  }
  
  # Combine results into one data.frame
  results_df <- do.call(rbind, results_list)
  return(results_df)
}


# Evaluate ensemble raster
gam_acc <- evaluate_sdm_thresholds(pred_raster = gam_pred_raster,
                                   test_pts = test_points,
                                   thresholds = thresholds,
                                   model_name= "GAM")

glm_acc <- evaluate_sdm_thresholds(pred_raster = glm_pred_raster,
                                   test_pts = test_points,
                                   thresholds = thresholds,
                                   model_name= "GLM")

rf_acc <- evaluate_sdm_thresholds(pred_raster = rf_pred_raster,
                                   test_pts = test_points,
                                   thresholds = thresholds,
                                  model_name= "RF")

brf_acc <- evaluate_sdm_thresholds(pred_raster = brt_pred_raster,
                                   test_pts = test_points,
                                   thresholds = thresholds,
                                   model_name= "BRF")

ens_acc <- evaluate_sdm_thresholds(pred_raster = ens_average_masked_all,
                                   test_pts = test_points,
                                   thresholds = thresholds,
                                   model_name= "ENS")

models_performance<- rbind(gam_acc, glm_acc, rf_acc, brf_acc, ens_acc)
# sensitivity = True Positive Rate (TPR), or recall

rownames(models_performance)<- NULL
#     ThresholdName Threshold       AUC Sensitivity Specificity       TSS Model
# 1  0.1%_Omission 0.2054824 0.8991285   0.9597070   0.6064982 0.5662052   GAM
# 2         MaxTSS 0.6940952 0.8991285   0.7802198   0.8267148 0.6069346   GAM
# 3  0.1%_Omission 0.2054824 0.8718874   0.9706960   0.5162455 0.4869415   GLM
# 4         MaxTSS 0.6940952 0.8718874   0.7179487   0.8447653 0.5627141   GLM
# 5  0.1%_Omission 0.2054824 0.9207495   0.9780220   0.5667870 0.5448090    RF
# 6         MaxTSS 0.6940952 0.9207495   0.7765568   0.8916968 0.6682535    RF
# 7  0.1%_Omission 0.2054824 0.9136351   0.9670330   0.6064982 0.5735312   BRF
# 8         MaxTSS 0.6940952 0.9136351   0.7802198   0.8628159 0.6430357   BRF
# 9  0.1%_Omission 0.2054824 0.9075993   0.9769231   0.5241935 0.5011166   ENS
# 10        MaxTSS 0.6940952 0.9075993   0.7692308   0.8306452 0.5998759   ENS


models_performance%>% filter(ThresholdName=="MaxTSS")
models_performance%>% filter(ThresholdName=="0.1%_Omission")
# ThresholdName Threshold       AUC Sensitivity Specificity       TSS Model
# 1 0.1%_Omission 0.2054824 0.8991285   0.9597070   0.6064982 0.5662052   GAM
# 2 0.1%_Omission 0.2054824 0.8718874   0.9706960   0.5162455 0.4869415   GLM
# 3 0.1%_Omission 0.2054824 0.9207495   0.9780220   0.5667870 0.5448090    RF
# 4 0.1%_Omission 0.2054824 0.9136351   0.9670330   0.6064982 0.5735312   BRF
# 5 0.1%_Omission 0.2054824 0.9075993   0.9769231   0.5241935 0.5011166   ENS

ggplot(models_performance, aes(x = Model, y = AUC, fill = Model)) +
  geom_col(position = "dodge", width = 0.7) +
  # facet_wrap(~ ThresholdName) +
  labs(
    title = "Model performance by threshold type",
    y = "True Skill Statistic (TSS)",
    x = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )


ggplot(models_performance, aes(x = Model, y = TSS, fill = Model)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ ThresholdName) +
  labs(
    title = "Model performance by threshold type",
    y = "True Skill Statistic (TSS)",
    x = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )






# Compute confusion-matrix-derived metrics given binary classification and truth
.cm_metrics <- function(pred_class, truth) {
  tp <- sum(pred_class == 1 & truth == 1, na.rm = TRUE)
  tn <- sum(pred_class == 0 & truth == 0, na.rm = TRUE)
  fp <- sum(pred_class == 1 & truth == 0, na.rm = TRUE)
  fn <- sum(pred_class == 0 & truth == 1, na.rm = TRUE)
  N  <- tp + tn + fp + fn
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall      <- sensitivity
  f1          <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
  tss         <- if (!is.na(sensitivity) && !is.na(specificity)) sensitivity + specificity - 1 else NA_real_
  # Kappa calculation:
  if (N > 0) {
    # Observed accuracy
    Po <- (tp + tn) / N
    # Expected accuracy under independence:
    row1 <- (tp + fp) / N  # predicted positive proportion
    row2 <- (fn + tn) / N  # predicted negative proportion
    col1 <- (tp + fn) / N  # true positive proportion
    col2 <- (fp + tn) / N  # true negative proportion
    Pe <- row1 * col1 + row2 * col2
    kappa <- if ((1 - Pe) != 0) (Po - Pe) / (1 - Pe) else NA_real_
  } else {
    kappa <- NA_real_
  }
  c(tp = tp, tn = tn, fp = fp, fn = fn,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    TSS = tss,
    Kappa = kappa)
}

# Main function
evaluate_sdm_metrics <- function(pred_raster, test_pts, truth_col, thresholds, model_name) {
  # Extract predictions at testing points
  pred_vals <- terra::extract(pred_raster, test_pts)[,2] |> as.vector()
  truth <- test_pts[[truth_col]]$kelp 
  
  # Handle NAs (your choice: replace or remove)
  pred_vals[is.na(pred_vals)] <- 0.02
  truth <- factor(truth, levels = c(0,1))
  
  # Now run ROC
  roc_obj <- pROC::roc(response = truth, predictor = pred_vals)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Remove NAs
  # valid_idx <- which(!is.na(pred_vals) )
  # pred_vals <- pred_vals[valid_idx]
  # truth     <- truth[valid_idx,]
  pred_vals[which(is.na(pred_vals))]<- 0.02
  
  # Initialize results
  results <- data.frame()
  
  # Compute AUC once (independent of threshold)
  roc_obj <- pROC::roc(response = truth, predictor = pred_vals)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Loop through thresholds
  for (thr_name in names(thresholds)) {
    thr <- thresholds[[thr_name]]
    
    pred_class <- as.integer(pred_vals >= thr)
    
    tp <- sum(pred_class == 1 & truth == 1)
    tn <- sum(pred_class == 0 & truth == 0)
    fp <- sum(pred_class == 1 & truth == 0)
    fn <- sum(pred_class == 0 & truth == 1)
    
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
    precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA
    f1          <- if (!is.na(precision) && !is.na(sensitivity) && 
                       (precision + sensitivity) > 0) {
      2 * precision * sensitivity / (precision + sensitivity)
    } else NA
    tss <- if (!is.na(sensitivity) && !is.na(specificity)) {
      sensitivity + specificity - 1
    } else NA
    kappa <- vcd::Kappa(table(factor(truth, levels=c(0,1)),
                              factor(pred_class, levels=c(0,1))))$Unweighted[1]
    
    results <- rbind(results, data.frame(
      ThresholdName = thr_name,
      Threshold     = thr,
      AUC           = auc_val,
      Sensitivity   = sensitivity,
      Specificity   = specificity,
      Precision     = precision,
      F1            = f1,
      TSS           = tss,
      Kappa         = kappa,
      model_name    =model_name
    ))
  }
  
  return(results)
}



# Suppose ensemble_raster is your predicted SpatRaster
# and test_points is a SpatVector with column 'kelp' (0/1).
# And you have training-derived thresholds:
thr_train<- thresholds

res <- evaluate_sdm_metrics(pred_raster = ens_average_masked_all,
                            test_pts = test_points,
                            truth_col = "kelp",
                            thresholds = thr_train,
                            model_name= "GAM")
print(res)

# Evaluate ensemble raster
gam_acc <- evaluate_sdm_metrics(pred_raster = gam_pred_raster,
                                   test_pts = test_points,
                                   truth_col = "kelp",
                                   thresholds = thresholds,
                                   model_name= "GAM")

glm_acc <- evaluate_sdm_metrics(pred_raster = glm_pred_raster,
                                   test_pts = test_points,
                                truth_col = "kelp",
                                   thresholds = thresholds,
                                   model_name= "GLM")

rf_acc <- evaluate_sdm_metrics(pred_raster = rf_pred_raster,
                                  test_pts = test_points,
                               truth_col = "kelp",
                                  thresholds = thresholds,
                                  model_name= "RF")

brf_acc <- evaluate_sdm_metrics(pred_raster = brt_pred_raster,
                                   test_pts = test_points,
                                truth_col = "kelp",
                                   thresholds = thresholds,
                                   model_name= "BRF")

ens_acc <- evaluate_sdm_metrics(pred_raster = ens_average_masked_all,
                                   test_pts = test_points,
                                   thresholds = thresholds,
                                   truth_col = "kelp",
                                   model_name= "ENS")

models_performance<- rbind(gam_acc, glm_acc, rf_acc, brf_acc, ens_acc)
rownames(models_performance)<- NULL
models_performance%>%
  filter(ThresholdName == "MaxTSS")
# ThresholdName Threshold       AUC Sensitivity Specificity Precision        F1       TSS     Kappa model_name
# 1        MaxTSS 0.6940952 0.8895854   0.7689531   0.8267148 0.8160920 0.7918216 0.5956679 0.5956679        GAM
# 2        MaxTSS 0.6940952 0.8612780   0.7075812   0.8447653 0.8200837 0.7596899 0.5523466 0.5523466        GLM
# 3        MaxTSS 0.6940952 0.9094866   0.7653430   0.8916968 0.8760331 0.8169557 0.6570397 0.6570397         RF
# 4        MaxTSS 0.6940952 0.9011717   0.7689531   0.8628159 0.8486056 0.8068182 0.6317690 0.6317690        BRF
# 5        MaxTSS 0.6940952 0.9187465   0.7220217   0.9241877 0.9049774 0.8032129 0.6462094 0.6462094        ENS


# 1. Discrimination ability (overall predictive power)
# These metrics assess how well the model separates presences from absences across the full range of thresholds.
# AUC (Area Under the ROC Curve)
# TSS (True Skill Statistic)
# Kappa (Cohen’s Kappa, sometimes grouped with agreement metrics)
# 
# 2. Classification accuracy (at a chosen threshold)
# These look at the confusion matrix after applying a threshold.
# Sensitivity (a.k.a. recall, true positive rate) → ability to detect presences
# Specificity (true negative rate) → ability to reject absences
# Precision (positive predictive value) → how reliable predicted presences are
#
# F1 score (harmonic mean of precision and recall) → balance between detecting presences and minimizing false positives
# 3. Agreement metrics (chance-corrected performance)
# These account for how much agreement is beyond what would be expected by chance.
# Kappa (if you keep it separate)

# Plots
df <- models_performance

# Reshape data to long format for plotting metrics
df$model_name <- factor(df$model_name, levels = c("GLM", "GAM", "RF", "BRF", "ENS"))

df_long <- df %>%
  pivot_longer(cols = c(Sensitivity, Specificity, Precision, F1, TSS, Kappa),
               names_to = "Metric", values_to = "Value")

# Plot
ggplot(df_long, aes(x = model_name, y = Value, fill = ThresholdName)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Performance metrics of SDMs under different thresholds",
       y = "Metric value", x = "Model") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/models_testing_metrics_plot.png", width = 17, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/models_testing_metrics_plot.pdf", width = 17, height = 12, dpi= 300, units="cm")

