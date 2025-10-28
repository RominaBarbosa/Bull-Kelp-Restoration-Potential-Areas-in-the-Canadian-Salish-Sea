###==================================================================
### Species Distribution models        SDMs          ################
### Model Predictions - thresholds and accuracy test ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Date: 24-Sep-2025                                ################
### Last edition: 24-Sep-2025                        ################
###==================================================================
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf) 

model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
postblob_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob"

### Load model predictions =====================================================
ens_average<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif", sep="/"))
ens_average_t2<- rast(paste(postblob_path, "M7/ens_ave_postblob_M7.tif", sep="/"))
names(ens_average)<- "ens_average"
names(ens_average_t2)<- "ens_average_postblob"

# #===============================================================================
# ### Create Substrate mask            ===========================================
# substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/SOG_substrate_20m.tif")
# substrate_west<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/WCVI_substrate_20m.tif")
# substrate_north<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/QCS_substrate_20m.tif")
# 
# substrate<- merge(substrate, substrate_north, substrate_west)
# plot(substrate)
# # The predicted raster files are classified as follows:
# # 1) Rock,
# # 2) Mixed,
# # 3) Sand,
# # 4) Mud
# 
# 
# # Mask substrate to study area
# substrate<- crop(substrate, ens_average)
# 
# # Align substrate layer to model prediction layer
# substrate_aligned <- terra::rast(ens_average)
# substrate_aligned <- terra::resample(substrate, substrate_aligned, method = "near")
# 
# # Convert categories of substrate in 1 (hard) and 2 (soft) substrate
# # substrate_aligned[substrate_aligned == 1]<- 1 # rocky is already 1
# substrate_aligned[substrate_aligned == 2]<- 1 # 2 was mixed substrate
# substrate_aligned[substrate_aligned == 3]<- 2 # 3 and 4 were mud and sand (I guess)
# substrate_aligned[substrate_aligned == 4]<- 2 # 3 and 4 were mud and sand (I guess)
# 
# # Save substrate, 1= hard substrate; 2= soft substrate
# # writeRaster(substrate_aligned, "substrate_SOG_aligned.tif", overwrite=T)


# ==============================================================================
#  Mask ensemble model predictions by substrate
# ==============================================================================
substrate_aligned<- rast(paste(model_results_path, "substrate_SOG_aligned.tif", sep="/"))
substrate_aligned[substrate_aligned == 2]<- NA
ens_averaget1_masked <- terra::mask(ens_average, substrate_aligned)
names(ens_averaget1_masked)<- "ens_average_t1"
ens_averaget2_masked <- terra::mask(ens_average_t2, substrate_aligned)
names(ens_averaget2_masked)<- "ens_average_t2"
# writeRaster(ens_averaget1_masked, "ens_average_masked_M6.tif", overwrite=T)
# writeRaster(ens_averaget2_masked, "ens_average_masked_M6.tif", overwrite=T)

# ==============================================================================
#  Identify thresholds based on training dataset and ensemble model (blob)
# ==============================================================================
# Load training and testing dataset 
traintest<- read.csv(paste(model_results_path,"training_testing_datasets_blob_M7.csv", sep="/"))
# train<- read.csv(paste(model_results_path,"Sep2025_M6_weightedPres/train_selected_table_FINALMODELS.csv", sep="/"))
traintest<- traintest[,-1]
colnames(traintest)
cols_selected<- c("kelp", "x", "y", "period", "set")
traintest<- traintest[,colnames(traintest)%in%cols_selected]
# traintest<- #traintest[,c(1,2,3,17,15,4)]

train<- traintest%>%
  filter(set=="train")

train_points<- vect(train, c("x","y"))
  
train$ensemble_t1pred <- terra::extract(ens_average, train_points)[,2]
train$ensemble_t2pred <- terra::extract(ens_average_t2, train_points)[,2]


### Get thresholds  based on training dataset =====
roc_train <- pROC::roc(response = train$kelp,
                         predictor = train$ensemble_t1pred)
plot(roc_train, col = "green", main = "ROC curve: Ensemble model")

# Find Max TSS threshold
thr_maxTSS <- roc_train$thresholds[which.max(roc_train$sensitivities + roc_train$specificities - 1)]# Using “drop point” from cumulative coverage or density
# [1] 0.5748859
# M7: 0.4733616

# Predicted suitability for presences in training data
train_ave_ens<- train

# Density estimates
pres_pred <- train_ave_ens$ensemble_t1pred[train_ave_ens$kelp == 1]
abs_pred  <- train_ave_ens$ensemble_t1pred[train_ave_ens$kelp == 0]

# Compute density estimates
dens_pres <- density(pres_pred, na.rm = TRUE)
dens_abs  <- density(abs_pred,  na.rm = TRUE)

# Build data frame for ggplot
df_plot <- data.frame(
  x = c(dens_abs$x, dens_pres$x),
  y = c(dens_abs$y, dens_pres$y),
  group = rep(c("Absences", "Presences"), 
              times = c(length(dens_abs$x), length(dens_pres$x)))
)

# Thresholds
threshold_drop <- quantile(pres_pred, 0.001, na.rm = TRUE)  # 0.1% presence omission
# 0.1% 
# 0.2785257 

#M7
# 0.1% 
# 0.1079078

threshold_drop <- quantile(pres_pred, 0.01, na.rm = TRUE)  # 0.1% presence omission
# M7: 1% ==> 0.2414677 

# threshold_drop <- quantile(pres_pred, 0.009, na.rm = TRUE)  # 0.9% presence omission = 0.2724581 


### Plot thresholds  ====
# Define label positions (x = middle of each region, y = top 90% of max density)
y_max <- 3.7#max(df_plot$y, na.rm = TRUE)
labels <- data.frame(
  label = c("Unsuitable", "Adequate", "Optimal"),
  x = c(threshold_drop / 2,
        (threshold_drop + thr_maxTSS) / 2,
        thr_maxTSS + (max(df_plot$x) - thr_maxTSS)/2),
  y = rep(0.9 * y_max, 3)
)
# labels[1,2]<- 

ggplot(df_plot, aes(x = x, y = y, color = group, fill = group)) +
  
  # Shaded suitability areas
  # geom_rect(aes(xmin = -Inf, xmax = threshold_drop, ymin = -Inf, ymax = Inf),
  #           fill = "blue", alpha = 0.01, inherit.aes = FALSE) +  # unsuitable
  # geom_rect(aes(xmin = threshold_drop, xmax = thr_maxTSS, ymin = -Inf, ymax = Inf),
  #           fill = "orange", alpha = 0.01, inherit.aes = FALSE) +      # moderate
  # geom_rect(aes(xmin = thr_maxTSS, xmax = Inf, ymin = -Inf, ymax = Inf),
  #           fill = "red", alpha = 0.01, inherit.aes = FALSE) +       # optimal
  geom_line(size = 1.2) +
  # geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.2, color = NA) + # adds semi-transparent fill
  geom_vline(xintercept = threshold_drop, color = "orange", linetype = "dashed", size = 1) +
  geom_vline(xintercept = thr_maxTSS, color = "red", linetype = "dashed", size = 1) +
  # geom_vline(xintercept = thr_maxF1, color = "darkorange", linetype = "dashed", size = 1) +
  # Labels for shaded regions
  lims(y=c(0, 3.5), x=c(0, 1))+
  geom_text(data = labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, color = "black", fontface = "bold") +
  labs(
    x = "Predicted Kelp Habitat Suitability",
    y = "Density",
    color = "",
    fill = ""
  ) +
  scale_color_manual(values = c("black", "grey")) +  # absences = blue, presences = red
  scale_fill_manual(values = c("black", "grey")) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.8, 0.7),
    legend.background = element_rect(fill = alpha("white", 0.1), color = NA),
    legend.key = element_blank()
  ) +
  annotate("text", x = threshold_drop, y = max(df_plot$y)*0.85,
           label = "0.1% Omission", color = "orange", angle = 90, vjust = -0.5) +
  annotate("text", x = thr_maxTSS, y = max(df_plot$y)*0.85,
           label = "Max TSS", color = "red", angle = 90, vjust = -0.5)

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Threshold_selection_plot_M7.pdf", width = 12, height = 8, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Threshold_selection_plot_M7.png", width = 12, height = 8, dpi= 300, units="cm")






source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")


tab_ensemble_thresholds <- get_thresh_table(
  roc_obj = roc_ave_ens,
  pred = train[which(!is.na(train$ensemble_t1pred)),"ensemble_t1pred"],
  truth = train[which(!is.na(train$ensemble_t1pred)),"kelp"],
  model_name = "Ensemble_ave"
)

#   Model        Criterion     Threshold Sensitivity Specificity     TSS   AUC
# 1 Ensemble_ave Max TSS         0.467         0.965     0.829   0.794   0.913
# 2 Ensemble_ave No omission  -Inf             1         0       0       0.913
# 3 Ensemble_ave 10% omission    0.00730       1         0.00724 0.00724 0.913 # it was calculated based on ROC curve and changed to quantile based for M7

#M7
# Model        Criterion   Threshold Sensitivity Specificity   TSS   AUC
# 1 Ensemble_ave Max TSS         0.480       0.918       0.927 0.845 0.930
# 2 Ensemble_ave No omission  -Inf           1           0     0     0.930
# 3 Ensemble_ave 1% omission     0.241       0.989       0.736 0.726 0.930

write.csv(tab_ensemble_thresholds, "Sep2025_M7_weightedPres/tab_ensemble_thresholds_M7.csv")


# ==============================================================================
#      Convert predictions in categories based on thresholds and 
#     evaluate its accuracy in both periods: 2014-2019 and 2020-2022
# ==============================================================================
### Extract ensemble predictions at testing points =============================
# traintest<- read.csv(paste(model_results_path,"training_testing_datasets_blob_M7.csv", sep="/"))
# traintest<- traintest[,-1]
# traintest<- traintest[,c(1,2,3,17,15,4)]

test<- traintest%>%
  filter(set=="test")

test_points<- vect(test, c("x", "y"))

test$ensemble_t1pred <- terra::extract(ens_average, test_points)[,2]
test$ensemble_t2pred <- terra::extract(ens_average_t2, test_points)[,2]
test$ensemble_t1masked <- terra::extract(ens_averaget1_masked, test_points)[,2]

# Testing period 2
# test_t2<- read.csv()
# test_t2_points<- vect(test_t2, c("x", "y"))
# 
# test_t2$ensemble_t1pred <- terra::extract(ens_average, test_t2_points)[,2]
# test_t2$ensemble_t2pred <- terra::extract(ens_average_t2, test_t2_points)[,2]

# Compute ROC for ensemble in both predicted periods based on testing points (30% blob and all points postblob)
roc_ave_ens <- pROC::roc(response = test$kelp,
                         predictor = test$ensemble_t1pred)
plot(roc_ave_ens, col = "green", main = "ROC curve: Ensemble model")

# roc_pca_ens2 <- pROC::roc(response = test_t2$kelp,
#                           predictor = test_t2$ensemble_t2pred)
# plot(roc_pca_ens2, col = "blue", main = "", add=T)



# --- Function to calculate evaluation metrics (threshold + AUC + TSS) ---
eval_metrics <- function(obs, pred, threshold) {
  # Convert predictions to binary using threshold
  pred_bin <- ifelse(pred >= threshold, 1, 0)
  
  # Confusion matrix elements
  TP <- sum(obs == 1 & pred_bin == 1)
  TN <- sum(obs == 0 & pred_bin == 0)
  FP <- sum(obs == 0 & pred_bin == 1)
  FN <- sum(obs == 1 & pred_bin == 0)
  
  # Metrics
  sensitivity <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))  # True Positive Rate
  specificity <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))  # True Negative Rate
  accuracy    <- (TP + TN) / (TP + TN + FP + FN)
  
  # Cohen's Kappa
  exp_acc <- (((TP+FP)*(TP+FN)) + ((FN+TN)*(FP+TN))) / ((TP+TN+FP+FN)^2)
  kappa   <- (accuracy - exp_acc) / (1 - exp_acc)
  
  # True Skill Statistic (TSS)
  TSS <- sensitivity + specificity - 1
  
  # AUC (threshold-independent)
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(obs, pred)))
  }, error = function(e) NA)
  
  return(data.frame(
    Threshold   = threshold,
    Accuracy    = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Kappa       = kappa,
    TSS         = TSS,
    AUC         = auc_val
  ))
}

# --- Evaluate each model ---
str(test)
test$kelp<- as.factor(test$kelp)
test<- na.exclude(test)

ens_metrics_thr1 <- eval_metrics(obs= test$kelp, pred= test$ensemble_t1pred, threshold= thr_maxTSS)
ens_metrics_thr2 <- eval_metrics(obs= test$kelp, pred= test$ensemble_t1pred, threshold= threshold_drop)

rbind(ens_metrics_thr1, ens_metrics_thr2)

ensMasked_metrics_thr1 <- eval_metrics(obs= test$kelp, pred= test$ensemble_t1masked, threshold= thr_maxTSS)
ensMasked_metrics_thr2 <- eval_metrics(obs= test$kelp, pred= test$ensemble_t1masked, threshold= threshold_drop)



# # --- Combine results into one table ---
# all_metrics <- rbind(
#   GLM = glm_metrics,
#   GAM = gam_metrics,
#   RF  = rf_metrics,
#   BRT = brt_metrics,
#   Ensemble = ens_metrics
# )

print(all_metrics)













# ### Evaluate and select threshold for binary results =====
# # Get metrics of model performance, AUC, Sensitivity and specificity
# results_selected <- bind_rows(
#   get_metrics_optimized2(model= glm_mod_s, test_data=test_sel, scale_params= scaling_params,
#                          model_name = "glm", threshold_type = "youden"),
#   get_metrics_optimized2(gam_mod_s, test_sel, scale_params= scaling_params_2, "gam",
#                          threshold_type = "youden"),#10pct_omission
#   get_metrics_optimized2(rf_mod_s, test_sel, "rf", threshold_type = "youden"),
#   get_metrics_optimized2(brt_mod_s, test_sel, "brt", threshold_type = "youden")
# )
# 
# 
# 
# 
# # Model 1 (includes SD_summer_salinity)
# # Model ThresholdType Threshold   AUC Sensitivity Specificity   TSS
# # 1 glm   youden            0.609 0.878       0.809       0.801 0.610
# # 2 gam   youden            0.645 0.910       0.798       0.872 0.670
# # 3 rf    youden            0.521 1           1           1     1    
# # 4 brt   youden            0.621 0.966       0.863       0.934 0.797
# 
# # Model 2 (exncludes SD_summer_salinity, changed Nitrate_summer_mean by Nitrate_summer_minimum (correlated)
# # Model ThresholdType Threshold   AUC Sensitivity Specificity   TSS
# # 1 glm   youden            0.621 0.888       0.826       0.812 0.638
# # 2 gam   youden            0.546 0.915       0.854       0.827 0.681
# # 3 rf    youden            0.508 1           1           1     1    
# # 4 brt   youden            0.522 0.959       0.901       0.894 0.795
# 
# # Model 4 (exncludes SD_summer_salinity, changed Nitrate_summer_mean by Nitrate_summer_minimum (correlated)
# # Removed presence and abcence points at soft substrate
# #   Model ThresholdType T hreshold  AUC Sensitivity Specificity   TSS
# # 1 glm   youden            0.637 0.900       0.825       0.848 0.673
# # 2 gam   youden            0.502 0.918       0.881       0.820 0.701
# # 3 rf    youden            0.5   1           1           1     1    
# # 4 brt   youden            0.560 0.979       0.932       0.931 0.864
# 
# 
# # Pivot the results for plotting
# results_long_selected <- results_selected %>%
#   pivot_longer(cols = c("AUC", "Sensitivity", "Specificity", "TSS"),
#                names_to = "Metric", values_to = "Value")
# 
# thresholds_table<- results_selected%>%
#   select(Model, Threshold)
# # Model Threshold
# # 1 glm       0.609
# # 2 gam       0.645
# # 3 rf        0.521
# # 4 brt       0.621
# 
# # M3:
# # Model Threshold
# # 1 glm       0.637
# # 2 gam       0.502
# # 3 rf        0.5  
# # 4 brt       0.560


