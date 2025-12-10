###==================================================================
### Species Distribution models        SDMs          ################
###     Ensemble model  - Blob conditions            ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Last version: 11-Aug-2025                        ################
###==================================================================
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf)

source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")
source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/plot_response_curves.R")

# This chunk was run and its outputs are stored in a RDS file
# # Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_FINAL2.tif")
names(raster_stack_20m)
terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)
terrain_vars<-  rast(tif_files[c(12)]) #,5,12,13,15
names(terrain_vars)
summary(terrain_vars)

# Merge rasters of all selected variables including terrain and NEMO
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)

# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask15_40<- bathy20m
bathy20m_mask15_40[!(bathy20m_mask15_40 >= -15 & bathy20m_mask15_40 <= 40)] <- NA # negative values are in land
plot(bathy20m_mask15_40)

# Mask all layers
raster_stack_predict<- raster_stack_20m_all[[vars_selected]] # subset selected variables
names(raster_stack_predict)
plot(raster_stack_predict[[1]])
raster_stack_predict<- mask(raster_stack_predict, bathy20m_mask15_40)
plot(raster_stack_predict[[1]])

saveRDS(raster_stack_predict, "raster_stack_predict.rds")


source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/plot_response_curves.R")

# ================================================================
# STEP 8: Predict in the entire area
# ================================================================
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025"
setwd(model_results_path)
dir()

glm_mod_s<- readRDS("glm_mod_s.rds")
gam_mod_s<- readRDS("gam_mod_s.rds")
rf_mod_s<-  readRDS("rf_mod_s.rds")
brt_mod_s<- readRDS("brt_mod_s.rds")

models <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")

# vars_selected<- c( "slope_7x7", "ammonium_spring_SD", "ammonium_winter_mean","PAR_summer_mean", "temperature_summer_mean",
#  "turbidity_summer_mean","salinity_summer_SD", "nitrate_summer_mean", "ammonium_summer_mean", "ammonium_spring_mean")   

# Model 1:
# vars_selected<- c( "temperature_summer_mean",  "slope_5x5" ,  "turbidity_summer_mean",
#                    "ammonium_spring_SD", "PAR_summer_mean", "currentSpeed_summer_mean", "salinity_summer_SD")      

# Model 2: added nitrate_winter_mean and nitrate_summer_minimum; replaced ammonium_spring_SD by ammonium_spring_mean; removed salinity_summer_SD
vars_selected<- c("slope_5x5",  "nitrate_summer_minimum", "turbidity_summer_mean", "nitrate_winter_mean", "temperature_summer_mean", 
                  "ammonium_spring_mean", "currentSpeed_summer_mean", "PAR_summer_mean")


train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_FINAL2.csv")
train_test_dataset<- (train_test_dataset[,-1])
names(train_test_dataset)

train<- train_test_dataset %>%
  filter(set == "train")

test<- train_test_dataset %>%
  filter(set == "test")


# Refit models using only selected vars
train_sel <- train %>% select(all_of(c("kelp", vars_selected)))
test_sel <- test %>% select(all_of(c("kelp", vars_selected)))

train_sel$kelp<- as.factor(train_sel$kelp)
test_sel$kelp<- as.factor(test_sel$kelp)

scaling_params_2 <- train_sel %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd), na.rm = TRUE))

train_sel_scaled <- scale_with_params(train_sel, scaling_params_2)
test_sel_scaled  <- scale_with_params(test_sel,  scaling_params_2)


# Loads raster inputs to predict the model
# raster_stack_predict<- readRDS("raster_stack_predict.rds")
# 
# # We need to scale the variables for predictions by considering the scaling parameters used when callibrating the model
# train_scaled_b <- scale(train_sel[, vars_selected])
# train_means <- attr(train_scaled_b, "scaled:center")
# train_sds   <- attr(train_scaled_b, "scaled:scale")
# 
# raster_stack_predict_scaled<- raster_stack_predict  # original raster stack
# 
# for (v in vars_selected) {
#   raster_stack_predict_scaled[[v]] <- (raster_stack_predict[[v]] - train_means[v]) / train_sds[v]
# }
# 
# names(raster_stack_predict_scaled)#<- names(raster_stack_predict)
# saveRDS(raster_stack_predict_scaled, "raster_stack_predict_scaled.rds")


raster_stack_predict<- readRDS("raster_stack_predict.rds")
raster_stack_predict_scaled<- readRDS("raster_stack_predict_scaled.rds")
# names(raster_stack) <- c("temperature_summer_mean",
#                          "slope_5x5",
#                          "turbidity_summer_mean",
#                          "ammonium_spring_SD",
#                          "PAR_summer_mean",
#                          "currentSpeed_summer_mean",
#                          "salinity_summer_SD")


# Predict probabilities
terra::predict(raster_stack_predict_scaled, #glm took a lot of memory so I saved directly to disk
                                  glm_mod_s,
                                  type = "response",
                                  na.rm = TRUE,
                                  filename = "tifs_Sep2025/glm2_pred_raster.tif",  # output written directly to disk
                           overwrite = TRUE)
glm_pred_raster <- rast("tifs_Sep2025/glm2_pred_raster.tif")
names(glm_pred_raster)<- "glm"

terra::predict(raster_stack_predict_scaled, #glm took a lot of memory so I saved directly to disk
               gam_mod_s,
               type = "response",
               na.rm = TRUE,
               filename = "tifs_Sep2025/gam2_pred_raster.tif",  # output written directly to disk
               overwrite = TRUE)

gam_pred_raster <- rast("tifs_Sep2025/gam2_pred_raster.tif")
names(gam_pred_raster)<- "gam"
rf_pred_raster  <- predict_raster_prob(rf_mod_s,  raster_stack_predict, "rf")
names(rf_pred_raster)<- "rf"
brt_pred_raster <- predict_raster_prob(brt_mod_s, raster_stack_predict, "brt")
names(brt_pred_raster)<- "brt"

### SAVE RESULT RASTERS =======================================================
# writeRaster(glm_pred_raster, "tifs_Sep2025/glm_pred_raster.tif", overwrite = TRUE)
# writeRaster(gam_pred_raster, "tifs_Sep2025/gam_pred_raster.tif", overwrite = TRUE)
# writeRaster(rf_pred_raster,  "tifs_Sep2025/rf2_pred_raster.tif",  overwrite = TRUE)
# writeRaster(brt_pred_raster, "tifs_Sep2025/brt2_pred_raster.tif", overwrite = TRUE)

setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025")
# glm_pred_raster<- rast("glm_pred_raster.tif")
# gam_pred_raster<- rast("gam_pred_raster.tif")
# rf_pred_raster<- rast("rf_pred_raster.tif")
# brt_pred_raster<- rast("brt_pred_raster.tif")

# ==============================================================================
# STEP 9: Binary results & Evaluation
# ==============================================================================
# Get metrics of model performance, AUC, Sensitivity and specificity

results_selected <- bind_rows(
  get_metrics_optimized2(model= glm_mod_s, test_data=train_sel, scale_params= scaling_params_2,
                         model_name = "glm", threshold_type = "youden"),
  get_metrics_optimized2(gam_mod_s, train_sel, scale_params= scaling_params_2, "gam",
                         threshold_type = "youden"),#10pct_omission
  get_metrics_optimized2(rf_mod_s, train_sel, "rf", threshold_type = "youden"),
  get_metrics_optimized2(brt_mod_s, train_sel, "brt", threshold_type = "youden")
)

# Model 1 (includes SD_summer_salinity)
# Model ThresholdType Threshold   AUC Sensitivity Specificity   TSS
# 1 glm   youden            0.609 0.878       0.809       0.801 0.610
# 2 gam   youden            0.645 0.910       0.798       0.872 0.670
# 3 rf    youden            0.521 1           1           1     1    
# 4 brt   youden            0.621 0.966       0.863       0.934 0.797

# Model 2 (exncludes SD_summer_salinity, changed Nitrate_summer_mean by Nitrate_summer_minimum (correlated)
# Model ThresholdType Threshold   AUC Sensitivity Specificity   TSS
# 1 glm   youden            0.621 0.888       0.826       0.812 0.638
# 2 gam   youden            0.546 0.915       0.854       0.827 0.681
# 3 rf    youden            0.508 1           1           1     1    
# 4 brt   youden            0.522 0.959       0.901       0.894 0.795

# Pivot the results for plotting
results_long_selected <- results_selected %>%
  pivot_longer(cols = c("AUC", "Sensitivity", "Specificity", "TSS"),
               names_to = "Metric", values_to = "Value")

thresholds_table<- results_selected%>%
  select(Model, Threshold)
# Model Threshold
# 1 glm       0.609
# 2 gam       0.645
# 3 rf        0.521
# 4 brt       0.621

# Convert continuous predictions into binary based on max sensitivity & specificity threshold:
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025")

# glm_bin <- binarize_raster(glm_pred_raster, "glm", thresholds_table)
# gam_bin <- binarize_raster(gam_pred_raster, "gam", thresholds_table)
# rf_bin  <- binarize_raster(rf_pred_raster, "rf", thresholds_table)
# brt_bin <- binarize_raster(brt_pred_raster, "brt", thresholds_table)

## SAVE RESULT RASTERS 
# writeRaster(glm_bin, "tifs/glm_bin.tif", overwrite = TRUE)
# writeRaster(gam_bin, "tifs/gam_bin.tif", overwrite = TRUE)
# writeRaster(rf_bin,  "tifs/rf_bin.tif",  overwrite = TRUE)
# writeRaster(brt_bin, "tifs/brt_bin.tif", overwrite = TRUE)

# glm_bin<- rast("tifs/glm_bin.tif")
# gam_bin<- rast("tifs/gam_bin.tif")
# rf_bin<- rast("tifs/rf_bin.tif")
# brt_bin<- rast("tifs/brt_bin.tif")

# # =================================================================
# # STEP 10: Plot predictions in the entire area
# # =================================================================
# library(sf)
# library(purrr)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(patchwork)
# 
# # Get country boundaries and coastline
# world <- ne_countries(scale = 10, returnclass = "sf")
# coastline <- ne_coastline(scale = 10, returnclass = "sf")
# crs_target <- "EPSG:3005"
# crs(glm_pred_raster) <- crs_target
# study_extent <- c(xmin = 984960, xmax = 1267240, ymin = 332520, ymax = 666640)
# world_proj <- st_transform(world, crs(glm_pred_raster))
# coastline_proj <- st_transform(coastline, crs(glm_pred_raster))
# 
# # Load coastline and land basemap (if not already loaded)
# coastline <- ne_coastline(scale = "large", returnclass = "sf")
# coastline <- st_transform(coastline, st_crs(glm_pred_raster))
# land_highres <- ne_download(scale = "large", type = "land", category = "physical", returnclass = "sf")
# land_highres <- st_transform(land_highres, st_crs(glm_pred_raster))
# 
# # plots
# glm_binary_map<- plot_raster_gg(glm_bin, title= "GLM Binary Prediction")
# gam_binary_map<- plot_raster_gg(gam_bin, title= "GAM Binary Prediction")
# rf_binary_map<- plot_raster_gg(rf_bin,  title= "RF Binary Prediction")
# brt_binary_map<- plot_raster_gg(brt_bin, title= "BRT Binary Prediction")
# 
# cowplot::plot_grid(glm_binary_map, gam_binary_map, rf_binary_map, brt_binary_map, nrow = 2,
#                    labels = c("A)", "B)", "C)", "D)"))
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/PredictionMaps_binary.pdf", width = 17, height = 17, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/PredictionMap_binary.png", width = 17, height = 17, dpi= 300, units="cm")



# ==============================================================================
# STEP 11: Calculate Ensemble Model in the entire area
# ==============================================================================
pred_rasters_list <- list(
  glm = glm2_pred_raster,
  gam = gam2_pred_raster,
  rf  = rf2_pred_raster,
  brt = brt2_pred_raster
)

# Ensemble model is an average of the 4 model approaches weighted by its TSS =====
# Convert to testing points into SpatVector
# test$kelp<- as.factor(test$kelp)
# test_points <- vect(test, geom = c("x", "y"), crs = "EPSG:3005")
# plot(test_points)

train$kelp<- as.factor(train$kelp)
train_points <- vect(train, geom = c("x", "y"), crs = "EPSG:3005")
plot(train_points)


# Calculate TSS + weights
tss_weights <- results_selected %>%
  select(Model, TSS) %>% 
      select(Model, TSS) %>%
      mutate(
        TSS = pmax(TSS, 0),              # Avoid negatives
        Weight = TSS / sum(TSS)          # Normalized weights
      )
  
  
print(tss_weights)
# Model   TSS Weight
# 1 glm   0.703  0.208
# 2 gam   0.859  0.254
# 3 rf    0.928  0.274
# 4 brt   0.892  0.264

# Model   TSS Weight  --> Last version Sep 1st
# 1 glm   0.610  0.198
# 2 gam   0.670  0.218
# 3 rf    1      0.325
# 4 brt   0.797  0.259

# tss_weights_1<- tss_weights
# tss_weights_1$Weight<- 1

tss_weights_1 <- results_selected %>%
  mutate(TSS=1)%>%
  select(Model, TSS) %>% 
  mutate(
    TSS = pmax(TSS, 0),              # Avoid negatives
    Weight = TSS / sum(TSS)          # Normalized weights
  )


# Build ensemble raster
ens_out <- ensemble_raster(pred_rasters_list, tss_weights)
plot(ens_out$raster)
# saveRDS(ensemble_model_rast, "ensemble_model_wAge_&_wSD_rast.rds")  

ens_out_normalAverage <- ensemble_raster(pred_rasters_list, tss_weights_1)
plot(ens_out_normalAverage$raster)

ensemble_SDweithed<- uncertainty_raster(pred_rasters_list, weights_df = tss_weights[,c(1,3)]) 
ensemble_SD<- uncertainty_raster(pred_rasters_list, weights_df = tss_weights_1[,c(1,3)]) 


# Save rasters
# terra::writeRaster(ens_out$raster, "ensemble_wAverage_suitability_blob_M2.tif", overwrite=TRUE)
# terra::writeRaster(ens_out_normalAverage$raster, "ensemble_Average_suitability_blob_M2.tif", overwrite=TRUE)
# writeRaster(ensemble_SD$max, "ensemble_Ave_model_max_M2.tif")
# writeRaster(ensemble_SD$min, "ensemble_Ave_model_min_M2.tif")
# writeRaster(ensemble_SD$sd, "ensemble_Ave_SD_blob_M2.tif")
# writeRaster(ensemble_SDweithed$max, "ensemble_wAve_model_max_M2.tif")
# writeRaster(ensemble_SDweithed$min, "ensemble_wAve_model_min_M2.tif")
# writeRaster(ensemble_SDweithed$sd, "ensemble_wAve_SDweithed_blob_M2.tif")


ens_average<- terra::rast("ensemble_Average_suitability_blob_M2.tif")
ens_wAverage<- terra::rast("ensemble_wAverage_suitability_blob_M2.tif") 
  
# Save weights + TSS table
write.csv(ensemble_SD$weights, "blob_model_tss_weights_M2.csv")



### Mark ensemble model by Substrate 
substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/SOG_substrate_20m.tif")
substrate_west<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/WCVI_substrate_20m.tif")
substrate_north<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/QCS_substrate_20m.tif")

substrate<- merge(substrate, substrate_north, substrate_west)
plot(substrate)
# The predicted raster files are classified as follows: 
# 1) Rock, 
# 2) Mixed, 
# 3) Sand, 
# 4) Mud


# Mask model predictions by substrate 
substrate<- crop(substrate, ens_average)

# Mask model predictions by substrate 
substrate_aligned <- terra::rast(ens_average)
substrate_aligned <- terra::resample(substrate, substrate_aligned, method = "near")

# substrate_aligned[substrate_aligned == 1]<- 1
substrate_aligned[substrate_aligned == 2]<- 1
substrate_aligned[substrate_aligned == 3]<- 2
substrate_aligned[substrate_aligned == 4]<- 2

# Save substrate, 1= hard substrate; 2= soft substrate
# writeRaster(substrate_aligned, "substrate_SOG_aligned.tif", overwrite=T)


substrate_aligned[substrate_aligned == 2]<- NA
ens_average_masked <- terra::mask(ens_average, substrate_aligned)
ens_wAverage_masked <- terra::mask(ens_wAverage, substrate_aligned)
# writeRaster(ens_average_masked, "ens_average_masked_M2.tif", overwrite=T)
# writeRaster(ens_wAverage_masked, "ens_wAverage_masked_M2.tif", overwrite=T)


ens_average_masked<- rast("ens_average_masked_2.tif")
ens_wAverage_masked<- rast("ens_wAverage_masked_2.tif")


# Calculate threshold for ENSEMBLE MODEL
# Extract ensemble predictions
train$ensemble_ave_pred <- terra::extract(ens_average_masked, train_points)[,2]
train$ensemble_ave_pred
train$ensemble_wAve_pred <- terra::extract(ens_wAverage_masked, train_points)[,2]


# Compute thrshold
# Compute ROC for ensemble
# train_masked<- train[,c(3,5,6,22,23)]
train_masked<- train%>%select(c(kelp, x, y, ensemble_ave_pred, ensemble_wAve_pred))
train_masked<- na.exclude(train_masked)

roc_ave_ens <- pROC::roc(response = train_masked$kelp,
                     predictor = train_masked$ensemble_ave_pred)
plot(roc_ave_ens, col = "green", main = "ROC curve: Ensemble model")

roc_wAve_ens <- pROC::roc(response = train_masked$kelp,
                     predictor = train_masked$ensemble_wAve_pred)
plot(roc_wAve_ens, col = "blue", main = "", add=T)



tab_ens <- get_thresh_table(
  roc_obj = roc_wAve_ens,
  pred = train[-which(is.na(train$ensemble_wAve_pred)),"ensemble_wAve_pred"],
  truth = train[-which(is.na(train$ensemble_wAve_pred)),"kelp"],
  model_name = "Ensemble_wAve"
)

tab_ens2 <- get_thresh_table(
  roc_obj = roc_ave_ens,
  pred = train[-which(is.na(train$ensemble_ave_pred)),"ensemble_ave_pred"],
  truth = train[-which(is.na(train$ensemble_ave_pred)),"kelp"],
  model_name = "Ensemble_Ave"
)

results_table <- dplyr::bind_rows(tab_ens, tab_ens2)
results_table

# With testing data: # Model         Criterion     Threshold Sensitivity Specificity     TSS   AUC
# <chr>         <chr>             <dbl>       <dbl>       <dbl>   <dbl> <dbl>
# 1 Ensemble_wAve Max TSS         0.588         0.846     0.798   0.645   0.910
# 2 Ensemble_wAve No omission  -Inf             1         0       0       0.910
# 3 Ensemble_wAve 10% omission    0.00893       1         0.00806 0.00806 0.910
# 4 Ensemble_Ave  Max TSS         0.535         0.869     0.774   0.643   0.908
# 5 Ensemble_Ave  No omission  -Inf             1         0       0       0.908
# 6 Ensemble_Ave  10% omission    0.00909       1         0.00806 0.00806 0.908

# With training data: 
# Model         Criterion     Threshold Sensitivity Specificity     TSS   AUC
# 1 Ensemble_wAve Max TSS         0.638         0.893     0.964   0.857   0.982
# 2 Ensemble_wAve No omission  -Inf             1         0       0       0.982
# 3 Ensemble_wAve 10% omission    0.00747       1         0.00331 0.00331 0.982
# 4 Ensemble_Ave  Max TSS         0.694         0.841     0.983   0.824   0.974
# 5 Ensemble_Ave  No omission  -Inf             1         0       0       0.974
# 6 Ensemble_Ave  10% omission    0.00707       1         0.00331 0.00331 0.974

# Model 2 : 
# Model         Criterion     Threshold Sensitivity Specificity     TSS   AUC
# 1 Ensemble_wAve Max TSS         0.603         0.903     0.935   0.838   0.981
# 2 Ensemble_wAve No omission  -Inf             1         0       0       0.981
# 3 Ensemble_wAve 10% omission    0.00945       1         0.00327 0.00327 0.981
# 4 Ensemble_Ave  Max TSS         0.542         0.924     0.889   0.813   0.973
# 5 Ensemble_Ave  No omission  -Inf             1         0       0       0.973
# 6 Ensemble_Ave  10% omission    0.00878       1         0.00327 0.00327 0.973


models_performance <- bind_rows(
  get_metrics_optimized2(model= glm_mod_s, test_data=train[-which(is.na(train$ensemble_ave_pred)),], scale_params= scaling_params_2,
                         model_name = "glm", threshold_type = "youden"),
  get_metrics_optimized2(gam_mod_s, train[-which(is.na(train$ensemble_ave_pred)),], scale_params= scaling_params_2, "gam",
                         threshold_type = "youden"),#10pct_omission
  get_metrics_optimized2(rf_mod_s, train[-which(is.na(train$ensemble_ave_pred)),], "rf", threshold_type = "youden"),
  get_metrics_optimized2(brt_mod_s, train[-which(is.na(train$ensemble_ave_pred)),], "brt", threshold_type = "youden")
)


ens_performance<- results_table%>%
  filter(Criterion=="Max TSS")%>%
  mutate(Criterion = "youden")

colnames(ens_performance)[2]<- "ThresholdType"

models_performance<- rbind(models_performance, ens_performance[,c(1:3,5:7,4)])
# write.csv(models_performance, "models_performance_calibration_Table_m2.csv")

# Model         ThresholdType Threshold   AUC Sensitivity Specificity   TSS
# 1 glm           youden            0.609 0.870       0.829       0.781 0.610
# 2 gam           youden            0.648 0.898       0.812       0.834 0.646
# 3 rf            youden            0.523 1           1           1     1    
# 4 brt           youden            0.621 0.963       0.878       0.907 0.785
# 5 Ensemble_wAve youden            0.638 0.982       0.893       0.964 0.857
# 6 Ensemble_Ave  youden            0.694 0.974       0.841       0.983 0.824


# Model 2
# Model         ThresholdType Threshold   AUC Sensitivity Specificity   TSS
# 1 glm           youden            0.634 0.873       0.825       0.807 0.633
# 2 gam           youden            0.546 0.906       0.868       0.791 0.659
# 3 rf            youden            0.527 1           1           1     1    
# 4 brt           youden            0.530 0.957       0.912       0.886 0.797
# 5 Ensemble_wAve youden            0.603 0.981       0.903       0.935 0.838
# 6 Ensemble_Ave  youden            0.542 0.973       0.924       0.889 0.813


### Evaluate and select threshold for binary results =====
# Extract threshold, sensitivity, specificity
thr_vec  <- roc_ave_ens$thresholds
sens_vec <- roc_ave_ens$sensitivities
spec_vec <- roc_ave_ens$specificities

# Keep only finite thresholds
finite_idx <- which(is.finite(thr_vec))
thr_vec_f  <- thr_vec[finite_idx]
sens_vec_f <- sens_vec[finite_idx]
spec_vec_f <- spec_vec[finite_idx]

# Compute TSS
tss_vec <- sens_vec_f + spec_vec_f - 1

# Create data frame for plotting
tss_df <- data.frame(
  Threshold = thr_vec_f,
  TSS = tss_vec
)

ggplot(tss_df, aes(x = Threshold, y = TSS)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(aes(x = Threshold[which.max(TSS)], y = max(TSS)),
             color = "red", size = 3) +
  labs(
    x = "Threshold",
    y = "TSS",
    title = "TSS vs Probability Threshold"
  ) +
  theme_minimal()


# Use predicted values for the training set
# train[-which(is.na(train$ensemble_ave_pred))
roc_train <- pROC::roc(response = train[-which(is.na(train$ensemble_ave_pred)),"kelp"],
                       predictor = train[-which(is.na(train$ensemble_ave_pred)),"ensemble_ave_pred"])

# Find Max TSS threshold
thr_maxTSS <- roc_train$thresholds[which.max(roc_train$sensitivities + roc_train$specificities - 1)]# Using “drop point” from cumulative coverage or density

# Predicted suitability for presences in training data
train_ave_ens<- train[-which(is.na(train$ensemble_ave_pred)),]

# Predicted suitability for presences and absences
pres_pred <- train_ave_ens$ensemble_ave_pred[train_ave_ens$kelp == 1]
abs_pred  <- train_ave_ens$ensemble_ave_pred[train_ave_ens$kelp == 0]

# Density estimates
dens_pres <- density(pres_pred, na.rm = TRUE)
dens_abs  <- density(abs_pred,  na.rm = TRUE)

# Plot both
library(ggplot2)

# Density estimates
dens_pres <- density(pres_pred, na.rm = TRUE)
dens_abs  <- density(abs_pred,  na.rm = TRUE)

# Build data frame for ggplot
df_plot <- data.frame(
  x = c(dens_pres$x, dens_abs$x),
  y = c(dens_pres$y, dens_abs$y),
  group = rep(c("Presences", "Absences"),
              times = c(length(dens_pres$x), length(dens_abs$x)))
)

# Thresholds
threshold_drop <- quantile(pres_pred, 0.001, na.rm = TRUE)  # 0.1% omission
thr_maxTSS <- thr_maxTSS  # make sure you have this value computed earlier

# Plot Density estimates
dens_pres <- density(pres_pred, na.rm = TRUE)
dens_abs  <- density(abs_pred,  na.rm = TRUE)

# Build data frame for ggplot
df_plot <- data.frame(
  x = c(dens_pres$x, dens_abs$x),
  y = c(dens_pres$y, dens_abs$y),
  group = rep(c("Absences", "Presences"),   # abs first, then pres
              times = c(length(dens_abs$x), length(dens_pres$x)))
)

# Thresholds
threshold_drop <- quantile(pres_pred, 0.001, na.rm = TRUE)  # 0.1% omission --> 0.2054824  --> M2= 0.200
thr_maxTSS <- thr_maxTSS  # 0.6940952  --> M2= 0.5420

# Define inputs
p <- train_ave_ens$ensemble_ave_pred   # predicted suitability from your model
y <- train_ave_ens$kelp                # true presence/absence (0/1)

# Extract thresholds
thr_vec <- roc_train$thresholds
thr_vec <- thr_vec[is.finite(thr_vec)]  # keep only finite

# Function to compute F1 at a given threshold
f1_at_thresh <- function(th){
  pred <- as.integer(p >= th)
  tp <- sum(pred == 1 & y == 1)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  
  precision <- if((tp + fp) > 0) tp / (tp + fp) else NA
  recall    <- if((tp + fn) > 0) tp / (tp + fn) else NA
  f1        <- if(is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA
  return(f1)
}

# Compute F1 across thresholds
f1_vals <- sapply(thr_vec, f1_at_thresh)

# Pick threshold with max F1
i_bestF1   <- which.max(f1_vals)
thr_maxF1  <- thr_vec[i_bestF1]
f1_max_val <- f1_vals[i_bestF1]

cat("Threshold maximizing F1:", thr_maxF1, "\n")
# Threshold maximizing F1: 0.4689823   
cat("Max F1 value:", f1_max_val, "\n")
# Max F1 value: 0.9283333 



# Plot
ggplot(df_plot, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = threshold_drop, color = "darkgreen", linetype = "dashed", size = 1) +
  geom_vline(xintercept = thr_maxTSS, color = "red", linetype = "dashed", size = 1) +
  # geom_vline(xintercept = thr_maxF1, color = "darkorange", linetype = "dashed", size = 1) +
  labs(
    x = "Predicted Kelp Habitat Suitability",
    y = "Density",
    color = ""
  ) +
  scale_color_manual(values = c("blue", "red")) +  # abs = blue, pres = red
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.15, 0.84),   # inside, top-left (x%, y%)
    legend.background = element_rect(fill = alpha("white", 0.1), color = NULL),
    legend.direction = "vertical",   # stacked vertically
    legend.key = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  # annotate("text", x = threshold_drop, y = max(df_plot$y)*0.8,
  #          label = "0.1% Omission", color = "darkgreen", angle = 90, vjust = -0.5) +
  annotate("text", x = thr_maxTSS, y = max(df_plot$y)*0.8,
           label = "Max TSS", color = "grey40", angle = 90, vjust = -0.5)+
  annotate("text", x = thr_maxTSS, y = max(df_plot$y)*0.8,
           label = "Max TSS", color = "grey40", angle = 90, vjust = -0.5)

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Threshold_selection_plot.png", width = 10, height = 8, dpi= 300, units="cm")




























### Binary predictions at optimal threshold  ======================
pred_bin_ave_ens <- ifel(ens_average_masked_all >= results_table[which(results_table$Model== "Ensemble_Ave" & results_table$Criterion== "Max TSS"), "Threshold"][[1]], 1, 0)
pred_bin_wAve_ens <- ifel(ens_wAverage_masked_all >= results_table[which(results_table$Model== "Ensemble_wAve" & results_table$Criterion== "Max TSS"), "Threshold"][[1]], 1, 0)
# writeRaster(pred_bin_ave_ens, "tifs/ens_ave_bin_blob.tif")
# writeRaster(pred_bin_wAve_ens, "tifs/ens_wAve_bin_blob.tif")


# Manual metrics
pred_vals <- terra::extract(pred_class, test_points)  # extract raster preds at test points
# Combine truth with predictions (drop ID col)
eval_df <- data.frame(
  truth = test_points$kelp,
  pred  = pred_vals[,2]   # second column = raster values
)

# Remove NAs (some points may fall outside raster extent)
eval_df <- na.omit(eval_df)
summary(eval_df)

# Confusion matrix components
tp <- sum(eval_df$pred == 1 & eval_df$truth == 1)
tn <- sum(eval_df$pred == 0 & eval_df$truth == 0)
fp <- sum(eval_df$pred == 1 & eval_df$truth == 0)
fn <- sum(eval_df$pred == 0 & eval_df$truth == 1)

sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
auc_val     <- as.numeric(pROC::auc(roc_ens))
tss         <- sensitivity + specificity - 1


## Add values of ensemble model into table of results ========
results_selected[5,1:6]<- list("ens", opt_thresh, auc_val, sensitivity, specificity,tss)
write.csv(results_selected, "results_modelsBlob_Eval_Table.csv")


thresholds_table[5,1]<- "ens"
thresholds_table[5,2]<- opt_thresh





# ==============================================================================
# STEP 12: RE-scale Ensemble Model (unsuitable vs suitability)
# ==============================================================================
ens_model_rescaled<- threshold_rescale(ens_raster= ensemble_model_rast$mean, opt_thresh = results_selected[which(results_selected$Model=="ens"),"Threshold"]$ Threshold)
plot(ens_model_rescaled$rescaled)
writeRaster(ens_model_rescaled$rescaled, "ens_model_rescaled.tif")
# writeRaster(ens_model_rescaled$binary, "ens_model_bin.tif")



# ==============================================================================
# STEP 13: Selection of the threshold 
# ==============================================================================
# thresholds to evaluate
thresholds <- seq(0, 1, by = 0.05)

# presence points (SpatVector or sf with presence-only points)
# assume you have: presence_points

# Use predicted values for the training set
roc_train <- pROC::roc(response = train_masked$kelp,
                       predictor = train_sel$ensemble_ave_pred)

# Find Max TSS threshold
thr_maxTSS <- roc_train$thresholds[which.max(roc_train$sensitivities + roc_train$specificities - 1)]
# Using “drop point” from cumulative coverage or density
# Predicted suitability for presences in training data
pres_pred <- train$ensemble_ave_pred[train$kelp == 1]

# Density or cumulative coverage
dens <- density(pres_pred)
plot(dens, main="Density of predicted presences (training)", xlab="Predicted suitability")

# e.g., 5th percentile of presence predictions
threshold_drop <- quantile(pres_pred, 0.05)
abline(v = threshold_drop, col="red", lty=2)


# Extract predictions for each model at presence locations
pred_rasters_list_2 <- list(
  glm = glm_pred_raster,
  gam = gam_pred_raster,
  rf  = rf_pred_raster,
  brt = brt_pred_raster,
  ens = ens_average
)

### Calculate omission error per threshold per model  
# test_points <- vect(test, geom = c("x", "y"), crs = "EPSG:3005")
train_points <- vect(train, geom = c("x", "y"), crs = "EPSG:3005")

model_preds <- lapply(pred_rasters_list_2, function(r) {
  terra::extract(r, train_points)[,2]  # column 2 = values
})

results_omission <- lapply(names(model_preds), function(m) {
  preds <- model_preds[[m]]
  data.frame(
    Model = m,
    Threshold = thresholds,
    Omission = sapply(thresholds, function(t) mean(preds < t, na.rm = TRUE))
  )
}) %>% bind_rows()

# Plot
ggplot(results_omission, aes(x = Threshold, y = Omission, color = Model)) +
  geom_line(size = 1) +
  geom_vline(xintercept = opt_thresh$Threshold, linetype = "dashed", color = "red", size = 1) +
  theme_bw() +
  labs(x = "Threshold", y = "Omission Error",
       title = "Omission Error vs. Threshold for Models")

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/omissionError_thresholds_models.pdf", width = 12, height = 9, dpi= 300, units="cm")

