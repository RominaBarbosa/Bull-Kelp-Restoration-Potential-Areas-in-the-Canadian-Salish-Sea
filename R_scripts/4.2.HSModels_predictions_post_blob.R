###==================================================================
### Species Distribution models        SDMs          ################
###                                                  ################
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

source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")

path_postBlob<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M2"
setwd(path_postBlob)




# Step 1 was run and now you only need to upload the inputs
# =============================================================================
# STEP 1: Load input variables of Post Blob Period of the entire area
# =============================================================================
# Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblobM2.tif")
names(raster_stack_20m)


### Load terrain variables =====================================================
terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)

# Load and stack all tif files
# terrain_vars<- terra::rast(paste(variables_selection_path,"terrain_rasters_selected.tif", sep="/"))
terrain_vars<-  rast(tif_files[c(13)])
names(terrain_vars)<- "slope_5x5"
 # "slope_7x7"

# crs(raster_stack_20m)<- crs(raster_stack_20m)


### Merge rasters of all selected variables including terrain and NEMO =========
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)
names(raster_stack_20m_all)


# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask15_40<- bathy20m
bathy20m_mask15_40[!(bathy20m_mask15_40 >= -15 & bathy20m_mask15_40 <= 40)] <- NA # negative values are in land
plot(bathy20m_mask15_40)

### Mask all layers by bathymetry  ==============================================
plot(raster_stack_20m_all[[1]])
raster_stack_20m_all<- mask(raster_stack_20m_all, bathy20m_mask15_40)
plot(raster_stack_20m_all[[1]])

setwd(path_postBlob)
saveRDS(raster_stack_20m_all, "raster_stack_predict_postBlob.rds")


### Load models  ==============================================
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M4"

glm_mod_s<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))
gam_mod_s<- readRDS(paste(model_results_path,"gam_mod_s.rds", sep="/"))
rf_mod_s<-  readRDS(paste(model_results_path,"rf_mod_s.rds", sep="/"))
brt_mod_s<- readRDS(paste(model_results_path,"brt_mod_s.rds", sep="/"))

models <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")

vars_selected<- colnames(glm_mod_s$data)[2:length(colnames(glm_mod_s$data))]

#  c( "slope_7x7", "ammonium_spring_SD", "ammonium_winter_mean","PAR_summer_mean", "temperature_summer_mean",
# "turbidity_summer_mean","salinity_summer_SD", "nitrate_summer_mean", "ammonium_summer_mean", "ammonium_spring_mean")   



### Scale variables before predictions =========================================
# We need the original training data to scale with same parameters
# train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob.csv")
train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_M4.csv")
train_test_dataset<- (train_test_dataset[,-1])
names(train_test_dataset)

train<- train_test_dataset%>%  filter(set== "train")
train_sel <- train %>% select(all_of(c("kelp", vars_selected)))
train_sel$kelp<- as.factor(train_sel$kelp)

train_scaled_b <- scale(train_sel[, vars_selected])
train_means <- attr(train_scaled_b, "scaled:center")
train_sds   <- attr(train_scaled_b, "scaled:scale")


raster_stack_predict_scaled<- raster_stack_20m_all[[names(raster_stack_20m_all)%in% vars_selected]]  # original raster stack

# for (v in vars_selected) {
#   raster_stack_predict_scaled[[v]] <- (raster_stack_20m_all[[v]] - train_means[v]) / train_sds[v]
# }
# names(raster_stack_predict_scaled)#<- names(raster_stack_predict)
# 
# saveRDS(raster_stack_predict_scaled, "raster_stack_predict_scaled_postBlob.rds")



# ================================================================
# STEP 2: Predict in the entire area
# ================================================================
# Predict probabilities
setwd(path_postBlob)
raster_stack_predict<- readRDS("raster_stack_predict_postBlob.rds")
raster_stack_predict_scaled<- readRDS("raster_stack_predict_scaled_postBlob.rds")

# Predict GLM and GAM and save results (it doesn't run if not saving directly on disk)
terra::predict(raster_stack_predict_scaled, #glm took a lot of memory so I saved directly to disk
                                  glm_mod_s,
                                  type = "response",
                                  na.rm = TRUE,
                                  filename = "glm_postBlob_pred_raster.tif",  # output written directly to disk
                           overwrite = TRUE)


# Align rasters in the same order than in the GAM model
raster_stack_aligned <- align_stack_to_model(raster_stack_predict_scaled, gam_mod_s)
terra::predict(raster_stack_aligned, #glm took a lot of memory so I saved directly to disk
               gam_mod_s,
               type = "response",
               na.rm = TRUE,
               filename = "gam_postBlob_pred_raster.tif",  # output written directly to disk
               overwrite = TRUE)



rf_pred_raster  <- predict_raster_prob(rf_mod_s,  raster_stack_predict, "rf")
brt_pred_raster <- predict_raster_prob(brt_mod_s, raster_stack_predict, "brt")

# SAVE RESULT RASTERS from RF and BRT
# writeRaster(rf_pred_raster,  "rf_postBlob_pred_raster.tif",  overwrite = TRUE)
# writeRaster(brt_pred_raster, "brt_postBlob_pred_raster.tif", overwrite = TRUE)

glm_pred_raster <- rast("glm_postBlob_pred_raster.tif")
gam_pred_raster <- rast("gam_postBlob_pred_raster.tif")
rf_pred_raster<- rast("rf_postBlob_pred_raster.tif")
brt_pred_raster<- rast("brt_postBlob_pred_raster.tif")


### Calculate ensemble model ================
# Ensemble model is an average of the 4 model approaches weighted by its original TSS 
# tss_weights<- read.csv(paste(model_results_path, "blob_model_tss_weights.csv", sep="/"))

tss_weights_1 <- data.frame(
  Model  = c("glm", "gam", "rf", "brt"),
  TSS    = c(0.673, 0.701 , 1, 0.864),
  Weight = c(0.25, 0.25, 0.25, 0.25)
)


pred_rasters_list <- list(
  glm = glm_pred_raster,
  gam = gam_pred_raster,
  rf  = rf_pred_raster,
  brt = brt_pred_raster
)


# Build ensemble raster
start<- Sys.time()
ens_out_postblob <- ensemble_raster(pred_rasters_list, tss_weights_1)
end<- Sys.time()
time_total<- end-start # 2 min

plot(ens_out_postblob$raster)
# writeRaster(ens_out_postblob$raster, "ens_output_postblob.tif")

start<- Sys.time()
ensemble_model_postBlob<- uncertainty_raster(pred_rasters_list, weights_df = tss_weights_1[,c(1,3)]) 
end<- Sys.time()
time_total<- end-start
time_total # 30 min ; M4 55.37769 mins

# saveRDS(ensemble_model_postBlob, "ensemble_model_postBlob.rds")

plot(ensemble_model_rast$max)


# Save rasters
writeRaster(ensemble_model_postBlob$mean, "ensemble_ave_postBlob.tif")
writeRaster(ensemble_model_postBlob$max, "ensemble_model_max_postBlob.tif")
writeRaster(ensemble_model_postBlob$min, "ensemble_model_min_postBlob.tif")
writeRaster(ensemble_model_postBlob$sd, "ensemble_model_wSD_postBlob.tif")


# ==============================================================================
# STEP 3: Binary results & Evaluation with presence/absence data from Post blob
# ==============================================================================
# Get metrics of model performance, AUC, Sensitivity and specificity

# thresholds_table<- read.csv(paste(model_results_path, "results_modelsBlob_Eval_Table.csv", sep="/"))
# ens_model_trheshold<- ens_model_trheshold[which(ens_model_trheshold$Model == "ens"), "Threshold"]

threshold_drop <-  0.1957
thr_maxTSS <- 0.6977754


thresholds_table<- thresholds_table%>%
  select(Model, Threshold)

# # Convert continuous predictions into binary based on max sensitivity & specificity threshold:
# glm_bin <- binarize_raster(glm_pred_raster, "glm", thresholds_table)
# gam_bin <- binarize_raster(gam_pred_raster, "gam", thresholds_table)
# rf_bin  <- binarize_raster(rf_pred_raster, "rf", thresholds_table)
# brt_bin <- binarize_raster(brt_pred_raster, "brt", thresholds_table)
# ens_bin <- binarize_raster(ensemble_model_postBlob$mean, "ens", thresholds_table)
# 
# # Save binary rasters if needed
# writeRaster(glm_bin, "glm_bin_postBlob.tif", overwrite = TRUE)
# writeRaster(gam_bin, "gam_bin_postBlob.tif", overwrite = TRUE)
# writeRaster(rf_bin,  "rf_bin_postBlob.tif",  overwrite = TRUE)
# writeRaster(brt_bin, "brt_bin_postBlob.tif", overwrite = TRUE)
# writeRaster(ens_bin, "ens_bin_postBlob.tif", overwrite = TRUE)
# 
# glm_bin<- rast("glm_bin_postBlob.tif")
# gam_bin<- rast("gam_bin_postBlob.tif")
# rf_bin<- rast("rf_bin_postBlob.tif")
# brt_bin<- rast("brt_bin_postBlob.tif")
# ens_bin<- rast("ens_bin_postBlob.tif")



results_selected <- bind_rows(
  get_metrics_optimized(model= glm_mod_s, test_data=test_sel_scaled, model_name = "glm"), # kelp must be a factor
  get_metrics_optimized(gam_mod_s, test_sel_scaled, "gam"), # kelp must be a factor
  get_metrics_optimized(rf_mod_s, test_sel, "rf"),
  get_metrics_optimized(brt_mod_s, test_sel, "brt")
)


# Compute ROC for ensemble
roc_ens <- pROC::roc(response = test_data$kelp,
                     predictor = test_data$ensemble_pred)

# Youden's J statistic to find optimal threshold
sens_vec <- roc_ens$sensitivities
spec_vec <- roc_ens$specificities
thr_vec  <- roc_ens$thresholds
J <- sens_vec + spec_vec - 1
i_best <- which.max(J)
opt_thresh <- thr_vec[i_best]


### Binary predictions at optimal threshold  ======================
pred_class <- ifel(ensemble_model_rast$mean >= opt_thresh, 1, 0)
# writeRaster(pred_class, "ens_bin_blob.tif")

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
# 
# # source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")
# 
# glm_binary_map<- plot_raster_gg(glm_bin, title= "GLM Binary Prediction")
# gam_binary_map<- plot_raster_gg(gam_bin, title= "GAM Binary Prediction")
# rf_binary_map<- plot_raster_gg(rf_bin,  title= "RF Binary Prediction")
# brt_binary_map<- plot_raster_gg(brt_bin, title= "BRT Binary Prediction")
# 
# cowplot::plot_grid(glm_binary_map, gam_binary_map, rf_binary_map, brt_binary_map, nrow = 2,
#                    labels = c("A)", "B)", "C)", "D)"))
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/PredictionMap_GLM.pdf", width = 14, height = 15, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/PredictionMap_GLM.png", width = 14, height = 15, dpi= 300, units="cm")
# 
# 
# 
# 
# # ==============================================================================
# # STEP 12: RE-scale Ensemble Model (unsuitable vs suitability)
# # ==============================================================================
# ens_model_rescaled<- threshold_rescale(ens_raster= ensemble_model_rast$mean, opt_thresh = results_selected[which(results_selected$Model=="ens"),"Threshold"]$ Threshold)
# plot(ens_model_rescaled$rescaled)
# writeRaster(ens_model_rescaled$rescaled, "ens_model_rescaled.tif")
# # writeRaster(ens_model_rescaled$binary, "ens_model_bin.tif")
# 
# 
# 
# # ==============================================================================
# # STEP 13: Evaluate the effect of thrsholding 
# # ==============================================================================
# # thresholds to evaluate
# thresholds <- seq(0, 1, by = 0.05)
# 
# # presence points (SpatVector or sf with presence-only points)
# # assume you have: presence_points
# 
# # Extract predictions for each model at presence locations
# pred_rasters_list_2 <- list(
#   glm = glm_pred_raster,
#   gam = gam_pred_raster,
#   rf  = rf_pred_raster,
#   brt = brt_pred_raster,
#   ens = ensemble_model_rast$mean
# )
# 
# ### Calculate omission error per threshold per model  
# test_points <- vect(test, geom = c("x", "y"), crs = "EPSG:3005")
# 
# model_preds <- lapply(pred_rasters_list_2, function(r) {
#   terra::extract(r, test_points)[,2]  # column 2 = values
# })
# 
# results_omission <- lapply(names(model_preds), function(m) {
#   preds <- model_preds[[m]]
#   data.frame(
#     Model = m,
#     Threshold = thresholds,
#     Omission = sapply(thresholds, function(t) mean(preds < t, na.rm = TRUE))
#   )
# }) %>% bind_rows()
# 
# # Plot
# ggplot(results_omission, aes(x = Threshold, y = Omission, color = Model)) +
#   geom_line(size = 1) +
#   geom_vline(xintercept = opt_thresh$Threshold, linetype = "dashed", color = "red", size = 1) +
#   theme_bw() +
#   labs(x = "Threshold", y = "Omission Error",
#        title = "Omission Error vs. Threshold for Models")
# 
# # ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/omissionError_thresholds_models.pdf", width = 12, height = 9, dpi= 300, units="cm")

