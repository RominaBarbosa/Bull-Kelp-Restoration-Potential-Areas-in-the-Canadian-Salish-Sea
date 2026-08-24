###==================================================================
### Species Distribution models        SDMs          ################
###     Ensemble model  - Blob conditions            ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Date: 11-Aug-2025                                ################
### Last edition: 23-Sep-2025                        ################
###==================================================================
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf)
library(raster)

source("/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")
source("/R_scripts/functions/plot_response_curves.R")

model_results_path<- "/SDM/SDM_results/Sep2025_M7_weightedPres"
terrain_path<- ("/SDM/environmental_layers/Topographic_Variables")


# This chunk was run and its outputs are stored in a RDS file
# # Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_blob.tif")
names(raster_stack_20m)
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)
terrain_vars<-  rast(tif_files[c(12)]) #,5,12,13,15
names(terrain_vars)
summary(terrain_vars)

# Merge rasters of all selected variables including terrain and NEMO
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)

# Create Mask from bathymetry
bathy20m<- rast("/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask10_30<- bathy20m
bathy20m_mask10_30[!(bathy20m_mask10_30 >= -10 & bathy20m_mask10_30 <= 30)] <- NA # negative values are in land
plot(bathy20m_mask10_30)

# Mask all layers
setwd(model_results_path)
glm_mod_s<- readRDS("glm_mod_s.rds")
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)

raster_stack_predict<- raster_stack_20m_all[[vars_selected]] # subset selected variables
names(raster_stack_predict)
plot(raster_stack_predict[[1]])
raster_stack_predict<- mask(raster_stack_predict, bathy20m_mask10_30)
plot(raster_stack_predict[[1]])

setwd(model_results_path)
saveRDS(raster_stack_predict, "raster_stack_predict.rds")
writeRaster(raster_stack_predict, "raster_stack_predict.tif")



# ================================================================
# STEP 8: Predict in the entire area
# ================================================================
setwd(model_results_path)
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
# vars_selected<- c("slope_5x5", "turbidity_summer_mean", "nitrate_winter_mean", "temperature_summer_mean",  #"nitrate_summer_minimum",
#                   "ammonium_spring_mean", "currentSpeed_summer_mean", "PAR_summer_mean")


# Load training and testing dataset 
train<- read.csv(paste(model_results_path,"train_selected_table_FINALMODELS_M7.csv", sep="/"))
train<- train[,-1]
names(train)

train_sel <- train %>% dplyr::select(all_of(c("kelp", vars_selected)))
train_sel$kelp<- as.factor(train_sel$kelp)

train_scaled_b <- scale(train_sel[, vars_selected])
train_means <- attr(train_scaled_b, "scaled:center")
train_sds   <- attr(train_scaled_b, "scaled:scale")

raster_stack_predict_scaled<- raster_stack_predict[[names(raster_stack_predict)%in% vars_selected]]  # original raster stack

for (v in vars_selected) {
  raster_stack_predict_scaled[[v]] <- (raster_stack_predict[[v]] - train_means[v]) / train_sds[v]
}
names(raster_stack_predict_scaled)#<- names(raster_stack_predict)


raster_stack_predict_scaled<- mask(raster_stack_predict_scaled, bathy20m_mask10_30)
# saveRDS(raster_stack_predict_scaled, "raster_stack_predict_scaled_Blob.rds")
# writeRaster(raster_stack_predict_scaled, "raster_stack_predict_scaled_Blob.tif")



# Loads raster inputs to predict the model
raster_stack_predict<- readRDS("raster_stack_predict.rds")
raster_stack_predict_scaled<- readRDS("raster_stack_predict_scaled_Blob.rds")
# names(raster_stack_predict) 

# ==============================================================================
# Predict probabilities  
# ==============================================================================
terra::predict(raster_stack_predict_scaled, #glm took a lot of memory so I saved directly to disk
                                  glm_mod_s,
                                  type = "response",
                                  na.rm = TRUE,
                                  filename = "tifs/glm_pred_raster.tif",  # output written directly to disk
                           overwrite = TRUE)

glm_pred_raster <- rast("tifs/glm_pred_raster.tif")
names(glm_pred_raster)<- "glm"


terra::predict(raster_stack_predict_scaled, #glm took a lot of memory so I saved directly to disk
               gam_mod_s,
               type = "response",
               na.rm = TRUE,
               filename = "tifs/gam_pred_raster.tif",  # output written directly to disk
               overwrite = TRUE)

gam_pred_raster <- rast("tifs/gam_pred_raster.tif")
names(gam_pred_raster)<- "gam"
rf_pred_raster  <- predict_raster_prob(rf_mod_s,  raster_stack_predict, "rf")
names(rf_pred_raster)<- "rf"
brt_pred_raster <- predict_raster_prob(brt_mod_s, raster_stack_predict, "brt")
names(brt_pred_raster)<- "brt"

### SAVE RESULT RASTERS =======================================================
writeRaster(glm_pred_raster, "tifs/glm_pred_raster.tif", overwrite = TRUE)
writeRaster(gam_pred_raster, "tifs/gam_pred_raster.tif", overwrite = TRUE)
writeRaster(rf_pred_raster,  "tifs/rf_pred_raster.tif",  overwrite = TRUE)
writeRaster(brt_pred_raster, "tifs/brt_pred_raster.tif", overwrite = TRUE)

# setwd("/SDM/SDM_results/Sep2025")
# glm_pred_raster<- rast("glm_pred_raster.tif")
# gam_pred_raster<- rast("gam_pred_raster.tif")
# rf_pred_raster<- rast("rf_pred_raster.tif")
# brt_pred_raster<- rast("brt_pred_raster.tif")



# ==============================================================================
# STEP 9: Calculate Ensemble Model in the entire area
# ==============================================================================
pred_rasters_list <- list(
  glm = glm_pred_raster,
  gam = gam_pred_raster,
  rf  = rf_pred_raster,
  brt = brt_pred_raster
)

# Ensemble model is an average of the 4 model approaches weighted by its TSS =====
# Convert to testing points into SpatVector
# test$kelp<- as.factor(test$kelp)
# test_points <- vect(test, geom = c("x", "y"), crs = "EPSG:3005")
# plot(test_points)

train$kelp<- as.factor(train$kelp)
train_points <- vect(train, geom = c("x", "y"), crs = "EPSG:3005")
plot(train_points)


### Calculate ensemble model ================
# Ensemble model is an average of the 4 model approaches weighted by its original TSS 
# tss_weights<- read.csv(paste(model_results_path, "blob_model_tss_weights.csv", sep="/"))
tss_weights_1 <- data.frame(
  Model  = c("glm", "gam", "rf", "brt"),
  # TSS    = c(0.673, 0.701 , 1, 0.864),
  Weight = c(0.25, 0.25, 0.25, 0.25)
)


pred_rasters_list <- list(
  glm = glm_pred_raster,
  gam = gam_pred_raster,
  rf  = rf_pred_raster,
  brt = brt_pred_raster
)


# Build ensemble raster
ens_ave <- ensemble_raster(pred_rasters_list, tss_weights_1) # 2 min
plot(ens_ave$raster)

# Uncertainty raster
# start<- Sys.time()
ensemble_model_Blob<- uncertainty_raster(pred_rasters_list, weights_df = tss_weights_1)
# end<- Sys.time()
# time_total<- end-start
# time_total # 30 min  ; M4 55.37769 mins; M6 46.9 mins
saveRDS(ensemble_model_Blob, "ensemble_model_Blob.rds")

# Save rasters
writeRaster(ens_ave$raster, "tifs/ensemble_ave_blob_M7.tif", overwrite=T)
writeRaster(ensemble_model_Blob$max, "ensemble_model_max_postBlob.tif")
writeRaster(ensemble_model_Blob$min, "ensemble_model_min_postBlob.tif")
writeRaster(ensemble_model_Blob$sd, "ensemble_model_wSD_postBlob.tif")



### PLOT RESULT MAPS ===========================================================
# Load high-resolution coastline (Natural Earth)
library(rnaturalearth)
library(rnaturalearthdata)
coast <- ne_coastline(scale = "large", returnclass = "sf")
coast_vect <- vect(coast)
# Reproject to EPSG:3005
# coast_vect <- project(coast_vect, "EPSG:3005")

coast_sf <- st_read("/Canada_Provincesgpr_000b11a_e/gpr_000b11a_e.shp")
coast_sf <- st_transform(coast_sf, "EPSG:3005")
coast_vect <- vect(coast_sf)
r_extent <- ext(ens_ave$raster)  # your raster extent
coast_crop <- crop(coast_vect, r_extent)

my_palette <- colorRampPalette(c("blue", "green", "yellow", "red"))

png("Ave_ensemble_blob_M7.png", width = 2000, height = 1500, res = 300)
terra::plot(ens_ave$raster, main = "Ensemble Model (average)", col = my_palette(100))
plot(coast_crop, col = "lightgrey", border = "darkgrey", lwd = 0.02, add = TRUE)
dev.off()

# png("PC2_ensemble_blob_M7.png", width = 2000, height = 1500, res = 300)
# terra::plot(ensemble_raster_pc2, main = "Model Disagreement (PC2)")
# plot(coast_crop, col = "lightgrey", border = "darkgrey", lwd = 0.02, add = TRUE)
# dev.off()
