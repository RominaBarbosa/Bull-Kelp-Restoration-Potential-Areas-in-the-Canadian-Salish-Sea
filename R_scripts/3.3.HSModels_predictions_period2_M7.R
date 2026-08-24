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
library(randomForest)
library(gbm)


source("/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")

path_postBlob<- "/SDM/SDM_predict_postBlob/M7"
setwd(path_postBlob)

# Step 1 was run and now you only need to upload the inputs
# =============================================================================
# STEP 1: Load input variables of Post Blob Period of the entire area
# =============================================================================
# Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblob.tif")
names(raster_stack_20m)


### Load terrain variables =====================================================
terrain_path<- ("/SDM/environmental_layers/Topographic_Variables")
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)

# Load and stack all tif files
# terrain_vars<- terra::rast(paste(variables_selection_path,"terrain_rasters_selected.tif", sep="/"))
terrain_vars<-  rast(tif_files[c(13)])
names(terrain_vars)<- "slope_5x5"
 # "slope_7x7"

# crs(raster_stack_20m)<- crs(raster_stack_20m)


### Merge rasters of all selected variables including terrain and NEMO =========
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)
# names(raster_stack_20m_all)[11]<- "slope_5x5"


# Create Mask from bathymetry
bathy20m<- rast("/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask10_30<- bathy20m
bathy20m_mask10_30[!(bathy20m_mask10_30 >= -10 & bathy20m_mask10_30 <= 30)] <- NA # negative values are in land

### Mask all layers by bathymetry  ==============================================
plot(raster_stack_20m_all[[1]])
raster_stack_20m_all<- mask(raster_stack_20m_all, bathy20m_mask10_30)
plot(raster_stack_20m_all[[1]])

setwd(path_postBlob)
# saveRDS(raster_stack_20m_all, "raster_stack_predict_postBlob.rds")
writeRaster(raster_stack_20m_all, "raster_stack_predict_postBlob.tif", overwrite=TRUE)


### Load models  ==============================================
model_results_path<- "/SDM/SDM_results/Sep2025_M7_weightedPres"

glm_mod_s<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))
gam_mod_s<- readRDS(paste(model_results_path,"gam_mod_s.rds", sep="/"))
rf_mod_s<-  readRDS(paste(model_results_path,"rf_mod_s.rds", sep="/"))
brt_mod_s<- readRDS(paste(model_results_path,"brt_mod_s.rds", sep="/"))

models <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")


# Get variables used in the model
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)


#  c( "slope_7x7", "ammonium_spring_SD", "ammonium_winter_mean","PAR_summer_mean", "temperature_summer_mean",
# "turbidity_summer_mean","salinity_summer_SD", "nitrate_summer_mean", "ammonium_summer_mean", "ammonium_spring_mean")   


### Scale variables before predictions =========================================
# We need the original training data to scale with same parameters
# train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob.csv")
train<- read.csv(paste(model_results_path,"train_selected_table_FINALMODELS_M7.csv", sep="/"))
train<- train[,-1]
names(train)

train_sel <- train %>% dplyr::select(all_of(c("kelp", vars_selected)))
train_sel$kelp<- as.factor(train_sel$kelp)

train_scaled_b <- scale(train_sel[, vars_selected])
train_means <- attr(train_scaled_b, "scaled:center")
train_sds   <- attr(train_scaled_b, "scaled:scale")

raster_stack_predict_scaled<- raster_stack_20m_all[[names(raster_stack_20m_all)%in% vars_selected]]  # original raster stack

for (v in vars_selected) {
  raster_stack_predict_scaled[[v]] <- (raster_stack_20m_all[[v]] - train_means[v]) / train_sds[v]
}
names(raster_stack_predict_scaled)#<- names(raster_stack_predict)

# saveRDS(raster_stack_predict_scaled, "raster_stack_predict_scaled_postBlob.rds")
writeRaster(raster_stack_predict_scaled, "raster_stack_predict_scaled_postBlob.tif", overwrite=TRUE)
writeRaster(raster_stack_20m_all, "raster_stack_predict_postBlob.tif", overwrite=TRUE)



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
predict_gam <- function(model, data) {
  mgcv::predict.gam(model, newdata = data, type = "response")
}

terra::predict(
  raster_stack_predict_scaled,
  gam_mod_s,                    # supply model explicitly
  fun = predict_gam,            # must accept (model, data)
  na.rm = TRUE,
  filename = "gam_postBlob_pred_raster.tif",
  overwrite = TRUE
)

# terra::predict(raster_stack_predict_scaled, #glm took a lot of memory so I saved directly to disk
#                gam_mod_s,
#                type = "response",
#                na.rm = TRUE,
#                filename = "gam_postBlob_pred_raster.tif",  # output written directly to disk
#                overwrite = TRUE)


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
ens_ave_postblob <- ensemble_raster(pred_rasters_list, tss_weights_1) # 2 min
plot(ens_ave_postblob$raster)


# # Uncertainty raster
# start<- Sys.time()
ensemble_model_postBlob<- uncertainty_raster(pred_rasters_list, weights_df = tss_weights_1[,c(1,2)])
# end<- Sys.time()
# time_total<- end-start
# time_total # 30 min ; M4 55.37769 mins
# 
# saveRDS(ensemble_model_postBlob, "ensemble_model_postBlob.rds")


# Save rasters
writeRaster(ens_ave_postblob$raster, "ens_ave_postblob_M7.tif", overwrite=T)
writeRaster(ensemble_model_postBlob$max, "ensemble_model_max_postBlob_M7.tif")
writeRaster(ensemble_model_postBlob$min, "ensemble_model_min_postBlob_M7.tif")
writeRaster(ensemble_model_postBlob$sd, "ensemble_model_wSD_postBlob_M7.tif")


# ==============================================================================
# STEP 3: Binary results & Evaluation with presence/absence data from Post blob
# ==============================================================================
# Load prediction rasters
ensemble_raster_PCA<- rast("ensemble_PCA_postblob_M7.tif")
ensemble_raster_ave<- rast("ens_ave_postblob_M7.tif")

# Load testing points (kelp centroids for the period 2020-2022 plus absenses, 1 per cell of 500x500m)
test_data_df<- read.csv("/SDM/Presence_absences_kelp/MoraSoto_postblob/Testing_Presence_absences_kelp_2020_2022_filtered.csv")
test_data_sat_pts<- vect(test_data_df, geom = c("x", "y"), crs = "EPSG:3005")

threshold_drop <-  0.1957
thr_maxTSS <- 0.6977754
thresholds_table<- tibble(Model= c("ave_postblob", "ave_postblob"), Threshold = c("threshold_drop", "thr_maxTSS"), threshold_value= c(threshold_drop, thr_maxTSS))
thresholds_table

# function_testing_table<- function(prediction_name= "post_blob",
#          testing_dataset_name= "satellite_dataset",
#          threshold= "thr_maxTSS", thresh_value=0.6977754,
#          pred_raster= ensemble_raster_ave, test_points= test_data_sat_pts){
#   # extract predicted values at testing points
#   pred_vals<- extract(ensemble_raster_ave, test_data_sat_pts)
#   test_data$ensemble_pred<- pred_vals[,-1]
#   
#   ### Get metrics of model performance, AUC, Sensitivity and specificity
#   # Compute ROC for ensemble
#   roc_ens <- pROC::roc(response = test_data$kelp,
#                        predictor = test_data$ensemble_pred)
#   
#   # Youden's J statistic to find optimal threshold
#   sens_vec <- roc_ens$sensitivities
#   spec_vec <- roc_ens$specificities
#   thr_vec  <- roc_ens$thresholds
#   J <- sens_vec + spec_vec - 1
#   
#   
#   ### Binary predictions at optimal threshold  ======================
#   if(threshold == "maxTSS"){
#     pred_class <- ifel(ensemble_raster_ave$lyr1 >= thresh_value, 1, 0)
#   }
#   if(threshold == "0.1_Omission"){
#     pred_class <- ifel(ensemble_raster_ave$lyr1 >= thresh_value, 1, 0)
#   }
#   
#   # Manual metrics
#   pred_vals <- terra::extract(pred_class, test_data_sat_pts)  # extract raster preds at test points
#   # Combine truth with predictions (drop ID col)
#   eval_df <- data.frame(
#     truth = test_data_sat_pts$kelp,
#     pred  = pred_vals[,2]   # second column = raster values
#   )
#   
#   # Remove NAs (some points may fall outside raster extent)
#   eval_df <- na.omit(eval_df)
#   summary(eval_df)
#   
#   # Confusion matrix components
#   tp <- sum(eval_df$pred == 1 & eval_df$truth == 1)
#   tn <- sum(eval_df$pred == 0 & eval_df$truth == 0)
#   fp <- sum(eval_df$pred == 1 & eval_df$truth == 0)
#   fn <- sum(eval_df$pred == 0 & eval_df$truth == 1)
#   
#   sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
#   specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
#   auc_val     <- as.numeric(pROC::auc(roc_ens))
#   tss         <- sensitivity + specificity - 1
#   
#   return(tibble(metrics= c("sensitivity","specificity","auc_val", "tss"), 
#                 metrics_vals= c(sensitivity,specificity,auc_val, tss),
#                 prediction= rep(paste(prediction_name),4),
#                 testing_dataset= rep(paste(testing_dataset_name),4),
#                 threshold= rep(paste(threshold),4)))
# 
# }

function_testing_table_balanced <- function(prediction_name = "post_blob",
                                            testing_dataset_name = "satellite_dataset",
                                            threshold = "maxTSS", 
                                            thresh_value = 0.6977754,
                                            pred_raster = ensemble_raster_ave, 
                                            test_points = test_data_sat_pts,
                                            n_reps = 100,        
                                            balanced = TRUE){
  # Extract predicted values at testing points
  pred_vals <- terra::extract(pred_raster, test_points)
  
  # Convert to data.frame (kelp must be in test_points!)
  test_data <- as.data.frame(test_points)
  test_data$ensemble_pred <- pred_vals[, -1]
  
  # Remove NAs (points falling outside raster)
  test_data <- na.omit(test_data)
  
  # If not balancing → run once -----------------
  if(!balanced){
    roc_ens <- pROC::roc(response = test_data$kelp,
                         predictor = test_data$ensemble_pred,
                         quiet = TRUE)
    pred_class <- ifelse(test_data$ensemble_pred >= thresh_value, 1, 0)
    
    tp <- sum(pred_class == 1 & test_data$kelp == 1)
    tn <- sum(pred_class == 0 & test_data$kelp == 0)
    fp <- sum(pred_class == 1 & test_data$kelp == 0)
    fn <- sum(pred_class == 0 & test_data$kelp == 1)
    
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    auc_val     <- as.numeric(pROC::auc(roc_ens))
    tss         <- sensitivity + specificity - 1
    
    return(tibble(metrics= c("sensitivity","specificity","auc_val","tss"), 
                  mean= c(sensitivity,specificity,auc_val,tss),
                  sd  = c(NA,NA,NA,NA),
                  prediction= prediction_name,
                  testing_dataset= testing_dataset_name,
                  threshold= threshold))
  }
  
  # If balancing → repeated subsampling ----------
  results <- replicate(n_reps, {
    pres <- dplyr::filter(test_data, kelp == 1)
    pres <- dplyr::sample_n(pres, 500)
    
    # abs <- dplyr::filter(test_data, kelp == 0)
    abs  <- dplyr::filter(test_data, kelp == 0) %>%
      dplyr::slice_sample(n = nrow(pres))
    samp <- dplyr::bind_rows(pres, abs)
    
    roc_ens <- pROC::roc(response = samp$kelp,
                         predictor = samp$ensemble_pred,
                         quiet = F)
    pred_class <- ifelse(samp$ensemble_pred >= thresh_value, 1, 0)
    
    tp <- sum(pred_class == 1 & samp$kelp == 1)
    tn <- sum(pred_class == 0 & samp$kelp == 0)
    fp <- sum(pred_class == 1 & samp$kelp == 0)
    fn <- sum(pred_class == 0 & samp$kelp == 1)
    
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    auc_val     <- as.numeric(pROC::auc(roc_ens))
    tss         <- sensitivity + specificity - 1
    
    c(sensitivity, specificity, auc_val, tss)
  })
  
  results <- t(results)
  colnames(results) <- c("sensitivity","specificity","auc_val","tss")
  
  # Summarize mean ± SD
  res_summary <- data.frame(
    metrics = colnames(results),
    mean    = apply(results, 2, mean, na.rm = TRUE),
    sd      = apply(results, 2, sd, na.rm = TRUE),
    prediction = prediction_name,
    testing_dataset = testing_dataset_name,
    threshold = threshold
  )
  
  return(res_summary)
}



test_data%>%
  group_by(kelp)%>%
  summarize(n=length(kelp))


post_blob_sat_thres1<- function_testing_table_balanced(prediction_name= "post_blob",
                                  testing_dataset_name= "satellite_dataset",
                                  threshold= "maxTSS",
                                  thresh_value= 0.6977754,
                                  balanced=T,
                                  n_reps=75,
                                  pred_raster= ensemble_raster_ave, test_points= test_data_sat_pts)

post_blob_sat_thres2<- function_testing_table_balanced(prediction_name= "post_blob",
                                              testing_dataset_name= "satellite_dataset",
                                              threshold= "0.1_Omission",
                                              thresh_value= threshold_drop,
                                              balanced=T,
                                              n_reps=75,
                                              pred_raster= ensemble_raster_ave, test_points= test_data_sat_pts)

rbind(post_blob_sat_thres2, post_blob_sat_thres2_balanced)


### Test predictions with substrate mask 
test_data_df$substrate<- (as.factor(test_data_sat_pts_substrate$substrate))

test_data_pts<-  vect(test_data_df, geom = c("x", "y"), crs = "EPSG:3005")
pred_vals <- terra::extract(pred_raster, test_data_pts)

# Convert to data.frame (kelp must be in test_points!)
test_data_df$ensemble_pred <- pred_vals[, -1]
test_data_df$ensemble_pred <- as.numeric(test_data_df$ensemble_pred )

test_data_df%>%
  group_by(kelp, substrate)%>%
  summarize(n= length(kelp))
#    kelp substrate     n
# 1     0 1          1795 
# 2     0 2          1900
# 3     1 1           863
# 4     1 2            46 # to be excluded form testing analysis

# Exclude presences in soft substrate --> potential bias in substrate model 
# test_data_df<- test_data_df[-which(test_data_df$substrate =="2" & test_data_df$kelp == 1),]

# Apply mask of substrate and depth to the testing points by comverting predicted values to insuitable (zero/0)
test_data_df[which(test_data_df$substrate!= 1), "ensemble_pred"] <- "0"
test_data_df[which(test_data_df$depth>= 10), "ensemble_pred"] <- "0"


test_data_pts_masked<-  vect(test_data_df, geom = c("x", "y"), crs = "EPSG:3005")

post_blob_sat_thres1_masked<- function_testing_table_balanced(prediction_name= "post_blob_masked",
                                                       testing_dataset_name= "satellite_dataset",
                                                       threshold= "maxTSS",
                                                       thresh_value= 0.6977754,
                                                       balanced=T,
                                                       n_reps=75,
                                                       pred_raster= ensemble_raster_ave, test_points= test_data_pts_masked)

post_blob_sat_thres2_masked<- function_testing_table_balanced(prediction_name= "post_blob_masked",
                                                                testing_dataset_name= "satellite_dataset",
                                                                threshold= "0.1_Omission",
                                                                thresh_value= threshold_drop,
                                                                balanced=T,
                                                                n_reps=75,
                                                                pred_raster= ensemble_raster_ave, test_points= test_data_pts_masked)

testing_table_result<- rbind(post_blob_sat_thres1, post_blob_sat_thres1_masked, post_blob_sat_thres2, post_blob_sat_thres2_masked)


# testing_table_result<- rbind(post_blob_sat_thres1, post_blob_sat_thres2, post_blob_sat, post_blob_sat_masked)


## Add values of ensemble model into table of results ========
# results_selected[5,1:6]<- list("ens", auc_val, sensitivity, specificity,tss)
write.csv(testing_table_result, "results_PostBlob_testing_Table_M7.csv")
