

# sensitivity / gradient mapping approach with an ensemble (average) of models (GLM, GAM, RF, BRT). It:
# Computes the ensemble prediction for the warm baseline (once).
# For each predictor, creates a warm stack where only that predictor is replaced by its cold-period layer.
# Predicts with each model (streaming to disk to avoid memory overload), averages model predictions to produce the ensemble prediction for that “single-variable-changed” scenario.
# Computes per-variable contribution rasters = (ensemble_with_only_var_changed − ensemble_warm).
# Writes contribution rasters to disk and computes a dominant-driver raster (which variable has the largest contribution at each pixel).

# Memory-efficient sensitivity mapping with terra + ensemble averaging
# Requirements: terra, pbapply (for progress), optionally foreach/doParallel if you want parallel
# Assumptions:
#  - env_warm and env_cold are SpatRasters (same extent/resolution/CRS and same layer order)
#  - names(env_warm) give predictor names
#  - mod_list is a list of fitted models. Each element can be:
#       * a model object (glm, gam, randomForest, gbm, etc.) OR
#       * a list with elements: $model (the model object) and $predict_args (a named list passed to predict)
#  - terra::predict works for your model objects (it calls base predict(model, newdata)). If one model needs special args
#    (e.g. gbm: n.trees), put them in predict_args.
#
# Example mod_list item:
#   mod_list <- list(glm_model, gam_model, list(model = rf_model, predict_args = list(type="response")),
#                    list(model = brt_model, predict_args = list(n.trees=1500)))
#
# Edit these paths and names:
out_dir <- "variables_contrib_change_outputs"   # directory to store temporary rasters and results
dir.create(out_dir, showWarnings = FALSE)
terraOptions(tempdir = out_dir)     # make terra use this directory for temp files

library(terra)
library(pbapply)

# ---------------------------
# Load or prepare your data
# ---------------------------
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres"

# Ensemble predictions for both periods (blob=env_warm and postblob=env_cold)
# env_warm <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif")   # SpatRaster with predictors as layers
# env_cold <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/ens_ave_postblob_M7.tif")

env_warm <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/raster_stack_predict.tif")   # SpatRaster with predictors as layers
env_cold <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/raster_stack_predict_postBlob.tif")

# Models 
glm_mod_s<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))
gam_mod_s<- readRDS(paste(model_results_path,"gam_mod_s.rds", sep="/"))
rf_mod_s<-  readRDS(paste(model_results_path,"rf_mod_s.rds", sep="/"))
brt_mod_s<- readRDS(paste(model_results_path,"brt_mod_s.rds", sep="/"))

mod_list <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)


# For the script to run, env_warm/env_cold/mod_list must exist in your environment.
# Check they align:
stopifnot(exists("env_warm"), exists("env_cold"), exists("mod_list"))
stopifnot(nlyr(env_warm) == nlyr(env_cold))
stopifnot(all(names(env_warm) == names(env_cold)))

pred_names <- names(env_warm)
nvars <- nlyr(env_warm)
nmodels <- length(mod_list)

# Optional: model weights (defaults to equal weighting)
# Provide a numeric vector of length nmodels if you want performance-weighted ensemble
if(!exists("model_weights")){
  model_weights <- rep(1, nmodels)
} else {
  stopifnot(length(model_weights) == nmodels)
}
model_weights <- model_weights / sum(model_weights)

# ---------------------------
# Helper: safe predict wrapper
# ---------------------------
# Accepts: model_entry (either model object or list with $model and $predict_args)
# Returns: SpatRaster (prediction) written to temp file (filename argument)
safe_predict_model <- function(env_stack, model_entry, filename, overwrite=TRUE){
  # Unpack model and predict_args
  if(is.list(model_entry) && !inherits(model_entry, "model")){
    model_obj <- model_entry$model
    predict_args <- model_entry$predict_args %||% list()
  } else {
    model_obj <- model_entry
    predict_args <- list()
  }
  # Use terra::predict which streams to disk if filename provided
  # Predict signature: predict(x=SpatRaster, model, ... , filename = "path")
  # merge predict_args with filename
  args <- c(list(x = env_stack, model = model_obj, filename = filename, overwrite = overwrite), predict_args)
  # call predict with do.call
  do.call(terra::predict, args)
}

# "null-coalesce" helper
`%||%` <- function(a, b) if(!is.null(a)) a else b

# ---------------------------
# Step 0: warm-period ensemble (baseline)
# ---------------------------
# We compute ensemble prediction for env_warm ONCE and store to disk.
warm_ensemble_file <- file.path(out_dir, "ensemble_warm.tif")
if(!file.exists(warm_ensemble_file)){
  message("Computing ensemble prediction for the warm baseline (streaming to disk)...")
  # For each model, predict and save a temporary raster
  tmp_files <- pbapply::pblapply(seq_len(nmodels), function(i){
    tf <- file.path(out_dir, paste0("pred_warm_model", i, ".tif"))
    safe_predict_model(env_warm, mod_list[[i]], filename = tf)
    return(tf)
  })
  # Read predicted rasters and compute weighted average (streamed)
  preds_rasters <- lapply(tmp_files, rast)
  # Weighted sum: multiply each raster by weight and sum; write result to disk
  # We'll compute sum(weights * raster) using terra::app and stack
  preds_stack <- rast(preds_rasters)
  # multiply each layer by its weight
  for(i in seq_len(nlayers(preds_stack))){
    preds_stack[[i]] <- preds_stack[[i]] * model_weights[i]
  }
  # Sum layers into ensemble_warm
  ensemble_warm <- app(preds_stack, sum, filename = warm_ensemble_file, overwrite=TRUE)
  # optionally remove per-model warm predictions to save disk (keep only ensemble)
  file.remove(unlist(tmp_files))
} else {
  message("Found existing warm ensemble raster; loading.")
  ensemble_warm <- rast(warm_ensemble_file)
}

# ---------------------------
# Step 1: Loop variables (streaming predictions)
# ---------------------------
message("Computing per-variable contributions (this may take a while).")
contrib_files <- character(nvars)

# We will process each variable sequentially; for each:
#  - create env_temp where only variable j is replaced by env_cold[[j]]
#  - predict with each model (stream to disk), compute weighted ensemble, subtract ensemble_warm
#  - write contribution raster to disk

for(j in seq_len(nvars)){
  varname <- pred_names[j]
  message(sprintf("[%d/%d] Processing variable: %s", j, nvars, varname))
  # Create env_temp as copy of env_warm, replace layer j
  env_temp <- env_warm
  env_temp[[j]] <- env_cold[[j]]
  # Predict per-model and write temp files
  tmp_files_j <- pbapply::pblapply(seq_len(nmodels), function(i){
    tf <- file.path(out_dir, sprintf("pred_var%s_model%d.tif", j, i))
    safe_predict_model(env_temp, mod_list[[i]], filename = tf)
    return(tf)
  })
  # Weighted ensemble for this variable
  preds_stack_j <- rast(unlist(tmp_files_j))
  for(i in seq_len(nlayers(preds_stack_j))){
    preds_stack_j[[i]] <- preds_stack_j[[i]] * model_weights[i]
  }
  ensemble_var_file <- file.path(out_dir, paste0("ensemble_var_", j, ".tif"))
  ensemble_var <- app(preds_stack_j, sum, filename = ensemble_var_file, overwrite=TRUE)
  # contribution = ensemble_var - ensemble_warm
  contrib_file <- file.path(out_dir, paste0("contrib_", j, "_", varname, ".tif"))
  # Use terra::writeRaster with calculation on the fly to avoid loading everything:
  contrib_rast <- ensemble_var - ensemble_warm
  writeRaster(contrib_rast, filename = contrib_file, overwrite = TRUE)
  contrib_files[j] <- contrib_file
  # Clean up per-model temp files for this variable and ensemble_var intermediate
  file.remove(unlist(tmp_files_j))
  file.remove(ensemble_var_file)
  # remove objects to free memory
  rm(env_temp, preds_stack_j, ensemble_var, contrib_rast); gc()
}

# ---------------------------
# Step 2: Stack contributions and compute dominant driver
# ---------------------------
message("Stacking contributions and computing dominant driver (which var had max contribution at each pixel).")
contrib_stack <- rast(contrib_files)
names(contrib_stack) <- pred_names

# Write stacked contributions to disk
contrib_stack_file <- file.path(out_dir, "contrib_stack.tif")
writeRaster(contrib_stack, filename = contrib_stack_file, overwrite = TRUE)

# Compute the dominant driver: index of layer with largest contribution in absolute value or positive contribution?
# Choose logic:
#  - If you want the variable with the largest POSITIVE contribution (increased suitability), use which.max(contrib_stack)
#  - If you want the variable with the largest ABSOLUTE contribution (biggest magnitude change, sign ignored), compute abs then which.max
# Here we'll compute BOTH.
dominant_positive_file <- file.path(out_dir, "dominant_positive.tif")
dominant_abs_file      <- file.path(out_dir, "dominant_abs.tif")

# dominant positive (largest positive contribution)
dominant_positive <- which.max(contrib_stack, filename = dominant_positive_file, overwrite = TRUE)

# dominant by absolute magnitude
contrib_abs_stack <- abs(contrib_stack)
dominant_abs      <- which.max(contrib_abs_stack, filename = dominant_abs_file, overwrite = TRUE)

# Save also the magnitude of the winning variable (absolute) for visualization
winning_magnitude_file <- file.path(out_dir, "winning_magnitude.tif")
winning_magnitude <- app(contrib_abs_stack, max, filename = winning_magnitude_file, overwrite = TRUE)

# ---------------------------
# Step 3: Optional summaries & maps
# ---------------------------
# Compute mean contribution per variable (global summary)
mean_contribs <- sapply(seq_len(nlayers(contrib_stack)), function(i){
  global(contrib_stack[[i]], fun = "mean", na.rm = TRUE)
})
mean_contribs_df <- data.frame(variable = pred_names, mean_contribution = as.numeric(mean_contribs))
mean_contribs_df <- mean_contribs_df[order(-abs(mean_contribs_df$mean_contribution)), ]
print("Mean contribution across landscape (ranked):")
print(mean_contribs_df)

# Also compute area (count of pixels) where each variable is the dominant positive driver
dom_pos_counts <- table(values(dominant_positive))
dom_pos_df <- data.frame(layer_index = as.integer(names(dom_pos_counts)), pixels = as.integer(dom_pos_counts))
if(nrow(dom_pos_df) > 0){
  dom_pos_df$variable <- pred_names[dom_pos_df$layer_index]
  dom_pos_df <- dom_pos_df[order(-dom_pos_df$pixels), ]
  print("Pixel counts where variable is dominant positive driver:")
  print(dom_pos_df)
}



# ---------------------------
# Step 4: Percent contributions per pixel
# ---------------------------

# contrib_stack already exists from Step 2:
#   a SpatRaster where each layer = contribution raster for one variable

# 1) Take absolute contributions (magnitude)
contrib_abs_stack <- abs(contrib_stack)

# 2) Compute total contribution magnitude at each pixel
total_abs <- app(contrib_abs_stack, sum, filename = file.path(out_dir, "total_abs_contrib.tif"), overwrite = TRUE)

# 3) Compute percent contribution for each variable
#    (|contribution_i| / total_abs) * 100
percent_stack <- contrib_abs_stack / total_abs * 100

# Write to disk
percent_stack_file <- file.path(out_dir, "percent_contrib_stack.tif")
writeRaster(percent_stack, filename = percent_stack_file, overwrite = TRUE)

# 4) Optionally, compute max percent and dominant variable (redundant with which.max)
max_percent_file <- file.path(out_dir, "max_percent.tif")
max_percent <- app(percent_stack, max, filename = max_percent_file, overwrite = TRUE)

dominant_percent_file <- file.path(out_dir, "dominant_percent.tif")
dominant_percent <- which.max(percent_stack, filename = dominant_percent_file, overwrite = TRUE)

# ---------------------------
# Step 5: Summaries
# ---------------------------
# Mean percent contribution of each variable across the landscape
mean_percents <- sapply(seq_len(nlayers(percent_stack)), function(i){
  global(percent_stack[[i]], fun = "mean", na.rm = TRUE)
})
mean_percents_df <- data.frame(variable = names(percent_stack), mean_percent = as.numeric(mean_percents))
mean_percents_df <- mean_percents_df[order(-mean_percents_df$mean_percent), ]
print("Mean percent contribution across landscape (ranked):")
print(mean_percents_df)

# ---------------------------
# Key outputs:
# ---------------------------
message("Percent contribution rasters written to: ", percent_stack_file)
message("Dominant variable (by percent) raster: ", dominant_percent_file)
message("Maximum percent contribution raster: ", max_percent_file)

# Example plotting:
# plot(rast(percent_stack_file)[[1]])   # percent contribution map for first variable
# plot(rast(dominant_percent_file))     # map showing index of dominant variable by percent













# ---------------------------
# Done
# ---------------------------
message("Finished. Outputs saved to: ", normalizePath(out_dir))
message("Key files:")
message(" - ensemble_warm: ", warm_ensemble_file)
message(" - contrib stack: ", contrib_stack_file)
message(" - dominant positive raster: ", dominant_positive_file)
message(" - dominant absolute raster: ", dominant_abs_file)
message(" - winning magnitude: ", winning_magnitude_file)

# Examples to load results for plotting:
# contrib_stack <- rast(file.path(out_dir, "contrib_stack.tif"))
# plot(contrib_stack[[1]])  # plot contribution of first variable
# plot(rast(dominant_positive_file)) # categorical map showing index of dominant positive variable
