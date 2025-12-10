


library(terra)
library(pROC)
library(dismo)
library(ecospat)
library(dplyr)

#---------------------------------------------------
# Helper to compute metrics for one model, threshold, and period
#---------------------------------------------------
calc_metrics <- function(true_vals, pred_vals, threshold, pred_raster = NULL, presence_pts = NULL) {
  
  # Remove NA predictions
  ok <- !is.na(pred_vals)
  true_vals <- true_vals[ok]
  pred_vals <- pred_vals[ok]
  
  # --- Threshold-independent metrics ---
  auc_val <- tryCatch(as.numeric(roc(true_vals, pred_vals)$auc), error = function(e) NA)
  
  # --- Point-based Boyce (presence/background only) ---
  boyce_point <- tryCatch({
    pres_scores <- pred_vals[true_vals == 1]
    bg_scores   <- pred_vals[true_vals == 0]
    ecospat::ecospat.boyce(fit = c(bg_scores, pres_scores),
                  obs = pres_scores,
                  nclass = 10, res = 100)$cor
  }, error = function(e) NA)
  
  # --- Raster-based Boyce (optional, if raster and presence points are provided) ---
  boyce_rast <- tryCatch({
    if (!is.null(pred_raster) && !is.null(presence_pts)) {
      pres_vals <- terra::extract(pred_raster, presence_pts)[, 2]
      pres_vals <- pres_vals[!is.na(pres_vals)]
      all_vals <- values(pred_raster)
      all_vals <- all_vals[!is.na(all_vals)]
      ecospat.boyce(fit = all_vals, obs = pres_vals, nclass = 10, res = 100)$cor
    } else NA
  }, error = function(e) NA)
  
  # --- Threshold-dependent metrics ---
  binary_pred <- as.numeric(pred_vals >= threshold)
  TP <- sum(binary_pred == 1 & true_vals == 1)
  TN <- sum(binary_pred == 0 & true_vals == 0)
  FP <- sum(binary_pred == 1 & true_vals == 0)
  FN <- sum(binary_pred == 0 & true_vals == 1)
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy    <- (TP + TN) / (TP + TN + FP + FN)
  TSS         <- sensitivity + specificity - 1
  
  data.frame(
    AUC = auc_val,
    Boyce_point = boyce_point,
    Boyce_raster = boyce_rast,
    TSS = TSS,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Accuracy = accuracy,
    N_eval = length(true_vals)
  )
}





# Model evaluation and calibration procedure
# -----------------------------
# Probability calibration
# -----------------------------
# This function rescales the model-predicted probabilities (pred_vals) to account for differences in prevalence between:
# the training dataset (prev_train), and
# the testing or projection dataset (prev_test)
# This correction is based on the relationship between prevalence and probability calibration described in Guillera-Arroita et al. (2015):
# “Predicted probabilities from presence–background models are typically conditional on the prevalence in the training sample and must be recalibrated if the prevalence in the target population differs.”
calibrate_prob <- function(pred_vals, prev_train, prev_test) {
  # Prevalence-based calibration following Guillera-Arroita et al. (2015)
  odds <- pred_vals / (1 - pred_vals)
  odds_adj <- odds * ((prev_test / (1 - prev_test)) / (prev_train / (1 - prev_train)))
  p_cal <- odds_adj / (1 + odds_adj)
  return(p_cal)
}


# Model evaluation for each period was conducted using a custom R function (process_period()) that automated prediction extraction, calibration, and performance assessment. This function integrated raster-based model predictions with independent testing datasets and computed a suite of evaluation metrics.
# Specifically, the function
#   Loaded model prediction rasters (GeoTIFF format) corresponding to each modelling algorithm (GLM, GAM, BRT, RF, and ensemble). Rasters were read using the terra package, and model names were assigned based on file names.
# Converted the testing dataset (containing coordinates, presence–absence labels, and associated environmental predictors) into a spatial object (SpatVector) matching the coordinate reference system of the raster layers.
# Extracted predicted suitability values from each model raster at the testing point locations using terra::extract(). These values were then appended to the testing dataset to create a combined evaluation table containing both observed (kelp) and predicted values.
# Estimated prevalence in the testing dataset (proportion of presence records) and applied a probability calibration step to each model’s predicted suitability values using the function calibrate_prob(). This correction adjusted predicted probabilities to account for potential differences in prevalence between the training and testing datasets, following the approach described by Platt (1999) and Guillera-Arroita et al. (2015).
# Computed model performance metrics for each model and evaluation threshold (1% omission and maximum TSS). For each threshold, the function called calc_metrics() to derive:
#   Threshold-independent metrics: Area Under the Receiver Operating Characteristic Curve (AUC) and continuous Boyce index (both raster- and point-based versions).
# Threshold-dependent metrics: True Skill Statistic (TSS), sensitivity, specificity, and accuracy.
# Metrics were computed by comparing observed presences and absences with binary predictions derived using each threshold.
# Compiled results into a summary table containing all metrics, thresholds, and model identifiers for the evaluated period.

# raster_folder_p1, test_t1, prev_train, thresholds, model_names, prevalence=prev_PeriodsCombined


process_period <- function(raster_folder=raster_folder_p1, test_data= test_t1, prev_train, thresholds, model_names, prevalence) {

  # Load rasters
  raster_files <- list.files(raster_folder, pattern = "\\.tif$", full.names = TRUE)
  rasters <- rast(raster_files)
  names(rasters)<- model_names

  # Convert test points to SpatVector
  test_vect <- vect(test_data, geom = c("x", "y"), crs = crs(rasters))

  # Extract raster values
  pred_vals <- terra::extract(rasters, test_vect)
  test_eval <- cbind(test_data, pred_vals[,-1])

  # # Prevalence in test dataset
  # prevalence <- mean(test_eval$kelp)

  # # Apply calibration
  test_eval_cal<- test_eval
  # for (m in model_names) {
  #   test_eval_cal[[paste0(m, "_cal")]] <- calibrate_prob(test_eval_cal[[m]], prev_train, prevalence)
  # }

  # Calculate metrics for all models and thresholds
  results <- data.frame()
  for (m in model_names) {
    r_path <- list.files(raster_folder, pattern = paste0(m, ".*\\.tif$"), full.names = TRUE)[1]
    r <- rast(r_path)
    presence_pts <- vect(test_data[test_data$kelp == 1, ], geom = c("x", "y"), crs = crs(r))

    for (th_name in names(thresholds)) {
      th <- thresholds[[th_name]]
      res <- calc_metrics(
        true_vals = test_eval$kelp,
        pred_vals = test_eval[[paste0(m)]],
        threshold = th,
        pred_raster = r,
        presence_pts = presence_pts
      )
      res$Model <- m
      res$Threshold <- th_name
      results <- bind_rows(results, res)
    }
  }

  return(results)
}



# NEW VERSION OF PROCESS_PERIOD TO ADD AS OUTPUT THE CLASSIFICATION OF EACH TESTING POINT INTO TP, TN, FP OR FN
process_period_v2 <- function(raster_folder=raster_folder_p1, test_data= test_t1, prev_train, thresholds, model_names, prevalence) {
  
  # Load rasters
  raster_files <- list.files(raster_folder, pattern = "\\.tif$", full.names = TRUE)
  rasters <- rast(raster_files)
  names(rasters)<- model_names
  
  # Convert test points to SpatVector
  test_vect <- vect(test_data, geom = c("x", "y"), crs = crs(rasters))
  
  # Extract raster values
  pred_vals <- terra::extract(rasters, test_vect)
  test_eval <- cbind(test_data, pred_vals[,-1])
  
  # Calculate metrics for all models and thresholds
  all_classifications <- list()
  results<- data.frame()
  
  for (m in model_names) {
    
    r_path <- list.files(raster_folder, pattern = paste0(m, ".*\\.tif$"), full.names = TRUE)[1]
    r <- rast(r_path)
    presence_pts <- vect(test_data[test_data$kelp == 1, ], geom = c("x", "y"), crs = crs(r))
    
    for (th_name in names(thresholds)) {
      
      th <- thresholds[[th_name]]
      
      # ----- metrics (your previous code) -----
      res_raw <- calc_metrics(
        true_vals = test_eval$kelp,
        pred_vals = test_eval[[m]],
        threshold = th
      )
      
      # res_cal <- calc_metrics(
      #   true_vals = test_eval$kelp,
      #   pred_vals = test_eval_cal[[paste0(m, "_cal")]],
      #   threshold = th,
      #   pred_raster = r,
      #   presence_pts = presence_pts
      # )[ , c("AUC", "Boyce_point", "Boyce_raster")]
      
      res <- cbind(Model = m, Threshold = th_name, 
                   # res_cal, 
                   res_raw[ , c("TSS", "Sensitivity", "Specificity", "Accuracy")])
      
      results <- bind_rows(results, res)
      
      # ----- NEW: classification table -----
      cls <- classify_points(
        test_data = test_eval,
        pred_vals = test_eval[[m]],   # raw predictions for binary classification
        threshold = th,
        model_name = m
      )
      
      cls$ThresholdName <- th_name   # keep track of which threshold rule was used
      all_classifications[[length(all_classifications) + 1]] <- cls
    }
  }
  
  return(list(
    Metrics = results,
    Classification = bind_rows(all_classifications)
  ))
}





classify_points <- function(test_data, pred_vals, threshold, model_name) {
  
  # Remove NAs
  ok <- !is.na(pred_vals)
  df <- test_data[ok, ]
  preds <- pred_vals[ok]
  obs <- df$kelp
  
  # Binary classification
  binary <- as.numeric(preds >= threshold)
  
  # Determine class
  class <- ifelse(binary == 1 & obs == 1, "TP",
                  ifelse(binary == 0 & obs == 0, "TN",
                         ifelse(binary == 1 & obs == 0, "FP",
                                "FN")))
  
  out <- data.frame(
    Model = model_name,
    Threshold = threshold,
    x = df$x,
    y = df$y,
    Observed = obs,
    Predicted = preds,
    Binary = binary,
    Classification = class
  )
  
  return(out)
}








### FUnction to create buffer area ===========
# create_buffers_from_xy <- function(df, x_col = "x", y_col = "y",
#                                    crs_epsg = 3005, radius = 30) {
#   # df: dataframe with x, y, and kelp or ID or model predictions
#   # radius: buffer radius in map units (meters for EPSG:3005)
#   
#   # Convert to SpatVector points
#   pts <- terra::vect(df, geom = c(x_col, y_col), crs = paste0("EPSG:", crs_epsg))
#   
#   # Generate buffers
#   buf <- terra::buffer(pts, width = radius)
#   
#   # Attach attributes back to polygons
#   buf_df <- as.data.frame(pts)
#   buf <- cbind(buf, buf_df)
#   
#   # Add ID if needed
#   if (!"ID" %in% names(buf)) {
#     buf$ID <- 1:nrow(buf)
#   }
#   
#   return(buf)
# }


build_buffer_pred_df <- function(
    xy_df,
    model_names = model_names,
    raster_folder,
    x_col = "x",
    y_col = "y",
    crs_epsg = 3005,
    buffer_radius = 180,
    threshold = thresholds$drop1[1],
    id_col = "ID"
){
  library(terra)
  library(dplyr)
  
  #---------------------------
  # 1. Create buffer polygons
  #---------------------------
  pts <- vect(xy_df, geom = c(x_col, y_col), crs = paste0("EPSG:", crs_epsg))
  buf <- terra::buffer(pts, width = buffer_radius)
  
  # Add ID if missing
  if (!id_col %in% names(xy_df)) {
    xy_df[[id_col]] <- seq_len(nrow(xy_df))
  }
  buf[[id_col]] <- xy_df[[id_col]]
  
  #---------------------------
  # 2. Load all prediction rasters
  #---------------------------
  raster_files <- list.files(raster_folder, pattern = "\\.tif$", full.names = TRUE)
  rasters <- rast(raster_files)
  names(rasters) <- model_names  # e.g. BRT, RF, SVM, XGB, MLP
  
  # Prepare output df
  out <- xy_df
  
  for (m in model_names) {
    out[[paste0(m, "_avg")]]  <- NA
    out[[paste0(m, "_n")]]    <- NA
  }
  
  #---------------------------
  # 3. Loop through models
  #---------------------------
  for (m in model_names) {
    r <- rasters[[m]]
    
    cat("Processing model:", m, "\n")
    
    # Extract raster values for ALL polygons, binding buffer attributes
    ext <- terra::extract(r, buf, bind = TRUE)
    # ext now includes *all* columns from buf, including id_col
    
    #---------------------------
    # Iterate over each buffer polygon
    #---------------------------
    for (i in seq_len(nrow(buf))) {
      
      # Correctly retrieve buffer ID value
      current_id <- buf[i, id_col][[1]]
      current_id<- as.numeric(current_id)
      
      # Select extracted rows belonging to this polygon
      row_sel <- ext[[id_col]] == current_id
      
      cell_vals <- ext[row_sel, m]
      
      # Clean
      cell_vals <- cell_vals[!is.na(cell_vals)]
      
      # Keep cells ≥ threshold
      good_vals <- cell_vals[cell_vals >= threshold]
      
      if (length(good_vals) == 0) {
        
        # --- NEW: use 90th percentile of all cell values ---
        p90 <- quantile(cell_vals, 0.90, na.rm = TRUE)
        
        top_vals <- cell_vals[cell_vals >= p90]
        
        out[i, paste0(m, "_avg")]  <- mean(top_vals)
        out[i, paste0(m, "_n")]    <- length(top_vals)
        
        
      } else {
        # Standard case: cells ≥ threshold exist
        out[i, paste0(m, "_avg")]  <- mean(good_vals)
        out[i, paste0(m, "_n")]    <- length(good_vals)
      }
    }
  }
  
  return(out)
}



### FUnction to test predictions using buffer area ===========
process_buffer_predictions <- function(
    buffer_pred_df,     # dataframe with: ID, x, y, kelp, model1, model2, ... model5 (averaged)
    thresholds,         # list of thresholds
    model_names         # vector of model names = column names in buffer_pred_df
) {
  
  # test_data identical to buffer-level predictions:
  test_data <- buffer_pred_df
  
  # Ensure classification functions have correct fields
  test_data$x <- buffer_pred_df$x
  test_data$y <- buffer_pred_df$y
  
  # Prepare output containers
  all_metrics <- list()
  all_classifications <- list()
  
  # Loop models
  for (m in model_names) {
    
    # Extract predicted values for this model (already averaged)
    pred_vals <- test_data[[paste(m, "avg", sep="_")]]
    
    # True observed kelp pres/abs
    true_vals <- test_data$kelp
    
    # Loop thresholds
    for (th_name in names(thresholds)) {
      
      th <- thresholds[[th_name]]
      
      # -------------------------
      # 1. Metrics (raw predictions)
      # -------------------------
      res_raw <- calc_metrics(
        true_vals = true_vals,
        pred_vals = pred_vals,
        threshold = th
      )
      
      # Keep only metrics that do not require raster
      res <- data.frame(
        Model = m,
        Threshold = th_name,
        AUC = res_raw$AUC,
        Boyce_point = res_raw$Boyce_point,
        Boyce_raster = NA,   # cannot be estimated without raster
        TSS = res_raw$TSS,
        Sensitivity = res_raw$Sensitivity,
        Specificity = res_raw$Specificity,
        Accuracy = res_raw$Accuracy,
        N_eval = res_raw$N_eval
      )
      
      all_metrics[[length(all_metrics) + 1]] <- res
      
      # -------------------------
      # 2. Classification table
      # -------------------------
      cls <- classify_points(
        test_data = test_data,
        pred_vals = pred_vals,
        threshold = th,
        model_name = m
      )
      
      cls$ThresholdName <- th_name
      all_classifications[[length(all_classifications) + 1]] <- cls
    }
  }
  
  # Return results in same structure as your original process function
  return(list(
    Metrics = bind_rows(all_metrics),
    Classification = bind_rows(all_classifications)
  ))
}

