
# Function to predict ensemble model with new environemtnal conditions (rasters or dataframe)
predict_model_sdm <- function(model, model_name, newdata, 
                              scale_params = NULL, 
                              quad_vars = NULL,
                              type = "response") {
  
  # If newdata is a SpatRaster -> extract values first
  is_raster <- inherits(newdata, "SpatRaster")
  if (is_raster) {
    newdata_df <- as.data.frame(newdata, na.rm = FALSE)
  } else {
    newdata_df <- newdata
  }
  
  # --- 1. Scaling step (only if scale_params provided) ---
  if (!is.null(scale_params) && model_name %in% c("glm", "gam")) {
    
    mean_vals <- scale_params[1, grep("_mean$", names(scale_params))]
    sd_vals   <- scale_params[1, grep("_sd$", names(scale_params))]
    names(mean_vals) <- sub("_mean$", "", names(mean_vals))
    names(sd_vals)   <- sub("_sd$", "", names(sd_vals))
    
    newdata_scaled <- as.data.frame(lapply(names(newdata_df), function(nm) {
      if (nm %in% names(mean_vals)) {
        (newdata_df[[nm]] - mean_vals[[nm]]) / sd_vals[[nm]]
      } else {
        newdata_df[[nm]]
      }
    }))
    names(newdata_scaled) <- names(newdata_df)
  } else {
    newdata_scaled <- newdata_df
  }
  
  # --- 2. Add quadratic terms for GLM only ---
  if (model_name == "glm" && !is.null(quad_vars)) {
    for (v in quad_vars) {
      if (v %in% names(newdata_scaled)) {
        newdata_scaled[[paste0("I(", v, "^2)")]] <- newdata_scaled[[v]]^2
      }
    }
  }
  
  # --- 3. Model-specific prediction logic ---
  if (model_name == "glm") {
    pred <- predict(model, newdata = newdata_scaled, type = "response")
  } else if (model_name == "gam") {
    pred <- predict(model, newdata = newdata_scaled, type = "response")
  } else if (model_name == "rf") {
    # Random forest: return probability of presence
      pred <- predict(model, newdata = newdata_df, type = "prob")[, 2]
  } else if (model_name == "brt") {
    # GBM / BRT models — need n.trees argument
    ntrees <- if (!is.null(model$n.trees)) {
      model$n.trees
    } else if (!is.null(model$gbm.call$best.trees)) {
      model$gbm.call$best.trees
    } else {
      ntrees <- 4650 # from the saived trained model info
    }
    # Make sure we use gbm::predict.gbm directly (not generic)
    # if ("gbm" %in% loadedNamespaces()) {
      pred <- gbm::predict.gbm(model, newdata = newdata_df, n.trees = ntrees, type = "response")
    # } else {
    #   pred <- predict(model, newdata = newdata_df, n.trees = ntrees, type = "response")
    # }
  } else {
    stop("Unknown model type: ", model_name)
  }
  
  # --- 4. Reattach to raster (if input was a raster) ---
  if (is_raster) {
    pred_rast <- newdata[[1]]
    values(pred_rast) <- pred
    return(pred_rast)
  } else {
    return(pred)
  }
}

