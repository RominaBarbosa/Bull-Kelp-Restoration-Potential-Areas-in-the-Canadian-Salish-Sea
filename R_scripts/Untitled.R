

library(raster)
library(mgcv)
library(randomForest)

#' Compute limiting variable per cell following Elith et al., 2010
#'
#' @param models List of models. Each model should be a named list: list(model=..., type="GLM/GAM/RF/BRF")
#' @param raster_stack RasterStack or RasterBrick of predictor variables
#' @param scale_params Optional: list of scale parameters (mean, sd) for each variable, needed for GLM/GAM
#' @param n_samples Optional: number of raster cells to sample for computation
#' @param buffer Optional: minimum distance (in cells) between sampled cells
#' @param var_range_points Number of points to sample along each variable range (default 50)
#' @return RasterLayer with limiting variable index (1 = first variable, etc.)
compute_limiting_variable <- function(models, raster_stack, scale_params=NULL, 
                                      n_samples=NULL, buffer=NULL, var_range_points=50) {
  
  # Get raster values
  vals <- as.data.frame(getValues(raster_stack))
  n_cells <- nrow(vals)
  n_vars <- ncol(vals)
  var_names <- names(raster_stack)
  
  # Subsample cells if requested
  sample_idx <- 1:n_cells
  if(!is.null(n_samples)) {
    if(!is.null(buffer)) {
      # Simple spatial sampling with buffer (distance in cell indices, approximate)
      # Here we assume regular raster and ignore exact coordinates for simplicity
      sampled <- c()
      remaining <- 1:n_cells
      while(length(sampled) < n_samples & length(remaining) > 0) {
        idx <- remaining[1]
        sampled <- c(sampled, idx)
        remaining <- remaining[abs(remaining - idx) > buffer]
      }
      sample_idx <- sampled
    } else {
      sample_idx <- sample(1:n_cells, n_samples)
    }
  }
  
  # Precompute variable ranges
  var_ranges <- lapply(vals, function(x) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=var_range_points))
  
  limiting_var_idx <- rep(NA, n_cells)
  
  for(cell in sample_idx) {
    obs_vals <- vals[cell, ]
    
    # Scale if needed
    if(!is.null(scale_params)) {
      for(i in 1:n_vars) {
        obs_vals[[i]] <- (obs_vals[[i]] - scale_params[[i]]$mean) / scale_params[[i]]$sd
      }
    }
    
    S_obs <- sapply(models, function(m) {
      if(m$type %in% c("GLM", "GAM")) {
        # For GLM with quadratic terms, square variables if needed
        pred_vals <- obs_vals
        if(m$type=="GLM" & any(grepl("\\^2", names(m$model$coefficients)))) {
          pred_vals2 <- pred_vals^2
          names(pred_vals2) <- paste0(names(pred_vals), "^2")
          pred_vals <- c(pred_vals, pred_vals2)
        }
        predict(m$model, newdata=as.data.frame(t(pred_vals)), type="response")
      } else if(m$type %in% c("RF","BRF")) {
        predict(m$model, newdata=as.data.frame(t(obs_vals)))
      } else stop("Unknown model type")
    })
    S_obs <- mean(S_obs) # ensemble average
    
    # Compute limiting variable
    L_i <- numeric(n_vars)
    for(i in 1:n_vars) {
      # For this variable, replace by sequence along range, others fixed
      pred_grid <- lapply(var_ranges[[i]], function(val) {
        new_vals <- obs_vals
        new_vals[[i]] <- val
        S_new <- sapply(models, function(m) {
          if(m$type %in% c("GLM", "GAM")) {
            pred_vals <- new_vals
            if(m$type=="GLM" & any(grepl("\\^2", names(m$model$coefficients)))) {
              pred_vals2 <- pred_vals^2
              names(pred_vals2) <- paste0(names(pred_vals), "^2")
              pred_vals <- c(pred_vals, pred_vals2)
            }
            predict(m$model, newdata=as.data.frame(t(pred_vals)), type="response")
          } else if(m$type %in% c("RF","BRF")) {
            predict(m$model, newdata=as.data.frame(t(new_vals)))
          }
        })
        mean(S_new)
      })
      L_i[i] <- max(unlist(pred_grid)) - S_obs
    }
    
    limiting_var_idx[cell] <- which.max(L_i)
  }
  
  # Create limiting variable raster
  limiting_raster <- raster_stack[[1]]
  limiting_raster[] <- limiting_var_idx
  return(limiting_raster)
}






vars_selected<- c( "temperature_summer_mean",  "slope_5x5" ,  "turbidity_summer_mean",
                   "ammonium_spring_SD", "PAR_summer_mean", "currentSpeed_summer_mean", "salinity_summer_SD")      

train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_FINAL.csv")
train_test_dataset<- (train_test_dataset[,-1])
names(train_test_dataset)

train<- train_test_dataset %>%
  filter(set == "train")

test<- train_test_dataset %>%
  filter(set == "test")

train_sel <- train %>% select(all_of(c("kelp", vars_selected)))
test_sel <- test %>% select(all_of(c("kelp", vars_selected)))

train_sel$kelp<- as.factor(train_sel$kelp)
test_sel$kelp<- as.factor(test_sel$kelp)

scaling_params_2 <- train_sel %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd), na.rm = TRUE))


