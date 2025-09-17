
# Loops over all predictors.
# Creates one-variable-at-a-time grids.
# Scale grids for GLM/GAM using the training mean and SD.
# Use the scaled grid for GLM/GAM predictions.
# Use raw grid for RF/BRT predictions.
# Combine all into a long-format dataset ready for ggplot.



library(ggplot2)
library(pdp)
library(visreg)
library(mgcv)    # for GAMs
library(gbm)     # for BRT

#-----------------------------
# 0. Define scaled variables
#-----------------------------
# vars <- attr(terms(glm_mod_quad), "term.labels")
# vars <- vars[!grepl("^I\\(", vars)]   # remove quadratic terms

# Save scaling info (mean and SD) used for GLM/GAM
# scale_params <- list(
#   mean = sapply(train[, vars], mean),
#   sd   = sapply(train[, vars], sd)
# )

#-----------------------------
# 1. Function to create one-variable-at-a-time grid
#-----------------------------
create_var_grid <- function(var, data, n_points = 50) {
  newdata <- data.frame(matrix(ncol = ncol(data), nrow = n_points))
  names(newdata) <- names(data)
  
  for (v in names(data)) {
    if (v == var) {
      if (is.numeric(data[[v]])) {
        newdata[[v]] <- seq(min(data[[v]], na.rm = TRUE),
                            max(data[[v]], na.rm = TRUE),
                            length.out = n_points)
      } else {
        newdata[[v]] <- levels(data[[v]])[1]  # reference level for factors
      }
    } else {
      if (is.numeric(data[[v]])) {
        newdata[[v]] <- mean(data[[v]], na.rm = TRUE)
      } else {
        newdata[[v]] <- levels(data[[v]])[1]
      }
    }
  }
  return(newdata)
}

#-----------------------------
# 2. Function to scale grid for GLM/GAM
#-----------------------------
scale_grid <- function(grid, scale_params) {
  for (v in names(grid)) {
    mean_name <- paste0(v, "_mean")
    sd_name   <- paste0(v, "_sd")
    if (mean_name %in% names(scale_params) && sd_name %in% names(scale_params)) {
      grid[[v]] <- (grid[[v]] - scale_params[[mean_name]]) / scale_params[[sd_name]]
    }
  }
  return(grid)
}

#-----------------------------
# 3. Function to get predicted curve for one model & variable
#-----------------------------
# grid <- create_var_grid("temperature_summer_mean", train)
# grid$kelp <- NULL  # remove response variable

get_curve <- function(model=glm_mod_quad, var= "temperature_summer_mean", grid, model_name="GAM", scale_params = NULL) {
  cls <- class(model)[1]
  
  if (cls %in% c("glm", "gam")) {
    if (is.null(scale_params)) stop("GLM/GAM models require scale_params.")
    
    # Extract original x for plotting
    x_unscaled <- grid[[var]]
    
    # Scale grid for prediction
    grid_scaled <- as.data.frame(scale_grid(grid, scale_params))
    
    # Predict
    pred <- predict(model, newdata = grid_scaled, type = "response")
    
    # Build df
    df <- data.frame(
      x = x_unscaled,
      fit = pred,
      model = model_name,
      var = var
    )
    
  } else if (cls %in% c("randomForest", "randomForest.formula")) {
    
    pd <- pdp::partial(
      object = model,
      pred.var = var,
      train = grid,
      prob = TRUE
    )
    df <- data.frame(
      x = pd[[var]],
      fit = 1 - pd$yhat,
      model = model_name,
      var = var
    )
    
  } else if (cls == "gbm") {
    
    pd <- pdp::partial(
      object = model,
      pred.var = var,
      train = grid,
      pred.fun = function(object, newdata) {
        predict(object, newdata, n.trees = model$n.trees, type = "response")
      },
      recursive = FALSE
    )
    
    df <- data.frame(
      x = pd[[var]],
      fit = pd$yhat,
      model = model_name,
      var = var
    )
    
  } else {
    stop("Unsupported model class")
  }
  
  return(df)
}



#-----------------------------
# 4. Loop over variables & models, collect curves
#-----------------------------
# all_curves <- do.call(rbind, lapply(vars, function(v) {
#   grid <- create_var_grid(v, train)   # one-variable-at-a-time grid
#   rbind(
#     get_curve(glm_mod_quad, v, grid, "GLM", scale_params),
#     get_curve(gam_mod,      v, grid, "GAM", scale_params),
#     get_curve(rf_mod,       v, grid, "randomForest"),
#     get_curve(brt_mod,      v, grid, "BRT")
#   )
# }))
# 
# #-----------------------------
# # 5. Plot all curves in facets
# #-----------------------------
# ggplot(all_curves, aes(x = x, y = fit, color = model)) +
#   geom_line(size = 1.2) +
#   facet_wrap(~var, scales = "free_x") +
#   theme_minimal() +
#   labs(
#     y = "Predicted response",
#     x = NULL,
#     color = "Model",
#     title = "Predicted response curves for all variables across models"
#   )





