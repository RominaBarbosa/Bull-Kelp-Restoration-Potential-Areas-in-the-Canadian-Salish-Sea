

library(raster)
library(mgcv)
library(randomForest)

#' Compute limiting variable per cell following Elith et al., 2010
#'
#' @param models List of models. Each model should be a named list: list(model=..., type="GLM/GAM/RF/BRF")
#' @param raster_stack RasterStack or RasterBrick of predictor variables
#' @param scale_params Optional: list of scale parameters (mean, sd) for each variable, needed for GLM/GAM
#' @param n_samples Optional: number of raster cells to be randomly sampled for computation
#' @param buffer Optional: minimum distance (in cells) between sampled cells
#' @param var_range_points Number of points to sample along each variable range (default 50)
#' @return RasterLayer with limiting variable index (1 = first variable, etc.)
#' 
#' If you set a buffer value, the function tries to enforce minimum distance between sampled cells (approximate, in cell units).
 # That way you get more spatially spread-out samples rather than purely random clustering.

n_samples= 6

compute_limiting_variable <- function(models, raster_stack, scale_params=NULL, 
                                      n_samples=NULL, buffer=NULL, var_range_points=50, model_type= "GLM") {
  
  n_cells <- terra::ncell(raster_stack)
  n_vars <- terra::nlyr(raster_stack)
  var_names <- names(raster_stack)
  
  # --- Subsample cells if requested ---
  sample_idx <- 1:n_cells
  if(!is.null(n_samples)) {
    if(!is.null(buffer)) {
      sampled <- c()
      remaining <- 1:n_cells
      while(length(sampled) < n_samples & length(remaining) > 0) {
        idx <- remaining[1]
        sampled <- c(sampled, idx)
        remaining <- remaining[abs(remaining - idx) > buffer]
      }
      sample_idx <- sampled
    } else {
      valid_cells <- which(!is.na(terra::values(raster_stack[[1]])))
      sample_idx <- sample(valid_cells, n_samples)
    }
  }
  
  # --- Precompute variable ranges ---
  # var_ranges <- lapply(1:n_vars, function(i) {
  #   vals <- terra::values(raster_stack[[i]], mat = FALSE)
  #   seq(min(vals, na.rm=TRUE), max(vals, na.rm=TRUE), length.out=var_range_points)
  # })
  # names(var_ranges) <- var_names
  # --- Precompute variable ranges using min and max (faster) ---
  var_ranges <- lapply(1:n_vars, function(i) {
    rng <- terra::minmax(raster_stack[[i]])
    from <- as.numeric(rng["min", 1])
    to   <- as.numeric(rng["max", 1])
    seq(from, to, length.out = var_range_points)
  })
  names(var_ranges) <- var_names
  
  
  
  # Initialize output for sampled cells
  limiting_var_idx <- rep(NA, n_cells)
  limiting_var <- rep(NA, n_cells)
  
  # Function to get values for a single cell
  get_cell_values <- function(cell_idx) {
    terra::extract(raster_stack, cell_idx)[1, ]  # single-row data.frame
  }
  
  # --- Loop over sampled cells ---
  for(cell in sample_idx) {
    obs_vals <- get_cell_values(cell)
    
    # Predict full suitability
    if(model_type %in% c("GLM", "GAM")) {
      # Scale if needed
      if(!is.null(scale_params)) {
        for(i in 1:n_vars) {
          var <- var_names[i]
          obs_vals[[var]] <- (obs_vals[[var]] - scale_list[[var]]$mean) / scale_list[[var]]$sd
        }
      }
      S_obs <- predict(model, newdata = obs_vals, type = "response")
    } else if(model_type %in% c("RF","BRF")) {
      S_obs <- predict(model, newdata= obs_vals)
    } else stop("Unknown model type")
    
    # Compute limiting variable per predictor
    L_i <- numeric(n_vars)
    
    for(i in 1:n_vars) {
      pred_grid <- sapply(var_ranges[[i]], function(val) {
        new_vals <- obs_vals
        
        # Scale the variable if needed (GLM/GAM)
        if(!is.null(scale_params) && model_type %in% c("GLM", "GAM")) {
          mean_val <- scale_list[[var_names[i]]]$mean
          sd_val   <- scale_list[[var_names[i]]]$sd
          val <- (val - mean_val) / sd_val
        }
        
        new_vals[[i]] <- val  # replace the variable with the current grid value
        
        # Predict using the appropriate model
        if(model_type %in% c("GLM", "GAM")) {
          predict(model, newdata = new_vals, type = "response")
        } else if(model_type %in% c("RF", "BRF")) {
          predict(model, newdata = new_vals)
        } else stop("Unknown model type")
      })  # end sapply
      
      # Difference between maximum suitability along the gradient and observed suitability
      L_i[i] <- max(pred_grid, na.rm = TRUE) - S_obs
    }
    
    # Store results
    limiting_var_idx[cell] <- which.max(L_i)           # index of limiting variable
    limiting_var[cell] <- var_names[which.max(L_i)]   # name of limiting variable
  }
  
  # Create output SpatRaster
  limiting_raster <- raster_stack[[1]]
  # limiting_raster2 <- raster_stack[[1]]
  terra::values(limiting_raster) <- limiting_var_idx
  names(limiting_raster)<- "var_index_limiting"
  # terra::values(limiting_raster2) <- limiting_var
  # limiting_results<- c(limiting_raster, limiting_raster2)
  # names(limiting_results)<- c("index_limiting", "variable_limiting")
  
  # Extract coordinates and raster values as a data.frame
  df <- terra::as.data.frame(limiting_raster, xy = TRUE, na.rm = TRUE)
  # df$limiting_variable<- as.factor(df[,3])
  
  return(list(limiting_raster, df))
}





## Load models  =====
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
setwd(model_results_path)
dir()

glm_mod_s<- readRDS("glm_mod_s.rds")
gam_mod_s<- readRDS("gam_mod_s.rds")
rf_mod_s<-  readRDS("rf_mod_s.rds")
brt_mod_s<- readRDS("brt_mod_s.rds")

models <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")


## Load raster stack  =====
raster_stack_predict<- readRDS("raster_stack_predict.rds")
raster_stack_r <- raster::stack(raster_stack_predict)

raster_stack_predict_scaled<- readRDS("raster_stack_predict_scaled.rds")

## Load the training dataset and define scaling parameters =====
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

# Convert into list as required for the function
var_names <- unique(gsub("_mean$|_sd$", "", names(scaling_params_2)))

# Build the list
scale_list <- list()

for(v in var_names) {
  mean_col <- paste0(v, "_mean")
  sd_col <- paste0(v, "_sd")
  
  scale_list[[v]] <- list(
    mean = scaling_params_2[[mean_col]][1],
    sd   = scaling_params_2[[sd_col]][1]
  )
}


# Your raster layer names
var_names <- names(obs_vals)

# Build scale list
scale_list <- list()
for(v in var_names) {
  # Find matching columns in scale_params
  mean_col <- grep(paste0("^", v, "_.*mean$"), names(scale_params), value = TRUE)
  sd_col   <- grep(paste0("^", v, "_.*sd$"), names(scale_params), value = TRUE)
  
  if(length(mean_col)==0 | length(sd_col)==0) stop(paste("Cannot find scale columns for", v))
  
  scale_list[[v]] <- list(
    mean = scale_params[[mean_col]][1],
    sd   = scale_params[[sd_col]][1]
  )
}


### Compute variables limiting =================================================
models_types<- c("GLM", "GAM", "RF", "BRT")
class(models[[3]])
  
for (i in 2:length(models)) {
  
  model_type= models_types[i]
  model <- models[[i]]
  
  glm_limiting<- compute_limiting_variable (model_type= models_types[i], model= models[[i]], raster_stack=raster_stack_predict, scale_params=scaling_params_2, 
                                            n_samples=10000, buffer=NULL, var_range_points=50)
  
  
  # Variable names in same order as your raster stack
  var_names <- names(raster_stack)
  
  # Example mapping
  index_to_name <- setNames(var_names, 1:length(var_names))
  
  # Apply mapping to your dataframe
  rast_df <- glm_limiting[[2]]
  rast_df$var_name <- index_to_name[as.character(rast_df$var_limiting)]
  rast_df$var_name<- as.factor(rast_df$var_name)
  # summary(rast_df)
  
  write.csv(rast_df, paste("limiting",model_type,".csv", sep="_"))
  
}

end_time<- Sys.time()








# Example: categorical fill for limiting variable names
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Download high-resolution land polygons
land_highres <- ne_download(
  scale = 10,            # 10 = high resolution (1:10m)
  type = "land",
  category = "physical",
  returnclass = "sf"
)

# Optional: crop to your study extent
study_extent <- c(xmin = -126, xmax = -122, ymin = 47, ymax = 51)  # example
land_highres <- st_crop(land_highres,
                        xmin = study_extent["xmin"],
                        xmax = study_extent["xmax"],
                        ymin = study_extent["ymin"],
                        ymax = study_extent["ymax"])

coastline <- ne_download(scale = 10, type = "coastline",
                         category = "physical", returnclass = "sf")

land_highres <- st_transform(land_highres, crs = st_crs(raster_stack))
coastline <- st_transform(coastline, crs = st_crs(raster_stack))


ggplot() +
  geom_sf(data = land_highres, fill = "grey85", color = NA) +
  geom_tile(data = rast_df, aes(x = x, y = y, fill = as.factor(var_limiting))) +
  geom_sf(data = coastline, color = "grey40", size = 0.3) +
  coord_sf(xlim = c(study_extent["xmin"], study_extent["xmax"]),
           ylim = c(study_extent["ymin"], study_extent["ymax"]),
           expand = FALSE) +
  scale_x_continuous(breaks = seq(-126, -122, by = 1)) +
  scale_y_continuous(breaks = seq(47, 51, by = 1)) +
  scale_fill_viridis_d(name = "Limiting variable", na.value = "transparent") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        legend.position = c(0.86, 0.90)) +
  labs(title = "Limiting Variable Map",
       fill = "Limiting Variable",
       x = "",
       y = "")


