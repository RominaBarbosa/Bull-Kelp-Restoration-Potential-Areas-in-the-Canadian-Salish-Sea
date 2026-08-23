
library(terra)
library(dplyr)
library(dismo)
library(ggplot2)
library(raster)

# ------------------------------------------------------------------------------
#  MESS analysis - evaluation of areas with environmental conditions 
#   outside the range of calibration data
# ------------------------------------------------------------------------------
# Stack predictors for each period
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres"
terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")
output_path<- "/Volumes/Romina_PSF/PSF/SDM/MESS"

# This chunk was run and its outputs are stored in a RDS file
## Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m_t1<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_blob.tif")
raster_stack_20m_t2<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblob.tif")


# load terrain variables
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)
terrain_vars<-  rast(tif_files[c(12)]) #,5,12,13,15

# Merge rasters of all selected variables including terrain and NEMO
raster_stack_20m_t1<- c(raster_stack_20m_t1, terrain_vars)
raster_stack_20m_t2<- c(raster_stack_20m_t2, terrain_vars)


# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask10_40<- bathy20m
bathy20m_mask10_40[!(bathy20m_mask10_40 >= -10 & bathy20m_mask10_40 <= 30)] <- NA # negative values are in land

# Mask all layers
glm_mod_s<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)



### Subset selected variables
raster_stack_20m_t1<- raster_stack_20m_t1[[vars_selected]] 
raster_stack_20m_t1<- mask(raster_stack_20m_t1, bathy20m_mask10_40)

raster_stack_20m_t2<- raster_stack_20m_t2[[vars_selected]] 
raster_stack_20m_t2<- mask(raster_stack_20m_t2, bathy20m_mask10_40)

names(raster_stack_20m_t2)
names(raster_stack_20m_t1)

plot(raster_stack_20m_t1[[1]])
plot(raster_stack_20m_t2[[1]])


env_t1 <- raster::stack(raster_stack_20m_t1)  # calibration (training period)
env_t2 <- raster::stack(raster_stack_20m_t2)  # projection (new period)

# Sample calibration conditions from t1
# (faster than using full raster if very large)
# Load training and testing dataset 
traintest<- read.csv(paste("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_M7.csv", sep="/"))
traintest<- traintest[,-1]
traintest<- traintest[,c("substrate",  "depth",   "bathymetry_20",  "slope_5x5", "TPI_3x3", "x",  "y", "set")]

train<- traintest%>%
  filter(set=="train")

calib_points<- vect(train, c("x","y"))
# Extract predictor values at training locations
calib_points <- extract(raster_stack_20m_t1, calib_points)
calib_points <- calib_points[ , -1]   # remove ID column


start=Sys.time()
mess_map <- dismo::mess(env_t1, calib_points)
end=Sys.time() # Time difference of 1.433886 hours

plot(mess_map)
end-start

# Area outside range
out_range_P2<- mess_map
out_range_P2[out_range_P2>= 0]<- NA
plot(out_range_P2)

start=Sys.time()
mess_map_t2 <- dismo::mess(env_t2, calib_points)
end=Sys.time()

writeRaster(mess_map, "mess_map_2014_2019_M7.tif")
# writeRaster(mess_map_t2, "/Volumes/Romina_PSF/PSF/SDM/MESS/mess_map_2020_2022.tif")

print(end-start)

library(terra)
library(future)
library(future.apply)
library(ggplot2)
library(dplyr)

run_mess_analysis <- function(env_raster, calib_points, n_sample = 70000, seed = 42) {
  
  set.seed(seed)
  
  # -----------------------------
  # 1. Sample points
  # -----------------------------
  samp <- spatSample(env_raster, size = n_sample, method = "regular",
                     na.rm = TRUE, as.df = TRUE, xy = TRUE)
  
  coords <- samp[, c("x", "y")]
  vals   <- samp[, !(names(samp) %in% c("x", "y"))]
  
  # -----------------------------
  # 2. Training ranges
  # -----------------------------
  mins <- apply(calib_points, 2, min, na.rm = TRUE)
  maxs <- apply(calib_points, 2, max, na.rm = TRUE)
  
  # -----------------------------
  # 3. Define MESS function
  # -----------------------------
  calc_mess <- function(vals, mins, maxs) {
    scores <- numeric(length(vals))
    for (i in seq_along(vals)) {
      v <- vals[i]
      if (is.na(v)) return(NA_real_)
      if (v >= mins[i] && v <= maxs[i]) {
        scores[i] <- min((v - mins[i])/(maxs[i] - mins[i]),
                         (maxs[i] - v)/(maxs[i] - mins[i])) * 100
      } else {
        scores[i] <- ifelse(v < mins[i],
                            -100 * (mins[i] - v)/(maxs[i] - mins[i]),
                            -100 * (v - maxs[i])/(maxs[i] - mins[i]))
      }
    }
    list(mess = min(scores, na.rm = TRUE),
         driver = which.min(scores))
  }
  
  # -----------------------------
  # 4. Parallel MESS calculation
  # -----------------------------
  plan(multisession, workers = availableCores() - 1)
  mess_results <- future_apply(vals, 1, calc_mess, mins = mins, maxs = maxs, future.seed = TRUE)
  plan(sequential)
  
  # -----------------------------
  # 5. Assemble data frame
  # -----------------------------
  mess_df <- data.frame(
    coords,
    mess   = sapply(mess_results, `[[`, "mess"),
    driver = sapply(mess_results, `[[`, "driver")
  )
  
  # -----------------------------
  # 6. Binary inside/outside
  # -----------------------------
  mess_df$inside_range <- factor(ifelse(mess_df$mess >= 0, "Inside range", "Outside range"),
                                 levels = c("Outside range", "Inside range"))
  
  # -----------------------------
  # 7. Driver factor for points outside range
  # -----------------------------
  outside_df <- mess_df %>% filter(mess < 0)
  unique_drivers <- sort(unique(outside_df$driver))
  
  # Assign labels only to those
  outside_df$driver_factor <- factor(outside_df$driver,
                                     levels = unique_drivers,
                                     labels = names(mins)[unique_drivers])
  
  # -----------------------------
  # 8. Count variables outside range
  # -----------------------------
  mess_df$num_vars_outside <- apply(vals, 1, function(x) sum(x < mins | x > maxs))
  
  # -----------------------------
  # 9. Generate plots
  # -----------------------------
  plot_mess_values <- ggplot(mess_df, aes(x = x, y = y, color = mess)) +
    geom_point(size = 0.8) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(color = "MESS score", title = "MESS values (sampled points)")
  
  plot_inside_outside <- ggplot(mess_df, aes(x = x, y = y, color = inside_range)) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c("red", "green")) +
    theme_minimal() +
    labs(color = "Environmental novelty", title = "Inside vs Outside training range")
  
  plot_driver <- ggplot(outside_df, aes(x = x, y = y, color = driver_factor)) +
    geom_point(size = 0.8, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Driver of novelty (outside training range)", color = "Limiting variable") +
    scale_color_viridis_d(option = "turbo")
  
  plot_num_outside <- mess_df %>% 
    filter(num_vars_outside >= 1) %>%
    ggplot(aes(x = x, y = y, color = as.factor(num_vars_outside))) +
    geom_point(size = 0.8) +
    theme_minimal() +
    labs(title = "Number of variables outside training range", color = "Variables outside range")
  
  # -----------------------------
  # 10. Return results
  # -----------------------------
  return(list(
    mess_df = mess_df,
    outside_df = outside_df,
    plots = list(
      mess_values = plot_mess_values,
      inside_outside = plot_inside_outside,
      driver = plot_driver,
      num_vars_outside = plot_num_outside
    )
  ))
}




################################################################################
# Same for t2 - postblob conditions 
results <- run_mess_analysis(env_t1, calib_points, n_sample = 80000)
results_t2 <- run_mess_analysis(env_t2, calib_points, n_sample = 80000)

# Access data
head(results$mess_df)
head(results$outside_df)

# Show plots individually
results$plots$mess_values
results$plots$inside_outside
results$plots$driver
results$plots$num_vars_outside

cowplot::plot_grid(results$plots$mess_values,
                   results$plots$inside_outside,
                   results$plots$driver,
                   results$plots$num_vars_outside)

cowplot::plot_grid(results_t2$plots$mess_values,
                   results_t2$plots$inside_outside,
                   results_t2$plots$driver,
                   results_t2$plots$num_vars_outside)

cowplot::plot_grid(results$plots$inside_outside,
                   results_t2$plots$inside_outside)

cowplot::plot_grid(results$plots$driver,
                   results_t2$plots$driver)
                   
cowplot::plot_grid(results$plots$mess_values,
                   results_t2$plots$mess_values)







library(terra)
library(future)
library(future.apply)

# ------------------------------
# Optimized full-raster MESS function
# ------------------------------
# full_raster_inside_outside <- function(env_raster, calib_points, cores = 4,
#                                        filename = NULL, mask_na = TRUE) {
  # ------------------------------
  # 1. Training ranges
  # ------------------------------
  mins <- apply(calib_points, 2, min, na.rm = TRUE)
  maxs <- apply(calib_points, 2, max, na.rm = TRUE)
  
  # ------------------------------
  # 2. Optional mask for fully NA pixels
  # ------------------------------
  if (mask_na) {
    valid_mask <- !is.na(env_raster[[1]])
    if (nlyr(env_raster) > 1) {
      for (i in 2:nlyr(env_raster)) valid_mask <- valid_mask & !is.na(env_raster[[i]])
    }
    env_raster <- mask(env_raster, valid_mask)
  }
  
  # ------------------------------
  # 3. Pixel-wise inside/outside function
  # ------------------------------
  calc_inside_outside <- function(...) {
    vals <- c(...)
    if (any(is.na(vals))) return(NA_real_)
    if (all(vals >= mins & vals <= maxs)) return(1)  # inside range
    return(0)                                      # outside range
  }
  
  # ------------------------------
  # 4. Apply to raster stack
  # ------------------------------
  inside_outside_rast <- app(env_raster,
                             fun = calc_inside_outside,
                             cores = cores,
                             filename = filename,
                             overwrite = TRUE)
  names(inside_outside_rast) <- "inside_outside"
  
  return(inside_outside_rast)
}






# inside_rast <- full_raster_inside_outside(env_raster = env_t1,
#                                           calib_points = calib_points,
#                                           cores = 4,
#                                           filename = "MESS_inside_outside_2014_2019.tif")
# 
# # Plot
# plot(inside_rast, col = c("red", "green"),
#      legend = TRUE,
#      main = "Inside (green) vs Outside (red) training range")
# 
# 
# # POst Blob period
# inside_rast_t2 <- full_raster_inside_outside(env_raster = env_t2,
#                                           calib_points = calib_points,
#                                           cores = 4,
#                                           filename = "MESS_inside_outside_2020_2022.tif")
# 
# # Plot
# plot(inside_rast_t2, col = c("red", "green"),
#      legend = TRUE,
#      main = "Inside (green) vs Outside (red) training range")
# 
# 
# 
# 
# 
# 
# 
# 

# 
# 
# 
# 
# 
# 
# 
# 
# # Function for MESS analysis with terra (optimize memory) to get also variable responsable of change
# terra_mess_extended <- function(env, mins, maxs, cores = 4,
#                                 filename_mess = NULL,
#                                 filename_novelty = NULL,
#                                 filename_driver = NULL) {
#   stopifnot(nlyr(env) == length(mins))
#   
#   varnames <- names(env)
#   
#   # ---- 1) Create a mask of all non-NA cells ----
#   valid_mask <- !is.na(env[[1]])
#   for (i in 2:nlyr(env)) {
#     valid_mask <- valid_mask & !is.na(env[[i]])
#   }
#   
#   # ---- 2) Define MESS function (embed mins/maxs) ----
#   mess_fun <- function(..., mins_local = mins, maxs_local = maxs) {
#     vals <- c(...)
#     if (any(is.na(vals))) return(NA_real_)
#     
#     scores <- rep(NA_real_, length(vals))
#     for (i in seq_along(vals)) {
#       v <- vals[i]
#       if (v >= mins_local[i] && v <= maxs_local[i]) {
#         scores[i] <- min((v - mins_local[i]) / (maxs_local[i] - mins_local[i]),
#                          (maxs_local[i] - v) / (maxs_local[i] - mins_local[i])) * 100
#       } else {
#         scores[i] <- ifelse(v < mins_local[i],
#                             -100 * (mins_local[i] - v)/(maxs_local[i] - mins_local[i]),
#                             -100 * (v - maxs_local[i])/(maxs_local[i] - mins_local[i]))
#       }
#     }
#     return(min(scores, na.rm = TRUE))
#   }
#   
#   mess_rast <- terra::app(env, fun = mess_fun, cores = cores,
#                           filename = filename_mess, overwrite = TRUE)
#   
#   # # ---- 3) Novelty maps per variable ----
#   # novelty_list <- list()
#   # for (i in seq_along(varnames)) {
#   #   novelty_list[[i]] <- terra::app(env[[i]],
#   #                                   fun = function(v, minv = mins[i], maxv = maxs[i]) {
#   #                                     if (is.na(v)) return(NA)
#   #                                     as.integer(v < minv | v > maxv)
#   #                                   },
#   #                                   cores = cores,
#   #                                   filename = if (!is.null(filename_novelty)) paste0(filename_novelty, "_", varnames[i], ".tif") else NULL,
#   #                                   overwrite = TRUE
#   #   )
#   #   names(novelty_list[[i]]) <- paste0("novel_", varnames[i])
#   # }
#   # novelty_stack <- rast(novelty_list)
#   # 
#   # # ---- 4) Variable driving novelty (embed mins/maxs) ----
#   # which_fun <- function(..., mins_local = mins, maxs_local = maxs) {
#   #   vals <- c(...)
#   #   if (any(is.na(vals))) return(NA_integer_)
#   #   
#   #   scores <- numeric(length(vals))
#   #   for (i in seq_along(vals)) {
#   #     v <- vals[i]
#   #     if (v >= mins_local[i] && v <= maxs_local[i]) {
#   #       scores[i] <- min((v - mins_local[i]) / (maxs_local[i] - mins_local[i]),
#   #                        (maxs_local[i] - v) / (maxs_local[i] - mins_local[i])) * 100
#   #     } else {
#   #       scores[i] <- ifelse(v < mins_local[i],
#   #                           -100 * (mins_local[i] - v)/(maxs_local[i] - mins_local[i]),
#   #                           -100 * (v - maxs_local[i])/(maxs_local[i] - mins_local[i]))
#   #     }
#   #   }
#   #   return(which.min(scores))
#   # }
#   # 
#   # driver_rast <- terra::app(env, fun = which_fun, cores = cores,
#   #                           filename = filename_driver, overwrite = TRUE)
#   # names(driver_rast) <- "novelty_driver"
#   
#   return(list(
#     mess = mess_rast
#     # novelty = novelty_stack,
#     # driver = driver_rast,
#     # valid_mask = valid_mask
#   ))
# }
# 
# 
# 
# parallel::detectCores()
# ncell(env_t2)         # total number of pixels (nrow * ncol * layers)
# ncell(env_t2) * 8 / 1e6  # approximate memory in MB for one layer (8 bytes per double)
# ncell(env_t2) * nlyr(env_t2) # if it is very large (>100 million), use 4 cores
# 
# 
# # Training ranges
# mins <- apply(calib_points, 2, min, na.rm = TRUE)
# maxs <- apply(calib_points, 2, max, na.rm = TRUE)
# names(env_t1)
# names(mins)
# names(maxs)
# 
# 
# # Run MESS on t1 and t2
# setwd(output_path)
# init= Sys.time()
# results_t1 <- terra_mess_extended(env_t1, mins, maxs, cores = 4,
#                                   filename_mess = "mess_map_t1.tif",
#                                   filename_novelty = "novelty_var_t1",
#                                   filename_driver = "driver_map_t1.tif")
# fin= Sys.time()
# fin-init
# 
# results_t2 <- terra_mess_extended(env_t2, mins, maxs, cores = 4,
#                                   filename_mess = "mess_map_t2.tif",
#                                   filename_novelty = "novelty_var_t2",
#                                   filename_driver = "driver_map_t2.tif")
# 
# # Access outputs
# mess_map_t1     <- results_t1$mess
# novelty_maps_t1 <- results_t1$novelty
# driver_map_t1   <- results_t1$driver
# 
# mess_map_t2     <- results_t2$mess
# novelty_maps_t2 <- results_t2$novelty
# driver_map_t2   <- results_t2$driver
# 
# 
# # Quick plots
# plot(mess_map_t1, main = "MESS index")
# plot(novelty_maps_t1, main = "Per-variable novelty (0=inside, 1=outside)")
# plot(driver_map_t1, main = "Variable driving novelty")
# 
# plot(mess_map_t2, main = "MESS index")
# plot(novelty_maps_t2, main = "Per-variable novelty (0=inside, 1=outside)")
# plot(driver_map_t2, main = "Variable driving novelty")
# 
# # Save results
# writeRaster(mess_map_t1,     "mess_index_t1.tif", overwrite = TRUE)
# writeRaster(novelty_maps_t1, "novelty_per_variable_t1.tif", overwrite = TRUE)
# writeRaster(driver_map_t1,   "novelty_driver_t1.tif", overwrite = TRUE)
# 
# writeRaster(mess_map_t2,     "mess_index_t2.tif", overwrite = TRUE)
# writeRaster(novelty_maps_t2, "novelty_per_variable_t2.tif", overwrite = TRUE)
# writeRaster(driver_map_t2,   "novelty_driver_t2.tif", overwrite = TRUE)
# 
# 
# 
# 
# 
# 
# library(reshape2)
# library(class)
# 
# # Step 1. PCA of training points
# # PCA on calibration values
# calib_pca <- prcomp(calib_points, center = TRUE, scale. = TRUE)
# 
# # Keep first few axes (say 5)
# pca_scores <- calib_pca$x[, 1:5]
# 
# # Step 2. K-means clustering
# set.seed(123)  # reproducible
# clusters <- kmeans(pca_scores, centers = 5, nstart = 20)
# 
# # Assign cluster ID to each training point
# training_pts$cluster <- clusters$cluster
# 
# # Step 3. Predict clusters across study area
# # We need to map those environmental clusters across all pixels of env_t2. A good trick:
# #   Project raster values into PCA space.
# # Predict cluster membership with stats::predict() + class::knn() or nearest-centroid.
# 
# # Extract raster values (in chunks)
# vals_t2 <- terra::values(env_t2, na.rm = TRUE)
# 
# # Project into PCA space (use same scaling/rotation as training PCA)
# pca_vals <- predict(calib_pca, newdata = vals_t2)[, 1:5]
# 
# # Predict cluster for each pixel using kNN with training centroids
# centroids <- clusters$centers
# # Use nearest centroid assignment
# cluster_id <- apply(pca_vals, 1, function(x) {
#   which.min(colSums((t(centroids) - x)^2))
# })
# 
# # Rebuild as raster
# cluster_rast <- env_t2[[1]]
# values(cluster_rast) <- NA
# values(cluster_rast)[!is.na(vals_t2[,1])] <- cluster_id
# names(cluster_rast) <- "env_cluster"
# 
# plot(cluster_rast, main="Environmental clusters (k=5)")
# 
# # Step 4. Summarize novelty by cluster
# # Suppose you already ran the extended MESS analysis and have novelty_maps (stack of binary rasters).
# # Add cluster map as mask
# novelty_with_clusters <- c(novelty_maps, cluster_rast)
# 
# # Extract summary
# df <- as.data.frame(terra::extract(novelty_with_clusters, cluster_rast, na.rm=TRUE))
# # df has one row per pixel: [novel_var1, novel_var2, ..., env_cluster]
# 
# novelty_summary <- df %>%
#   group_by(env_cluster) %>%
#   summarise(across(starts_with("novel_"), mean, na.rm = TRUE))
# 
# print(novelty_summary)
# # This gives the fraction of pixels novel per variable per cluster.
# 
# # Step 5. Visualize
# df_long <- melt(novelty_summary, id.vars = "env_cluster")
# 
# ggplot(df_long, aes(x = variable, y = value, fill = factor(env_cluster))) +
#   geom_bar(stat="identity", position="dodge") +
#   ylab("Fraction novel") +
#   xlab("Variable") +
#   theme_minimal()
# 


