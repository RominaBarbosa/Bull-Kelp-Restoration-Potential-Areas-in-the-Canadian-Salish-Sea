##==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
### 1.2.4- Derive terrain features from bathymetry (20m res.)   ################
### Author: Romina Barbosa                                      ################
### Date last edition: 5-August2025                             ################
###==================================================================
# Load libraries

library(spatialEco)
library(raster)
library(RSAGA)
library(terra)

my_path_bathy<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres"
setwd(my_path_bathy)  
dir()
output_dir<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables"

# Load bathymetry raster (20 m resolution)
bathy <- rast("coastwide_20m.tif")  # Replace with your file path
bathy_r <- raster(bathy)

# bathy_11x11 <- rast("coastwide_smooth_11x11.tif") 
# bathy_15x15 <- rast("coastwide_smooth_15x15.tif") 

# Define window sizes (cells)
window_sizes <- c(3, 5, 7) #60, 100m, 220m, and 300m

# ---- FUNCTIONS ----

# Slope
get_slope <- function(r, size, filename, suffix="") {
  # Smooth the input raster
  smooth <- focal(r, w = matrix(1, size, size), fun = mean, na.rm = TRUE)
  # Compute slope (absolute values)
  slope_r <- terrain(smooth, v = "slope", unit = "degrees", neighbors = 8)
  # Assign name
  names(slope_r) <- paste0("slope", suffix)
  # Save to disk
  writeRaster(slope_r, filename, overwrite = TRUE)
  
  return(slope_r)
}


# Roughness (standard deviation)
get_roughness <- function(r, size, filename) {
  focal(
    r,
    w = matrix(1, size, size),
    fun = sd,
    na.rm = TRUE,
    filename = filename,
    overwrite = TRUE
  )
}

# TPI
get_TPI <- function(r=bathy, size=15, filename= paste(output_dir,paste0("TPI", suffix, ".tif"), sep="/")) {
  # Convert SpatRaster to RasterLayer if needed
  # if (inherits(r, "SpatRaster")) {
  #   r <- raster::raster(r)
  # }
  
  tpi_raster <- spatialEco::tpi(r, scale = size)
  layer_name <- paste0("TPI_", size, "x", size)
  names(tpi_raster) <- layer_name
  tpi_raster<- raster(tpi_raster)
  writeRaster(tpi_raster, filename= filename, overwrite = TRUE)
  # return(tpi_raster)
}


# TRI (use terra instead of spatialEco to avoid focal errors)
get_TRI <- function(r, size, filename) {
  focal(r, w = matrix(1, size, size), fun = function(x) {
    center <- x[ceiling(length(x)/2)]
    mean(abs(x - center), na.rm = TRUE)
  }, filename = filename, overwrite = TRUE)
}


# Northerness and Easterness
get_orientation <- function(bathy, size, suffix, output_dir) {
  library(terra)
  
  # Step 1: Smooth bathymetry
  w <- matrix(1, size, size)
  smooth_bathy <- focal(bathy, w = w, fun = mean, na.rm = TRUE)
  
  # Step 2: Calculate aspect from smoothed bathymetry (in radians)
  aspect <- terrain(smooth_bathy, v = "aspect", unit = "radians", neighbors = 8)
  
  # Step 3: Compute orientation
  northerness <- cos(aspect)
  easterness  <- sin(aspect)
  
  # Step 4: Rename layers
  names(northerness) <- paste0("northerness", suffix)
  names(easterness)  <- paste0("easterness", suffix)
  
  # Step 5: Save rasters
  writeRaster(northerness, file.path(output_dir, paste0("northerness", suffix, ".tif")), overwrite = TRUE)
  writeRaster(easterness,  file.path(output_dir, paste0("easterness", suffix,  ".tif")), overwrite = TRUE)
  
  # Optionally return the layers
  return(list(northerness = northerness, easterness = easterness))
}




# ---- LOOP THROUGH WINDOW SIZES ----

# for (w in window_sizes) { # R crash with the loop, so I run each at the time manually

  w= 5
  suffix <- paste0("_", w, "x", w)
  
  # Calculate each layer
  get_slope(bathy, w, filename= paste(output_dir, paste0("slope", suffix, ".tif"), sep="/"))
  
  get_roughness(bathy, w, filename= paste(output_dir,paste0("roughness", suffix, ".tif"), sep="/")) 
  
  get_TRI(bathy, w, paste(output_dir,paste0("TRI", suffix, ".tif"), sep="/"))
  
  get_TPI(bathy, w, paste(output_dir,paste0("TPI", suffix, ".tif"), sep="/"))
  
  get_orientation(bathy, w, output_dir=output_dir)
 
  cat("Saved layers for window size", w, "\n")
# }




# Final raster stack
terrain_stack <- stack(terrain_layers)
names(terrain_stack) <- names(terrain_layers)
terrain_stack<- rast(terrain_stack)

terrain_stack0<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/terrain_variables.tif")
terrain_stack1<- c(terrain_stack, terrain_stack0)
names(terrain_stack1)

# Save updated stack
writeRaster(terrain_stack1, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/terrain_variables.tif", overwrite = TRUE)

# Plot
plot(terrain_stack[["northerness_11x11"]])  # Example


