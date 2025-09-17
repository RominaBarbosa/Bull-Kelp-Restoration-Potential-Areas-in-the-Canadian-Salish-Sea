###=============================================================================
###     NEMO model - SalishSeaCast Variables INTERPOLATION      ################
###                                                             ################
### input data acquired with: "2.2.4.Variables_correlation_all" ################
### and "1.2.3.Climatology_metrics_Xperiod.R"                   ################
### input data: selected variable metrics (climatology)         ################
### output data:1- rasters of variables  (metrics)              ################
###             2- 20 m resolution multiraster (terra) of variables ############
### Author: Romina Barbosa                                      ################
### Date: 21-July-2025                                          ################
### Last edition: 15-Sep-2025   
###=============================================================================
# 19-Aug-2025 : updated selected variables based on GLM and GAM with scaled variables

library(sf)
library(terra)
library(dplyr)


# After running models GLM and GAM with scaled variables
# [1] "slope_7x7"               "ammonium_spring_SD"      "ammonium_winter_minimum"    "PAR_summer_mean"         "temperature_summer_mean"
# [6] "turbidity_summer_mean"   "salinity_summer_SD"      "nitrate_summer_mean"     "ammonium_summer_mean"    "ammonium_spring_mean"  

# selected_variables<- c("slope_7x7",   "ammonium_spring_SD",   "ammonium_summer_SD", 
#                        "ammonium_winter_mean",    "PAR_summer_mean", 
#                        "temperature_summer_mean",
#               "turbidity_summer_mean",   "salinity_summer_SD",      
#               "nitrate_summer_minimum",   "nitrate_winter_minimum", 
#               "ammonium_summer_minimum",  "ammonium_spring_minimum", 
#               "PAR_summer_maximum", "currentSpeed_summer_mean")

# selected_variables<- c(
#   "nitrate_summer_minimum",
#   "nitrate_winter_minimum", # correlated with   "salinity_summer_minimum",
#   "currentSpeed_summer_mean",
#   "salinity_summer_SD", # correlated with DIC_spring_SD
#   "ammonium_winter_minimum",
#   "ammonium_spring_SD",
#   "ammonium_summer_SD",
#   "temperature_summer_mean",
#   "PAR_summer_mean",# "", # correlated with "PAR_spring_mean",
#   "PAR_summer_maximum",
#   "ammonium_spring_mean",
#   "turbidity_summer_mean",# corr with turbidity_spring_mean
#   "ammonium_summer_minimum" # corr with ammonium_summer_mean and maximum
# )

selected_variables<- c("ammonium_spring_mean",
                       "ammonium_summer_minimum",
                       "ammonium_winter_mean",
                       "currentSpeed_summer_mean",
                       "nitrate_summer_minimum",
                       "nitrate_winter_mean",
                       "PAR_summer_mean",
                       "temperature_summer_mean",
                       "turbidity_summer_mean")
# [1]] "ammonium_spring_mean"     "ammonium_summer_minimum"  "ammonium_winter_mean"     "currentSpeed_summer_mean" "nitrate_summer_minimum"  
# [6] "nitrate_winter_mean"      "PAR_summer_mean"          "temperature_summer_mean"  "turbidity_summer_mean"    "longitude"               
# [11] "latitude" 


# Upload variables from the SalishSeaCast model and merge 
# my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_blob_")
my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_post_blob")
setwd(my_path)
dir()
files <- list.files(my_path, pattern = "\\.csv$", full.names = TRUE)
length(files)# 140


# Helper function to extract variable name from filename
extract_var_name <- function(file_path, index, period="post_blob") {
  filename <- basename(tools::file_path_sans_ext(file_path))
  parts <- unlist(strsplit(filename, "_"))
  
  if(period =="blob"){
    # Your previous logic for variable naming by index ranges:
    if (index %in% c(1:16, 65:128, 145:160)) {
      var_name <- paste(parts[c(1, 3, 4)], collapse = "_")
    } else if (index %in% c(17:32)) {
      var_name <- paste("currentDirection", paste(parts[4:5], collapse = "_"), sep = "_")
    } else if (index %in% c(33:48)) {
      var_name <- paste("currentSpeed", paste(parts[4:5], collapse = "_"), sep = "_")
    } else if (index %in% c(49:64)) {
      var_name <- paste("DIC", paste(parts[5:6], collapse = "_"), sep = "_")
    } else if (index %in% c(129:144)) {
      var_name <- paste("alkalinity", paste(parts[4:5], collapse = "_"), sep = "_")
    } else {
      var_name <- NA_character_
    }
  }
  
  
  if(period =="post_blob"){
    # Your previous logic for variable naming by index ranges:
    if (index %in% c(1:16, 45:108, 125:140)) {
      var_name <- paste(parts[c(1, 4, 5)], collapse = "_")
    # } else if (index %in% c(17:32)) {
    #   var_name <- paste("currentDirection", paste(parts[5:6], collapse = "_"), sep = "_")
    } else if (index %in% c(17:28)) {
      var_name <- paste("currentSpeed", paste(parts[5:6], collapse = "_"), sep = "_")
    } else if (index %in% c(29:44)) {
      var_name <- paste("DIC", paste(parts[6:7], collapse = "_"), sep = "_")
    } else if (index %in% c(109:124)) {
      var_name <- paste("alkalinity", paste(parts[5:6], collapse = "_"), sep = "_")
    } else {
      var_name <- NA_character_
    }
  }
  
  return(var_name)
}



# Initialize merged data with the first file that matches selected variables
df_all <- NULL
first_file_found <- FALSE

for (i in seq_along(files)) { #
  
  var_name <- extract_var_name(files[i], i, period ="post_blob")
  
  # Only process if variable is in selected_variables
  if (!is.na(var_name) && var_name %in% selected_variables) {
    df <- read.csv(files[i])
    df <- df[, c(5:6, ncol(df))]  # Keep relevant columns
    
    colnames(df)[ncol(df)] <- var_name
    
    if (!first_file_found) {
      df_all <- df
      first_file_found <- TRUE
      
      print(dim(df_all))
      head(df_all)
    } else {
      # Merge by spatial coordinates & bathymetry
      unmatched_rows <- anti_join(df_all, df, by = c("latitude", "longitude"))
      print(paste("Number of unmatched rows:", nrow(unmatched_rows), "file:", files[i]))
      df_all <- merge(df_all, df, by = c("latitude", "longitude"))
    }
  }
}

# df_all now contains only selected variables merged, aligned by spatial coords
head(df_all)
selected_variables %in% colnames(df_all) # all TRUE except the slopevariables

selected_variables_df<- df_all
colnames(selected_variables_df)
# colnames(selected_variables_df)
# [1] "latitude"                 "longitude"                "ammonium_spring_mean"     "ammonium_spring_SD"      
# [5] "ammonium_summer_minimum"  "ammonium_summer_SD"       "ammonium_winter_minimum"  "currentSpeed_summer_mean"
# [9] "nitrate_summer_minimum"   "nitrate_winter_minimum"   "PAR_summer_maximum"       "PAR_summer_mean"         
# [13] "salinity_summer_SD"       "temperature_summer_mean"  "turbidity_summer_mean"   

# [1] "latitude"                 "longitude"                "ammonium_spring_mean"     "ammonium_summer_minimum"  "ammonium_winter_mean"    
# [6] "currentSpeed_summer_mean" "nitrate_summer_minimum"   "nitrate_winter_mean"      "PAR_summer_mean"          "temperature_summer_mean" 
# [11] "turbidity_summer_mean"  


# M3:


###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
# write.csv(selected_variables_df, "/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_to_interpolations_SepFINAL.csv")
# write.csv(selected_variables_df, "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/selected_variables_to_interpolations_postBlob_M2.csv") # 19 August 2025



## Add temperature variables
my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_temp_tolerance_metrics")
setwd(my_path)
files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE) # 32 files
files

df_variable= read.csv(files[20])
df_variable<- df_variable[,c(5:6, ncol(df_variable))]
filename <- strsplit(tools::file_path_sans_ext(files[20]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(3:6)]
filename <- paste0(filename, collapse ="_")
colnames(df_variable)[ncol(df)]<- filename

df_variable2= read.csv(files[28])
df_variable2<- df_variable2[,c(5:6, ncol(df_variable2))]
filename <- strsplit(tools::file_path_sans_ext(files[28]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(3:7)]
filename <- paste0(filename, collapse ="_")
colnames(df_variable2)[ncol(df_variable2)]<- filename

variables_1_2= merge(df_variable, df_variable2, by= c("latitude", "longitude"))
head(variables_1_2)


####============================================================================
### INTERPOLATION of Variables 
library(ncdf4)
library(sp)
library(raster)
library(gstat)
library(sf)
library(dplyr)
library(terra)

# Load input data
bathy <- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
# output_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated_Blob"
# variables_df <- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_to_interpolations_SepFINAL.csv")
variables_df <- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/selected_variables_to_interpolations_postBlob_M2.csv")
output_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated_postBlob"
variables_df<- variables_df[,-1] 
colnames(variables_df)


# variables_interoplation_function<- function(bathy= bathy, variables_data= variables_df, output_path= output_path){
  
  # if(colnames(variables_df)[1] != "gridY" & colnames(variables_df)[1] != "gridX"){
  #   variables_df<- variables_df[,-1]
  # }
  
# variables_df<- variables_1_2
  
  for (i in 3:ncol(variables_df)) { # 
    variable_name<- colnames(variables_df)[i]
    variable<- variables_df[c("longitude", "latitude", paste(variable_name))]
    
    if(variable_name == "currentSpeed_summer_mean"){
      variable[which(variable$currentSpeed_summer_mean == 0),]<- NA
    }
    
    variable<- na.exclude(variable)
    
    # Replace 'longitude' and 'latitude' with your actual column names
    points_sf <- st_as_sf(variable, coords = c("longitude", "latitude"), crs = 4326)
    
    # 2. Reproject to match your bathymetry raster (assuming EPSG:3005)
    points_sf <- st_transform(points_sf, crs = 3005)
    
    # 3. Convert to SpatialPointsDataFrame
    points_sp <- as(points_sf, "Spatial")
    names(points_sp)<- "value"
    points_sp <- points_sp[!is.na(points_sp$value), ]
    
    # 4. Create raster grid 
    template <- rast(ext(bathy), resolution = 500, crs = crs(bathy))
    template_raster <- raster(template)  # Convert SpatRaster to RasterLayer
    template_sp <- as(template_raster, "SpatialPixelsDataFrame")
    
    # 5. Run IDW interpolation using sp::idw()
    idw_result <- gstat::idw(formula = value ~ 1, locations = points_sp, newdata = template_sp, idp = 3.5, nmax = 12)#inverse distance weighted interpolation
    idw_raster <- rast(idw_result)# Convert idw result back to terra raster
    
    raster<- idw_raster[[1]]
    
    if(variable_name %in% c("currentSpeed_summer_mean")){ # current speed was calculated from u and v velocity which are information from vectors positioned at the bottom left of the cell, so I moved the raster to correct this
      r<- idw_raster$var1.pred
      res_x <- xres(r)
      res_y <- yres(r)
      
      # Create a new raster with the same extent + resolution
      r_shifted <- r
      
      # Shift the raster values bottom-left by 1 cell
      r_shifted[] <- NA  # fill new raster with NAs
      nrows <- nrow(r)
      ncols <- ncol(r)
      
      # Paste original values starting at row 2 (shift down)
      # Column starts at 1 (no horizontal shift) for bottom-left
      r_shifted[2:nrows, 1:ncols] <- r[1:(nrows-1), 1:ncols]
      r_shifted[2:nrows, 1:(ncols-1)] <- r[1:(nrows-1), 2:ncols]
      
      raster<- r_shifted
    }
  
    # 6. Plot or save
    # plot(r_corrected, main= paste("interpolation:", variable_name))
    terra::writeRaster(raster, paste(output_path, paste(variable_name , "postblob_interpolated_output.tif", sep="_"), sep="/"), overwrite = TRUE)
    
  }  



# # Load your bathymetry raster (20m resolution)
# bathymetry_20m <- rast("/your/path/to/bathymetry_20m.tif")
# 
# # Resample 500m raster stack to 20m grid using bathymetry as template
# # method = "bilinear" or "near" (nearest neighbor, for categorical)
# raster_stack_20m <- resample(raster_stack_500m, bathymetry_20m, method = "bilinear")
# 
# # Check the result
# print(raster_stack_20m)
# plot(raster_stack_20m[[1]])  # plot first layer

# Downscale SalishSeaCast model variables to 20 m resolution ========================================
# List all tif files
rasters_path<- output_path
tif_files <- list.files(rasters_path, pattern = "\\.tif$", full.names = TRUE)


# Load and stack all tif files
raster_stack <- rast(tif_files)
# layer_names <- basename(tif_files) %>%
#   sub("_interpolated_output.*", "", .)  # Remove everything from '_interpolated_output' onward

layer_names <- basename(tif_files) %>%
  sub("_postblob_interpolated_output.*", "", .)  # Remove everything from '_interpolated_output' onward


# Assign names to raster layers
names(raster_stack) <- layer_names


# selected_variables<- c(
#   "nitrate_summer_minimum",
#   "nitrate_winter_minimum", # correlated with   "salinity_summer_minimum",
#   "currentSpeed_summer_mean",
#   "salinity_summer_SD", # correlated with DIC_spring_SD
#   "ammonium_winter_minimum",
#   "ammonium_spring_SD",
#   "ammonium_summer_SD",
#   "temperature_summer_mean",
#   "PAR_summer_mean",# "", # correlated with "PAR_spring_mean",
#   "PAR_summer_maximum",
#   "ammonium_spring_mean",
#   "turbidity_summer_mean",# corr with turbidity_spring_mean
#   "ammonium_summer_minimum" # corr with ammonium_summer_mean and maximum
# )

# selected_variables<- c("ammonium_spring_mean",
#                        "ammonium_summer_minimum",
#                        "ammonium_winter_mean",
#                        "currentSpeed_summer_mean",
#                        "nitrate_summer_minimum",
#                        "nitrate_winter_mean",
#                        "PAR_summer_mean",
#                        "temperature_summer_mean",
#                        "turbidity_summer_mean")

raster_stack <- raster_stack[[selected_variables]]

# Load your bathymetry raster (20m resolution)
bathymetry_20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
# Resample 500m raster stack to 20m grid using bathymetry as template
# method = "bilinear" or "near" (nearest neighbor, for categorical)

raster_stack_20m <- resample(raster_stack, bathymetry_20m, method = "bilinear")
# names(raster_stack_20m) <- layer_names

# Check the result
print(raster_stack_20m)
plot(raster_stack_20m[[1]])  # plot first layer

# Save raster stack of resampled variables at 20 m resolution
terra::writeRaster(raster_stack_20m, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblobM2.tif",
                   overwrite=TRUE)
# terra::writeRaster(raster_stack_20m, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_postBlob.tif",
#                    overwrite=TRUE)


names(raster_stack_20m)
# terra::writeRaster(raster_stack_20m$currentSpeed_summer_mean, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/currentSpeed_summer_mean_20m_resolution.tif",
#                    overwrite=TRUE)






###=============================================================================
### Interpolate tidal current from Foreman model ===============================
tidal_current_pts<- st_read("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Foreman RMS Tidal Model/Foreman_rms_tidal.shp")

# Reproject to match your bathymetry raster (assuming EPSG:3005)
points_sf <- st_transform(tidal_current_pts, crs = 3005)
points_sp <- as(points_sf, "Spatial")
# names(points_sp)<- "value"

# Create raster grid 
template <- rast(ext(bathy), resolution = 500, crs = crs(bathy))
template_raster <- raster(template)  # Convert SpatRaster to RasterLayer
template_sp <- as(template_raster, "SpatialPixelsDataFrame")

# Run IDW interpolation using sp::idw()
idw_result <- idw(formula = RMS_Tidal_ ~ 1, locations = points_sp, newdata = template_sp, idp = 3, nmax = 12)#inverse distance weighted interpolation
idw_raster <- rast(idw_result)# Convert idw result back to terra raster

# Plot or save
plot(idw_raster, main= paste("interpolation:", "RMS_Tidal_"))

output_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Foreman RMS Tidal Model"
# writeRaster(idw_raster[[1]], paste(output_path, paste("RMS_Tidal" , "interpolated.tif", sep="_"), sep="/"), overwrite = TRUE)




# ###===============================s==============================================
# # The entire area includes a lot of grids and we need to reduce the number of cells
# # We will restrict the interpolated area to areas with depth 0-30 m
# 
# crs(bathy)
# 
# # Stack all env. variables 
# # rasters.t<- stack_rasters_path_funtion(path= output_path)
# files_rasters <- list.files(path = output_path, pattern = "\\.tif$", full.names = TRUE)
# 
# rasters<- stack(files_rasters[1])
# filename <- strsplit(tools::file_path_sans_ext(files_rasters[1]), "/")[[1]][7]
# filename <- strsplit(filename, "_")[[1]][c(1:3)]
# filename <- paste0(c(filename, "blob"), collapse ="_")
# names(rasters)<- filename
# 
# 
# for (i in files_rasters[2:length(files_rasters)]) {
#   raster_i<- raster(i)
#   filename <- strsplit(tools::file_path_sans_ext(i), "/")[[1]][7]
#   filename <- strsplit(filename, "_")[[1]][c(1:3)]
#   filename <- paste0(c(filename, "blob"), collapse ="_")
#   names(raster_i)<- filename
#   rasters<- stack(rasters, raster_i)
# }
# 
# 
# # Stack all TOPOGRAPHIC variables
# files_rasters.t <- list.files(path = "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables", pattern = "\\.tif$", full.names = TRUE)
# 
# files_rasters.t<- files_rasters.t[-2]# - easterness
# files_rasters.t<- files_rasters.t[-3]# - bathy salishseacast model
# files_rasters.t<- files_rasters.t[-5]# - TRI
# 
# 
# rasters.t<- stack(files_rasters.t[1])
# 
# for (i in files_rasters.t[2:length(files_rasters.t)]) {
#   raster_i<- raster(i)
#   rasters.t<- stack(rasters.t, raster_i)
# }
# 
# 
# ### Mask rasters by bathymetry =================================================
# bathy <- raster("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
# summary(bathy)
# bathy_filtered <- bathy
# bathy_filtered[ bathy>= 150 ]<- NA
# bathy_filtered[ bathy_filtered<= -50 ]<- NA
# plot(bathy_filtered)
# 
# 
# ## Stack env. variables and topographyc variables 
# variables_rasters<- stack(rasters, rasters.t)
# rasters_masked <- raster::mask(variables_rasters, bathy_filtered)
# 
# # Plot rasters
# plot(rasters_masked, ncol=3)
# plot(rasters_masked[["PAR_summer_mean_blob"]])
# plot(rasters_masked[["turbidity_summer_mean_blob"]])
# plot(rasters_masked[["temperature_summer_mean_blob"]])
# plot(rasters_masked[["salinity_summer_mean_blob"]])
