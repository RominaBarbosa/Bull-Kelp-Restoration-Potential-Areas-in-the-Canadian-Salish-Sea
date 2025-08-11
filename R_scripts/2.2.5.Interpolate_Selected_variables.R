###=============================================================================
###     NEMO model - SalishSeaCast Variables INTERPOLATION      ################
###                                                             ################
### input data acquired with: "2.2.4.Variables_correlation_all" ################
### input data: selected variable metrics (dataframe)           ################
### output data: rasters of variables  (metrics)                ################
### Author: Romina Barbosa                                      ################
### Date: 21-July-2025                                          ################
### Last edition: 11-Aug-2025   
###=============================================================================
# 11-Aug-2025 : updated selected variables 

library(sf)
library(terra)
library(dplyr)



selected_variables<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_nemo_masked.csv")
colnames(selected_variables)
# [1] "latitude"                        "longitude"                       
# "ammonium_spring_mean"           
# [4] "ammonium_spring_SD"              "ammonium_summer_mean"            "ammonium_summer_SD"             
# [7] "ammonium_winter_mean"            "currentDirection_spring_min"     "currentDirection_summer_modal"  
# [10] "currentSpeed_summer_mean"        "nitrate_summer_mean"             "PAR_summer_maximum"             
# [13] "PAR_summer_mean"                 "salinity_summer_mean"            "salinity_summer_SD"             
# [16] "temperature_summer_mean"         "turbidity_summer_mean"           "summer_cumulated_degrees_18"    
# [19] "summer_hours_above_threshold_18"

selected_variables<- colnames(selected_variables)
selected_variables<- selected_variables[3:length(selected_variables)]

# Upload variables from the SalishSeaCast model and merge 
my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_blob_")
setwd(my_path)
dir()
files <- list.files(my_path, pattern = "\\.csv$", full.names = TRUE)
length(files)# 160


# Helper function to extract variable name from filename
extract_var_name <- function(file_path, index) {
  filename <- basename(tools::file_path_sans_ext(file_path))
  parts <- unlist(strsplit(filename, "_"))
  
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
  
  return(var_name)
}

# Initialize merged data with the first file that matches selected variables
df_all <- NULL
first_file_found <- FALSE

for (i in seq_along(files)) {
  var_name <- extract_var_name(files[i], i)
  
  # Only process if variable is in selected_variables
  if (!is.na(var_name) && var_name %in% selected_variables) {
    df <- read.csv(files[i])
    df <- df[, c(2:6, ncol(df))]  # Keep relevant columns
    
    colnames(df)[ncol(df)] <- var_name
    
    if (!first_file_found) {
      df_all <- df
      first_file_found <- TRUE
    } else {
      # Merge by spatial coordinates & bathymetry
      df_all <- merge(df_all, df, by = c("gridY", "gridX", "bathymetry", "latitude", "longitude"))
    }
  }
}

# df_all now contains only selected variables merged, aligned by spatial coords
head(df_all)

# # Load file 1 to start building the merged dataset
# df_all=read.csv(files[1])
# head(df_all)
# df_all<- df_all[,c(2:6,ncol(df_all))]
# filename <- strsplit(tools::file_path_sans_ext(files[1]), "/")[[1]][7]
# filename <- strsplit(filename, "_")[[1]][c(1,3,4)]
# filename <- paste0(filename, collapse ="_")
# colnames(df_all)[ncol(df_all)]<- filename
# 
# # Merge all files by adding a column of the metric to the df_all dataset
# for (i in 117:160) {
#   
#   df= read.csv(files[i])
#   df= df[,c(2:6,ncol(df))]
#   
#   if(i %in% c(1:16, 65:128, 145:160)){
#     filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
#     filename <- strsplit(filename, "_")[[1]][c(1,3,4)]
#     filename <- paste0(filename, collapse ="_")
#   }
#   
#   if(i %in% c(17:32)){
#     filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
#     filename <- strsplit(filename, "_")[[1]][c(4:5)]
#     filename <- paste("currentDirection", paste0(filename, collapse ="_"), sep="_")
#   }
#   
#   if(i %in% c(33:48)){
#     filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
#     filename <- strsplit(filename, "_")[[1]][c(4:5)]
#     filename <- paste("currentSpeed", paste0(filename, collapse ="_"), sep="_")
#   }
#   
#   if(i %in% c(49:64)){
#     filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
#     filename <- strsplit(filename, "_")[[1]][c(5:6)]
#     filename <- paste("DIC", paste0(filename, collapse ="_"), sep="_")
#   }
#   
#   if(i %in% c(129:144)){
#     filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
#     filename <- strsplit(filename, "_")[[1]][c(4:5)]
#     filename <- paste("alkalinity", paste0(filename, collapse ="_"), sep="_")
#   }
# 
#   colnames(df)[ncol(df)]<- filename
#   df_all= merge(df_all, df, by=c("gridY", "gridX", "bathymetry", "latitude" , "longitude"))
#   
# }
# 
# selected_variables_df<- df_all[,which(colnames(df_all) %in% colnames(selected_variables))]
# colnames(selected_variables_df)
# rm(df_all)


## Add temperature variables 
my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_temp_tolerance_metrics_")
files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE) # 32 files

df_variable= read.csv(files[24])
df_variable<- df_variable[,c(2:4)]
filename <- strsplit(tools::file_path_sans_ext(files[24]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(3:6)]
filename <- paste0(filename, collapse ="_")
colnames(df_variable)[3]<- filename

selected_variables_df= merge(selected_variables_df, df_variable, by=c("latitude" , "longitude"))

df_variable= read.csv(files[32])
df_variable<- df_variable[,c(2:4)]
filename <- strsplit(tools::file_path_sans_ext(files[32]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(3:7)]
filename <- paste0(filename, collapse ="_")
colnames(df_variable)[3]<- filename

selected_variables_df= merge(selected_variables_df, df_variable, by=c("latitude" , "longitude"))
colnames(selected_variables_df)
# [1] "latitude"                        "longitude"                       
# "ammonium_spring_mean"           
# [4] "ammonium_spring_SD"              "ammonium_summer_mean"            "ammonium_summer_SD"             
# [7] "ammonium_winter_mean"            "currentDirection_spring_min"     "currentDirection_summer_modal"  
# [10] "currentSpeed_summer_mean"        "nitrate_summer_mean"             "PAR_summer_maximum"             
# [13] "PAR_summer_mean"                 "salinity_summer_mean"            "salinity_summer_SD"             
# [16] "temperature_summer_mean"         "turbidity_summer_mean"           "summer_cumulated_degrees_18"    
# [19] "summer_hours_above_threshold_18"


###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
# write.csv(selected_variables, "/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_to_interpolations.csv")



library(ncdf4)

####============================================================================
### INTERPOLATION of Variables 
library(sp)
library(raster)
library(gstat)
library(sf)
library(dplyr)
library(terra)

# Load input data
bathy <- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
output_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated"
variables_df <- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_to_interpolations.csv")

# variables_interoplation_function<- function(bathy= bathy, variables_data= variables_df, output_path= output_path){
  
  if(colnames(variables_df)[1] != "latitude" & colnames(variables_df)[1] != "longitude"){
    variables_df<- variables_df[,-1]
  }
  
  
  for (i in 3:ncol(variables_df)) {
    variable_name<- colnames(variables_df)[i]
    variable<- variables_df[c("longitude", "latitude", paste(variable_name))]
    
    # Replace 'longitude' and 'latitude' with your actual column names
    points_sf <- st_as_sf(variable, coords = c("longitude", "latitude"), crs = 4326)
    
    # 2. Reproject to match your bathymetry raster (assuming EPSG:3005)
    points_sf <- st_transform(points_sf, crs = 3005)
    
    # 3. Convert to SpatialPointsDataFrame
    points_sp <- as(points_sf, "Spatial")
    names(points_sp)<- "value"
    
    # 4. Create raster grid 
    template <- rast(ext(bathy), resolution = 500, crs = crs(bathy))
    template_raster <- raster(template)  # Convert SpatRaster to RasterLayer
    template_sp <- as(template_raster, "SpatialPixelsDataFrame")
    
    # 5. Run IDW interpolation using sp::idw()
    idw_result <- idw(formula = value ~ 1, locations = points_sp, newdata = template_sp, idp = 3.5, nmax = 12)#inverse distance weighted interpolation
    idw_raster <- rast(idw_result)# Convert idw result back to terra raster
    
    # 6. Plot or save
    plot(idw_raster, main= paste("interpolation:", variable_name))
    terra::writeRaster(idw_raster[[1]], paste(output_path, paste(variable_name , "interpolated_output.tif", sep="_"), sep="/"), overwrite = TRUE)
    
  }  
# }



###===============================s==============================================
# The entire area includes a lot of grids and we need to reduce the number of cells
# We will restrict the interpolated area to areas with depth 0-30 m

crs(bathy)

# Stack all env. variables 
# rasters.t<- stack_rasters_path_funtion(path= output_path)
files_rasters <- list.files(path = output_path, pattern = "\\.tif$", full.names = TRUE)

rasters<- stack(files_rasters[1])
filename <- strsplit(tools::file_path_sans_ext(files_rasters[1]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(1:3)]
filename <- paste0(c(filename, "blob"), collapse ="_")
names(rasters)<- filename


for (i in files_rasters[2:length(files_rasters)]) {
  raster_i<- raster(i)
  filename <- strsplit(tools::file_path_sans_ext(i), "/")[[1]][7]
  filename <- strsplit(filename, "_")[[1]][c(1:3)]
  filename <- paste0(c(filename, "blob"), collapse ="_")
  names(raster_i)<- filename
  rasters<- stack(rasters, raster_i)
}


# Stack all TOPOGRAPHIC variables
files_rasters.t <- list.files(path = "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables", pattern = "\\.tif$", full.names = TRUE)

files_rasters.t<- files_rasters.t[-2]# - easterness
files_rasters.t<- files_rasters.t[-3]# - bathy salishseacast model
files_rasters.t<- files_rasters.t[-5]# - TRI


rasters.t<- stack(files_rasters.t[1])

for (i in files_rasters.t[2:length(files_rasters.t)]) {
  raster_i<- raster(i)
  rasters.t<- stack(rasters.t, raster_i)
}


### Mask rasters by bathymetry =================================================
bathy <- raster("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
summary(bathy)
bathy_filtered <- bathy
bathy_filtered[ bathy>= 150 ]<- NA
bathy_filtered[ bathy_filtered<= -50 ]<- NA
plot(bathy_filtered)


## Stack env. variables and topographyc variables 
variables_rasters<- stack(rasters, rasters.t)
rasters_masked <- raster::mask(variables_rasters, bathy_filtered)

# Plot rasters
plot(rasters_masked, ncol=3)
plot(rasters_masked[["PAR_summer_mean_blob"]])
plot(rasters_masked[["turbidity_summer_mean_blob"]])
plot(rasters_masked[["temperature_summer_mean_blob"]])
plot(rasters_masked[["salinity_summer_mean_blob"]])




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


