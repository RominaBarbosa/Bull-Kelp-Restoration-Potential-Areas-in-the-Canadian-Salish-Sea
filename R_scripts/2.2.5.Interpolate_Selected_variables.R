###=============================================================================
###     NEMO model - SalishSeaCast Variables INTERPOLATION      ################
###                                                             ################
### input data acquired with: "2.2.4.Variables_correlation_all" ################
### input data: selected variable metrics (dataframe)           ################
### output data: rasters of variables  (metrics)                ################
### Author: Romina Barbosa                                      ################
### Date: 21-July-2025                                          ################
### Last edition: 22-July-2025   
###=============================================================================
library(sf)
library(terra)
library(dplyr)



Selected_variables_names<- read.csv( "/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_WOut_bathymetricVars_FINAL.csv")
Selected_variables_names<- colnames(Selected_variables_names)[4:ncol(Selected_variables_names)]
# "ammonium_fall_mean"              "ammonium_spring_SD"              "ammonium_summer_mean"           
# [4] "PAR_summer_mean"                 "salinity_fall_SD"                "salinity_summer_mean"           
# [7] "salinity_summer_SD"              "temperature_summer_mean"         "turbidity_summer_mean"          
# [10] "uVelocity_summer_maximum"        "uVelocity_summer_mean"           "uVelocity_summer_minimum"       
# [13] "vVelocity_summer_maximum"        "vVelocity_summer_mean"           "vVelocity_summer_minimum"       
# [16] "ammonium_spring_mean.y"          "summer_hours_above_threshold_18"

# Add ".csv" to each
# desired_filenames <- paste0(Selected_variables_names, ".csv")


## SET PATCH AND UPLOAD FILES of env variables
### Upload variables from the SalishSeaCast model and merge ====================
my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics")
setwd(my_path)
dir()

# List all CSVs in your folder
files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE)


# Load file 1 to start building the merged dataset
df_all=read.csv(files[1])
head(df_all)
df_all<- df_all[,c(2:4)]
filename <- strsplit(tools::file_path_sans_ext(files[1]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(1,3,4)]
filename <- paste0(filename, collapse ="_")
colnames(df_all)[3]<- filename

# Merge all files by adding a column of the metric to the df_all dataset
for (i in 2:length(files)) {
  
  df= read.csv(files[i])
  df= df[,c(2:4)]
  
  if(i %in% c(1:16, 33:76, 93:140)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(1,3,4)]
    filename <- paste0(filename, collapse ="_")
  }
  
  if(i %in% c(17:32)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(5:6)]
    filename <- paste("DIC", paste0(filename, collapse ="_"), sep="_")
  }
  
  if(i %in% c(77:92)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(4:5)]
    filename <- paste("alkalinity", paste0(filename, collapse ="_"), sep="_")
  }
  
  colnames(df)[3]<- filename
  df_all= merge(df_all, df, by=c("latitude" , "longitude"))
  
}

SSCast_variables<- df_all



# Filter the files to include only those in desired_filenames
Selected_variables_names<- c( "latitude", "longitude", Selected_variables_names)
selected_variables<- SSCast_variables[, (names(SSCast_variables) %in% Selected_variables_names)]
head(selected_variables, 2)


###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
# write.csv(selected_variables, "/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_to_interpolations.csv")





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
output_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers"
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
    writeRaster(idw_raster[[1]], paste(output_path, paste(variable_name , "interpolated_output.tif", sep="_"), sep="/"), overwrite = TRUE)
    
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


