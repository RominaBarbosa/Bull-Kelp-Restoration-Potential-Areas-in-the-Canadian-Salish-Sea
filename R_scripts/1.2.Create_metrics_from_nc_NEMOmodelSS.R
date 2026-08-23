###=============================================================================
###     NEMO model - SalishSeaCast data process                 ################
###                                                             ################
### 1.2- Loop to merge monthly data into seasons                ################
###    - and calculate metrics
### https://salishsea.eos.ubc.ca/                               ################
### input data acquired with: "1.1.Download_nc_NEMOmodel.R"     ################
### input data: month data of hourly series                     ################
### output data: metric of a season                             ################
### Author: Romina Barbosa                                      ################
### Date: 15-Feb-2025                                           ################
### Last edition: 25-March-2025
### converted into metrics 
###=============================================================================

# Load packages
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(ggplot2) # package for plotting
library(lubridate)
library(chron)
library(mondate)
library(abind)

# memory.limit(size = 20000) # only works for windows


## Set Paths ===================================================================
input_path= "/modeled_variables_original/Monthly_nc"
setwd(input_path)
output_path=  "/modeled_variables_original/seasonal_metrics_0.5m"
path_bathy= "Volumes/Romina_PSF/modeled_variables_original"


### Different metrics are calculated:
## mean = average from hourly data from the entire season
## SD = Standard deviation of the hourly data from the entire season
## minimum = average value from the 0.1 percentile, i.e., mean of minimum values
## maximum = average value from the 0.9 percentile, i.e., mean of maximum values

source(paste("/Volumes/Romina_PSF/PSF/R_scripts/functions","RFunction_merge_nemo_seasonal_data.R", sep="/"))
source(paste("/Volumes/Romina_PSF/PSF/R_scripts/functions","RFunction_calculate_simple_metrics.R", sep="/"))

years<- c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")
# variables<- c("ammonium", "uVelocity", "total_alkalinity")# vVelocity, "u_wind","turbidity",  "salinity", "nitrate", 
variables<- c( "PAR") #"PAR", "u_wind",
season<- c( "fall", "spring", "winter", "summer") 


for (v in 1:length(variables)) {
  variable_v<- variables[v]
    
  for (s in 1:length(season)) {
    season_s<- season[s]
    
    for(year in years){
      merged_spring_data <- merge_nemo_seasonal_data(depth = "surface", year = year, variable = variable_v, 
                                                     input_path = input_path, 
                                                     season = season_s)
      
      calculate_simple_metrics(data = merged_spring_data, 
                               season= season_s,
                               year= year,
                               depth = "surface",
                               variable= variable_v,
                               output_path= output_path, 
                               path_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/"))
      
      time= Sys.time()
      print(paste("finished", year, ";", variable_v, time))
    }
  }
}






















# # Install and load necessary libraries
# install.packages("gstat")
# install.packages("raster")
# install.packages("sp")
# library(gstat)
# library(raster)
# library(sp)
# 
# # Check for missing values in the original data
# sum(is.na(merged_df$longitude))  # Check if longitude has NA values
# sum(is.na(merged_df$latitude))   # Check if latitude has NA values
# sum(is.na(merged_df$value))      # Check if value has NA values
# 
# # Remove rows with NA values
# merged_df <- na.omit(merged_df)
# 
# # Convert the data into a SpatialPointsDataFrame
# coordinates(merged_df) <- ~longitude + latitude
# proj4string(merged_df) <- CRS("+proj=longlat +datum=WGS84")  # Set CRS (WGS84)
# 
# # Define a regular grid (for example, resolution of 0.1 degrees)
# grid_lon <- seq(min(merged_df$longitude), max(merged_df$longitude), by = 0.1)
# grid_lat <- seq(min(merged_df$latitude), max(merged_df$latitude), by = 0.1)
# 
# # Create a grid of coordinates
# grid <- expand.grid(longitude = grid_lon, latitude = grid_lat)
# coordinates(grid) <- ~longitude + latitude
# proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")
# 
# # Perform ordinary kriging interpolation
# kriging_model <- gstat(id = "value", formula = value ~ 1, data = merged_df)
# interpolated_values <- tryCatch({
#   predict(kriging_model, newdata = grid)
# }, error = function(e) {
#   message("Error in kriging interpolation: ", e)
#   NULL
# })
# 
# # Check the result of the interpolation
# if (!is.null(interpolated_values)) {
#   # Convert to raster
#   raster_data <- rasterFromXYZ(cbind(interpolated_values$longitude, interpolated_values$latitude, interpolated_values$var1.pred))
#   
#   # Save the raster as a GeoTIFF file
#   writeRaster(raster_data, filename = "output_interpolated_raster.tif", format = "GTiff", overwrite = TRUE)
#   
#   # Plot the raster
#   plot(raster_data)
# } else {
#   message("Interpolation failed.")
# }




## ===================================================================
## Calculate cumulated degrees above temperature threshold metric
## ===================================================================

# Steps:
# Set the threshold (e.g., 18°C).
# Filter the temperature values above the threshold (18°C).
# Calculate the cumulated degrees: This is the sum of the temperatures that are above the threshold for each location.
# Count the number of hours: This is the number of time steps where the temperature exceeds the threshold for each location.


# Assuming 'temperature_combined' is the 3D array (time x lat x lon)
# With dimensions (time, latitude, longitude)

bathy_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/")
calculate_DegreeHours_Hours_OverThreshold<- function(data= merged_spring_data$merged_data, threshold_temp= 18, bathy_file=bathy_file){
  
  # Step 1: Set the threshold temperature (e.g., 15°C)
  # threshold_temp <- 18
  
  # Step 2: Initialize matrices to store the cumulated temperature and number of hours
  cumulated_degrees <- matrix(NA, nrow = (dim(data)[1]), ncol = (dim(data)[2]))
  num_hours_above_threshold <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
  
  # Step 3: Iterate over each latitude and longitude
  for (i in 1:dim(data)[1]) {
    for (j in 1:dim(data)[2]) {
      
      # Step 4: Extract all temperature values for this lat, lon
      location_values <- data[ i, j, ]  # Extract values for this lat, lon
      
      # Step 5: Identify the values above the threshold (15°C)
      above_threshold_values <- location_values[location_values > threshold_temp]
      
      which(!is.na(above_threshold_values))
      
      # Step is.na()# Step 6: Calculate the accumulated degrees (sum of temperatures above the threshold)
      cumulated_degrees[i, j] <- sum(above_threshold_values, na.rm = TRUE)
      
      # Step 7: Calculate the number of hours with temperatures above the threshold
      num_hours_above_threshold[i, j] <- length(above_threshold_values)
    }
    
    merged_array<- abind(cumulated_degrees, num_hours_above_threshold, along = 3, 
                         make.names = TRUE)
  }
  
  # Step 8: add lat and lon coordinates
  for (x in 1:dim(merged_array)[3]) {
    # Convert the matrix into a dataframe with row and column names
    value_df <- as.data.frame(as.table(merged_array[,,x]))
    
    # Rename columns to make it more intuitive
    colnames(value_df) <- c( "gridX", "gridY","value")
    length(unique(value_df$gridX))
    
    bathy_grid<-  read.csv(bathy_file)
    
    
    # Merge the latitude-longitude dataframe with the reshaped values dataframe
    merged_df <- merge(bathy_grid, value_df, by = c("gridY", "gridX"))
    merged_df <- merged_df[, c("latitude", "longitude", "value")]
    
    # Remove rows with missing values
    merged_df <- na.omit(merged_df)
    
    metric= dimnames(merged_array)[[3]][x]
    write.csv(merged_df, paste(output_path, paste("temperature", year, season, metric, threshold_temp, depth,  ".csv", sep="_"), sep="/"))
    
  }
  
  # Step 9: Optionally, convert the result matrices to data frames for easier visualizatio
  # # Check if there are cells with values above the threshold =====================
  # if(which(!is.na(cumulated_degrees_df$cumulated_degrees) & cumulated_degrees_df$cumulated_degrees >= 0) > 1){
  #   cumulated_degrees_df <- data.frame(latitude = bathy_grid$latitude, longitude =  bathy_grid$longitude, cumulated_degrees = as.vector(cumulated_degrees))
  # }else{
  #   print("there is no values over the threshold")
  # }
  # 
  if(which(!is.na(num_hours_above_threshold_df$num_hours_above_threshold) & num_hours_above_threshold_df$num_hours_above_threshold >= 0) > 1){
    num_hours_above_threshold_df <- data.frame(bathy_grid$latitude, longitude =  bathy_grid$longitude, num_hours_above_threshold = as.vector(num_hours_above_threshold))
    
    return(list(num_hours_above_threshold_df, cumulated_degrees_df))
    
  }else{
    print("there is no values over the threshold")
  }
}


# Assuming cumulated_degrees_df has columns: longitude, latitude, cumulated_degrees
# Ensure no missing values and that the data is in the correct format
# cumulated_degrees_df <- na.omit(cumulated_degrees_df)
# num_hours_above_threshold_df<- na.omit(num_hours_above_threshold_df)



# output_path= "/SDM"
# write.csv(num_hours_above_threshold_df, paste(output_path, "spring_2021_hours_above18.csv", sep="/"))



# Ensure valid coordinate ranges
if (any(cumulated_degrees_df$latitude < -90 | cumulated_degrees_df$latitude > 90)) {
  stop("Latitude values out of bounds!")
}
if (any(cumulated_degrees_df$longitude < -180 | cumulated_degrees_df$longitude > 180)) {
  stop("Longitude values out of bounds!")
}

## Convert to a SpatialPointsDataFrame
# coordinates(cumulated_degrees_df) <- ~longitude + latitude
# 
# # Check the coordinate reference system (CRS)
# # If needed, define the CRS (WGS84 in this example)
# crs(cumulated_degrees_df) <- CRS("+proj=longlat +datum=WGS84")
# 
# # Convert the SpatialPointsDataFrame to a raster
# raster_cumulated_degrees <- rasterFromXYZ(cumulated_degrees_df[, c("longitude", "latitude", "cumulated_degrees")])
# 
# # Plot the raster
# plot(raster_cumulated_degrees, main = "Cumulated Temperature Above 15°C", col = terrain.colors(100))
# 
# # Optionally, use mapview for interactive map
# library(mapview)
# mapview(raster_cumulated_degrees, zcol = "cumulated_degrees", 
#         layer.name = "Cumulated Temperature Above 15°C", legend = TRUE)





## ===================================================================
## Convert month metrics into seasonal metrics
## ===================================================================
### Add coordinates to the metrics =============================================
path_data<- "path_name"

bathy_grid<-  read.csv(paste(path_data, "modeled_variables_original/ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/"))
head(bathy_grid)
bathy_grid<- bathy_grid[-1,]
bathy_grid<- bathy_grid[,c("gridY", "gridX", "latitude", "longitude")]


# data<- read.csv(paste(path_data, "modeled_variables_original/Temp_Sal_hourly_2020_2025_ubcSSg3DPhysicsFields1moV21-11_65e9_c913_7128.csv", sep="/"))
# head(data)
# data<- data[-1,]
# colnames(data)
data<- data[,c("time", "gridY", "gridX", "salinity", "temperature")]
data<- merge(data, bathy_grid, by= c("gridY", "gridX"))
head(data)
rm(bathy_grid)
colnames(data)[3]<- "date"
colnames(data)[6:7]<- c("lat", "lon") 
data$lon<- as.numeric(data$lon)
data$lat<- as.numeric(data$lat)
data$temperature<- as.numeric(data$temperature)
data$salinity<- as.numeric(data$salinity)
str(data)
data<- data%>% filter(!is.na(salinity))


# ### Create an array with the csv table ===================================
# # Get unique latitudes, longitudes, and dates
# unique_lats <- sort(unique(data$lat))
# unique_lons <- sort(unique(data$lon))
# unique_dates <- sort(unique(data$date))
# 
# # Create an empty array to store the temperature data
# # The dimensions of the array will be (date, lat, lon)
# temperature_array <- array(NA, dim = c(length(unique_dates), length(unique_lats), length(unique_lons)))
# 
# # Map lat, lon, and temperature to the array for each date
# for(i in 1:nrow(data)) {
#   date_idx <- which(unique_dates == data$date[i])
#   lat_idx <- which(unique_lats == data$lat[i])
#   lon_idx <- which(unique_lons == data$lon[i])
#   
#   # Assign temperature to the corresponding position in the array
#   temperature_array[date_idx, lat_idx, lon_idx] <- data$temperature[i]
# }


# Convert the date column to Date object and extract month
data$date <- as.POSIXct(data$date, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
data$month <- month(data$date)
data$year <- year(data$date)

# Create a 'season' column based on the month
data$season <- case_when(
  data$month %in% c(3, 4, 5) ~ "Spring",
  data$month %in% c(6, 7, 8) ~ "Summer",
  data$month %in% c(9, 10, 11) ~ "Autumn",
  data$month %in% c(12, 1, 2) ~ "Winter",
  TRUE ~ "Unknown"
)

# Summarize temperature and salinity by season, lat, and lon
summary_data_spring <- data %>%
  filter(year != 2025)%>%
  filter(temperature != 0)%>%
  filter(season == "Spring")%>%
  group_by(lat, lon) %>%
  summarize(
    avg_Spring_SST = mean(temperature, na.rm = TRUE),
    min_Spring_SST = min(temperature, na.rm = TRUE),
    max_Spring_SST = max(temperature, na.rm = TRUE),
    sd_Spring_SST = sd(temperature, na.rm = TRUE),
    avg_Spring_Surf_Sal = mean(salinity, na.rm = TRUE),
    max_Spring_Surf_Sal = max(salinity, na.rm = TRUE),
    min_Spring_Surf_Sal = min(salinity, na.rm = TRUE),
    sd_Spring_Surf_Sal = sd(salinity, na.rm = TRUE)
  )

# Summarize temperature and salinity by season, lat, and lon
summary_data_summer <- data %>%
  filter(year != 2025)%>%
  filter(temperature != 0)%>%
  filter(season == "Summer")%>%
  group_by(lat, lon) %>%
  summarize(
    avg_Summer_SST = mean(temperature, na.rm = TRUE),
    min_Summer_SST = min(temperature, na.rm = TRUE),
    max_Summer_SST = max(temperature, na.rm = TRUE),
    sd_Summer_SST = sd(temperature, na.rm = TRUE),
    avg_Summer_Surf_Sal = mean(salinity, na.rm = TRUE),
    max_Summer_Surf_Sal = max(salinity, na.rm = TRUE),
    min_Summer_Surf_Sal = min(salinity, na.rm = TRUE),
    sd_Summer_Surf_Sal = sd(salinity, na.rm = TRUE)
  )

# Summarize temperature and salinity by season, lat, and lon
summary_data_fall <- data %>%
  filter(year != 2025)%>%
  filter(temperature != 0)%>%
  filter(season == "Autumn")%>%
  group_by(lat, lon) %>%
  summarize(
    avg_Fall_SST = mean(temperature, na.rm = TRUE),
    min_Fall_SST = min(temperature, na.rm = TRUE),
    max_Fall_SST = max(temperature, na.rm = TRUE),
    sd_Fall_SST = sd(temperature, na.rm = TRUE),
    avg_Fall_Surf_Sal = mean(salinity, na.rm = TRUE),
    max_Fall_Surf_Sal = max(salinity, na.rm = TRUE),
    min_Fall_Surf_Sal = min(salinity, na.rm = TRUE),
    sd_Fall_Surf_Sal = sd(salinity, na.rm = TRUE)
  )

# Summarize temperature and salinity by season, lat, and lon
summary_data_winter <- data %>%
  filter(year != 2025)%>%
  filter(temperature != 0)%>%
  filter(season == "Winter")%>%
  group_by(lat, lon) %>%
  summarize(
    avg_Winter_SST = mean(temperature, na.rm = TRUE),
    min_Winter_SST = min(temperature, na.rm = TRUE),
    max_Winter_SST = max(temperature, na.rm = TRUE),
    sd_Winter_SST = sd(temperature, na.rm = TRUE),
    avg_Winter_Surf_Sal = mean(salinity, na.rm = TRUE),
    max_Winter_Surf_Sal = max(salinity, na.rm = TRUE),
    min_Winter_Surf_Sal = min(salinity, na.rm = TRUE),
    sd_Winter_Surf_Sal = sd(salinity, na.rm = TRUE)
  )

# Print the summarized data
head(summary_data_winter)
# write.csv(summary_data_spring, paste(path_data, "/modeled_variables_original/seasonal_metrics/", "nemo_Spring_S&T_Jan2020_Jan2024.csv", sep="/"))
# write.csv(summary_data_summer, paste(path_data, "/modeled_variables_original/seasonal_metrics/", "nemo_Summer_S&T_Jan2020_Jan2024.csv", sep="/"))
# write.csv(summary_data_fall, paste(path_data, "/modeled_variables_original/seasonal_metrics/", "nemo_Fall_S&T_Jan2020_Jan2024.csv", sep="/"))
# write.csv(summary_data_winter, paste(path_data, "/modeled_variables_original/seasonal_metrics/", "nemo_Winter_S&T_Jan2020_Jan2024.csv", sep="/"))

# Print the summarized data
write.csv(summary_data, paste0(path_data, "/modeled_variables_original/seasonal_metrics/", paste(substr(unique(data_i$date), 1, 10), "nemo_surface_S&T_Jan2020_Jan2024.csv", sep="_")))


# Summarize temperature and salinity by season, lat, and lon
summary_data <- data %>%
  filter(year != 2025)%>%
  filter(temperature != 0)%>%
  group_by(lat, lon, season) %>%
  summarize(
    avg_temperature = mean(temperature, na.rm = TRUE),
    avg_salinity = mean(salinity, na.rm = TRUE)
  )

# Print the summarized data
head(summary_data)
# write.csv(summary_data, paste0(path_data, "/modeled_variables_original/seasonal_metrics/", paste(substr(unique(data_i$date), 1, 10), "nemo_surface_S&T_Jan2020_Jan2025.csv", sep="_")))




### PLOT MAPS with ggplot2 =====================================================
library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)

# Assuming you already have the summarized data `summary_data` with lat, lon, season, etc
# Convert summary data to an sf object
summary_sf <- st_as_sf(summary_data, coords = c("lon", "lat"), crs = 4326) #2369

# Download the world map
world_map <- ne_countries(scale = "medium", returnclass = "sf")

# Define area limits
min_lon = min(data$lon)
max_lon = max(data$lon)
min_lat = min(data$lat)
max_lat = max(data$lat)

# # Plot temperature or salinity on the map by season
# ggplot() +
#   geom_sf(data = world_map, fill = "lightgray") +  # Base world map
#   geom_sf(data = summary_sf, aes(color = avg_temperature), size = 1) +  # Plot temperature as points
#   scale_color_viridis_c(name = "Avg Temperature (°C)") +  # Add color scale
#   theme_minimal() +  # Clean theme
#   coord_sf(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat), expand = FALSE)+
#   labs(title = "Average Temperature by Location and Season",
#        subtitle = "Data points represent temperature for each location and season",
#        x = "Longitude", y = "Latitude") +
#   facet_wrap(~ season)  # Facet by season





### Rasterize the seasonal variable metrics ====================================

for (i in 2:length(i_time)) {
  data_i<- data[which(data$date == i_time[i]),]
  
  data_i<- data_i[,-c(1:2)]
  data_i<- data_i%>% filter(temperature != 0)
  # write.csv(data_i, paste0(path_data, "/modeled_variables_original/monthly_csvs/", paste(substr(unique(data_i$date), 1, 10), "nemo_surface_S&T.csv", sep="_")))
  
  
  # library(terra)
  # # Create a SpatVector for the points
  # points <- vect(data_i, geom = c("lon", "lat"))
  # 
  # # Define the raster
  # r <- rast(ext(min(data_i$lon), max(data_i$lon), min(data_i$lat), max(data_i$lat)), res = 0.005)
  # 
  # # Rasterize
  # rasterized <- rasterize(points, r, field = "temperature", fun = "mean")
  # 
  # # Plot the result
  # plot(rasterized)
  # 
  # raster_stack<- stack(raster_stack, rasterized)
  # 
  # # Save the raster
  # # writeRaster(rasterized, "raster_output.tif", filetype = "GTiff")
  
}

