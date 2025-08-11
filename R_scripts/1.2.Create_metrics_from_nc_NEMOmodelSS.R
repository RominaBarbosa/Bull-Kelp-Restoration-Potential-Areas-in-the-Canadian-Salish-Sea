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
# input_path= "Y:/Romina/Monthly_nc"
# input_path= "Y:/Romina/Montly_nc_midwater5m"
input_path= "D:/PSF/modeled_variables_original/Monthly_nc"

#input_path= "/Volumes/Romina_PSF/PSF/modeled_variables_original/Monthly_nc"
setwd(input_path)
dir()

# output_path=  "F:/PSF/modeled_variables_original/seasonal_metrics_0.5m"
output_path= "Volumes/Romina_PSF/modeled_variables_original/seasonal_metrics_0.5m_depth"

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
season<- c( "fall", "spring", "winter", "summer") #
  
# "finished 2014 ; dissolved_inorganic_carbon 2025-03-24 18:06:33.594026"
# years<- c("2022")
# variables<- c("temperature")
# season<- c("spring")


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


# To process: 
# seasonal metrics of SST 2022, 2023, 2024
# metrics for turbidity 2015
























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











# Steps:
# Set the threshold (e.g., 15°C).
# Filter the temperature values above the threshold (15°C).
# Calculate the cumulated degrees: This is the sum of the temperatures that are above the threshold for each location.
# Count the number of hours: This is the number of time steps where the temperature exceeds the threshold for each location.


# Assuming 'temperature_combined' is the 3D array (time x lat x lon)
# With dimensions (time, latitude, longitude)

bathy_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/")
calculate_DegreeHours_Hours_OverThreshold<- function(data= merged_spring_data$merged_data, threshold_temp= 12, bathy_file=bathy_file){
  
  # Step 1: Set the threshold temperature (e.g., 15°C)
  # threshold_temp <- 15
  
  
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





output_path= "C:/Users/romi_/OneDrive - University of Victoria/Romina_personal/PSF/SDM"
write.csv(cumulated_degrees_df, paste(output_path, "spring_2021_degreeHours15.csv", sep="/"))
write.csv(num_hours_above_threshold_df, paste(output_path, "spring_2021_hours_above15.csv", sep="/"))






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

































path_data<- "E:/Pacific_Salmon_Fundation"

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



# Load required libraries
library(dplyr)
library(lubridate)

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


### Separate in months to calculate monthly metrics ============================

start <- "202--01-15 00:00:00"
end <- "2025-01-15 00:00:00"
start.time.num <- as.numeric(as.chron(start))

# +1 means one month.  Use +12 if you want one year.
end.time.num <- as.numeric(as.chron(paste(mondate(start)+1, start)))

# 1/24 means one hour.  Change as needed.
hours <- as.chron(seq(start.time.num, end.time.num, 1/24))
as.chron(seq(from=start, by=interval*60, to=end))

set.seed(1);
R <- dim(no3.array)[1];
C <- dim(no3.array)[2];
Z <- dim(no3.array)[3];
N <- R*C*Z;

# end<-  ymd("2023-05-01")


# mean_metric <- tapply(no3.array,list(rep(1:R,C*Z),rep(1:C,each=R,times=Z),rep(strftime(dates,'%m'),each=R*C)),mean)

df<- data.frame(matrix(nrow = 16, ncol = 2))
colnames(df)<- c("month", "hours_count")
df$month<- c(seq(1, 12, 1), seq(1,4,1))
df$day<- rep(1, 16)
df$year<- c(rep(2022,12), rep(2023, 4))
df$start_time<- "00:00"
df$start<- ymd_hm(paste(df$year, df$month, df$day, df$start_time))
df$first_hour_num<- NA
df[1,"first_hour_num"]<- 1
df$last_hour_num<- NA
df[1,"first_hour_num"]<- 1

# no3_metrics<- array(data= NA, dim = c(dim(no3.array)[1], dim(no3.array)[2], 12))
# no3_metrics<- list()

start <- as.POSIXct("2022-01-01")
end <- as.POSIXct("2023-05-01")

time_hour<- seq(from=start, by=interval*60, to=end)

df$start_time <- c(which(time_hour == "2022-01-01 00:00:00 UTC"), which(time_hour == "2022-02-01 00:00:00 UTC"),
                   which(time_hour == "2022-03-01 00:00:00 UTC"), which(time_hour == "2022-04-01 00:00:00 UTC"),
                   which(time_hour == "2022-05-01 00:00:00 UTC"), which(time_hour == "2022-06-01 00:00:00 UTC"),
                   which(time_hour == "2022-07-01 00:00:00 UTC"), which(time_hour == "2022-08-01 00:00:00 UTC"),
                   which(time_hour == "2022-09-01 00:00:00 UTC"), which(time_hour == "2022-10-01 00:00:00 UTC"),
                   which(time_hour == "2022-11-01 00:00:00 UTC"), which(time_hour == "2022-12-01 00:00:00 UTC"),
                   which(time_hour == "2023-01-01 00:00:00 UTC"), which(time_hour == "2023-02-01 00:00:00 UTC"),
                   which(time_hour == "2023-03-01 00:00:00 UTC"), which(time_hour == "2023-03-01 00:00:00 UTC"))

df$end_time <- c(which(time_hour == "2022-02-01 00:00:00 UTC"),
                 which(time_hour == "2022-03-01 00:00:00 UTC"), which(time_hour == "2022-04-01 00:00:00 UTC"),
                 which(time_hour == "2022-05-01 00:00:00 UTC"), which(time_hour == "2022-06-01 00:00:00 UTC"),
                 which(time_hour == "2022-07-01 00:00:00 UTC"), which(time_hour == "2022-08-01 00:00:00 UTC"),
                 which(time_hour == "2022-09-01 00:00:00 UTC"), which(time_hour == "2022-10-01 00:00:00 UTC"),
                 which(time_hour == "2022-11-01 00:00:00 UTC"), which(time_hour == "2022-12-01 00:00:00 UTC"),
                 which(time_hour == "2023-01-01 00:00:00 UTC"), which(time_hour == "2023-02-01 00:00:00 UTC"),
                 which(time_hour == "2023-03-01 00:00:00 UTC"), which(time_hour == "2023-03-01 00:00:00 UTC"),
                 which(time_hour == "2023-05-01 00:00:00 UTC"))

calc_monthly_metrics<- function(data= no3.array, metric= "max", df= df){
  
  no3_month_1<- no3.array[,, 1:  df$end_time[1]]
  # save this data in a raster. Note that we provide the coordinate reference system “CRS” in the standard well-known text format. For this data set, it is the common WGS84 system.
  # We will need to transpose and flip to orient the data correctly. The best way to figure this out is through trial and error, but remember that most netCDF files record spatial data from the bottom left corner.
  no3_month_m1<- flip(raster(t(as.matrix(apply(no3_month_1,c(1,2),metric))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
  names(no3_month_m1)<- paste("NO3", df$month[1], df$year[1], metric, sep="_")
  # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
  # colnames(csv_file)<- c("NO3", "lat", "lon")
  # csv_file$NO3<- c(as.matrix(apply(no3_month_1,c(1,2),metric)))
  # # csv_file[is.na(csv_file$NO3), "NO3"]<- 9999
  # csv_file$lon<- c(lon)
  # csv_file$lat<- c(lat)  
  # csv_file$NO3<- as.numeric(csv_file$NO3)
  # csv_file<- csv_file%>%
  #   filter(!is.na(NO3))
  # name<- paste(names(no3_month_m1), ".csv", sep="")
  # write.csv(csv_file, paste(output_path, name, sep="/"))
  
  
  no3_month_2<- no3.array[,, 2:  df$end_time[2]]
  no3_month_m2<- flip(raster(t(as.matrix(apply(no3_month_2,c(1,2),metric))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
  names(no3_month_m2)<- paste("NO3", df$month[2], df$year[1], metric, sep="_")
  # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
  # colnames(csv_file)<- c("NO3", "lat", "lon")
  # csv_file$NO3<- c(as.matrix(apply(no3_month_2,c(1,2),metric)))
  # # csv_file[is.na(csv_file$NO3), "NO3"]<- 9999
  # csv_file$lon<- c(lon)
  # csv_file$lat<- c(lat)  
  # csv_file$NO3<- as.numeric(csv_file$NO3)
  # csv_file<- csv_file%>%
  #   filter(!is.na(NO3))
  # name<- paste(names(no3_month_m2), ".csv", sep="")
  # write.csv(csv_file, paste(output_path, name, sep="/"))
  
  no3_metrics<- stack(no3_month_m1, no3_month_m2)
  
  for (m in 3:16) {
    df_m<- df[m,]
    start<- df_m$start
    # end.time.num <- as.numeric(as.chron(paste(mondate(start)+1, start)))# +1 means one month. 
    end.time <- mondate(start)+1# +1 means one month. 
    n_hours<- length(seq(from=start, by=interval*60, to= as.POSIXct(end.time)-2))
    
    
    no3_month<- no3.array[,, df_m$start_time:  df_m$end_time]
    no3_month_m<- flip(raster(t(as.matrix(apply(no3_month,c(1,2),paste(metric)))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
    names(no3_month_m)<- paste("NO3", df$month[m], df$year[m], metric, sep="_")
    no3_metrics<- stack(no3_metrics, no3_month_m)
    
    ## no3<- as.data.frame(c(apply(no3_month,c(1,2),paste(metric))))
    # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
    # colnames(csv_file)<- c("NO3", "lat", "lon")
    # csv_file$NO3<- c(apply(no3_month, c(1, 2), paste(metric)))
    # csv_file$lon<- c(lon)
    # csv_file$lat<- c(lat)  
    # 
    # csv_file$NO3<- as.numeric(csv_file$NO3)
    # csv_file<- csv_file%>%
    #   filter(!is.na(NO3))
    # name<- paste(names(no3_month_m), ".csv", sep="")
    # write.csv(csv_file, paste(output_path, name, sep="/"))
    
  }
  return(no3_metrics)
}

plot(no3_metrics)

max_metrics<- calc_monthly_metrics(data= no3.array, metric= "max", df= df)
min_metrics<- calc_monthly_metrics(data= no3.array, metric= "min", df= df)
mean_metrics<- calc_monthly_metrics(data= no3.array, metric= "mean", df= df)
sd_metrics<- calc_monthly_metrics(data= no3.array, metric= "sd", df= df)

plot(max_metrics[[1:12]])#
plot(min_metrics[[1:12]])
plot(mean_metrics[[1:12]])
plot(sd_metrics[[1:12]])#NO3_sd_2022

plot(max_metrics[[13:16]])#
plot(min_metrics[[13:16]])
plot(mean_metrics[[13:16]])
plot(sd_metrics[[13:16]])#NO3_sd_2022


output_path<- "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/metrics_NO3"
setwd(paste(output_path))
# terra::writeRaster(max_metrics, "max_metrics.tif", overwrite=TRUE)

stack_NO3<- stack(max_metrics,min_metrics,mean_metrics,sd_metrics)
terra::writeRaster(stack_NO3, "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/SDMs/layers/processed_layers/NO3_stacklayers2.tif", overwrite=TRUE)


# # ============================================================================##
## Plot metric function ========================================================

plot_metrics<- function(metrics_stack= max_metrics){
  # obtain max and min values
  # maxv <- max(maxValue(metrics_stack), na.rm = T)+1
  # minv <- min(minValue(metrics_stack), na.rm = T)
  maxv <- 25
  minv <- 0
  
  # set the breaks between min and max values
  brks <- seq(minv,maxv,by=0.1)
  nbrks <- length(brks)-1
  r.range <- c(minv, maxv)
  
  # generate palette
  colfunc<-colorRampPalette(c("springgreen", "royalblue", "yellow", "red"))
  
  # par(mfrow=c(4,3))
  for(i in 1:12) {# seq_len(nlayers(metrics_stack))){
    tmp <- metrics_stack[[i]]
    plot(tmp, breaks=brks,col=colfunc(nbrks), legend = F, zlim=c(minv,maxv),
         main = names(tmp)) 
    plot(metrics_stack[[i]], legend.only=TRUE, col=colfunc(nbrks),
         legend.width=1, legend.shrink=0.75,
         legend.args=list(text='value', side=4, font=2, line=2.5, cex=0.8))
  }
  
  
}





################################################################################
## Process second timeseries to calculate monthly metrics 2023 =================
## barbosa0_surf_2023.05.10_2023.12.31
nc_data <- nc_open('C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/barbosa0_surf_2023.05.10_2023.12.31.nc')

lon <- ncvar_get(nc_data, "lon_rho")
lat <- ncvar_get(nc_data, "lat_rho", verbose = F)
t <- ncvar_get(nc_data, "ocean_time")
length(t)# [1] 5665


# calculate duration to check the time of the simulation dataset
stdate<- ymd("2023-05-10")

difftime(ymd_h("2024-01-01 01"), ymd_h("2023-05-10 00"), units = "hours")# Time difference of 5665 hours

start <- as.POSIXct("2023-05-10")
interval <- 60
end <- start + as.difftime(length(t), units="hours")

# Read in the data from the NO3 variable and verify the dimensions of the array.
no3.array23 <- ncvar_get(nc_data, "NO3") # store the data in a 3-dimensional array
dim(no3.array23) 

fillvalue <- ncatt_get(nc_data, "NO3", "_FillValue")
fillvalue


# All done reading in the data. We can close the netCDF file.
nc_close(nc_data)

# First, a little housekeeping. Let’s replace all those pesky fill values with the R-standard ‘NA’.
no3.array23[no3.array23 == fillvalue$value] <- NA






################################################################################
### Separate in months to calculate monthly metrics ============================
library(chron)
library(mondate)
start <- "2023-05-10 00:00:00"
end <- "2024-01-31 01:00:00"
# start.time.num <- as.numeric(as.chron(start))

# 1/24 means one hour.  Change as needed.
# hours <- as.chron(seq(start.time.num, end.time.num, 1/24))
# as.chron(seq(from=start, by=interval*60, to=end))

set.seed(1);
R <- dim(no3.array23)[1];
C <- dim(no3.array23)[2];
Z <- dim(no3.array23)[3];
N <- R*C*Z;

df<- data.frame(matrix(nrow = 7, ncol = 2))
colnames(df)<- c("month", "hours_count")
df$month<- c(seq(6, 12, 1))
df$day<- rep(1, 7)
df$year<- c(rep(2023,7))
df$start_time<- "00:00"
df$start<- ymd_hm(paste(df$year, df$month, df$day, df$start_time))

start <- as.POSIXct("2023-05-10")
end <- as.POSIXct("2024-01-01")

time_hour<- seq(from=start, by=interval*60, to=end)

df$start_time <- c(which(time_hour == "2023-06-01 00:00:00 UTC"),
                   which(time_hour == "2023-07-01 00:00:00 UTC"), which(time_hour == "2023-08-01 00:00:00 UTC"),
                   which(time_hour == "2023-09-01 00:00:00 UTC"), which(time_hour == "2023-10-01 00:00:00 UTC"),
                   which(time_hour == "2023-11-01 00:00:00 UTC"), which(time_hour == "2023-12-01 00:00:00 UTC"))

df$end_time <- c(which(time_hour == "2023-07-01 00:00:00 UTC"), which(time_hour == "2023-08-01 00:00:00 UTC"),
                 which(time_hour == "2023-09-01 00:00:00 UTC"), which(time_hour == "2023-10-01 00:00:00 UTC"),
                 which(time_hour == "2023-11-01 00:00:00 UTC"), which(time_hour == "2023-12-01 00:00:00 UTC"),
                 which(time_hour == "2023-12-31 00:00:00 UTC"))

df$end_date <- time_hour[df$end_time]


csv_outputs_path<- "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/metrics_NO3"

calc_monthly_metrics2023<- function(data= no3.array, metric= "max", df= df, csv_outputs_path= csv_outputs_path){
  
  no3_month_1<- no3.array[,, 1:  df$end_time[1]]
  # save this data in a raster. Note that we provide the coordinate reference system “CRS” in the standard well-known text format. For this data set, it is the common WGS84 system.
  # We will need to transpose and flip to orient the data correctly. The best way to figure this out is through trial and error, but remember that most netCDF files record spatial data from the bottom left corner.
  no3_month_m1<- flip(raster(t(as.matrix(apply(no3_month_1,c(1,2),metric))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
  names(no3_month_m1)<- paste("NO3", df$month[1], df$year[1], metric, sep="_")
  # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
  # colnames(csv_file)<- c("NO3", "lat", "lon")
  # csv_file$NO3<- c(as.matrix(apply(no3_month_1,c(1,2),metric)))
  # # csv_file[is.na(csv_file$NO3), "NO3"]<- 9999
  # csv_file$lon<- c(lon)
  # csv_file$lat<- c(lat)
  # csv_file$NO3<- as.numeric(csv_file$NO3)
  # csv_file<- csv_file%>%
  #   filter(!is.na(NO3))
  # name<- paste(names(no3_month_m1), ".csv", sep="")
  # write.csv(csv_file, paste(csv_outputs_path, name, sep="/"))
  
  
  no3_month_2<- no3.array[,, 2:  df$end_time[2]]
  no3_month_m2<- flip(raster(t(as.matrix(apply(no3_month_2,c(1,2),metric))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
  names(no3_month_m2)<- paste("NO3", df$month[2], df$year[1], metric, sep="_")
  # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
  # colnames(csv_file)<- c("NO3", "lat", "lon")
  # csv_file$NO3<- c(as.matrix(apply(no3_month_2,c(1,2),metric)))
  # # csv_file[is.na(csv_file$NO3), "NO3"]<- 9999
  # csv_file$lon<- c(lon)
  # csv_file$lat<- c(lat)
  # csv_file$NO3<- as.numeric(csv_file$NO3)
  # csv_file<- csv_file%>%
  #   filter(!is.na(NO3))
  # name<- paste(names(no3_month_m2), ".csv", sep="")
  # write.csv(csv_file, paste(csv_outputs_path, name, sep="/"))
  
  no3_metrics<- stack(no3_month_m1, no3_month_m2)
  
  if (length(df[,1]) >= 3){
    
    for (m in 3:7) {
      df_m<- df[m,]
      start<- df_m$start
      # end.time.num <- as.numeric(as.chron(paste(mondate(start)+1, start)))# +1 means one month. 
      end.time <- mondate(start)+1# +1 means one month. 
      n_hours<- length(seq(from=start, by=interval*60, to= as.POSIXct(end.time)-2))
      
      
      no3_month<- no3.array[,, df_m$start_time:  df_m$end_time]
      no3_month_m<- flip(raster(t(as.matrix(apply(no3_month,c(1,2),paste(metric)))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
      names(no3_month_m)<- paste("NO3", df$month[m], df$year[m], metric, sep="_")
      no3_metrics<- stack(no3_metrics, no3_month_m)
      
      # # no3<- as.data.frame(c(apply(no3_month,c(1,2),paste(metric))))
      # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
      # colnames(csv_file)<- c("NO3", "lat", "lon")
      # csv_file$NO3<- c(apply(no3_month, c(1, 2), paste(metric)))
      # csv_file$lon<- c(lon)
      # csv_file$lat<- c(lat)
      # 
      # csv_file$NO3<- as.numeric(csv_file$NO3)
      # csv_file<- csv_file%>%
      #   filter(!is.na(NO3))
      # name<- paste(names(no3_month_m), ".csv", sep="")
      # write.csv(csv_file, paste(csv_outputs_path, name, sep="/"))
    }
  }
  
  return(no3_metrics)
}



max_metrics23<- calc_monthly_metrics2023(data= no3.array23, metric= "max", df= df)
min_metrics23<- calc_monthly_metrics2023(data= no3.array23, metric= "min", df= df)
mean_metrics23<- calc_monthly_metrics2023(data= no3.array23, metric= "mean", df= df)
sd_metrics23<- calc_monthly_metrics2023(data= no3.array23, metric= "sd", df= df)

par(mfrow=c(3,4))
plot(stack(max_metrics[[13:16]],max_metrics23))#
plot(stack(min_metrics[[13:16]], min_metrics23))
plot(stack(mean_metrics[[13:16]], mean_metrics23))
plot(stack(sd_metrics[[13:16]], sd_metrics23))#NO3_sd_2022

plot(max_metrics23)#NO3_max_JuneDec2023
plot(min_metrics23)
plot(mean_metrics23)
plot(sd_metrics23)


################################################################################
### Separate in seasons and calculate seasonal metrics =========================

"metrics_NO3_csv_month"

df<- data.frame(matrix(nrow = 2, ncol = 1))
df$month<- c(seq(1, 2, 1))
df$day<- rep(1, 2)
df$year<- c(rep(2023,2))
df$start_time<- "00:00"
df$start<- ymd_hm(paste(df$year, df$month, df$day, df$start_time))

start <- as.POSIXct("2023-05-10")
end <- as.POSIXct("2024-01-01")

time_hour<- seq(from=start, by=interval*60, to=end)

df$start_time <- c(which(time_hour == "2023-07-01 00:00:00 UTC"), which(time_hour == "2023-10-01 00:00:00 UTC"))
df$end_time <- c(which(time_hour == "2023-09-01 00:00:00 UTC"), which(time_hour == "2023-12-31 00:00:00 UTC"))

df$end_date <- time_hour[df$end_time]
df$start_date <- time_hour[df$start_time]
df$month<- c("summer", "fall")

df.season<- df


csv_outputs_path<- "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready"

max_metrics23wintfall<- calc_monthly_metrics2023(data= no3.array, metric= "max", df= df, csv_outputs_path= paste(csv_outputs_path, "metrics_NO3_csv_seasons", sep="/"))
min_metrics23wintfall<- calc_monthly_metrics2023(data= no3.array, metric= "min", df= df, csv_outputs_path= paste(csv_outputs_path, "metrics_NO3_csv_seasons", sep="/"))
mean_metrics23wintfall<- calc_monthly_metrics2023(data= no3.array, metric= "mean", df= df, csv_outputs_path= paste(csv_outputs_path, "metrics_NO3_csv_seasons", sep="/"))
sd_metrics23wintfall<- calc_monthly_metrics2023(data= no3.array, metric= "sd", df= df, csv_outputs_path= paste(csv_outputs_path, "metrics_NO3_csv_seasons", sep="/"))



##############################################################################
### Winter and Spring Seasonal metrics from 2023 #############################
nc_data <- nc_open('C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/barbosa0_surf_2022.01.01_2023.05.10.nc')
lon <- ncvar_get(nc_data, "lon_rho")
lat <- ncvar_get(nc_data, "lat_rho", verbose = F)
t <- ncvar_get(nc_data, "ocean_time")

# Read in the data from the NO3 variable and verify the dimensions of the array.
no3.array2022 <- ncvar_get(nc_data, "NO3") # store the data in a 3-dimensional array

fillvalue <- ncatt_get(nc_data, "NO3", "_FillValue")
fillvalue

# All done reading in the data. We can close the netCDF file.
nc_close(nc_data)

# First, a little housekeeping. Let’s replace all those pesky fill values with the R-standard ‘NA’.
no3.array2022[no3.array2022 == fillvalue$value] <- NA

### 2023
nc_data <- nc_open('C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/barbosa0_surf_2023.05.10_2023.12.31.nc')
# lon <- ncvar_get(nc_data, "lon_rho")
# lat <- ncvar_get(nc_data, "lat_rho", verbose = F)
t <- ncvar_get(nc_data, "ocean_time")

# Read in the data from the NO3 variable and verify the dimensions of the array.
no3.array2023 <- ncvar_get(nc_data, "NO3") # store the data in a 3-dimensional array
fillvalue <- ncatt_get(nc_data, "NO3", "_FillValue")
fillvalue

# All done reading in the data. We can close the netCDF file.
nc_close(nc_data)

# First, a little housekeeping. Let’s replace all those pesky fill values with the R-standard ‘NA’.
no3.array2023[no3.array2023 == fillvalue$value] <- NA


# Df 
df<- data.frame(matrix(nrow = 2, ncol = 1))
df$month<- c(seq(1, 2, 1))
df$day<- rep(1, 2)
df$year<- c(rep(2023,2))
df$start_time<- "00:00"
df$start<- ymd_hm(paste(df$year, df$month, df$day, df$start_time))

time_hour22<- seq(from=as.POSIXct("2022-01-01"), by=interval*60, to= as.POSIXct("2023-06-01"))
time_hour23<- seq(from=as.POSIXct("2023-05-10"), by=interval*60, to= as.POSIXct("2023-12-31"))#

df$start_time <- c(which(time_hour22 == "2023-01-01 00:00:00 UTC"))
df$end_time <- c(which(time_hour22 == "2023-04-01 00:00:00 UTC"))




calc_metrics_csv<- function(data= no3.array, metric= "max", csv_outputs_path= csv_outputs_path, season_name= "NO3_winter2023"){
  no3_array= data
  no3_array1<- flip(raster(t(as.matrix(apply(data,c(1,2),metric))), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")), direction='y')
  names(no3_array1)<- paste(season_name, metric, sep="_")
  # csv_file<- as.data.frame(matrix(ncol=3, nrow= length(lat)))
  # colnames(csv_file)<- c("NO3", "lat", "lon")
  # csv_file$NO3<- c(as.matrix(apply(data,c(1,2),metric)))
  # csv_file$lon<- c(lon)
  # csv_file$lat<- c(lat)
  # csv_file$NO3<- as.numeric(csv_file$NO3)
  # csv_file<- csv_file%>%
  #   filter(!is.na(NO3))
  # name<- paste(names(no3_array1), ".csv", sep="")
  # write.csv(csv_file, paste(csv_outputs_path, "metrics_NO3_csv_seasons", name, sep="/"))
  return(no3_array1)
  # return(csv_file)
}

### spring 2023 
no3_spring1<- no3.array2022[,,1:  df$end_time[1]]
no3_spring2<- abind::abind(no3.array2022[,, which(time_hour22 == "2023-04-01 00:00:00 UTC"):  which(time_hour22 == "2023-05-11 00:00:00 UTC")], 
                           no3.array2023[,, which(time_hour23 == "2023-05-11 00:00:00 UTC"):  which(time_hour23 == "2023-07-01 00:00:00 UTC")])

no3_spring<- abind::abind(no3_spring1, no3_spring2)
dim(no3_spring)

spring2023mean<- calc_metrics_csv(data= no3_spring, metric= "mean", csv_outputs_path= csv_outputs_path, season_name= "NO3_spring2023")
plot(spring2023mean)
spring2023min<- calc_metrics_csv(data= no3_spring, metric= "min", csv_outputs_path= csv_outputs_path, season_name= "NO3_spring2023")
spring2023sd<- calc_metrics_csv(data= no3_spring, metric= "sd",  csv_outputs_path= csv_outputs_path, season_name= "NO3_spring2023")
spring2023max<- calc_metrics_csv(data= no3_spring, metric= "max", csv_outputs_path= csv_outputs_path, season_name= "NO3_spring2023")


### Winter 2023 
no3_winter<- no3.array2022[,, which(time_hour22 == "2023-01-01 00:00:00 UTC"):  which(time_hour22 == "2023-04-01 00:00:00 UTC")]
dim(no3_winter)[3]/24/3

winter2023mean<- calc_metrics_csv(data= no3_winter, metric= "mean", csv_outputs_path= csv_outputs_path, season_name= "NO3_winter2023")
winter2023min<- calc_metrics_csv(data= no3_winter, metric= "min", csv_outputs_path= csv_outputs_path, season_name= "NO3_winter2023")
winter2023sd<- calc_metrics_csv(data= no3_winter, metric= "sd", csv_outputs_path= csv_outputs_path, season_name= "NO3_winter2023")
winter2023max<- calc_metrics_csv(data= no3_winter, metric= "max", csv_outputs_path= csv_outputs_path, season_name= "NO3_winter2023")





### Output plots       =========================================================
output_path<- "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/plots_metrics_NO3"

par(mfrow=c(2,2))

# png(paste(output_path, "sd_surfaceNO3.png", sep="/"), units= "cm", height = 20, width = 19, res = 300)
plot_metrics(winter2023mean)
# dev.off()

# png(paste(output_path, "max_surfaceNO3.png", sep="/"), units= "cm",height = 20, width = 19, res = 300)
plot_metrics(winter2023min)
# dev.off()

# png(paste(output_path, "min_surfaceNO3.png", sep="/"), units= "cm",height = 20, width = 19, res = 300)
plot_metrics(winter2023sd)
# dev.off()

# png(paste(output_path, "mean_surfaceNO3.png", sep="/"), units= "cm", height = 20, width = 19, res = 300)
plot_metrics(winter2023max)
# dev.off()


par(mfrow=c(4,4))


plot_metrics(winter2023mean)
plot_metrics(winter2023min)
plot_metrics(winter2023max)
plot_metrics(winter2023sd)

plot_metrics(spring2023mean)
plot_metrics(spring2023min)
plot_metrics(spring2023max)
plot_metrics(spring2023sd)

plot_metrics(mean_metrics23wintfall[[1]])
plot_metrics(min_metrics23wintfall[[1]])
plot_metrics(max_metrics23wintfall[[1]])
plot_metrics(sd_metrics23wintfall[[1]])

plot_metrics(mean_metrics23wintfall[[2]])
plot_metrics(min_metrics23wintfall[[2]])
plot_metrics(max_metrics23wintfall[[2]])
plot_metrics(sd_metrics23wintfall[[2]])



# > save.image("C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/Environmental_cond/Washington_model_MacCready/NO3_metrics_prep.RData")