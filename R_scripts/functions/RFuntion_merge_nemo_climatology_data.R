### ============================================================================
### Function to ensemble the datasets (yearly seasonal metrics) per period
### ============================================================================
# the function works for different depth layers including: 
# surface = 0.5 m depth; 
# bottom_10 = 10 m depth
# There are different variables available including:

input_path= "/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth"
# Function to merge data by period of years (blob, etc.), keeping memory usage low
# merge_nemo_climatology_data <- function(depth = "surface", time_period = "blob", 
#                                         variable = "ammonium", 
#                                         metric= "mean",
#                                         input_path = "D:/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth", 
#                                         season = "winter") {
#   require(lubridate)
#   # Define years based on the time period
#   
#   if (time_period == "blob") {
#     years <- c( "2014","2015", "2016", "2017", "2018", "2019")#"2014", fro wind without data
#   } else if (time_period == "post_blob") {
#     months <- c("2020", "2021", "2022", "2023")
#   } else {
#     stop("Invalid time period; Choose 'blob' or 'post_blob'")
#   }
#   
#   # Initialize a list to hold merged data for the selected season
#   merged_data <- data.frame(matrix(ncol=3, nrow = 0))
#   colnames(merged_data)<- c("latitude", "longitude", "year")
#   
#   # Loop through the months in the selected season
#   for (year in years) {
#     
#     # ammonium_2014_fall_maximum_surface_
#     
#     # Define file name based on the variable
#     file_name <- paste(input_path, paste(variable, year, season, metric, depth,".csv", sep = "_"), sep = "/")
#     data<- read.csv(file_name)
#     
#     if(length(colnames(data))==4){
#       data<- data[,-1]
#     }
#     
#     data$year<- year
#     merged_data <- rbind(merged_data, data)
#     
#   }
#   
#   
#   colnames(merged_data)[3]<- "value"
#   # print("Metric value must be in column 4 of the input table")
#   
#   # Calculate the average of the metric for the entire period
#   merged_data_mean<- merged_data%>%
#     group_by(latitude, longitude)%>%
#     summarize(climatology_ave = mean(value))
#   
#   # write.csv(merged_data_mean, paste(output_path, paste(variable, time_period, season, metric, depth,".csv", sep = "_"), sep = "/"))
#   return(merged_data_mean)
#   
# }



#===============================================================================


merge_nemo_climatology_data <- function(
    depth = "surface",
    time_period = "blob",
    variable = "ammonium",
    metric = "mean",
    input_path = "D:/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth",
    season = "winter",
    grid_template  # df with latitude, longitude
) {
  library(dplyr)
  
  # Years by time period
  if (time_period == "blob") {
    years <- c("2014", "2015", "2016", "2017", "2018", "2019")
  } else if (time_period == "post_blob") {
    years <- c("2020", "2021", "2022") #, "2023"
  } else {
    stop("Invalid time period; Choose 'blob' or 'post_blob'")
  }
  
  # Prepare template
  merged_data <- grid_template %>%
    mutate(
      latitude  = round(as.numeric(latitude), 6),
      longitude = round(as.numeric(longitude), 6)
    )
  
  # Loop through years
  for (year in years) {
    file_name <- file.path(input_path, paste(variable, year, season, metric, depth, ".csv", sep = "_"))
    data <- read.csv(file_name)
    
    if (ncol(data) == 4) {
      data <- data[, -1]
    }
    
    data <- data %>%
      mutate(
        latitude  = round(as.numeric(latitude), 6),
        longitude = round(as.numeric(longitude), 6)
      )
    
    colnames(data)[3] <- paste0("value_", year)
    
    merged_data <- left_join(merged_data, data, by = c("latitude", "longitude"))
  }
  
  # Calculate climatology mean
  value_cols <- grep("^value_", names(merged_data))
  merged_data$climatology_ave <- rowMeans(merged_data[, value_cols], na.rm = TRUE)
  
  return(merged_data)
}

  
  

# How to use the function:
# # grid_template must have:
# # gridX, gridY, latitude, longitude for ALL cells
# grid_template <- read.csv("bathy_template.csv") %>%
#   mutate(gridX = as.integer(gridX), gridY = as.integer(gridY))
# 
# 
# climatology <- merge_nemo_climatology_data(
#   depth = "surface",
#   time_period = "blob",
#   variable = "ammonium",
#   metric = "mean",
#   input_path = "D:/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth",
#   season = "winter",
#   grid_template = grid_template
# )

# merged_sf <- st_as_sf(
#   merged_data %>% filter(!is.na(climatology_ave)),
#   coords = c("longitude", "latitude"),
#   crs = 4326 # WGS84
# )

# Quick plot
# ggplot(merged_sf) +
#   geom_sf(aes(color = climatology_ave), size = 0.5) +
#   scale_color_viridis_c() +
#   theme_minimal()
