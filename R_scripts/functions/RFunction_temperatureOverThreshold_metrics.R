# Steps:
# Set the threshold (e.g., 15°C).
# Filter the temperature values above the threshold (15°C).
# Calculate the cumulated degrees: This is the sum of the temperatures that are above the threshold for each location.
# Count the number of hours: This is the number of time steps where the temperature exceeds the threshold for each location.


## Create an exmple to text the function
# Define the number of matrices (slices)
num_matrices <- 3  

# Create values for the matrices (cycling through 11:15)
matrix_values <- rep(11:15, length.out = 20)  # 20 values per matrix

# Create a 3D array with dimensions (4 rows, 5 columns, 3 matrices)
my_array <- array(matrix_values, dim = c(4, 5, num_matrices))

# Print the array
print(my_array)



# Assuming 'temperature_combined' is the 3D array (time x lat x lon)
# With dimensions (time, latitude, longitude)

bathy_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/")

calculate_DegreeHours_Hours_OverThreshold<- function(data= merged_spring_data, threshold_temp= 12, 
                                                     bathy_file=bathy_file, depth= "surface", season= season_s){
  
  # Step 1: Set the threshold temperature (e.g., 15°C)
  # threshold_temp <- 15
  
  gridX<- data$gridX
  gridY<- data$gridY 
  data<- data$merged_data
  
  
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
    
    # Set dimension names
    dimnames(merged_array) <- list(as.character(gridX),  # Row names (Latitude)
                                   as.character(gridY),  # Column names (Longitude)
                                   c(paste("cumulated_degrees", threshold_temp, sep="_"), 
                                     paste("hours_above_threshold", threshold_temp, sep="_")))
    
  }
  
  # Step 8: add lat and lon coordinates
  for (x in 1:dim(merged_array)[3]) {
    # Convert the matrix into a dataframe with row and column names
    value_df <- as.data.frame(as.table(merged_array[,,x]))
    
    # Rename columns to make it more intuitive
    colnames(value_df) <- c( "gridX", "gridY","value")
    length(unique(value_df$gridX))
    
    bathy_grid<-  read.csv(bathy_file)
    bathy_grid<- bathy_grid[-1,]
    
    # Merge the latitude-longitude dataframe with the reshaped values dataframe
    merged_df <- merge(bathy_grid, value_df, by = c("gridY", "gridX"))
    
    # Remove rows outside the area (where bathymetry is NaN)
    merged_df<- merged_df%>% filter(!is.nan(bathymetry))
    merged_df <- merged_df[, c("latitude", "longitude", "value")]
    
    # Rename the column of metric with specific threshold used
    metric=   dimnames(merged_array)[[3]][x]
    colnames(merged_df)<- c("latitude", "longitude", metric)
    
    
    if(length(!is.na(merged_df[,metric])) >=1) {
      
      output_name= paste("temperature", year, season, metric, depth,  ".csv", sep="_")
      write.csv(merged_df, paste(output_path, output_name, sep="/"))
      
    }
  }
  
}
