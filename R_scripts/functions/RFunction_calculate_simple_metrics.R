

calculate_simple_metrics<- function(data = merged_spring_data, 
                                    season= "spring",
                                    year= "2014",
                                    depth = "bottom_10",
                                    variable= variable,
                                    output_path= output_path, 
                                    path_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/"),
                                    bathy_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/")){
  
  data_a<- data
  data<- data_a[[1]]
  # Assuming 'data' is the 3D array ( lat x lon x time )
  # With dimensions (time, latitude, longitude)
  
  # Step 1: Calculate the 0.9 and 0.1 percentile thresholds for each latitude and longitude
  percentile_90 <- apply(data, c(1, 2), function(x) quantile(x, 0.9, na.rm = TRUE))
  percentile_10 <- apply(data, c(1, 2), function(x) quantile(x, 0.1, na.rm = TRUE))
  
  # Step 2: Initialize matrices for the averages (above 0.9 percentile, below 0.1 percentile, and general average)
  avg_above_90 <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
  avg_below_10 <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
  general_avg  <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
  general_sd   <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
  
  # Step 3: Iterate over each latitude and longitude
  for (i in 1:dim(data)[1]) {
    for (j in 1:dim(data)[2]) {
      
      # Step 4: Extract all temperature values for this lat, lon
      location_values <- data[i, j, ]  # Extract values for this lat, lon
      
      # Filter out zero values
      location_values <- location_values[location_values != 0]
      
      
      # Step 5: Calculate the general average and SD (all data)
      general_avg[i, j] <- mean(location_values, na.rm = TRUE)
      general_sd[i, j] <- sd(location_values, na.rm = TRUE)
      
      # Get the 0.9 and 0.1 percentile values for this location
      threshold_90 <- percentile_90[i, j]
      threshold_10 <- percentile_10[i, j]
      
      # Step 6: Filter values above the 0.9 percentile for this location
      values_above_90 <- location_values[location_values >= threshold_90]
      
      # Step 7: Filter values below the 0.1 percentile for this location
      values_below_10 <- location_values[location_values <= threshold_10]
      
      # Step 8: Calculate the average above the 0.9 percentile
      if (length(values_above_90) > 0) {
        avg_above_90[i, j] <- mean(values_above_90, na.rm = TRUE)
      }
      
      # Step 9: Calculate the average below the 0.1 percentile
      if (length(values_below_10) > 0) {
        avg_below_10[i, j] <- mean(values_below_10, na.rm = TRUE)
      }
    }
  }
  
  # Step 4: merge matrix of metrics into array 
  merged_array<- abind(avg_above_90, avg_below_10, general_avg, general_sd,
                       along = 3, 
                       make.names = TRUE)
  
  # Set dimension names
  dimnames(merged_array) <- list(as.character(data_a$gridX),  # Row names (Latitude)
                                 as.character(data_a$gridY),  # Column names (Longitude)
                                 c("maximum", "minimum", "mean", "SD"))
  
  
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
    write.csv(merged_df, paste(output_path, paste(variable, year, season, metric, depth,  ".csv", sep="_"), sep="/"))
    
  }
  
}
