# NEMO model dataset extraction and merging 


### ============================================================================
### Function to automatically ensemble the datasets (hourly data) per year or season
### ============================================================================
# the function works for different depth layers including: surface = 0.5 m depth; bottom_10 = 10 m depth; bottom_15 = 15 m depth
# There are different variables available including:


# Function to merge data by season (spring, summer, etc.), keeping memory usage low
merge_nemo_seasonal_data <- function(depth = "surface", year = "2015", variable = "temperature", 
                                     input_path = "D:/PSF/modeled_variables_original/Monthly_nc", 
                                     season = "winter") {
  require(lubridate)
  # Convert depth based on input
  # if (depth == "surface") {
  #   depth = 0.5000003
  # }  
  # if (depth == "bottom_10") {
  #   depth = 10.5047655
  # } 
  # if (depth == "bottom_15") {
  #   depth = 15.634288
  # }
  # if (depth == "bottom_5"){
  #   depth = 5.634288
  # }
  
  # Define months based on the season
  if (season == "winter") {
    months <- c("12", "01", "02")
  } else if (season == "spring") {
    months <- c("03", "04", "05")
  } else if (season == "summer") {
    months <- c("06", "07", "08")
  } else if (season == "fall") {
    months <- c("09", "10", "11")
  } else {
    stop("Invalid season. Choose from 'winter', 'spring', 'summer', 'fall'.")
  }
  
  # Initialize a list to hold merged data for the selected season
  merged_data <- NULL
  
  # Loop through the months in the selected season
  for (month in months) {
    
    
    if(season == "winter" & month == "12"){
      # Set initial date 1st day of December from the previous year
      
      init_date <- paste(as.numeric(year)-1, month, "01", sep = "-")
       
      # Use lubridate to calculate the first day of the following month
      end_date <- lubridate::ceiling_date(as.Date(init_date), "month")
      end_date <- format(end_date, "%Y-%m-%d")  # Format as YYYY-MM-DD
      
      
    }else{
      # Set initial date (1st day of the month)
    init_date <- paste(year, month, "01", sep = "-")
    
    # Use lubridate to calculate the first day of the following month
    end_date <- lubridate::ceiling_date(as.Date(init_date), "month")
    end_date <- format(end_date, "%Y-%m-%d")  # Format as YYYY-MM-DD
    
    }
    
    # Define file name based on the variable
    if (depth == "bottom_10"){
      file_name <- paste(input_path, paste(variable, init_date, end_date, "hourly", "bottom_10.nc", sep = "_"), sep = "/")
    }else    {
      file_name <- paste(input_path, paste(variable, init_date, end_date, "hourly.nc", sep = "_"), sep = "/")
    }
    
    # Open the .nc file and read the variable data
    nc_data <- nc_open(file_name)
    data <- ncdf4::ncvar_get(nc_data, variable)  # Read the variable data
    
    
    # get lan and lon data to create array
    gridX <- ncvar_get(nc_data, "gridX")
    gridY <- ncvar_get(nc_data, "gridY", verbose = F)
    
    # Identify the indices for June and July
    t_data <- as.POSIXct(ncvar_get(nc_data, "time"), origin = "1970-01-01", tz = "UTC")
    length(t_data)
    
    # data.frame(lon, lat)
    
    # Check the dimensions of the data
    data_dims <- dim(data)
    print(paste("Data dimensions for month", month, ": ", paste(data_dims, collapse = "x")))
    
    # If data has 3 dimensions, remove the last time step (assuming the 3rd dimension is time)
    if (length(data_dims) == 3) {
      data <- data[, , -data_dims[3]]  # Remove the last time step in the 3rd dimension (time)
    } else {
      # If data has fewer dimensions, print a warning
      warning("Data for month ", month, " has fewer than 3 dimensions. Skipping last row removal.")
    }
    
    # Close the nc file
    nc_close(nc_data)
    rm(nc_data)
    
    # Merge data for the selected season (month by month)
    if (is.null(merged_data)) {
      merged_data <- data  # Initialize merged data with the first month's data
    } else {
      # Check if the dimensions match before merging
      merged_dims <- dim(merged_data)
      if (length(merged_dims) == length(data_dims) && all(merged_dims[1:2] == data_dims[1:2])) {
        # If the first 2 dimensions match, merge the data
        merged_data <- abind(merged_data, data, along = 3)  # Concatenate along the 3rd dimension (time/month)
      } else {
        # Handle dimension mismatch (pad shorter data with NA values)
        cat("Dimension mismatch for month", month, ". Padding data.\n")
        
        # Find the maximum length of the 3rd dimension
        max_time_length <- max(merged_dims[3], data_dims[3])
        
        # Pad the shorter array with NA values to match the maximum time length
        if (merged_dims[3] < max_time_length) {
          pad_length <- max_time_length - merged_dims[3]
          merged_data <- abind(merged_data, array(NA, dim = c(merged_dims[1:2], pad_length)), along = 3)
        }
        if (data_dims[3] < max_time_length) {
          pad_length <- max_time_length - data_dims[3]
          data <- abind(data, array(NA, dim = c(data_dims[1:2], pad_length)), along = 3)
        }
        
        # Now merge the data
        merged_data <- abind(merged_data, data, along = 3)
      }
    }
    
    # cat("Merged data for month:", month, "Season:", season, "\n")
  }
  
  return(list(merged_data= merged_data, time= t_data, gridX= gridX, gridY= gridY))
}
