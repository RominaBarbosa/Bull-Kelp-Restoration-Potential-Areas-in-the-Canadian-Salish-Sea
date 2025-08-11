###=============================================================================
###     NEMO model - SalishSeaCast data process                 ################
###                                                             ################
### 1.2- Loop to merge monthly data into seasons                ################
###    - and calculate metrics od current (speed and direction)
### https://salishsea.eos.ubc.ca/                               ################
### input data acquired with: "1.1.Download_nc_NEMOmodel.R"     ################
### input data: month data of hourly series                     ################
### output data: metric of a season                             ################
### Author: Romina Barbosa                                      ################
### Date: 31-July-2025                                           ################
### Last edition: 04-August-2025
###=============================================================================


library(ncdf4)
library(abind)
library(lubridate)

my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/Monthly_nc")
setwd(my_path)
dir()

# Helper to split indices evenly
split_indices <- function(n, parts) {
  splits <- floor(seq(0, n, length.out = parts + 1))
  splits <- splits + 1
  splits[parts + 1] <- n
  idx_list <- lapply(1:parts, function(i) splits[i]:(splits[i+1]))
  return(idx_list)
}

# Function that processes date ranges and chunks to compute metrics
calculate_seasonal_metrics_from_date_ranges <- function(date_ranges, season_name, spatial_splits = 3,
                                                        data_path = "/Volumes/Romina_PSF/PSF/modeled_variables_original/Monthly_nc") {
  
  chunk_results <- vector("list", spatial_splits)
  for (i in 1:spatial_splits) {
    chunk_results[[i]] <- NULL
  }
  # Variables to store lon/lat dims - will be set on first chunk
  lon <- NULL
  lat <- NULL
  x_dim <- NULL
  y_dim <- NULL
  y_splits <- NULL

  for (chunk_id in 1:spatial_splits) {
    cat("Processing spatial chunk", chunk_id, "of", spatial_splits, "\n")
    
    # Initialize empty lists to accumulate seasonal values
    seasonal_speed_vals <- list()
    seasonal_dir_vals <- list()
    
    
    
    for (month_ in 1:length(date_ranges)) {
      date_init <- as.Date(date_ranges[[month_]][[1]])
      date_end <- as.Date(date_ranges[[month_]][[2]])
      
      file_name_u <- file.path(data_path, paste0("uVelocity_", date_init, "_", date_end, "_hourly.nc"))
      file_name_v <- file.path(data_path, paste0("vVelocity_", date_init, "_", date_end, "_hourly.nc"))
      
      u_nc <- nc_open(file_name_u)
      v_nc <- nc_open(file_name_v)
      
      if (is.null(lon)) {
        lon <- ncvar_get(u_nc, "gridX")
        lat <- ncvar_get(u_nc, "gridY")
        x_dim <- length(lon)
        y_dim <- length(lat)
        y_splits <- split_indices(y_dim, spatial_splits)
      }
      
      
      
      fillvalue_u <- ncatt_get(u_nc, "uVelocity", "_FillValue")$value
      fillvalue_v <- ncatt_get(v_nc, "vVelocity", "_FillValue")$value
      
      
      if (is.null(lon)) {
        lon <- ncvar_get(u_nc, "gridX")
        lat <- ncvar_get(u_nc, "gridY")
        x_dim <- length(lon)
        y_dim <- length(lat)
        y_splits <- split_indices(y_dim, spatial_splits)
      }
      
      y_range <- y_splits[[chunk_id]]
     
       y_start <- y_range[1]
      y_count <- length(y_range)
      
      start_u <- c(1, y_start, 1, 1)
      count_u <- c(x_dim, y_count, 1, -1)
      u_chunk <- ncvar_get(u_nc, "uVelocity", start = start_u, count = count_u)
      
      start_v <- c(1, y_start, 1, 1)
      count_v <- c(x_dim, y_count, 1, -1)
      v_chunk <- ncvar_get(v_nc, "vVelocity", start = start_v, count = count_v)
      
      nc_close(u_nc)
      nc_close(v_nc)
      rm(u_nc, v_nc)
      
      u_chunk[u_chunk == fillvalue_u] <- NA
      v_chunk[v_chunk == fillvalue_v] <- NA
      gc()
      
      speed <- sqrt(u_chunk^2 + v_chunk^2)
      direction <- (atan2(v_chunk, u_chunk) * 180 / pi) %% 360
      
      # Accumulate seasonal time series
      seasonal_speed_vals[[month_]] <- speed
      seasonal_dir_vals[[month_]] <- direction
      
      rm(u_chunk, v_chunk, speed, direction)
      gc()
    }
    
    
    # Combine seasonal data across months: [x, y, time]
    seasonal_speed <- abind::abind(seasonal_speed_vals, along = 3)
    seasonal_dir <- abind::abind(seasonal_dir_vals, along = 3)
    
    # Allocate arrays for output metrics
    speed_metrics <- array(NA, dim = c(x_dim, length(y_range), 4),
                           dimnames = list(NULL, NULL, c("min", "max", "mean", "modal")))
    dir_metrics <- array(NA, dim = c(x_dim, length(y_range), 4),
                         dimnames = list(NULL, NULL, c("min", "max", "mean", "modal")))
    
    
    
    # Metrics function
    modal <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    calc_metrics <- function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return(c(min=NA, max=NA, mean=NA, modal=NA))
      q10 <- quantile(x, probs=0.1, na.rm=TRUE)
      q90 <- quantile(x, probs=0.9, na.rm=TRUE)
      vals_below_10 <- x[x <= q10]
      min_val <- if (length(vals_below_10) > 0) mean(vals_below_10) else q10
      vals_above_90 <- x[x >= q90]
      max_val <- if (length(vals_above_90) > 0) mean(vals_above_90) else q90
      mean_val <- mean(x)
      modal_val <- modal(x)
      return(c(min=min_val, max=max_val, mean=mean_val, modal=modal_val))
    }
    
    # Loop through pixels
    for (i in 1:x_dim) {
      for (j in 1:length(y_range)) {
        speed_metrics[i, j, ] <- calc_metrics(seasonal_speed[i, j, ])
        dir_metrics[i, j, ] <- calc_metrics(seasonal_dir[i, j, ])
      }
    }
    
    # Store in chunk list 
    chunk_results[[chunk_id]] <- list(speed = speed_metrics, direction = dir_metrics)
    
    # Clean monthly lists to free memory
    rm(speed_metrics, dir_metrics)
    gc()
  } # end chunk_id loop
  
  
  # Merge spatial chunks along y dimension
  merge_metric_chunks <- function(metric_name) {
    merged_array <- NULL
    for (chunk_id in 1:spatial_splits) {
      chunk_data <- chunk_results[[chunk_id]][[metric_name]]
      if (is.null(merged_array)) {
        merged_array <- chunk_data
      } else {
        merged_array <- abind::abind(merged_array, chunk_data, along = 2)
      }
    }
    return(merged_array)
  }
  
  speed_full <- merge_metric_chunks("speed")
  direction_full <- merge_metric_chunks("direction")
  
  results<- list(
    lon = lon,
    lat = lat,
    speed = speed_full,
    direction = direction_full,
    season = season_name
  )
  
  return(results)
}



# Wrapper function for user-friendly calling by year and season name
calculate_seasonal_metrics <- function(year, season_name, spatial_splits = 3,
                                       data_path = "/Volumes/Romina_PSF/PSF/modeled_variables_original/Monthly_nc") {
  # Define months per season
  season_months <- switch(season_name,
                          "summer" = c(6,7,8),
                          "winter" = c(12,1,2),
                          "spring" = c(3,4,5),
                          "fall" = c(9,10,11),
                          stop("Unknown season name"))
  
  date_ranges <- lapply(season_months, function(m) {
    actual_year <- ifelse(m < 3 & season_name == "winter", year + 1, year)
    start_date <- as.Date(sprintf("%d-%02d-01", actual_year, m))
    # Calculate end_date as the first day of next month
    end_date <- start_date %m+% months(1)
    list(start_date, end_date)
  })
  
  calculate_seasonal_metrics_from_date_ranges(date_ranges, season_name, spatial_splits, data_path)
}




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


save_metrics_from_final_df <- function(final_df, year=2014, season, variable_name = "current", output_path=output_path) {
  # Columns to save (all except lon, lat, gridID)
  metric_cols <- setdiff(colnames(final_df), c("gridX", "gridY"  ,"gridID", "longitude","latitude", "bathymetry"))
  
  for (metric_col in metric_cols) {
    # Create a subset dataframe with only lon, lat, and the metric
    df_sub <- final_df[, c("longitude","latitude", metric_col)]
    colnames(df_sub)[3] <- "value"  # rename metric column to 'value'
    
    # Construct filename:
    # Example: current_speed_2014_fall_mean_surface_.csv
    measurement= strsplit(metric_col, "_")[[1]][1]
    metric= strsplit(metric_col, "_")[[1]][2]
    
    fname <- paste0(variable_name,"_",measurement,"_", year, "_", season, "_", metric,"_surface_.csv")
    fname <- paste(output_path, fname, sep="/")
    
    # Save CSV
    write.csv(df_sub, fname, row.names = FALSE)
  }
}


# ✅ Overview of What to Do
# 1️⃣ Load your depth CSV → it should have columns like gridX, gridY, longitude, latitude.
# 2️⃣ From your results, get each metric array (e.g., results$speed[,, "mean"]), flatten them into vectors.
# 3️⃣ Create a dataframe where each row corresponds to a grid point (matching gridX, gridY).
# 4️⃣ Merge metrics and coordinates using gridX, gridY as keys.
# 5️⃣ Write the final combined dataframe to CSV.

# Load your depth CSV: adjust path and column names as needed
path_bathy= "/Volumes/Romina_PSF/PSF/modeled_variables_original"
bathy_df <- read.csv(paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/"))  # should have gridX, gridY, lon, lat
colnames(bathy_df)
output_path= "/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth"


seasons= c("summer", "fall", "spring", "winter") #


for (year in 2015:2022) {
  
  for (season in seasons) {
    
    results <- calculate_seasonal_metrics(year=year, season_name= season, spatial_splits=3)
    
    dim_speed <- dim(results$speed)  # Should be [x, y, 4]
    
    grid_df <- expand.grid(
      gridY = 1:dim_speed[2],  # 900
      gridX = 1:dim_speed[1]   # 398
    )
    grid_df$gridID <- 1:nrow(grid_df)
    
    
    
    # Step 2: Extract metrics from results
    extract_metric_df <- function(metric_array, metric_name) {
      metric_values <- as.vector(aperm(metric_array, c(2, 1)))  # Flatten [y, x]
      if (length(metric_values) != nrow(grid_df)) {
        stop("Mismatch between metric size and grid_df")
      }
      df <- grid_df
      df[[metric_name]] <- metric_values
      return(df[, c("gridID", "gridX", "gridY", metric_name)])
    }
    
    # Step 3: Extract each metric
    speed_mean_df  <- extract_metric_df(results$speed[, , "mean"],   "speed_mean")
    speed_max_df   <- extract_metric_df(results$speed[, , "max"],    "speed_max")
    speed_min_df   <- extract_metric_df(results$speed[, , "min"],    "speed_min")
    speed_modal_df <- extract_metric_df(results$speed[, , "modal"],  "speed_modal")
    
    dir_mean_df    <- extract_metric_df(results$direction[, , "mean"],   "dir_mean")
    dir_max_df     <- extract_metric_df(results$direction[, , "max"],    "dir_max")
    dir_min_df     <- extract_metric_df(results$direction[, , "min"],    "dir_min")
    dir_modal_df   <- extract_metric_df(results$direction[, , "modal"],  "dir_modal")
    
    # Step 4: Merge all metric dataframes by gridID, gridX, gridY
    metrics_df <- Reduce(function(x, y) merge(x, y, by = c("gridID", "gridX", "gridY"), all = TRUE), 
                         list(speed_mean_df, speed_max_df, speed_min_df, speed_modal_df,
                              dir_mean_df, dir_max_df, dir_min_df, dir_modal_df))
    
    # Step 5: Merge with bathymetry using gridX and gridY
    # Assume bathy_df has columns: gridX, gridY, depth or similar
    final_df <- merge(metrics_df, bathy_df, by = c("gridX", "gridY"), all.x = TRUE)
    
    
    save_metrics_from_final_df(final_df, year, season, variable_name = "current", output_path="/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth")
    
  }
  
}





