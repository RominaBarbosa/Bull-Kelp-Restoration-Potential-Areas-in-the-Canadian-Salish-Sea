###=============================================================================
###     NEMO model - SalishSeaCast data process                 ################
###                                                             ################
### - Sensitivity test to select the idp value for IDW interpolations ##########
### input data acquired with: "1.1.Download_nc_NEMOmodel.R"     ################
### input data: climatology metrics (dataframe)                 ################
### output data: table with accuracy metrics of interpolations  ################
### Author: Romina Barbosa                                      ################
### Date: 21-July-2025                                          ################
### Last edition: 22-July-2025   
###=============================================================================

# Load packages
library(raster) # package for raster manipulation
library(ggplot2) 
library(gstat)


### Varify the interpolations by Computing validation metrics  =======
# Test multiple IDP values to optimize interpolation
bathy <- rast("/Volumes/Romina_PSF/PSF/Env_Variables/bathy_coastwide_500m.tif")
output_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers"
variables_df <- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/selected_variables_to_interpolations.csv")


# eval_interpolation<- function(variables_data= variables_df, output_path= output_path){
  
  # Define the values of idp to test
  idp_values <- seq(1, 4, by = 0.5)
  
  # Initialize a data frame to store results
  results <- data.frame(
    idp = idp_values,
    RMSE = NA,
    MAE = NA,
    R2 = NA,
    variable=NA
  )
  
  for (i in 3:ncol(variables_df)) { #
    # Initialize a data frame to store results
    results_i <- data.frame(
      idp = idp_values,
      RMSE = NA,
      MAE = NA,
      R2 = NA
    )
    
    
    variable_name<- colnames(variables_df)[i]
    variable<- variables_df[c("longitude", "latitude", paste(variable_name))]
    
    # Replace 'longitude' and 'latitude' with your actual column names
    points_sf <- st_as_sf(variable, coords = c("longitude", "latitude"), crs = 4326)
    
    # 2. Reproject to match your bathymetry raster (assuming EPSG:3005)
    points_sf <- st_transform(points_sf, crs = 3005)
    
    # 3. Convert to SpatialPointsDataFrame
    points_sp <- as(points_sf, "Spatial")
    names(points_sp)<- "value"
    
    
    # Loop over idp values
    for (a in seq_along(idp_values)) {
      idp <- idp_values[a]
      
      # # Run leave-one-out cross-validation
      # cv <- gstat::krige.cv(
      #   value ~ 1,
      #   locations = points_sp,
      #   nfold = nrow(points_sp),  # LOOCV
      #   set = list(idp = idp)
      # )
      
      # 10-fold cross-validation (much faster)
      cv <- gstat::krige.cv(
        value ~ 1,
        locations = points_sp,
        nfold = 10,
        set = list(idp = idp)
      )
      
      # Calculate metrics
      rmse <- sqrt(mean((cv$observed - cv$var1.pred)^2, na.rm = TRUE))
      mae <- mean(abs(cv$observed - cv$var1.pred), na.rm = TRUE)
      r2 <- cor(cv$observed, cv$var1.pred, use = "complete.obs")^2
      
      # Store
      results_i$RMSE[a] <- rmse
      results_i$MAE[a] <- mae
      results_i$R2[a] <- r2
    }
    
    
    results_i$variable<- variable_name
    results<-  rbind(results, results_i)
    
    
  }
#   # View the results
#   print(results)
#   
#   return(results)
# }


# write.csv(results, paste(output_path, "results_idw.csv", sep="/"))

# interpolation_evaluation<- eval_interpolation(variables_data= variables_df[,1:5], output_path= output_path)


### Plot results 
results<- na.exclude(results)

results%>%
  filter(variable=="temperature_summer_mean")%>%
  ggplot( aes(x = idp)) +
  geom_line(aes(y = RMSE), color = "blue") +
  geom_point(aes(y = RMSE), color = "blue") +
  geom_line(aes(y = MAE), color = "orange") +
  geom_point(aes(y = MAE), color = "orange") +
  labs(
    title = "Cross-Validation Error vs IDW Power",
    y = "Error (RMSE / MAE)",
    x = "IDW Power (idp)"
  ) +
  facet_wrap(~variable, ncol=3)+
  theme_minimal()

# ggsave(paste(output_path, "results_idw_plot_temperature_summer_mean.png", sep="/"), width = 14, height = 15, units = "cm")
