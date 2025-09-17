###=============================================================================
###     NEMO model - SalishSeaCast data process                 ################
###                                                             ################
### 1.2.3- Transform yearly seasonal metrics into climatology metrics   ##########
### https://salishsea.eos.ubc.ca/                               ################
### input data acquired with: "1.1.Download_nc_NEMOmodel.R"     ################
### input data: seasonl metrics per year                        ################
### output data: climatology of a period, e.g., "blob"          ################
### Author: Romina Barbosa                                      ################
### Date: 5-July-2025                                           ################
### Last edition: 11-Aug-2025 -corrected function to include lat and long ######
###=============================================================================
# Load packages
library(dplyr)
library(sf)

source(paste("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions","RFuntion_merge_nemo_climatology_data.R", sep="/"))

variables<- c("current_speed", "current_dir")
variables<- c(  "nitrate", "ammonium", "dissolved_inorganic_carbon",
              "turbidity",  "salinity", "temperature", "total_alkalinity")
seasons<- c( "winter", "fall","spring","summer") #
# seasons<- c() # Winter of wind was calculated from 2016 to 2019 due to lck of data from 2014 
# seasons<- c( "fall","winter","spring","summer") #

metrics<- c("mean", "maximum", "minimum", "SD")
metrics<- c("mean", "max", "min", "modal")

# output_path<- "/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_blob_"
output_path<- "/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_post_blob"
time_period = "post_blob"
depth = "surface"


for (variable in variables) {
  
  for (season in seasons) {
    
    for (metric in metrics) {
      
      # merged_data_mean<- merge_nemo_climatology_data(depth = "surface", time_period = "post_blob", 
      #                                           variable = variable, 
      #                                           metric= metric,
      #                                           input_path = "/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth", 
      #                                           season = season)
      
      # gridX, gridY, latitude, longitude for ALL cells
      grid_template <- read.csv("/Volumes/Romina_PSF/PSF/modeled_variables_original/ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv") 
      grid_template<- grid_template[-1,]
      
      grid_template$latitude <- as.numeric(as.character(grid_template$latitude))
      grid_template$longitude <- as.numeric(as.character(grid_template$longitude))


      merged_data_mean <- merge_nemo_climatology_data(
        depth = "surface",
        time_period = time_period,
        variable = variable,
        metric = metric,
        input_path = "/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth",
        season = season,
        grid_template = grid_template
      )
      
      write.csv(merged_data_mean, paste(output_path, paste(variable, time_period, season, metric, depth,".csv", sep = "_"), sep = "/"))

    }
  }
  
}









#### Calculate metric for metrics of temperature above threshold ===============
output_path<- "/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_temp_tolerance_metrics"
setwd("/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth")
dir()
setwd(output_path)


# Set up what metrics you're working with in order to set the names for input files
metric_temperature = "cumulated_degrees" # this or the  metric type in the next row
metric_temperature = "hours_above_threshold"


if(metric_temperature == "cumulated_degrees"){
  metrics<- c()
  for (i in 13:18) {
  metric<- paste("cumulated_degrees", i, sep="_")
  metrics<- c(metrics, metric)
  }

 }else{
  metrics<- c()
  for (i in 11:18) {
    metric<- paste("hours_above_threshold", i, sep="_")
    metrics<- c(metrics, metric)
  }
 }


seasons<- c("spring","summer")

# gridX, gridY, latitude, longitude for ALL cells
grid_template <- read.csv("/Volumes/Romina_PSF/PSF/modeled_variables_original/ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv") 
grid_template<- grid_template[-1,]

grid_template$latitude <- as.numeric(as.character(grid_template$latitude))
grid_template$longitude <- as.numeric(as.character(grid_template$longitude))


for (season in seasons) {
    for (metric in metrics) {
      
      
      merged_data_mean <- merge_nemo_climatology_data(
        depth = "surface",
        time_period = "blob",
        variable = "temperature",
        metric = metric,
        input_path = "/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth",
        season = season,
        grid_template = grid_template
      )
      
       # merged_data_mean<- merge_nemo_climatology_data(depth = "surface", time_period = "blob", 
      #                                                variable = "temperature", 
      #                                                metric= metric,
      #                                                input_path = "/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth", 
      #                                                season = season)
      
      write.csv(merged_data_mean, paste(output_path, paste(variable, time_period, season, metric, depth,".csv", sep = "_"), sep = "/"))
    }
}
