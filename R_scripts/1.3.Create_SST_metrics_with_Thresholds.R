###=============================================================================
###     NEMO model - SalishSeaCast data process                 ################
###                                                             ################
### 1.3- Loop to merge monthly data into seasons                ################
###    - and calculate metrics
### https://salishsea.eos.ubc.ca/                               ################
### input data acquired with: "1.1.Download_nc_NEMOmodel.R"     ################
### Author: Romina Barbosa                                      ################
### Date: 15-Feb-2025                                           ################
### Last edition: 18-March-2025
### converted into metrics 
###=============================================================================

# Load packages
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ncdf4)
library(ggplot2) # package for plotting
library(lubridate)
library(chron)
library(mondate)
library(abind)

## Set Paths ===================================================================
input_path= "/modeled_variables_original/Monthly_nc"
setwd(input_path)  

output_path=  "/modeled_variables_original/seasonal_metrics_0.5m_depth"
path_bathy= "/modeled_variables_original"   # used bathymetry layer to extract coordinates: ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv

source(paste("/R_scripts/functions","RFunction_merge_nemo_seasonal_data.R", sep="/"))
source(paste("/R_scripts/functions","RFunction_calculate_simple_metrics.R", sep="/"))
source(paste("/R_scripts/functions","RFunction_temperatureOverThreshold_metrics.R", sep="/"))


years<- c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")
seasons<- c("spring", "fall", "winter", "summer")  # select seasons 

for (s in 1:2) {
  require(abind)
  
  season_s<- seasons[s]
  
  for(year in years){
    
    merged_season_data <- merge_nemo_seasonal_data(depth = "surface", year = year, variable = "temperature", 
                                                   input_path = input_path, 
                                                   season = season_s)
    
    for (t in 11:18) {
      
      calculate_DegreeHours_Hours_OverThreshold(data= merged_season_data, threshold_temp= t, 
                                                bathy_file= paste(path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/"),
                                                depth= "surface", season= season_s)
    }
   
  time= Sys.time()
  print(paste("finished", year, ";", season_s, time))
  }
  
}




