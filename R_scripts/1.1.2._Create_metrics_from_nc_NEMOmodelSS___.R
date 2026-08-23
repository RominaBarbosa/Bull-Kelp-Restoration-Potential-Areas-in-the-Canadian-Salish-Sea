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





## ===================================================================
### Calculate temperature Over Threshold metrics
## ===================================================================
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
