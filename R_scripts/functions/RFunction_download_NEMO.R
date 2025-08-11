###=============================================================================
###     NEMO model - SalishSeaCast data download                ################
###                                                             ################
### 1.1- Loop to Download the data from the website             ################
### https://salishsea.eos.ubc.ca/erddap/griddap/index.html?page=1&itemsPerPage=1000
### https://salishsea.eos.ubc.ca/                               ################
### Author: Romina Barbosa                                      ################
### Date: 15-Feb-2025                                           ################
### Last edition: 28-March-2025    
###=============================================================================


## NOTE: WInd conditions from 2023 and 2024 are not available in the website

# Load packages
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(lubridate)
library(chron)
library(mondate)

## Download files from the web
### ============================================================================
## To download products from: Green, Salish Sea, 3d Physics Fields, Hourly, v21-11
# The available variables are:
# salinity (Reference Salinity, g kg-1)
# temperature (Conservative Temperature, degC)
# sigma_theta (Potential Density (sigma_theta), kg m-3)

#To download SST and Salinity:
#"https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DPhysicsFields1hV21-11.nc?salinity[(2021-08-01T00:30:00Z):1:(2021-09-01T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)],temperature[(2021-08-01T00:30:00Z):1:(2021-09-01T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]"

# To download only SST:
# download.file(url="https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DPhysicsFields1hV21-11.nc?temperature[(2021-08-20T00:30:00Z):1:(2021-09-01T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]",
#               destfile=paste(path,"2021_August_temp.nc", sep="/"))


### ============================================================================
## Example for HRDPS, SalishSeaCast, Atmospheric Forcing Fields, Hourly, v23-02
# The available variables are:
# precip (Precipitation Flux, kg m-2 s-1)
# solar (Downward Short-Wave (Solar) Radiation Flux, W m-2)
# tair (Air Temperature at 2m, K)
# u_wind (U-Component of Wind at 10m, m s-1)
# v_wind (V-Component of Wind at 10m, m s-1)


### ============================================================================
### To download products from the surface atmosphere fields: # HRDPS, SalishSeaCast, Atmospheric Forcing Fields, Hourly, v23-02
#"https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSaSurfaceAtmosphereFieldsV23-02.nc?u_wind[(2025-02-20T23:00:00Z):1:(2025-02-20T23:00:00Z)][(0.0):1:(229.0)][(0.0):1:(189.0)],v_wind[(2025-02-20T23:00:00Z):1:(2025-02-20T23:00:00Z)][(0.0):1:(229.0)][(0.0):1:(189.0)]"
# The available variables are:
#   atmpres (Air Pressure at MSL, Pa)
#   precip (Precipitation Flux, kg m-2 s-1)
#   qair (Specific Humidity at 2m, kg kg-1)
#   rhair (2 metre relative humidity, %)
#   solar (Downward Short-Wave (Solar) Radiation Flux, W m-2)
#   tair (Air Temperature at 2m, K)
#   therm_rad (Downward Long-Wave (Thermal) Radiation Flux, W m-2)
#   lhtfl (Latent heat net flux, W m-2)
#   u_wind (U-Component of Wind at 10m, m s-1)
#   v_wind (V-Component of Wind at 10m, m s-1)


### ============================================================================
### To download current speed from Green, Salish Sea, 3d v Grid Variable Fields, Hourly, v21-11
# The available variable is: vVelocity (Ocean Current Along y-axis, m s-1)
# "https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DvGridFields1hV21-11.nc?vVelocity[(2020-01-01T00:30:00Z):1:(2020-01-31T23:30:00Z)][(15.634288):1:(0.5000003)][(0.0):1:(897)][(0.0):1:(397.0)]"


### ============================================================================
### To download current speed from Green, Salish Sea, 3d u Grid Variable Fields, Hourly, v21-11
#  The available variable is: uVelocity (Ocean Current Along x-axis, m s-1)
# "https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DuGridFields1hV21-11.nc?uVelocity[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]"


### ============================================================================
### To download products from Green, Salish Sea, 3d PAR and Turbidity Fields, Hourly, v21-11
#   The available variable is:
#   PAR (Photosynthetically Available Radiation, W m-2)
#   turbidity (Fraser River Turbidity, NTU)
# https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DLightFields1hV21-11.nc?PAR[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)],turbidity[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]


### ============================================================================
### To download products from Green, Salish Sea, 3d Chemistry Fields, Hourly, v21-11
# The available variables are:
#   dissolved_inorganic_carbon (mmol m-3)
#   total_alkalinity (Total Alkalinity Concentration, mmol m-3)
#   dissolved_oxygen (Dissolved Oxygen Concentration, mmol m-3)
# https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DChemistryFields1hV21-11.nc?dissolved_inorganic_carbon[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)],total_alkalinity[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)],dissolved_oxygen[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]


### ============================================================================
### To download products from Green, Salish Sea, 3d Biology Fields, Hourly, v21-11
# https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1hV21-11.nc?nitrate[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(441.4661)][(0.0):1:(897.0)][(0.0):1:(397.0)]
# The available variables are:
#   ammonium (Ammonium Concentration, mmol m-3)
#   biogenic_silicon (Biogenic Silicon Concentration, mmol m-3)
#   diatoms (Diatoms Concentration, mmol m-3)
#   dissolved_organic_nitrogen (mmol m-3)
#   flagellates (Flagellates Concentration, mmol m-3)
#   z2_zooplankton (Z2 Zooplankton Concentration, mmol m-3)
#   z1_zooplankton (Z1 Zooplankton Concentration, mmol m-3)
#   nitrate (Nitrate Concentration, mmol m-3)
#   particulate_organic_nitrogen (mmol m-3)
#   silicon (Silicon Concentration, mmol m-3)


### ============================================================================
### Function to automatically download the datasets per month (hourly data)
### ============================================================================

# the function works for different depth layers including: surface = 0.5 m depth; bottom_10 = 10 m depth; bottom_15 = 15 m depth
# There are different variables available including:

download_NEMO<- function(depth= "surface", years_all = c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022"), variable= "temperature",
                         output_path= "D:/PSF/modeled_variables_original/Monthly_nc",
                         months= c( "06", "07", "08", "09", "10", "11", "12")){
  
  require(curl)
  # months<- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  
  if (depth == "surface"){
    depth = 0.5000003
  }
  if (depth == "bottom_5"){
    depth = 5.5001507
  }
  if (depth == "bottom_10"){
    depth = 10.5047655
  }
  if (depth == "bottom_15"){
    depth = 15.634288
  }
  
  for (y in 1:(length(years_all)-1)) {
    year= years_all[y]
    
    for (i in 1:(length(months))) {
      month = months[i]
      day = "01"
      
      if(month == "12"){
        # Set initial and end date
        end_year = years_all[y+1]
        init_date = paste(year, month, day, sep="-")
        end_date = paste(end_year, "01", day, sep="-")
        if(is.na(end_year)){
          print("if you are downloading data of DECEMBER, you must add an extra year in years_all")
        }
      }else{
        
        end_month = months[i+1]
        # Set initial and end date
        init_date = paste(year, month, day, sep="-")
        end_date = paste(year, end_month, day, sep="-")
        
        if(is.na(end_month)){
          print("you need to give an extra month after the end month you need, if you need Jan and Feb them give months from Jan to April")
        }
      }
      
      if (variable == "temperature" | variable == "salinity" ){
        source= "ubcSSg3DPhysicsFields1hV21"
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, "-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
      }
      
      if (variable == "PAR" | variable == "turbidity" ){
        source= "ubcSSg3DLightFields1hV21"
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, "-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
      }
      
      if (variable == "dissolved_inorganic_carbon" | variable == "total_alkalinity" | variable == "dissolved_oxygen"){
        # https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DChemistryFields1hV21-11.nc?dissolved_inorganic_carbon[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)],total_alkalinity[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)],dissolved_oxygen[(2025-02-20T23:30:00Z):1:(2025-02-20T23:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]
        source= "ubcSSg3DChemistryFields1hV21"
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, "-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
      }
      
      if (variable == "nitrate" | variable == "ammonium" ){
        source= "ubcSSg3DBiologyFields1hV21"
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, "-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
      }
      
      if (variable == "u_wind" | variable == "v_wind" ){
        source= "ubcSSaSurfaceAtmosphereFieldsV1"
        # https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSaSurfaceAtmosphereFieldsV1.nc?u_wind[(2023-02-22T23:00:00Z):1:(2023-02-22T23:00:00Z)][(0.0):1:(662500.0)][(0.0):1:(637500.0)]
        #file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, ".nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date, "T00:30:00Z)][(0.0):1:(662500.0)][(0.0):1:(637500.0)]", sep="")
        
        ##Second download _ 2014
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, ".nc?",variable,"%5B(", init_date, "T00:30:00Z):1:(", end_date, "T00:30:00Z)%5D%5B(0.0):1:(662500.0)%5D%5B(0.0):1:(637500.0)%5D", sep="")
      }
      
      
      if (variable == "uVelocity"){
        source= "ubcSSg3DuGridFields1hV21"
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, "-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
      }
      
      if (variable == "vVelocity"){
        source= "ubcSSg3DvGridFields1hV21"
        file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/", source, "-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
      }
      
      
      
      
      if(depth == 0.5000003){
        
        save_file_name<- paste(variable, init_date, end_date, "hourly", sep= "_")
      }else{
        
        save_file_name<- paste(variable, init_date, end_date, round(depth, 1), "hourly", sep= "_")
      }
      
      # Download data from website
      curl_download(file_name,
                    destfile = paste(output_path, paste(save_file_name, ".nc", sep=""), sep="/"))
      
    }
  }
}


### ============================================================================
### Download the data ==========================================================
### ============================================================================
# Increase the time to download to avoid error when the dataset is too big
options(timeout = 900)  # Set timeout to 300 seconds

## Download Temperature ====
init_time= Sys.time()


download_NEMO(depth= "surface", years_all= c("2022"), 
              variable= "total_alkalinity", 
              output_path= "D:/PSF/modeled_variables_original/Monthly_nc",  #"G:/Pacific_Salmon_Fundation//Monthly_nc",
              months= c( "11"))




init_time2= Sys.time() # started at 4:25 pm 27-02-2025





###=============================================================================
### PAR and Turbidity ==========================================================
download_NEMO(depth= "surface", years_all= c("2022", "2023"), 
              variable= "PAR", 
              output_path= "Y:/Romina/Monthly_nc",  #"G:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
              months= c("12"))

###=============================================================================
### nitrate   ==========================================================
download_NEMO(depth= "surface", years_all= c("2023", "2024"),
              variable= "v_wind", 
              output_path= "D:/PSF/modeled_variables_original/Monthly_nc",  #"G:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
              months= c("01", "02", "03", "04","05", "06", "07", "08", "09", "10", "11", "12"))


## Download salinity ====
download_NEMO(depth= "surface", years_all= c("2022", "2023", "2024"), 
              variable= "vVelocity", 
              output_path= "Y:/Romina/Monthly_nc", #"G:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
              months= c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"))


###=============================================================================
# variables= c("dissolved_inorganic_carbon", "PAR", "u_wind", "ammonium", "nitrate")

variables= c( "u_wind", "v_wind")

for (a in 2:length(variables)) {
  variable_a= variables[a]
  download_NEMO(depth= "surface", years_all= c("2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"),
                variable= variable_a, 
                output_path= "PSF:/modeled_variables_original/Monthly_nc",  #"G:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
                months= c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"))
}


# for (a in 1:length(variables)) {
#   variable_a= variables[a]
#   download_NEMO(depth= "surface", years_all= c("2023", "2024", "2025"),
#                 variable= variable_a, 
#                 output_path= "Y:/Romina/Montly_nc_midwater5m",  #"G:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
#                 months= c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"))
# }



