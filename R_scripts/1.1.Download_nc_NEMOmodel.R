###=============================================================================
###     NEMO model - SalishSeaCast data download                ################
###                                                             ################
### 1.1- Loop to Download the data from the website             ################
### https://salishsea.eos.ubc.ca/erddap/griddap/index.html?page=1&itemsPerPage=1000
### https://salishsea.eos.ubc.ca/                               ################
### Author: Romina Barbosa                                      ################
### Date: 15-Feb-2025                                           ################
### Last edition: 21-Feb-2025     
###=============================================================================


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

download_NEMO<- function(depth= "surface", year = "2014", variable= "u_wind",
                         output_path= "E:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
                         months= c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){
  
  require(curl)
  # months<- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  
  if (depth == "surface"){
    depth = 0.5000003
  }
  if (depth == "bottom_10"){
    depth = 10.5047655
  }
  if (depth == "bottom_15"){
    depth = 15.634288
  }
  
  for (i in 1:length(months)) {
    month = months[i]
    day = "01"
    end_month = months[i+1]
    
    if(is.na(end_month)){
      print("you need to give an extra month after the end month you need, if you need Jan and Feb them give months from Jan to April")
    }
    
    # Set initial and end date
    init_date = paste(year, month, day, sep="-")
    end_date = paste(year, end_month, day, sep="-")
    
    # Set input file name to download data and output file name to save them
    file_name<- paste("https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DPhysicsFields1hV21-11.nc?",variable,"[(", init_date, "T00:30:00Z):1:(", end_date,"T00:30:00Z)][(0.5000003):1:(0.5000003)][(0.0):1:(897.0)][(0.0):1:(397.0)]", sep="")
    save_file_name<- paste(variable, init_date, end_date, "hourly", sep= "_")
    
    # Download data from website
    curl_download(file_name,
                  destfile = paste(output_path, paste(save_file_name, ".nc", sep=""), sep="/"))
    
    }
}


### ============================================================================
### Download the data ==========================================================
### ============================================================================
# Increase the time to download to avoid error when the dataset is too big
options(timeout = 500)  # Set timeout to 300 seconds

## Download Temperature ====
init_time= Sys.time()

download_NEMO(depth= "surface", year = "2015", variable= "v_wind", output_path= "D:/PSF/modeled_variables_original/Monthly_nc",
              months= c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"))

end_time= Sys.time()
end_time - init_time
# Time difference of 10.7525 mins


## Download salinity ====
init_time= Sys.time()

download_NEMO(depth= "surface", year = "2020", variable= "salinity", output_path= "E:/Pacific_Salmon_Fundation/modeled_variables_original/Monthly_nc",
              months= c("01", "02", "03", "04", "05", "06", "07", "08"))

end_time= Sys.time()
end_time - init_time



