###=============================================================================
###     NEMO model - SalishSeaCast Variables Selection          ################
###                                                             ################
### input data acquired with: "1.4.Climatology_metrics_Xperiod.R" ############
### input data: climatology metrics (dataframe)                 ################
### output data: table with selected variables  (metrics)       ################
### Author: Romina Barbosa                                      ################
### Date: 17-July-2025                                          ################
### Last edition: 01-August-2025   
# added nitrate variables and calculated speed (vector)and direction of water current
###=============================================================================



### Filter the data to exclude values below 70 m depth (I selected this depth 
### in order to keep enough values (res. 500 m) for a further interpolation, if needed)

library(data.table)
library(purrr)
library(readr)
library(dplyr)
library(tidyverse)
library(stringr)
library(rasterVis)

###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++####
cor_plot_function<- function(data= selected_variables, 
                             path_plot= "/SDM/Variables_selection/plots",
                             plot_name= "CorrVars_clustered",
                             plot_type= "pdf", 
                             pdf_height= 25,
                             pdf_width= 30){
  
  COL2 <- function(name, n) {
    colorRampPalette(RColorBrewer::brewer.pal(11, name))(n)
  }
  setwd(path_plot)
  cor_matrix.d<- cor(data, use = "pairwise.complete.obs")
  
  if(plot_type == "png"){
    
    plot_name1<- paste(path_plot, paste(plot_name,"matrix.png", sep="_"), sep="/")
    png(plot_name1,
        height = 18, width = 18, units="cm", res=300)
    
    corrplot::corrplot(cor_matrix.d, type = 'lower', tl.col = 'black',#order = 'hclust',
                       cl.ratio = 0.15, tl.srt = 45, col = COL2('PuOr', 30), cl.cex = 0.8,  tl.cex = 0.6) #, tl.pos = "n"
    
  }
  
  if(plot_type == "pdf"){
    
    plot_name1<- paste(path_plot, paste(plot_name,"matrix.pdf", sep="_"), sep="/")
    
    pdf(plot_name1, 
        height = pdf_height, width = pdf_width)
    
    corrplot::corrplot(cor_matrix.d, type = 'lower', tl.col = 'black',#order = 'hclust',
                       cl.ratio = 0.15, tl.srt = 45, col = COL2('PuOr', 30), cl.cex = 2,  tl.cex = 2) #, tl.pos = "n"
    
  }
  
  
  dev.off()
  
  
  cor_variables<- as.matrix(round(cor_matrix.d, 3))
  dissimilarity = 1 - cor_variables
  distance = as.dist(dissimilarity) 
  
  hClust<- hclust(distance, method="average")
  # clust <- rect.hclust(hClust, h=0.7, border=0)
  # plot(hClust, main = "Clustering of Correlated Variables", xlab = "", sub = "")
  # abline(h = 0.3, col = "red", lty = 2) 
  # rect.hclust(hClust, h = 0.3, border = 2:6) #Draw cluster rectangles at same height
  
  
  ## text labels rotated 45 degrees and  wider color legend with numbers right aligned
  
  if(plot_type == "png"){
    plot_name2<- paste(path_plot, paste(plot_name,"clusters.png", sep="_"), sep="/")
    
    png(plot_name2,
      height = 14, width = 17, units="cm", res=300)
    
    plot(hClust, main = "Clustering of Correlated Variables", xlab = "", sub = "", cex = 0.6)
    abline(h = 0.3, col = "red", lty = 0.3) 
    
    rect.hclust(hClust, h = 0.3, border = 2:6) #Draw cluster rectangles at same height
    
  }
  
  if(plot_type == "pdf"){
    plot_name2<- paste(path_plot, paste(plot_name,"clusters.pdf", sep="_"), sep="/")
    
    pdf(plot_name2, 
        height = pdf_height,
        width = pdf_width)
    
    plot(hClust, main = "Clustering of Correlated Variables", xlab = "", sub = "", cex = 2)
    abline(h = 0.3, col = "red", lty = 8) 
    
    rect.hclust(hClust, h = 0.3, border = 2:6) #Draw cluster rectangles at same height
    
  }
  
   dev.off()
  
  return(hClust) 
}


###=============================================================================
### Prepare bathymetry data to reduce size of datasets by excluding areas depper than 70 m =================
depth_raster <- raster("/opographic_variables_20mres/coastwide_20m.tif")
# depth_raster <- raster("SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
crs(depth_raster)

bathy_filtered <- depth_raster
bathy_filtered[ bathy_filtered>= 70 ]<- NA
bathy_filtered[ bathy_filtered<= -50 ]<- NA

# Visualize the mask 
plot(depth_raster, main = "Depth Raster")
plot(bathy_filtered) #, add = TRUE, col = "red", pch = 20)
names(bathy_filtered)<- "depth_m"

crs(bathy_filtered)<- CRS("+init=epsg:3005") 


## SET PATCH AND UPLOAD FILES of env variables
### Upload variables from the SalishSeaCast model and merge ====================
my_path<-("/modeled_variables_original/climatology_metrics")
setwd(my_path)
dir()
files <- list.files(my_path, pattern = "\\.csv$", full.names = TRUE)
length(files)# 188


# Load file 1 to start building the merged dataset
df_all=read.csv(files[1])
head(df_all)
df_all<- df_all[,c(2:4)]
filename <- strsplit(tools::file_path_sans_ext(files[1]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(1,3,4)]
filename <- paste0(filename, collapse ="_")
colnames(df_all)[3]<- filename

# Merge all files by adding a column of the metric to the df_all dataset
for (i in 2:156) {
  
  df= read.csv(files[i])
  df= df[,c(2:4)]
  
  if(i %in% c(1:16, 65:124, 141:156)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(1,3,4)]
    filename <- paste0(filename, collapse ="_")
  }
  
  if(i %in% c(17:32)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(4:5)]
    filename <- paste("currentDirection", paste0(filename, collapse ="_"), sep="_")
  }
  
  if(i %in% c(33:48)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(4:5)]
    filename <- paste("currentSpeed", paste0(filename, collapse ="_"), sep="_")
  }
  
  if(i %in% c(49:64)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(5:6)]
    filename <- paste("DIC", paste0(filename, collapse ="_"), sep="_")
  }
  
  if(i %in% c(125:140)){
    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(4:5)]
    filename <- paste("alkalinity", paste0(filename, collapse ="_"), sep="_")
  }
  
  # if(i %in% c(125:140)){
  #   filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
  #   filename <- strsplit(filename, "_")[[1]][c(4:5)]
  #   filename <- paste("wind", paste0(filename, collapse ="_"), sep="_")
  # }
  
  colnames(df)[3]<- filename
  df_all= merge(df_all, df, by=c("latitude" , "longitude"))
  
}

colnames(df_all)

SSCast_variables<- df_all


# ### Make sure the coordinates are correct, x and y (crs=3005)
# nemo=read.csv(files[1])
# nemo_path_bathy<-"/Volumes/Romina_PSF/PSF/modeled_variables_original"
# bathy_nemo= read.csv(paste(nemo_path_bathy, "ubcSSnBathymetryV21-08_d0bf_13e6_3a75.csv", sep="/"))
# bathy_nemo[length(bathy_nemo$gridY), "gridY"] #897
# bathy_nemo[length(bathy_nemo$gridY), "gridX"] #397
# length(bathy_nemo$gridY)

# length(unique(SSCast_variables$latitude))
# length(unique(SSCast_variables$longitude))
# length(SSCast_variables$longitude)
# 
# merged_df <- merge(bathy_grid, value_df, by = c("gridY", "gridX"))
# merged_df <- merged_df[, c("latitude", "longitude", "value")]



###=============================================================================
### Add depth value to variables data by location===============================
# Transform your spatial object 
coordinates(SSCast_variables) <- ~longitude + latitude  # Convert to spatial object
proj4string(SSCast_variables) <- CRS("+proj=longlat +datum=WGS84")  # Set CRS (WGS84)

# Convert to sf with original long/lat and WGS84 CRS first
sf_points <- st_as_sf(SSCast_variables, coords = c("longitude", "latitude"), crs = 4326)

bc_albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
# bc_albers <- CRS("+init=epsg:3005")  # or full string:
SSCast_bc <- spTransform(SSCast_variables, bc_albers)

# Step 2: Convert to sf object (for mapview)
sf_points <- st_as_sf(SSCast_bc)
head(st_coordinates(sf_points))
plot(sf_points[1], cex=0.04)


# CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Note: If you get a warning about +init=epsg:3005 being deprecated, use the full string instead.

merged_df <- spTransform(SSCast_variables, bc_albers)
values_depth <- raster::extract(bathy_filtered, merged_df)
summary(values_depth)
df_with_values <- cbind(as.data.frame(merged_df), depth = values_depth)
df_with_values<- df_with_values[,c(157:158, 1:156, 159)]
names(df_with_values)[159]<- "depth_m"
names(df_with_values)[1:2]<- c("longitude", "latitude")

# Exclude grids outside of bathymetry thresholds
SSCast_variables<- as.data.frame(df_with_values)%>%
  filter(!is.na(depth_m))


sf_points <- st_as_sf(SSCast_variables, coords = c("longitude", "latitude"), crs = bc_albers)
plot(sf_points[1])
mapview::mapview(sf_points[1])
mapview(bathy_filtered)
# # Save bathy filtered csv file (  70 <= depth <= -50)
# write.csv(df, "filtered_raster_data.csv", row.names = FALSE)



################################################################################
####============================================================================
### Correlation among variables from Salish Sea Cast Model =====================
library(corrplot)
length(SSCast_variables[which(is.na(SSCast_variables)),])
SSCast_variables<- na.exclude(SSCast_variables)
names(SSCast_variables)
vars_to_remove<- c( "ammonium_spring_SD", "ammonium_winter_SD", "ammonium_fall_SD", "ammonium_summer_SD", 
                    "nitrate_spring_SD", "nitrate_winter_SD", "nitrate_fall_SD", "nitrate_summer_SD", 
                    "nitrate_spring_SD", "nitrate_winter_SD", 
                    "turbidity_spring_SD", "turbidity_winter_SD", "turbidity_fall_SD","turbidity_summer_SD",
                    "PAR_spring_SD", "PAR_winter_SD", "PAR_fall_SD", "PAR_summer_SD", 
                    "DIC_spring_SD", "DIC_summer_SD", "DIC_fall_SD", "DIC_winter_SD",
                    "temperature_spring_SD", "temperature_summer_SD", "temperature_fall_SD", "temperature_winter_SD",
                    "alkalinity_spring_SD", "alkalinity_summer_SD", "alkalinity_fall_SD", "alkalinity_winter_SD",
                    "salinity_spring_SD", "salinity_summer_SD", "salinity_fall_SD", "salinity_winter_SD",
                    "currentDirection_spring_modal", "currentDirection_summer_min",   "currentSpeed_fall_modal",
                    "currentDirection_fall_max" ,    "currentDirection_fall_mean"   ,
                    "currentDirection_fall_min",     "currentDirection_fall_modal" ,  "currentDirection_spring_max" ,  "currentDirection_spring_mean" ,
                     "currentDirection_spring_min",   "currentDirection_summer_max" ,  "currentDirection_summer_mean",  "currentDirection_summer_modal",
                    "currentDirection_winter_max",   "currentDirection_winter_mean" , "currentDirection_winter_min" ,  "currentDirection_winter_modal",
                    "currentSpeed_summer_modal", "currentSpeed_spring_modal", "currentSpeed_winter_modal"
                    ) 
tokeep <- setdiff(colnames(SSCast_variables), vars_to_remove)
SSCast_variables <- SSCast_variables[, tokeep]


cor_matrix<- cor(SSCast_variables[, c(3:ncol(SSCast_variables))])

# Get absolute correlation matrix (ignore sign)
abs_cor <- abs(cor_matrix)

# Set diagonal to 0 to ignore self-correlation
diag(abs_cor) <- 0

# Find columns where all correlations < 0.7
keep_vars <- apply(abs_cor, 2, function(x) all(x < 0.7))

# Subset the dataframe by removing correlated variables
#  remove variables with cor > 0.7

library(caret)
# cor_matrix.a is your correlation matrix
vars_to_remove <- findCorrelation(cor_matrix, cutoff = 0.7, names = TRUE)


# Subset the original dataframe, keeping only the good ones
# vars_to_remove <- remove_highly_correlated(cor_matrix, cutoff = 0.7)
tokeep <- setdiff(colnames(SSCast_variables), vars_to_remove)
nonCorr_variables <- SSCast_variables[, tokeep]

head(nonCorr_variables)
colnames(nonCorr_variables)
# 4 September 2025 --> removed SD metrics (Low Salinity SD was being identified as limiting factor in the Juan de Fuca Strait and norhtern of Salish Sea, but this does not make sense)
# 1] "longitude"                "latitude"                 "ammonium_fall_minimum"    "ammonium_spring_minimum"  "currentSpeed_fall_max"   
# [6] "nitrate_spring_maximum"   "PAR_fall_maximum"         "turbidity_winter_minimum" "depth_m"  


cor_matrix.a<- cor(nonCorr_variables)
corrplot::corrplot(cor_matrix.a, type = 'lower', tl.col = 'black',#order = 'hclust',
                   cl.ratio = 0.15, tl.srt = 45, col = COL2('PuOr', 30), cl.cex = 0.8,  tl.cex = 0.6,
                   addCoef.col = 'black',  # Add this line to show the correlation values
                   number.cex = 0.5) #, tl.pos = "n"


nonCorr_variables<- cbind(nonCorr_variables, SSCast_variables[, 1:2])

# There are some of these still highly correlated, I removed then manually
vars_to_remove<- c("latitude", "longitude", "depth_m") #"turbidity_winter_mean", "turbidity_summer_SD", "direction_summer_SD", 
tokeep <- setdiff(colnames(nonCorr_variables), vars_to_remove)
nonCorr_variables <- nonCorr_variables[, tokeep]

cor_matrix.a<- cor(nonCorr_variables)
corrplot::corrplot(cor_matrix.a, type = 'lower', tl.col = 'black',#order = 'hclust',
                   cl.ratio = 0.15, tl.srt = 45, col = COL2('PuOr', 30), cl.cex = 0.8,  tl.cex = 0.6,
                   addCoef.col = 'black',  # Add this line to show the correlation values
                   number.cex = 0.5) #, tl.pos = "n"

colnames(nonCorr_variables)
# 4 Sep 2025
# [1] "ammonium_fall_minimum"    "ammonium_spring_minimum"  "currentSpeed_fall_max"    "nitrate_spring_maximum"   "PAR_fall_maximum"        
# [6] "turbidity_winter_minimum"

nonCorr_variables<- cbind(SSCast_variables[,1:2], nonCorr_variables)


### Select one variable from each group of correlated variables ================
corr_variables <- setdiff(colnames(SSCast_variables), c(colnames(nonCorr_variables), "depth_m"))#"latitude", "longitude", 
# corr_variables<- corr_variables[c(1:13,28:length(corr_variables))]
Corr_variables <- SSCast_variables[, corr_variables]
colnames(Corr_variables) #143 # 129 after removing current direction variables # 101 after removing SD metrics of all variables

cor_nemo<- cor_plot_function(data= Corr_variables, 
                                    path_plot= "/SDM/Variables_selection/plots",
                                    plot_name= "correlated_salishseacast_variables",
                             plot_type = "pdf",
                             pdf_height= 25, 
                             pdf_width= 43)


cor_matrix.d<- cor(Corr_variables)

cor_variables<- as.matrix(round(cor_matrix.d, 3))
dissimilarity = 1 - cor_variables
distance = as.dist(dissimilarity) 

hClust<- hclust(distance, method="average")

path_plot= "PSF/SDM/Variables_selection/plots"
plot_name2<- paste(path_plot, paste("cor_salishseacast_vars_q","clusters_v2.png", sep="_"), sep="/")

png(plot_name2,
      height = 13, width = 30, units="cm", res=300)
  
plot(hClust, main = "Clustering of Correlated Variables", xlab = "", sub = "", cex = 0.6)
abline(h = 0.3, col = "red", lty = 0.2) 
  
rect.hclust(hClust, h = 0.3, border = 2:6) #Draw cluster rectangles at same height
dev.off()


### Select one variable per group of correlated variables ======================
variables<- c(
  "currentSpeed_summer_mean",
  "nitrate_summer_minimum",
  "nitrate_winter_mean", # correlated with   "salinity_summer_minimum", and nitrate winter mean
  "temperature_summer_mean",
  "ammonium_spring_mean", 
  "PAR_summer_mean",# "", # correlated with "PAR_spring_mean",
  "ammonium_spring_mean",
  "ammonium_winter_mean",
  "PAR_summer_maximum", 
  "turbidity_summer_mean",# corr with turbidity_spring_mean
  "ammonium_summer_minimum" )# corr with ammonium_summer_mean and maximum


Corr_variables_selected<- Corr_variables[, (names(Corr_variables) %in% variables)]
colnames(Corr_variables_selected)  

setdiff(variables, names(Corr_variables_selected)) # check for any missing variable

Corr_variables_selected<- cbind(Corr_variables_selected, SSCast_variables[,1:2])

# write_csv(Corr_variables_selected, "SDM/Variables_selection/selected_variables_to_interpolations_v4.csv")


### ============================================================================
### Upload files of temperature variables and merge ============================
my_path<-("/modeled_variables_original/climatology_temp_tolerance_metrics")
setwd(my_path)
dir()

# List all CSVs in your folder
files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE) # 32 files


# Load file 1 to start building the merged dataset
df_all=read.csv(files[1])
head(df_all)
df_all<- df_all[,c(5:6,ncol(df))]
filename <- strsplit(tools::file_path_sans_ext(files[1]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(3:6)]
filename <- paste0(filename, collapse ="_")
colnames(df_all)[3]<- filename

# Merge all files by adding a column of the metric to the df_all dataset
for (i in 2:length(files)) {
  
  df= read.csv(files[i])
  df= df[,c(5:6,ncol(df))]
  
  filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
  
  if(length(strsplit(filename, "_")[[1]]) == 8){
    
    filename <- strsplit(filename, "_")[[1]][c(3:7)]
    filename <- paste0(filename, collapse ="_")
  }
  
  if(length(strsplit(filename, "_")[[1]]) == 7){
    
    filename <- paste0(strsplit(filename, "_")[[1]][c(3:6)], collapse ="_")
  }
  
  
  colnames(df)[3]<- filename
  df_all= merge(df_all, df, by=c("latitude" , "longitude"))
  
}


temperature_variables<- df_all
head(temperature_variables)

###=============================================================================
### Add depth value to variables data by location===============================
# Transform your spatial object 
coordinates(temperature_variables) <- ~longitude + latitude  # Convert to spatial object
proj4string(temperature_variables) <- CRS("+proj=longlat +datum=WGS84")  # Set CRS (WGS84)
bc_albers <- CRS("+init=epsg:3005")  # or full string:
# CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Note: If you get a warning about +init=epsg:3005 being deprecated, use the full string instead.

temperature_variables_sp <- spTransform(temperature_variables, bc_albers)
values_depth <- raster::extract(bathy_filtered, temperature_variables_sp)
temperature_variables <- cbind(as.data.frame(temperature_variables_sp), depth = values_depth)
temperature_variables<- temperature_variables[,c(33:34, 1:32, 35)]
names(temperature_variables)[ncol(temperature_variables)]<- "depth_m"
names(temperature_variables)[1:2]<- c("latitude", "longitude")

temperature_variables<- as.data.frame(temperature_variables)%>%
  filter(!is.na(depth_m))

summary(temperature_variables$depth_m)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -49.999 -15.338   3.558   7.848  32.726  69.990 


### Correlation among temperature tolerance metrics ============================
colnames(temperature_variables)
cor_temp_variables<- cor_plot_function(data= temperature_variables[,c(3:34)], 
                                       path_plot= "/SDM/Variables_selection/plots",
                                       plot_name= "cor_temp_variables")

colnames(temperature_variables)
### Select one variable per group of correlated variables 
variables<- c("summer_cumulated_degrees_18", 
              "summer_hours_above_threshold_18",
              "longitude", "latitude" 
)

temperature_variables_selected<- temperature_variables[, (names(temperature_variables) %in% variables)]
colnames(temperature_variables_selected)  

Corr_variables_selected2<- merge(Corr_variables_selected, temperature_variables_selected, by= c("latitude", "longitude"))

# Merge with other variables to intoerpolate 
# write_csv(Corr_variables_selected2, "/SDM/Variables_selection/selected_variables_to_interpolations_v2.csv")


###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
### ============================================================================
### Upload variables derived from Bathymetry and merge =========================
# (20 m res original bathymetry layer and terrain derived features)
# bathymetry was degraded to 500 m res using QGis: Raster/Projection/Warp(with average) )

my_path_bathy<- "/SDM/environmental_layers/Topographic_Variables"
setwd(my_path_bathy)  
dir()

# List all files in your folder
files <- list.files(path = my_path_bathy, pattern = "\\.tif$", full.names = TRUE) # 7 files
files<- files[-1]

rasters_bathy<- stack(paste(files[1], sep="/"))
bc_albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
crs(rasters_bathy)<- bc_albers


for (i in 2:length(files)) {
  raster_i<- raster(paste(files[i], sep="/"))
  rasters_bathy<- stack(rasters_bathy, raster_i)
}

names(rasters_bathy)
# [1] "easterness_3x3"                     "easterness_5x5"                     "easterness_7x7"                     "northerness_3x3"                   
# [5] "northerness_5x5"                    "northerness_7x7"                    "focal_sd.1"                         "roughness_5x5"                     
# [9] "focal_sd.2"                         "slope.1"                            "slope_5x5"                          "slope.2"                           
# [13] "terrain_rasters_masked_multiband_1" "TPI_3x3"                            "TPI_5x5"                            "TPI_7x7"                           
# [17] "coastwide_20m.1"                    "TRI_5x5"                            "coastwide_20m.2"       

names_files<- basename(files)
names_files<- tools::file_path_sans_ext(names_files)
names(rasters_bathy)<- names_files
rasters_bathy<- rasters_bathy[[-13]] # remove "terrain_rasters_masked_multiband"
# [1] "easterness_3x3"  "easterness_5x5"  "easterness_7x7"  "northerness_3x3" "northerness_5x5"
# [6] "northerness_7x7" "roughness_3x3"   "roughness_5x5"   "roughness_7x7"   "slope_3x3"      
# [11] "slope_5x5"       "slope_7x7"       "TPI_3x3"         "TPI_5x5"         "TPI_7x7"        
# [16] "TRI_3x3"         "TRI_5x5"         "TRI_7x7" 


# Mask to exclude deep areas 
bathy_20m_2<- raster(paste(my_path_bathy,"topographic_variables_20mres/coastwide_20m.tif", sep="/"))
bathy_20m_2[bathy_20m_2 >= 90]<- NA
bathy_20m_2[bathy_20m_2 <= 50]<- NA
plot(bathy_20m_2)

rasters_bathy_masked<- mask(rasters_bathy, bathy_20m_2)
plot(rasters_bathy_masked[[15:18]])

# Save the masked rasters
# terra::writeRaster(rasters_bathy_masked, "/SDM/Variables_selection/terrain_rasters_masked_multiband2.tif",  overwrite = TRUE)

rasters_bathy_masked<- terra::rast("/SDM/Variables_selection/terrain_rasters_masked_multiband2.tif")
names(rasters_bathy_masked)


# # Number of random points you want
n <- 10000
# Extract random sample (all layers at same pixels)
set.seed(42)  # reproducibility
sample_df <- sampleRandom(rasters_bathy_masked, size = n, na.rm = TRUE, xy = TRUE)
# Check the result
head(sample_df)
# write.csv(sample_df, "/SDM/Variables_selection/XY_topographyc_variables_sample_table.csv")


# Correlation matrix
cor_matrix <- cor(sample_df, use = "pairwise.complete.obs", method = "pearson")

cor_topographic<- cor_plot_function(data= sample_df[,c(3:ncol(sample_df))],
                                 path_plot= "/SDM/Variables_selection/plots",
                                 plot_name= "cor_topographic_vars")

# Exclude "Easterness" (correlated with Northerness) and TRI (correlated with Slope)
exclude_topo<- c("easterness_5x5", "easterness_7x7",
                 "northerness_5x5", "northerness_7x7",
                 "roughness_5x5", "roughness_7x7", "roughness_3x3",
                 "slope_3x3", "slope_7x7",
                 "TRI_3x3", "TRI_7x7",
                 "TPI_5x5", "TPI_7x7",
                 "TRI_5x5", "TRI_7x7")

keep_layers <- which(!(names(rasters_bathy_masked) %in% exclude_topo))
topo_selected <- rasters_bathy_masked[[keep_layers]]
names(topo_selected)

# [1] "easterness_3x3"  "northerness_3x3" "slope_5x5"  "TPI_3x3"

topo_selected2 <- rasters_bathy_masked[[keep_layers]]
names(topo_selected2)

# Save the masked rasters
# terra::writeRaster(topo_selected2, "/SDM/Variables_selection/terrain_rasters_selected.tif", overwrite=T)
topo_selected<- terra::rast("SDM/Variables_selection/terrain_rasters_selected.tif")
names(topo_selected)<- c( "easterness_3x3",  "northerness_3x3", "slope_5x5" ,      "TPI_3x3"   )


# Check the correlation among the selected topographic variables
# topo_selected_df<- as.data.frame(topo_selected, xy = TRUE)
set.seed(123)  # for reproducibility
topo_selected<- stack(topo_selected)
vals <- sampleRandom(topo_selected, size = 10000, na.rm = TRUE, sp = FALSE, xy=T)

# Convert to data frame
vals_df <- as.data.frame(vals)
head(vals_df)

# Correlation matrix (Pearson by default)
library(corrplot)
cor_matrix.d<- cor(vals_df)
corrplot::corrplot(cor_matrix.d, type = 'lower', tl.col = 'black',#order = 'hclust',
                   cl.ratio = 0.15, tl.srt = 45, col = COL2('PuOr', 30), cl.cex = 2,  tl.cex = 2) #, tl.pos = "n"



# ### ============================================================================
# ### Merge all selected variables and correlate them ==========
# ### ============================================================================
# select_variables<- merge(Corr_variables_selected, nonCorr_variables, by= c("longitude", "latitude"))
# select_variables<- merge(select_variables, temperature_variables_selected, by= c("longitude", "latitude"))
# # select_variables<- merge(select_variables, topo_selected, by= c("longitude", "latitude"))
# colnames(select_variables)
# 
# # We need to get values of terrain variables at coordinates of Nemo model variables
# terrain_vals<- raster::extract(topo_selected2, select_variables[,c("longitude", "latitude")])
# head(terrain_vals)
# 
# 
# test_point <- select_variables[1, 1:2]
# raster::extract(topo_selected, test_point)
# plot(topo_selected[[1]])
# points(select_variables[1:10, 1:2], col = "red", pch = 16)
# 
# # Convert dataframe to sf object, specifying coordinate columns and CRS
# sf_points <- st_as_sf(select_variables, coords = c("longitude", "latitude"), 
#                       crs = projection(topo_selected))
# 
# # Plot raster
# plot(topo_selected[[1]], main = "Raster with points")
# 
# # Add points on top
# plot(sf_points$geometry, col = "blue", pch = 16, add = TRUE)
# 
# 
# # Combine with Nemo data
# select_variables <- cbind(select_variables, terrain_vals)
# summary(select_variables)



### Correlation among all NEMO metrics ============================
names(Corr_variables_selected)
cor_selected_variables<- cor_plot_function(data= Corr_variables_selected[,-c((ncol(Corr_variables_selected)-1),ncol(Corr_variables_selected))], 
                                       path_plot= "/SDM/Variables_selection/plots",
                                       plot_name= "cor_selected_variables_noTerrain_final2",
                                       plot_type = "png")

exclude.2 = c(
              "PAR_summer_maximum"
)

select_variables<- Corr_variables_selected
selected_variables2<- select_variables[, !(names(select_variables) %in% exclude.2)]
colnames(selected_variables2)# 25 metrics

correlation_table<- cor(selected_variables2[1:9], use = "pairwise.complete.obs")
correlation_table<-as.data.frame(correlation_table)
write.csv(correlation_table, "/SDM/Variables_selection/correlation_table_nemoSelectedVariables_Sep2025.csv")

cor_selected_final<- cor_plot_function(data= selected_variables2[1:9], 
                                       # pdf_height= 25,
                                       # pdf_width= 30,
                                 path_plot= "/PSF/SDM/Variables_selection/plots",
                                 plot_name= "cor_selected_variables_noTerrain_final_v2", plot_type="png")

cor_selected_final<- cor_plot_function(data= selected_variables2[1:9], 
                                       pdf_height= 15,
                                       pdf_width= 15,
                                       path_plot= "/SDM/Variables_selection/plots",
                                       plot_name= "cor_selected_variables_noTerrain_final_v2", plot_type="pdf")


colnames(selected_variables2)
# [1] "longitude"                  "latitude"                   "ammonium_fall_mean"         "ammonium_summer_mean"      
# [5] "nitrate_summer_mean"        "PAR_spring_maximum"         "PAR_summer_mean"            "salinity_summer_mean"      
# [9] "salinity_summer_SD"         "temperature_summer_mean"    "turbidity_summer_mean"      "speedVector_summer_maximum"
# [13] "speedVector_summer_mean"    "ammonium_fall_SD"           "ammonium_spring_mean"       "nitrate_spring_SD"         
# [17] "nitrate_summer_SD"          "nitrate_winter_SD"          "direction_fall_mean"        "direction_summer_mean"     
# [21] "direction_summer_minimum"   "direction_winter_maximum"   "direction_winter_mean"   

# 4 Sep
# [1]] "ammonium_spring_mean"     "ammonium_summer_minimum"  "ammonium_winter_mean"     "currentSpeed_summer_mean" "nitrate_summer_minimum"  
# [6] "nitrate_winter_mean"      "PAR_summer_mean"          "temperature_summer_mean"  "turbidity_summer_mean"    "longitude"               
# [11] "latitude" 

###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
# write.csv(selected_variables2, "/SDM/Variables_selection/selected_variables_WOut_bathymetricVars_4Sept2025.csv")
# save.image("/SDM/Variables_selection/variables_selection_1August2025.R")



