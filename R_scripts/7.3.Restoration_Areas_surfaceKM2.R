###====================================================================================
###     Species Distribution models    SDMs                                      ######
###                                                                              ######
### 7.3- Assessing the Area for Restoration, with masking by depth and substrate ######
### Content:                                                                     ######
### Export tables of summary area (km2) of each restoration potential category   ######
### #
### Author: Romina Barbosa                                                       ######
### Date last edition: 03-Dec-2025                                               ######
###====================================================================================
# Load packages
library(terra)
library(dplyr)
library(dismo)
library(ggplot2)
library(raster)
library(tidyr)


# Paths
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
setwd(model_results_path)
postblob_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob"
out_dir <- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7"

### Load substrate layer =====================================================
substrate_aligned<- rast(paste(model_results_path, "substrate_SOG_aligned.tif", sep="/"))
substrate_aligned[substrate_aligned == 2]<- NA
substrate_aligned <- terra::mask(substrate_aligned, ens_average)

# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask10_30<- bathy20m
bathy20m_mask10_30[!(bathy20m_mask10_30 >= -10 & bathy20m_mask10_30 <= 30)] <- NA # negative values are in land

bathy20m_mask10_15<- bathy20m
bathy20m_mask10_15[!(bathy20m_mask10_15 >= -10 & bathy20m_mask10_15 <= 15)] <- NA # negative values are in land


### Summarize environmental conditions during each periods at each Restoration category ====
# Load raster stack of variables resampled at 20 m resoltuion
r<- rast( paste(out_dir, "stability_layer_categoric.tif", sep="/"))
r<- terra::mask(r, bathy20m_mask10_15)
r_df <- as.data.frame(r, xy = TRUE)  # adds x and y columns
colnames(r_df)[3] <- "category"
r_df$category <- as.factor(r_df$category)


raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_blob.tif")
names(raster_stack_20m)

cummulated_18SST<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated_Blob/summer_cumulated_degrees_18_interpolated_output.tif")
hourAbove_18SST<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated_Blob/summer_hours_above_threshold_18_interpolated_output.tif")


raster_stack_20m_masked <- terra::mask(raster_stack_20m, bathy20m_mask10_30)
plot(raster_stack_20m_masked[[1]])
# writeRaster(raster_stack_20m_masked, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/raster_stack_blob_masked.tif", overwrite=T)

# raster_stack_postblob<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblob.tif")
names(raster_stack_postblob)
raster_stack_postblob_masked <- terra::mask(raster_stack_postblob, bathy20m_mask10_30)
plot(raster_stack_postblob_masked[[1]])
# writeRaster(raster_stack_postblob_masked, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/raster_stack_postblob_masked.tif", overwrite=T)


# Extract env conditions at each cell and merge with Restoration category
# Select the period, period 1 or period 2 (commented lines)
df_env <- cbind(
  r_df,
  terra::extract(raster_stack_20m, r_df[, c("x", "y")])
)

# df_env <- cbind(
#   r_df,
#   terra::extract(raster_stack_postblob, r_df[, c("x", "y")])
# )

# Remove NA classes
df_env <- df_env[!is.na(df_env$category), ]
df_env$class <- as.factor(df_env$category)


### Add depth and susbstrate to calculate change in area (surface in m2) ======
r_df_depth <- cbind(
  r_df,
  terra::extract(bathy20m_mask10_15, r_df[, c("x", "y")])
)

substrate<- rast(paste(model_results_path, "substrate_SOG_aligned_.tif", sep="/"))
substrate <- terra::mask(substrate, bathy20m_mask10_15)


r_df_substrate <- cbind(
  r_df,
  terra::extract(substrate, r_df[, c("x", "y")])
)

levels(r_df_2$substrate_aligned) 
summary(r_df_2$substrate_aligned)

r_df_2<- merge(r_df_substrate, df_env, by= c("x", "y", "category"))
head(r_df_2)

r_df_2<- merge(r_df_2, r_df_depth, by= c("x", "y", "category"))

# create categories of susbtrate and depth features
r_df_2$depth_categ<- r_df_2$coastwide_20m
r_df_2[which(r_df_2$coastwide_20m <=15 & r_df_2$coastwide_20m >= -3 ), "depth_categ"]<- "suitable_depth"
r_df_2[which(r_df_2$coastwide_20m >15), "depth_categ"]<- "unsuitable_depth"
r_df_2[which(r_df_2$coastwide_20m <3), "depth_categ"]<- "unsuitable_depth"


r_df_2_ras <- rast(
  r_df_2[,c(1,2,3,17,18,5)],
  type = "xyz",
  crs  = "EPSG:3005"     # your CRS here
)

plot(r_df_2_ras)



area_summary<- r_df_2[,c(1,2,3,17,18,5)]%>%
  group_by(category)%>%
  summarize(area_km2= length(x)*0.02*0.02)

r_df_2$depth_categ<- as.factor(r_df_2$depth_categ)
area_summary_depth<- r_df_2[,c(1,2,3,17,18,5)]%>%
  filter(depth_categ == "suitable_depth")%>%
  group_by(category)%>%
  summarize(area_km2_maskDepth= length(x)*0.02*0.02)

area_summary_depthSubstrate<- r_df_2[,c(1,2,3,17,18,5)]%>%
  filter(depth_categ == "suitable_depth")%>%
  filter(substrate_aligned == "1")%>%
  group_by(category)%>%
  summarize(area_km2_maskDepthSubstrate= length(x)*0.02*0.02)

final_summary <- area_summary %>%
  left_join(area_summary_depth, by = "category") %>%
  left_join(area_summary_depthSubstrate, by = "category")

final_summary

# Add percentage relative to the total area before masking
final_summary <- final_summary %>%
  mutate(
    perc_depth = area_km2_maskDepth / area_km2 * 100,
    perc_depthSubstrate = area_km2_maskDepthSubstrate / area_km2 * 100,
    perc_depthSubstrate_ofSuitableDepth = area_km2_maskDepthSubstrate / area_km2_maskDepth * 100
  )

levels(final_summary$category)<- c("Non-viable", "Low-Potential", "Medium-Potential", "High-Potential" )


print(final_summary)

total_area<- sum(final_summary$area_km2)
final_summary$area_km2/total_area
mean(final_summary$perc_depthSubstrate_ofSuitableDepth)


