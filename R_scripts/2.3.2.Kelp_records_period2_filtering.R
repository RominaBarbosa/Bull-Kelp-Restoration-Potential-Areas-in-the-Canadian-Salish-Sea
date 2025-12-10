###==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
###  Presence and speudo-absence data Preparation               ################
### Author: Romina Barbosa                                      ################
### Date last edition: 1-Oct-2025                                           ################
###==================================================================
library(raster)
library(dplyr)
library(terra)
library(sf)
library(flexsdm)
library(ggplot2)
library(tools)
options(digits=13)

# layers_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/clipped_variables"
# SDM_path<- "/Volumes/Romina_PSF/PSF/SDM"
# 
# plotspath<- "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Plots"
# mypath<- "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp"


###==============================================================================
## Load  data  =============
absences_sat<- read.csv( paste( "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/MoraSoto_postblob/absence_points_2020_2022_500m_.csv", sep="/"))
colnames(absences_sat)
summary(absences_sat$bathy1)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
# -9.943054199  2.465224028 13.693238258 14.652496671 26.409198761 39.993675232 

summary(as.factor(absences_sat$layer))
# Cover_2020     Cover_2021 Cover_SVI_2020 Cover_SVI_2021 Cover_SVI_2022 
# 1588            668           2214           1538           2106 

absences_sat<- absences_sat%>%
  select(x,y,Year,Date)

absences_sat$kelp<- 0

precense_sat<- read.csv( paste( "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/MoraSoto_postblob/Kelp_postblob_centroids_merged.csv", sep="/"))
colnames(precense_sat)
precense_sat<- precense_sat%>%
  dplyr::select(x_,y,Year,Date)

precense_sat$kelp<- 1
summary(precense_sat)

colnames(precense_sat)[1]<- "x"

# Merge presences and abcenses
occ_data_sdm<- rbind(absences_sat, precense_sat)
summary(as.factor(occ_data_sdm$kelp))

substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/substrate_SOG_aligned_.tif")
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
names(bathy20m)<- "bathy20m"


kelp_points <- vect(occ_data_sdm, geom = c("x", "y"))
crs(kelp_points) <- "EPSG:3005"
plot(kelp_points, cex=0.2)

extracted_values <- terra::extract(substrate, kelp_points)
extracted_values2 <- terra::extract(bathy20m, kelp_points)


# Combine extracted values with kelp coordinates (optional)
occ_data_sdm <- cbind(occ_data_sdm, extracted_values[, -1])  # remove ID column from extract
occ_data_sdm <- cbind(occ_data_sdm, extracted_values2[, -1])  # remove ID column from extract
colnames(occ_data_sdm)[6:7]<- c("substrate", "depth")

occ_data_sdm%>%
  group_by(kelp, substrate)%>% summarize(n_records= length(kelp))

#     kelp substrate n_records
# 1     0 1              4316
# 2     0 2              3618
# 3     0 NA              180
# 4     1 1              7266
# 5     1 2               270


# Explore number of points with NAs (Presences and absences)
occ_data_sdm$substrate<- as.factor(occ_data_sdm$substrate)
summary(occ_data_sdm)
substrate_NAs<- occ_data_sdm[which(is.na(occ_data_sdm$substrate)),]
substrate_NAs%>%
  group_by(kelp, substrate)%>% summarize(n_records= length(kelp))# only absences without substrate info
#    kelp substrate n_records
# 1     0 NA             4454
# 2     1 NA              310


# Nov 2025 (seems tthat Substrate layer only had hard substrate, I corrected it)
#      kelp substrate n_records
#   1     0 NA              180

occ_data_sdm<- occ_data_sdm%>%
  filter(!is.na(substrate)) # remove presnce and absences without substrate info

kelp_points <- vect(occ_data_sdm, geom = c("x", "y"))
crs(kelp_points) <- "EPSG:3005"
plot(kelp_points, cex=0.2, col="red", add=T)


# ### Extracted cluster values by proximity  =====================================
# library(sf)
# library(nngeo) 
# clusters<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/clusters_MoraSoto2024_3005.csv")
# clusters <- st_as_sf(clusters[,c(31:33)], coords = c("lon_3005", "lat_3005"), crs = 3005)  # or your spatial data frame
# 
# # # Find nearest cluster for each kelp point --> I did this in QGis Sep 2025
# # nearest_ix <- st_nn(kelp, clusters, k = 1, returnDist = F ) # maxdist = 1000   ---Ensure both are in same CRS (important!)
# # # Attach cluster data to kelp
# # kelp$env_cluster <- clusters$Cluster[unlist(nearest_ix)]


## %######################################################%##
#                                                          #
####          Filtering occurrences                     ####
## to keep only one per cell/grid (500 m resolution layers) ====================================
kelp_vect <- vect(occ_data_sdm, geom = c("x", "y"))
crs(kelp_vect) <- "EPSG:3005"# Convert sf to SpatVector (terra uses this)
mapview::mapview(kelp_vect)

# Create a raster of cell IDs
grid_500m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
grid_id_raster <- setValues(grid_500m, 1:ncell(grid_500m))
names(grid_id_raster) <- "grid_ID"

# Extract grid_ID values at point locations
grid_ids <- extract(grid_id_raster, kelp_vect)

# Combine with point data
points_with_grid <- cbind(kelp_vect, grid_ids["grid_ID"])
plot(kelp_vect, cex=0.2, col="blue")
head(points_with_grid)

# Extract coordinates
coords <- crds(points_with_grid)

# Combine with attribute data
points_with_grid <- cbind(as.data.frame(points_with_grid), coords)


## Keep only one point per cell, with preference for kelp presence, i.e., if there is any presence record it will be keept to represent the cell
length(points_with_grid$grid_ID) #[1] 15470 (with absences)
length(points_with_grid[which(points_with_grid$kelp==1),"kelp"]) #7536 presences
length(unique(points_with_grid[which(points_with_grid$kelp==1),"grid_ID"])) # 909 unique cells with presences

df_filtered <- as.data.frame(points_with_grid) %>%
  group_by(grid_ID) %>%
  # arrange so that kelp == 1 comes first
  arrange(desc(kelp), .by_group = TRUE) %>%
  # keep only the first row per group
  slice(1) %>%
  ungroup()


summary(as.factor(df_filtered$kelp)) # total 4604
# 0    1 
# 3695  909 

df_filtered%>%
  group_by(kelp, substrate)%>%
  summarize(n_records= length(substrate))
#    kelp substrate n_records
# 1     0 1              1795
# 2     0 2              1900
# 3     1 1               863
# 4     1 2                46

summary(as.factor(df_filtered$substrate))
# 1    2 
# 2658 1946 

length(unique(df_filtered$grid_ID))# [1] 4604 
# write.csv(df_filtered, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/MoraSoto_postblob/Presence_absences_kelp_2020_2022_filtered.csv")



# Initial dataset 
df_filtered.st<- st_as_sf(df_filtered, coords = c("x", "y"), crs = 3005)
df.filtered.ll <- st_transform(df_filtered.st, 4326)

# plot(df.filtered.ll[3])
plot(df_filtered.st[3], cex=0.1, col="red", add=T)


df_filtered%>%
  filter(kelp=="1")%>%
  filter(depth<=-10)%>%
  summarize(n= length(depth))

df_filtered%>%
  filter(kelp=="1")%>%
  ggplot(aes(x=kelp, y=depth))+ 
  geom_boxplot()

plot(df_filtered$kelp,df_filtered$depth)




