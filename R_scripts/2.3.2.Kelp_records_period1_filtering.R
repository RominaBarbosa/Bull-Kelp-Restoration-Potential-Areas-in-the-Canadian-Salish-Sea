###==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
### 2- Presence, speudo-absence and background data Preparation ################
### Author: Romina Barbosa                                      ################
### Date last edition: 8-Sep-2025                                           ################
### Code edited from: https://sjevelazco.github.io/flexsdm/articles/v01_pre_modeling.html#data-species-occurrence-and-background-data ####
###==================================================================
library(raster)
library(dplyr)
library(terra)
library(sf)
library(flexsdm)
library(ggplot2)
library(tools)
options(digits=13)

layers_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/clipped_variables"
SDM_path<- "/Volumes/Romina_PSF/PSF/SDM"

plotspath<- "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Plots"
mypath<- "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp"


###==============================================================================
## Load occurrence data and exclude absences from shoreline surveys =============
###==============================================================================
# occ_data_sdm<- read.csv( paste(SDM_path, "Presence_absences_kelp/presence_2014_to_2019_filtered_plus_AbsencesFINAL.csv", sep="/"))
occ_data_sdm<- read.csv( paste(SDM_path, "Presence_absences_kelp/presence_2014_to_2019_filtered_plus_AbsencesFINAL_cluster.csv", sep="/"))
# columns<- c("latitude", "longitude", "Class_name", "Region", "Year", "Area_m2", "layer")
# occ_data_sdm<- occ_data_sdm[,columns]
head(occ_data_sdm)
summary(occ_data_sdm)

occ_data_sdm%>%group_by(Cluster)%>%
  summarize(n_Cluster= length(Cluster))


substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/substrate_SOG_aligned_.tif")
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
names(bathy20m)<- "bathy20m"

kelp_points <- vect(occ_data_sdm, geom = c("x", "y"))
crs(kelp_points) <- "EPSG:3005"

extracted_values <- terra::extract(substrate, kelp_points)
unique(as.numeric(extracted_values$substrate_aligned))
depth_values <- terra::extract(bathy20m, kelp_points)

# Combine extracted values with kelp coordinates (optional)
occ_data_sdm <- cbind(occ_data_sdm, extracted_values[, -1])  # remove ID column from extract
colnames(occ_data_sdm)[73]<- "substrate"
occ_data_sdm <- cbind(occ_data_sdm, depth_values[, -1])  # remove ID column from extract
colnames(occ_data_sdm)[74]<- "depth"


# Explore number of points with NAs (Presences and absences)
occ_data_sdm$substrate<- as.factor(occ_data_sdm$substrate)
summary(occ_data_sdm)
substrate_NAs<- occ_data_sdm[which(is.na(occ_data_sdm$substrate)),]

substrate_NAs%>%group_by(kelp, Cluster)%>%
  summarize(n=length(kelp))
#    kelp Cluster     n
# 1     0       1    68
# 2     0       2   118
# 3     0       3   252
# 4     0       4   311
# 5     0       5     6
# 6     0      NA     1
# 7     1       1     1
# 8     1       2     3
# 9     1       3     2


length(substrate_NAs$field_1)
# [1] 762


# There were few presences without substrate info (6 records of presence and many of absence) --> removed
occ_data_sdm<- occ_data_sdm[which(!is.na(occ_data_sdm$substrate)),]

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
# 
# 
as.data.frame(occ_data_sdm)%>%
  group_by(Cluster, kelp)%>%
  summarize(n_occ= length(Year))

#    Cluster  kelp n_occ
# 1       1     0  1062
# 2       1     1   408
# 3       2     0  1906
# 4       2     1   665
# 5       3     0  2876
# 6       3     1   155
# 7       4     0  2030
# 8       4     1    53
# 9       5     0   130

#    kelp n_occ
# 1     0  8004
# 2     1  1281

## %######################################################%##
#                                                          #
####          Filtering occurrences                     ####
#               by distance and depth                      #
## %######################################################%##
options(digits=13)
# ### Load env layers 
# source("/Volumes/Romina_PSF/PSF/R_scripts/stack_rasters_path_function.R")
# # stack_vars<- stack_rasters_path_funtion(layers_path, terra_class= "Y")
# 
# raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_FINAL2.tif")
# names(raster_stack_20m)
# terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")
# tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)
# terrain_vars<-  rast(tif_files[c(12,15)]) #,5,12,13,15
# names(terrain_vars)
# summary(terrain_vars)
# 
# # Merge rasters of all selected variables including terrain and NEMO
# raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)
# 
# # Create Mask from bathymetry
# bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
# bathy20m_mask15_40<- bathy20m
# 
# 
# bathy20m_mask15_40[!(bathy20m_mask15_40 >= -15 & bathy20m_mask15_40 <= 40)] <- NA # negative values are in land
# plot(bathy20m_mask15_40)
# 
# bathymetry_20m<- crop(bathymetry_20m, stack_vars)
# bathy20m



### Filter and Determine the depth of the 99% of kelp presence records ====================
summary(occ_data_sdm)

filtered_occ <- occ_data_sdm%>%
  filter(depth >= -15 & depth <= 40)

as.data.frame(filtered_occ)%>%
  group_by(Cluster, kelp)%>%
  summarize(n_occ= length(Year))


# Cluster  kelp n_occ
# 1       1     0   523
# 2       1     1   407
# 3       2     0  1301
# 4       2     1   657
# 5       3     0  1983
# 6       3     1   153
# 7       4     0  1592
# 8       4     1    53
# 9       5     0   123

# kelp n_occ
# 1     0  5522
# 2     1  1270


# Nov 2025
# Cluster  kelp n_occ
# 1       1     0   493
# 2       1     1   407
# 3       2     0  1247
# 4       2     1   657
# 5       3     0  1842
# 6       3     1   153
# 7       4     0  1405
# 8       4     1    53
# 9       5     0   117


Q99_kelp_depth<- quantile(filtered_occ[which(filtered_occ$kelp==1),"depth"], probs =0.99, na.rm = TRUE)
# 95% =  9.97939825058 
# 99% =  21.51443481445 
# 1% = -6.653466176987 
# 5% = -3.00764799118 

# September 
# 90% = 8.472120285034
# 95% = 11.71061553955
# 99% = 17.63191188812 
# 98% = 15.86030252457

# Nov 
# 99% = 17.65239244461 
# 98% = 15.86030252457 
# 90% = 8.640222263336 


Q99_kelp_rocky_depth<- quantile(filtered_occ[which(filtered_occ$kelp==1 & filtered_occ$substrate==1),"depth"], probs =0.99, na.rm = TRUE)
# 17.48803371429 
# 98% 
# 15.75731815338 


as.data.frame(filtered_occ)%>%
  group_by(Cluster, kelp, substrate)%>%
  summarize(n_occ= length(Year))
#   Cluster  kelp substrate n_occ
# 1       1     0 1           345
# 2       1     0 2           148
# 3       1     1 1           391
# 4       1     1 2            16
# 5       2     0 1           497
# 6       2     0 2           750
# 7       2     1 1           609
# 8       2     1 2            48
# 9       3     0 1          1042
# 10       3     0 2           800
# 11       3     1 1           148
# 12       3     1 2             5
# 13       4     0 1           682
# 14       4     0 2           723
# 15       4     1 1            43
# 16       4     1 2            10
# 17       5     0 1            41
# 18       5     0 2            76
summary(filtered_occ)

##========================================================================================================================
### Select records with soft substrate and filter per distance (500 m layers' res) for validation of models  ==============
filtered_occ_soft<- filtered_occ%>%
  filter(substrate == 2)

filtered_occ_soft%>%group_by(kelp, substrate)%>%
  summarize(n_records= length(kelp))
# kelp substrate n_records
# 1     0 2              2497
# 2     1 2                79

kelp_vect_soft <- vect(filtered_occ_soft, geom = c("x", "y"))
crs(kelp_vect_soft) <- "EPSG:3005"# Convert sf to SpatVector (terra uses this)

# Create a raster of cell IDs
grid_500m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/bathy_coastwide_500m.tif")
grid_id_raster <- setValues(grid_500m, 1:ncell(grid_500m))
names(grid_id_raster) <- "grid_ID"

# Extract grid_ID values at point locations
grid_ids_all <- terra::extract( grid_id_raster, kelp_vect_soft)
length(grid_ids_all$grid_ID)
length(kelp_vect_soft$kelp)
length(filtered_occ_soft$field_1)
length(grid_ids_all["grid_ID"])

# Combine with point data
points_with_grid_s <- cbind(filtered_occ_soft, grid_ids_all[,"grid_ID"])
plot(kelp_vect_soft)

# Extract coordinates
coords_soft <- crds(points_with_grid_s)

# Combine with attribute data
points_with_grid_s <- cbind(as.data.frame(points_with_grid_s), coords_soft)
colnames(points_with_grid_s)[75]<- "grid_ID"

## Keep only one point per cell, with preference for kelp presence, i.e., if there is any presence record it will be keept to represent the cell
length(points_with_grid_s$grid_ID) #[1] 2576
length(unique(points_with_grid_s$grid_ID)) # 2576

points_with_grid_s[which(points_with_grid_s$kelp==1),]

df_filtered_soft <- as.data.frame(points_with_grid_s) %>%
  group_by(grid_ID) %>%
  # arrange so that kelp == 1 comes first
  arrange(desc(kelp), .by_group = TRUE) %>%
  # keep only the first row per group
  slice(1) %>%
  ungroup()


length(df_filtered_soft$kelp) #2576
length(unique(df_filtered_soft$grid_ID))# 2576

# write.csv(df_filtered_soft, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&softSubstrate.csv")

# Plot dataset 2025-2019
df_filtered_soft.st<- st_as_sf(df_filtered_soft, coords = c("x", "y"), crs = 3005)

plot(df_filtered_soft.st[11], cex=0.1) # presence/absence
plot(df_filtered_soft.st[71], cex=0.1) # Substrate

df_filtered_soft%>%group_by(kelp, substrate)%>%
  summarize(n_records= length(kelp))
#    kelp substrate n_records
# 1     0 2              2497
# 2     1 2                79

##==============================================================================
### Select records with hard substrate and filter per distance for calibration of models  ==============
filtered_occ_rocky<- filtered_occ%>%
  filter(substrate == 1)

as.data.frame(filtered_occ_rocky)%>%
  group_by(Cluster, kelp)%>%
  summarize(n_occ= length(Year))

#  Cluster  kelp n_occ
# 1       1     0   775
# 2       1     1   391
# 3       2     0   698
# 4       2     1   617
# 5       3     0  1298
# 6       3     1   150
# 7       4     0   910
# 8       4     1    43
# 9       5     0    43

# kelp n_occ
# 1     0  3724
# 2     1  1201


# Nov 2025
# Cluster  kelp n_occ
# 1       1     0   345
# 2       1     1   391
# 3       2     0   497
# 4       2     1   609
# 5       3     0  1042
# 6       3     1   148
# 7       4     0   682
# 8       4     1    43
# 9       5     0    41

# kelp n_occ
# 1     0  2607
# 2     1  1191


## Filtering kelp presences to keep only one per cell/grid (500 m resolution layers) ===========
kelp_vect <- vect(filtered_occ_rocky, geom = c("x", "y"))
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
plot(kelp_vect)
head(points_with_grid)

# Extract coordinates
coords <- crds(points_with_grid)

# Combine with attribute data
points_with_grid <- cbind(as.data.frame(points_with_grid), coords)


## Keep only one point per cell, with preference for kelp presence, i.e., if there is any presence record it will be keept to represent the cell
length(points_with_grid$grid_ID) #[1] 1287  -- 3794 (with absences)  # Nov 25 [1] 3798
length(unique(points_with_grid$grid_ID)) # [1] 3798

points_with_grid[which(points_with_grid$kelp==1),]

df_filtered <- as.data.frame(points_with_grid) %>%
  group_by(grid_ID) %>%
  # arrange so that kelp == 1 comes first
  arrange(desc(kelp), .by_group = TRUE) %>%
  # keep only the first row per group
  slice(1) %>%
  ungroup()


length(df_filtered$kelp) #3794  # Nov [1] 3798
length(unique(df_filtered$grid_ID))# [1] 1287  -- 3794 (with absences)  # Nov 25 [1] 3798

as.data.frame(df_filtered)%>%
  group_by(Cluster, kelp)%>%
  summarize(n_occ= length(Year))
# 1       1     0   345
# 2       1     1   391
# 3       2     0   497
# 4       2     1   609
# 5       3     0  1042
# 6       3     1   148
# 7       4     0   682
# 8       4     1    43
# 9       5     0    41

as.data.frame(df_filtered)%>%
     group_by(kelp)%>%
     summarize(n_occ= length(Year))
# kelp n_occ
# 1     0  2607
# 2     1  1191


# Plot dataset 2025-2019
df_filtered.st<- st_as_sf(df_filtered, coords = c("x", "y"), crs = 3005)
df.filtered.ll <- st_transform(df_filtered.st, 4326)

plot(df.filtered.ll[11], cex=0.1) # presence/absence
plot(df.filtered.ll[71], cex=0.1) # Substrate


# There are no duplicated cell ID between soft and hard substrate --> because the absences were created at 500 m resolution so there is no overlap
df_filtered_entire<- rbind(df_filtered, df_filtered_soft)
df_filtered_entire<- df_filtered_entire%>%
  dplyr::select(cell_id, Year, kelp, substrate, x, y)

length((df_filtered_entire$cell_id))
length(unique(df_filtered_entire$cell_id))
df_filtered_entire[which(duplicated(df_filtered_entire$cell_id)), "cell_id"]

# # Calculate the proportion of area occupied by each environmental cluster ====
# clusters<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/clusters_MoraSoto2024_3005.csv")
# # clusters_pts <- st_as_sf(clusters[,c(31:33)], coords = c("lon_3005", "lat_3005"), crs = 3005)  # or your spatial data frame
# 
# cl_rel_area<- as.data.frame(clusters) %>%
#   group_by(Cluster) %>% 
#   summarise(n_segments = length(Cluster)) # same than using n()
# 
# # Cluster    n_segments
# # 1       1        430 
# # 2       2        864
# # 3       3       1137
# # 4       4       1753
# # 5       5        124
# # 6      NA        419
# 
# 
# ## Add information about the cluster area/proportion of total study area =======
# cl_rel_area<- na.exclude(cl_rel_area)
# total_segments<- sum(as.numeric(cl_rel_area$n_segments)) 
# cl_rel_area$prop_area<- cl_rel_area$n_segments / total_segments * 100
# #   Cluster n_segments prop_area
# # 1       1        430      9.98
# # 2       2        864     20.1 
# # 3       3       1137     26.4 
# # 4       4       1753     40.7 
# # 5       5        124      2.88
# 
# 
# colnames(cl_rel_area)[1]<- "Cluster"
# 
# # Merge kelp data with information of cluster area 
# kelp_unique<- df_filtered%>%select(kelp, Tidal_clas, Fetch_clas, SpTSM_clas, SumTSM_cla, Spring_mea,  Summer_mea, Cluster,
#                                    lat_3005,       lon_3005, substrate,   depth, grid_ID,       x,      y)
# 
# occ_data_cl<- merge(kelp_unique, cl_rel_area, by="Cluster") 
# colnames(occ_data_cl)
# 
# 
# ### Summarize number of kelp record per Environmental Cluster ==================
# occ_data_cl_summary<- occ_data_cl %>%
#   group_by( kelp,Cluster, prop_area) %>% 
#   summarise(n_occ = n())
# # kelp Cluster prop_area n_occ
# # 1     0       1      9.98   345
# # 2     0       2     20.1    494
# # 3     0       3     26.4   1042
# # 4     0       4     40.7    682
# # 5     0       5      2.88    41
# # 6     1       1      9.98   390
# # 7     1       2     20.1    609
# # 8     1       3     26.4    148
# # 9     1       4     40.7     43
# 
# 
# ggplot(occ_data_cl_summary, aes(x=(prop_area), y= n_occ, color= as.factor(Cluster)))+
#   geom_point()+
#   facet_wrap(~ kelp)
# 
# # ### Compare entire kelp records dataset and the filtered one in the map ========
# # # Final dataset
# # occ_data_cl_sf <- occ_data_cl %>%
# #   st_as_sf(coords = c("x", "y"), crs = 3005)
# # mapview::mapview(occ_data_cl_sf, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check
# # 
# # # Initial dataset 
# # mapview::mapview(kelp, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check
# 
# 
# 
# ### Plot the relationship of number of kelp presences by Cluster
# all_pres_plot<- occ_data_cl%>%
#   group_by(Cluster, kelp, prop_area) %>% 
#   summarise(n_occ = n())%>%
#   ggplot(aes(x= prop_area, y= n_occ, color=as.factor(Cluster)))+
#   geom_point(size=3)+
#   facet_wrap(~kelp)+
#   # ylim(c(0, 700))+
#   # geom_abline(intercept = 0, slope = 1, size = 0.5, color="grey") +
#   labs(x="Proportion of area (%)", y="Number of records", color="Cluster")+
#   theme_bw()#+ theme(legend.position="NULL")
# 
# all_pres_plot

# write.csv(occ_data_cl, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&rocky.csv")

occ_data_cl0<- df_filtered
## %######################################################%##
#                                                          #
####          Modeling with weighted occurrences        ####
#           by environmental cluster                       #
## %######################################################%##
# occ_data_cl<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&rocky.csv")
# summary(occ_data_cl)
# 
# press<- occ_data_cl%>%
#   filter(kelp==1)
# 
# abs<- occ_data_cl%>%
#   filter(kelp==0)%>%
#   sample(size = length(press$kelp))
#   
# # frequency per cluster
# freq <- table(occ_data_cl$Cluster)
# occ_data_cl$raw_weight <- 1 / as.numeric(freq[as.character(occ_data_cl$Cluster)])  # inverse freq
# 
# # Normalize weights: make mean weight = 1 (or sum = nrow)
# occ_data_cl$weight <- occ_data_cl$raw_weight * (length(occ_data_cl$raw_weight) / sum(occ_data_cl$raw_weight))
# 
# 
# 
# plot(occ_data_cl$weight, occ_data_cl$incl_prob)
# 
# min_per_cluster <- 43  # target samples per cluster per class = lower n in cluster 4
# 
# train_data_equal <- occ_data_cl %>%
#   group_by(kelp, Cluster) %>%
#   sample_n(min_per_cluster, replace = TRUE) %>%
#   ungroup()
# 
# train_data_equal %>%
#   group_by(kelp, Cluster) %>%
#   summarize(length(kelp))
# 
# 
# plot(vect(train_data_equal, geom= c("x","y")))
# 
# # write.csv(train_data_equal, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_filteredCluster.csv")




