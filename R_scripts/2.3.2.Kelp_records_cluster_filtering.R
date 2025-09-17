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

substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/substrate_SOG_aligned.tif")
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
names(bathy20m)<- "bathy20m"

kelp_points <- vect(occ_data_sdm, geom = c("x", "y"))
crs(kelp_points) <- "EPSG:3005"

extracted_values <- terra::extract(substrate, kelp_points)
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

# There were few presences without substrate info (6 records) --> removed
substrate_NAs%>%group_by(kelp, Cluster)%>%
  summarize(n=length(kelp))
# kelp Cluster     n
# <int>   <int> <int>
#   1     0       1    68
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

# 
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
# 1       1     0  1130
# 2       1     1   409
# 3       2     0  2024
# 4       2     1   668
# 5       3     0  3128
# 6       3     1   157
# 7       4     0  2341
# 8       4     1    53
# 9       5     0   136
# 10      NA     0     1


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

Q99_kelp_depth<- quantile(filtered_occ[which(filtered_occ$kelp==1),"depth"], probs =0.99, na.rm = TRUE)
# 95% =  9.97939825058 
# 99% =  21.51443481445 
# 1% = -6.653466176987 
# 5% = -3.00764799118 

# September 
# 90% = 8.472120285034
# 95% = 11.71061553955
# 99% = 17.63191188812 

Q99_kelp_rocky_depth<- quantile(filtered_occ[which(filtered_occ$kelp==1 & filtered_occ$substrate==1),"depth"], probs =0.99, na.rm = TRUE)
summary(filtered_occ)

filtered_occ_rocky<- filtered_occ%>%
  filter(substrate == 1)

as.data.frame(filtered_occ_rocky)%>%
  group_by(Cluster, kelp)%>%
  summarize(n_occ= length(Year))

# Cluster  kelp n_occ
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



## Filtering kelp presences to keep only one per cell/grid (500 m resolution layers) ====================================
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
length(points_with_grid$grid_ID) #[1] 1287  -- 3794 (with absences)
points_with_grid[which(points_with_grid$kelp==1),]

df_filtered <- as.data.frame(points_with_grid) %>%
  group_by(grid_ID) %>%
  # arrange so that kelp == 1 comes first
  arrange(desc(kelp), .by_group = TRUE) %>%
  # keep only the first row per group
  slice(1) %>%
  ungroup()


length(df_filtered$kelp) #3794
length(unique(df_filtered$grid_ID))# [1] 1287  -- 3794 (with absences)



# Initial dataset 
df_filtered.st<- st_as_sf(df_filtered, coords = c("x", "y"), crs = 3005)
df.filtered.ll <- st_transform(df_filtered.st, 4326)

mapview::mapview(df.filtered.ll, zcol= "Cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check



# Calculate the proportion of area occupied by each environmental cluster ====
clusters<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/clusters_MoraSoto2024_3005.csv")
# clusters_pts <- st_as_sf(clusters[,c(31:33)], coords = c("lon_3005", "lat_3005"), crs = 3005)  # or your spatial data frame

cl_rel_area<- as.data.frame(clusters) %>%
  group_by(Cluster) %>% 
  summarise(n_segments = length(Cluster)) # same than using n()

# Cluster    n_segments
# 1       1        430 
# 2       2        864
# 3       3       1137
# 4       4       1753
# 5       5        124
# 6      NA        419


## Add information about the cluster area/proportion of total study area =======
cl_rel_area<- na.exclude(cl_rel_area)
total_segments<- sum(as.numeric(cl_rel_area$n_segments)) 
cl_rel_area$prop_area<- cl_rel_area$n_segments / total_segments * 100
#   Cluster n_segments prop_area
# 1       1        430      9.98
# 2       2        864     20.1 
# 3       3       1137     26.4 
# 4       4       1753     40.7 
# 5       5        124      2.88


colnames(cl_rel_area)[1]<- "Cluster"

# Merge kelp data with information of cluster area 
kelp_unique<- df_filtered%>%select(kelp, Tidal_clas, Fetch_clas, SpTSM_clas, SumTSM_cla, Spring_mea,  Summer_mea, Cluster,
                                   lat_3005,       lon_3005, substrate,   depth, grid_ID,       x,      y)

occ_data_cl<- merge(kelp_unique, cl_rel_area, by="Cluster") 
colnames(occ_data_cl)


### Summarize number of kelp record per Environmental Cluster ==================
occ_data_cl_summary<- occ_data_cl %>%
  group_by( kelp,Cluster, prop_area) %>% 
  summarise(n_occ = n())
# kelp Cluster prop_area n_occ
# 1     0       1      9.98   345
# 2     0       2     20.1    494
# 3     0       3     26.4   1042
# 4     0       4     40.7    682
# 5     0       5      2.88    41
# 6     1       1      9.98   390
# 7     1       2     20.1    609
# 8     1       3     26.4    148
# 9     1       4     40.7     43


ggplot(occ_data_cl_summary, aes(x=(prop_area), y= n_occ, color= as.factor(Cluster)))+
  geom_point()+
  facet_wrap(~ kelp)

# ### Compare entire kelp records dataset and the filtered one in the map ========
# # Final dataset
# occ_data_cl_sf <- occ_data_cl %>%
#   st_as_sf(coords = c("x", "y"), crs = 3005)
# mapview::mapview(occ_data_cl_sf, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check
# 
# # Initial dataset 
# mapview::mapview(kelp, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check



### Plot the relationship of number of kelp presences by Cluster
all_pres_plot<- occ_data_cl%>%
  group_by(Cluster, kelp, prop_area) %>% 
  summarise(n_occ = n())%>%
  ggplot(aes(x= prop_area, y= n_occ, color=as.factor(Cluster)))+
  geom_point(size=3)+
  facet_wrap(~kelp)+
  # ylim(c(0, 700))+
  # geom_abline(intercept = 0, slope = 1, size = 0.5, color="grey") +
  labs(x="Proportion of area (%)", y="Number of records", color="Cluster")+
  theme_bw()#+ theme(legend.position="NULL")

all_pres_plot

# write.csv(occ_data_cl, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&rocky.csv")


# ## %######################################################%##
# #                                                          #
# ####          Filtering occurrences                     ####
# #           by environmental cluster                       #
# ## %######################################################%##
# # 1. Define exponent for weighting
# a <- 0.5  # adjust based on how strong the correction should be; 1 gives equal weight to all clusters
# 
# # 2. Compute area-based weights per cluster
# cluster_weights <- (1 / unique(occ_data_cl$prop_area))^a
# names(cluster_weights) <- unique(occ_data_cl$Cluster)
# 
# # 3. Rescale weights to 0–1
# min_wt <- min(cluster_weights)
# max_wt <- max(cluster_weights)
# cluster_probs <- (cluster_weights - min_wt) / (max_wt - min_wt)
# cluster_probs <- cluster_probs + 0.05
# 
# # 4. Assign inclusion probability to each row
# occ_data_cl$incl_prob <- cluster_probs[as.character(occ_data_cl$Cluster)]
# 
# # 5. Sample presences and absences separately
# set.seed(42)
# 
# # Define desired sample sizes
# N_pres <- 100*5#sum(occ_data_cl$kelp == 1) # 1201
# N_abs <- N_pres #sum(occ_data_cl$kelp == 0)  # 3724
# 
# # Or choose fixed sizes
# # Sample presence rows
# sel_pres <- occ_data_cl[occ_data_cl$kelp == 1, ]
# sel_pres <- sel_pres[sample(
#   nrow(sel_pres), size = N_pres, prob = sel_pres$incl_prob, replace = F
# ), ]
# 
# # Sample absence rows
# sel_abs <- occ_data_cl[occ_data_cl$kelp == 0, ]
# sel_abs <- sel_abs[sample(
#   nrow(sel_abs), size = N_abs, prob = sel_abs$incl_prob, replace = F
# ), ]
# 
# # Combine into balanced training set
# training_data <- rbind(sel_pres, sel_abs)
# 
# training_data%>%group_by(kelp, Cluster)%>%
#   summarize(length(kelp))
# # kelp `length(kelp)`
# # 1     0            215
# # 2     1            215
# 
# # With alpha= 3
# #     kelp Cluster `length(kelp)`
# # 1     0       1             55
# # 2     0       2              7
# # 3     0       3              8
# # 4     0       5            145
# # 5     1       1            179
# # 6     1       2             33
# # 7     1       3              3
# 
# head(training_data)
# 
# # With alpha= 1
# #    kelp Cluster `length(kelp)`
# # 1     0       1            125
# # 2     0       2             36
# # 3     0       3             37
# # 4     0       5             17
# # 5     1       1            134
# # 6     1       2             75
# # 7     1       3              6
# 
# 
# training_data %>%
#   group_by(Cluster) %>%
#   summarize(mean_prob = mean(incl_prob, na.rm = TRUE))
# 
# 
# # Cluster 4 gets low probability and then there are no records selected for the training,
# # so we force to have a minimum of records per cluster
# 
# min_per_cluster <- 45  # minimum records per cluster for each class (presence/absence)
# target_total  <- 200 # total desired training records per class
# 
# # Prepare data subsets
# pres <- occ_data_cl %>% filter(kelp == 1)
# absn <- occ_data_cl %>% filter(kelp == 0)
# 
# sample_cluster <- function(subset_data=pres) {
#   subset_data %>%
#     group_by(Cluster) %>%
#     mutate(incl_prob = incl_prob / sum(incl_prob)) %>%  # normalize probabilities
#     summarise(
#       # Basic sampling
#       base_sample = min_per_cluster,
#       # Additional proportional allocation
#       extra_sample = round(incl_prob[1] * (target_total - n_distinct(Cluster) * min_per_cluster)),
#       total_sample = base_sample + extra_sample
#     ) -> cluster_counts
#   
#   # Now sample accordingly
#   sampled <- subset_data %>%
#     group_by(Cluster) %>%
#     do(sample_n(., size = cluster_counts$total_sample[cluster_counts$Cluster == unique(.$Cluster)], replace = TRUE))
#   
#   return(sampled)
# }
# 
# sel_pres <- sample_cluster(pres)
# sel_abs  <- sample_cluster(absn)
# training_data <- bind_rows(sel_pres, sel_abs)




## %######################################################%##
#                                                          #
####          Modeling with weighted occurrences        ####
#           by environmental cluster                       #
## %######################################################%##
occ_data_cl<- read.csv("presence_data_2014_to_2019_filtered_depth&rocky.csv")

press<- occ_data_cl%>%
  filter(kelp==1)

abs<- occ_data_cl%>%
  filter(kelp==0)%>%
  sample(size = length(press$kelp))
  
# frequency per cluster
freq <- table(occ_data_cl$Cluster)
occ_data_cl$raw_weight <- 1 / as.numeric(freq[as.character(occ_data_cl$Cluster)])  # inverse freq

# Normalize weights: make mean weight = 1 (or sum = nrow)
occ_data_cl$weight <- occ_data_cl$raw_weight * (length(occ_data_cl$raw_weight) / sum(occ_data_cl$raw_weight))



plot(occ_data_cl$weight, occ_data_cl$incl_prob)

min_per_cluster <- 43  # target samples per cluster per class = lower n in cluster 4

train_data_equal <- occ_data_cl %>%
  group_by(kelp, Cluster) %>%
  sample_n(min_per_cluster, replace = TRUE) %>%
  ungroup()

train_data_equal %>%
  group_by(kelp, Cluster) %>%
  summarize(length(kelp))


plot(vect(train_data_equal, geom= c("x","y")))

# write.csv(train_data_equal, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_filteredCluster.csv")
































###=============================================================================
### Create dataset to run models ===============================================
rdom_occs_df<- data.frame(matrix(ncol=5, nrow = 0))
colnames(rdom_occs_df)<- c("cluster6","occ", "x","y", "rdomreplicate")
  
set.seed(1)
for (a in 1:50) {
  for (i in 2:6) {

  cluster_df<- filt_0.1km[which(filt_0.1km$cluster6 == i),]
  n<- c(n_occs_perCluster[which(n_occs_perCluster$cluster6 == i), "ideal_n"])
  rdom_occs_cl<- sample_n( cluster_df, size= as.numeric(n), replace=F)
  rdom_occs_cl$rdomreplicate<- a
  # rdom_occs_cl<- rdom_occs_cl[,c("cluster6","occ", "lon","lat", "rdomreplicate")]
  
  rdom_occs_df<- rbind(rdom_occs_df, rdom_occs_cl)
  }
  
}

summary(as.factor(rdom_occs_df$cluster6))/50
head(rdom_occs_df)
unique(rdom_occs_df$rdomreplicate)

# write.csv(rdom_occs_df, paste(SDM_path, "csvs_to_SDMs/datset_occ_cl_replicates/50_rdom_newoccs_df_28May2024.csv", sep="/"))




################################################################################
### Make same process for absence records
filtabs_0.1km <- occfilt_geo(
  data = abs,
  x = "x",
  y = "y",
  env_layer = dem,
  method = c("defined", d = "0.1"),
  prj = crs(dem)
)

# Extracting values from raster ... 
# Number of unfiltered records: 1999
# Distance threshold(km): 0.1
# Number of filtered records: 942

# Extracting values from raster ... 28 May 2024
# Number of unfiltered records: 4484
# Distance threshold(km): 0.1
# Number of filtered records: 3356


# pLOT MAP WITH ABSENCES
filtabs_0.1km_sf <- filtabs_0.1km %>%
  st_as_sf(coords = c("x", "y"), crs = 3005)
mapview::mapview(filtabs_0.1km_sf, cex=2, zcol= "origin", col.regions= c("red", "blue", "green", "orange"))# a visual check

Focal_sites_sf<- Focal_sites %>%
  st_as_sf(coords = c("x", "y"), crs = 3005)
mapview::mapview(Focal_sites_sf, cex=2, zcol= "origin", col.regions= c("red", "blue", "green", "orange"))# a visual check
Focal_sites

filtabs_0.1km %>%
  group_by(cluster6, prop_area) %>% 
  summarise(n_occ = n())

# cluster6 prop_area n_occ
# 1        2     14.3     91
# 2        3     14.8    105
# 3        4     35.7    405
# 4        5     23.2    269
# 5        6      8.02    72

# cluster6 prop_area n_occ --- 28 May 2024
# 1        2     14.3    435
# 2        3     14.8    539
# 3        4     35.7   1630
# 4        5     23.2    337
# 5        6      8.02   415
n_occs_perCluster

plot_abs_all<- as.data.frame(filtabs_0.1km)%>%
  group_by(cluster6, prop_area) %>% 
  summarise(n_occ = n())%>%
  mutate(index= n_occ/prop_area)%>%
  ggplot(aes(x= prop_area, y= n_occ, color=as.factor(cluster6)))+#, size=index
  geom_point(size=3)+
  ylim(c(0, 700))+
  # geom_abline(intercept = 0, slope = 1, size = 0.5, color="grey") +
  # ylim(c(0, 45))+xlim(c(0, 45))+
  labs(x="Proportion of area (%)", y="Number of absence records", color="Cluster")+
  theme_bw()#+ theme(legend.position="NULL")

# Determine the number of absences from each cluster 
as.data.frame(filtabs_0.1km)%>%
  group_by(cluster6, prop_area) %>% 
  summarise(n_occ = length(occ))%>%
  group_by(cluster6)%>%
  mutate(index= n_occ/prop_area)

# cluster6 prop_area n_occ index
# 1        2     14.3    435  30.5
# 2        3     14.8    539  36.4
# 3        4     35.7   1630  45.6
# 4        5     23.2    337  14.5
# 5        6      8.02   415  51.7

# I select the cluster with the lowest index (Cluster 4; with lowest number of absences in relation to its surface)
# then I use divide this index per 3 and calculate N for each cluster
selected_index_abs<- round((337 / 14.5)/2, 0) # 3 index from cluster 2

filtabs_0.1km$occ<- 1
filtabs_0.1km<- as.data.frame(filtabs_0.1km)

n_abs_perCluster<- filtabs_0.1km%>%
  group_by(cluster6, prop_area) %>% 
  summarise(n_occ = length(occ))%>%
  group_by(cluster6)%>%
  mutate(index= n_occ/prop_area, ideal_n= round((selected_index_abs*prop_area),0))

#   cluster6 prop_area n_occ index ideal_n
# 1        2     14.3    435  30.5     171
# 2        3     14.8    539  36.4     178
# 3        4     35.7   1630  45.6     429
# 4        5     23.2    337  14.5     279
# 5        6      8.02   415  51.7      96



# Filter one presence per pixel 
df_raster<- as.data.frame(demp, xy = TRUE)
df_raster$pixel_ID<- seq(1, length(df_raster[,1]), 1)
head(df_raster)
demp_a<- rasterFromXYZ(df_raster) #[,c(1:2,4)]
crs(demp_a)<- crs(demp)
demp_a<- rast(demp_a)


filtabs_0.1km<- as.data.frame(filtabs_0.1km)

abs2<- sdm_extract(
  data = filtabs_0.1km,
  x = "x",
  y = "y",
  env_layer = dem_a, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = TRUE
)

head(abs2)
which(duplicated(abs2$pixel_ID))# there are no duplicates!


# Exclude duplicated records in same cell by selecting presences over absences
abs2$occ<- 0
pres_abs<- rbind(press2, abs2)
pres_abs<- as.data.frame(pres_abs)
duplicates_exclude<- which(duplicated(pres_abs$pixel_ID) & pres_abs$occ== 0)
pres_abs[which(duplicated(pres_abs$pixel_ID)),]
pres_abs[which(pres_abs$pixel_ID == "10494322"), c("x", "y")]

length(pres_abs[duplicates_exclude,"pixel_ID"] ) # 7 absences were removed 

pres_abs<- pres_abs[-duplicates_exclude,] 

filtabs_0.1km<- pres_abs[which(pres_abs$occ==0),]

plot_abs_filtered<- filtabs_0.1km%>%
  group_by(cluster6, prop_area) %>% 
  summarise(n_occ = n())%>%
  group_by(cluster6)%>%
  mutate(index= n_occ/prop_area, ideal_n=  round((selected_index_abs*prop_area),0))%>%
  # ggplot(aes(x= prop_area, y= n_occ, color=as.factor(cluster6)))+
  ggplot(aes(x= prop_area, y= ideal_n, color=as.factor(cluster6)))+ #, size=index
  geom_point(size=3)+
  ylim(c(0, 700))+
  # geom_abline(intercept = 0, slope = 1, size = 0.5, color="grey") +
  labs(x="Proportion of area (%)", y="Number of absence records", color="Cluster")+
  theme_bw()#+ theme(legend.position="NULL")


cowplot::plot_grid(all_pres_plot, filtered_pres_plot, plot_abs_all, plot_abs_filtered)
# ggsave("")



###=============================================================================
### Create dataset to run models ===============================================
rdom_abs_df<- data.frame(matrix(ncol=5, nrow = 0))
colnames(rdom_abs_df)<- c("cluster6","occ", "x","y", "rdomreplicate")

for (a in 1:50) {
  for (i in 2:6) {
    
    cluster_df<- filtabs_0.1km[which(filtabs_0.1km$cluster6 == i),]
    # n<- c(n_abs_perCluster[which(n_abs_perCluster$cluster6 == i), "ideal_n"])
    n<- c(n_abs_perCluster[which(n_occs_perCluster$cluster6 == i), "ideal_n"]) ## same number of absences than presences
    rdom_occs_cl<- sample_n(cluster_df, size= as.numeric(n), replace=F)
    rdom_occs_cl$rdomreplicate<- a
    # rdom_occs_cl<- rdom_occs_cl[,c("cluster6","occ", "lon","lat", "rdomreplicate")]
    
    rdom_abs_df<- rbind(rdom_abs_df, rdom_occs_cl)
  }
}


summary(as.factor(rdom_abs_df$cluster6))/50
n_abs_perCluster

colnames(rdom_abs_df)
colnames(rdom_occs_df)
unique(rdom_abs_df$rdomreplicate)



### MERGE PRESeNces AND ABSENCES datasets ========================================
colnames(rdom_occs_df)== colnames(rdom_abs_df[c(1:66,69)])
rdom_pres_abs_df<- rbind(rdom_occs_df, rdom_abs_df[c(1:66, 69)])
colnames(rdom_pres_abs_df)
head(rdom_pres_abs_df)


# Add focal sites dataset for all replicates 
Focal_sites
# Focal_sites<- occ_data_cl%>%
#   filter(origin =="focal_sites")

pres_focal<- Focal_sites[which(Focal_sites$occ==1),]
abs_focal<- Focal_sites[which(Focal_sites$occ==0),]

filt_focal_0.04km <- occfilt_geo(
  data = pres_focal,
  x = "x",
  y = "y",
  env_layer = dem,
  method = c("defined", d = "0.04"),
  prj = crs(dem)
)
# Extracting values from raster ... 
# Number of unfiltered records: 100
# Distance threshold(km): 0.1
# Number of filtered records: 35

# Extracting values from raster ... 
# Number of unfiltered records: 100
# Distance threshold(km): 0.02
# Number of filtered records: 52

# Extracting values from raster ... 
# Number of unfiltered records: 100
# Distance threshold(km): 0.04
# Number of filtered records: 41


# Filter absences
filt_absfocal_0.04km <- occfilt_geo(
  data = abs_focal,
  x = "x",
  y = "y",
  env_layer = dem,
  method = c("defined", d = "0.04"),
  prj = crs(dem)
)
# Extracting values from raster ... 
# Number of unfiltered records: 127
# Distance threshold(km): 0.04
# Number of filtered records: 60




# Exclude duplicated records in same cell by selecting presences over absences

filt_focal_0.04km<- sdm_extract(
  data = filt_focal_0.04km,  x = "x",  y = "y",
  env_layer = dem_a, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = TRUE
)

filt_absfocal_0.04km<- sdm_extract(
  data = filt_absfocal_0.04km,  x = "x",  y = "y",
  env_layer = dem_a, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = TRUE
)

filt_absfocal_0.04km$occ<- 0
pres_abs_focal<- rbind(filt_focal_0.04km, filt_absfocal_0.04km)
pres_abs_focal<- as.data.frame(pres_abs_focal)
duplicates_exclude<- which(duplicated(pres_abs_focal$pixel_ID) & pres_abs_focal$occ== 0)
pres_abs_focal[which(duplicated(pres_abs_focal$pixel_ID)),]
pres_abs_focal[which(pres_abs_focal$pixel_ID == "11766514"), c("x", "y")]

length(pres_abs_focal[duplicates_exclude,"pixel_ID"] ) # 5 absences were removed 

pres_abs_focal<- pres_abs_focal[-duplicates_exclude,] 

filtabs_focal_0.04km<- pres_abs_focal[which(pres_abs_focal$occ==0),]


filt_presabs_focal_0.04km<- rbind(filt_focal_0.04km, filtabs_focal_0.04km)
filt_presabs_focal_0.04km<- as.data.frame(filt_presabs_focal_0.04km)


# pLOT MAP WITH ABSENCES
filt_focal_sf <- filt_presabs_focal_0.04km %>%
  st_as_sf(coords = c("x", "y"), crs = 3005)
mapview::mapview(filt_focal_sf, cex=3, zcol= "occ", col.regions= c("red", "blue"))# a visual check



### Merge all datasets to create the final table for SDMs ======================
rdom_pres_abs_df
data_focal<- filt_presabs_focal_0.04km
data_focal$rdomreplicate<- 1
i=2
for (i in 2:50) {
  df_i<- filt_presabs_focal_0.04km
  df_i$rdomreplicate<- i
  data_focal<- rbind(data_focal, df_i )
}


colnames(data_focal[,c(1:66, 69)])
colnames(rdom_pres_abs_df)



FINAL_OCC_DATASET<- rbind(rdom_pres_abs_df, data_focal[,c(1:66, 69)])

# write.csv(FINAL_OCC_DATASET, paste(SDM_path,
#                                   "csvs_to_SDMs/datset_occ_cl_replicates/50_rdom_occs_abs_28May2024_FINAL.csv",
#                                   sep="/"), row.names = F)

n_abs_perCluster$occ<- "absences"
n_occs_perCluster$occ<- "presences"
n_abs_perCluster$prop_Used<- round(n_abs_perCluster$ideal_n/ n_abs_perCluster$n_occ *100, 2)
n_occs_perCluster$prop_Used<- round(n_occs_perCluster$ideal_n/ n_occs_perCluster$n_occ *100, 2)

table<- rbind(n_abs_perCluster,n_occs_perCluster)
# write.csv(table, "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/SDMs/records_per_cluster_filtering.csv")


summary_occs<- FINAL_OCC_DATASET%>%
  group_by( rdomreplicate, occ, cluster6)%>%
  summarize(n())

as.data.frame(summary_occs)


###=============================================================================
### Plot diverse replicate datasets  ===========================================
rdom_pres_abs_df$cluster6
rdom_pres_abs_df$cluster6<- as.factor(rdom_pres_abs_df$cluster6)
rep1<- FINAL_OCC_DATASET[which(FINAL_OCC_DATASET
                          $rdomreplicate ==1),]
rep2<- FINAL_OCC_DATASET[which(FINAL_OCC_DATASET$rdomreplicate ==2),]
rep6<- FINAL_OCC_DATASET[which(FINAL_OCC_DATASET$rdomreplicate ==1),]

summary((rep1$cluster6))
summary((rep6$cluster6))


par(mfrow=c(2,2))
newColors <- c("springgreen", "blue", "magenta", "orange", "darkgrey", "black")


plot(stack_vars$bathymetry, main= "Presences, rep1", col="grey", legend = FALSE)
points(rep1[which(rep1$occ ==1),c("x", "y")], 
       col= (newColors)[rep1$cluster6], cex=0.75, pch=16)

legend("bottomright",   # location of legend
       legend = levels(rep1$cluster6), # categories or elements to render in
       # the legend
       fill = newColors) # color palette to use to fill objects in legend.

plot(stack_vars$bathymetry, main= "Presences, rep2", col="grey", legend = FALSE)
points(rep2[which(rep2$occ ==1),c("x", "y")], 
       col= (newColors)[rep2$cluster6], cex=0.75, pch=16)


plot(stack_vars$bathymetry, main= "Absences, rep2", col="grey", legend = FALSE)
points(rep2[which(rep2$occ ==0),c("x", "y")], 
       col= (newColors)[rep2$cluster6], cex=0.75, pch=16)

# add a legend to your map
legend("bottomright",   # location of legend
       legend = levels(rep2$cluster6), # categories or elements to render in
       # the legend
       fill = newColors) # color palette to use to fill objects in legend.


i=1
par(mfrow=c(4,4))
for (i in 1:8) {
  plot(stack_vars$bathymetry, main= paste("Presences", i), col= "grey")
  points(FINAL_OCC_DATASET[which(FINAL_OCC_DATASET$rdomreplicate ==i & FINAL_OCC_DATASET$occ == 1),c("x", "y")],
         cex=0.75, pch=16)#col= FINAL_OCC_DATASET$cluster6, 
  
  plot(stack_vars$bathymetry, main= paste("Absences", i), col= "grey")
  points(FINAL_OCC_DATASET[which(FINAL_OCC_DATASET$rdomreplicate ==i & FINAL_OCC_DATASET$occ == 0),c("x", "y")],
         cex=0.75, pch=16) #col= FINAL_OCC_DATASET$cluster6,
  
}

legend("bottomright",   # location of legend
       legend = levels(FINAL_OCC_DATASET$cluster6), # categories or elements to render in
       # the legend
       fill = newColors) # color palette to use to fill objects in legend.




# par(mfrow=c(3,4))
# for (i in 1:10) {
#   plot(clusters, main= paste("Absences", i), col= "grey")
#   points(rdom_pres_abs_df[which(rdom_pres_abs_df$rdomreplicate ==i & rdom_pres_abs_df$occ == 0),c("x", "y")], 
#          col= rdom_pres_abs_df$cluster6, cex=0.75, pch=16)
# }







