###==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
### 2- Presence, speudo-absence and background data Preparation ################
### Author: Romina Barbosa                                      ################
### Date last edition: 28-May-2024                                           ################
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
occ_data_sdm<- read.csv( paste(SDM_path, "Presence_absences_kelp/merged_kelp_pts_2014_2019.csv", sep="/"))
columns<- c("latitude", "longitude", "Class_name", "Region", "Year", "Area_m2", "layer")
occ_data_sdm<- occ_data_sdm[,columns]
head(occ_data_sdm)
summary(occ_data_sdm)


### Load env layers =============================================================
source("/Volumes/Romina_PSF/PSF/R_scripts/stack_rasters_path_function.R")
stack_vars<- stack_rasters_path_funtion(layers_path, terra_class= "Y")
unique(names(stack_vars))
names(stack_vars)[4]<- "bathymetry"

bathymetry_20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathymetry_20m<- crop(bathymetry_20m, stack_vars)

### Extract env. conditions at the kelp presence locations  ====================
occ_data_vars<- sdm_extract(
  data = occ_data_sdm,
  x = "longitude",
  y = "latitude",
  env_layer = stack_vars, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = F  #NA's   :661, but they should be 315 based on the corrected bathymetry points 
)

occ_data_vars<- sdm_extract(
  data = occ_data_vars,
  x = "longitude",
  y = "latitude",
  env_layer = bathymetry_20m, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = F  #NA's   :661, but they should be 315 based on the corrected bathymetry points 
)

colnames(occ_data_vars)[ncol(occ_data_vars)]<- "bathymetry_20m"
length(occ_data_vars$Year)# [1] 17724
occ_data_vars<- as.data.frame(occ_data_vars)


### Determine the depth of the 99% of kelp presence records ====================
Q99_kelp_depth<- quantile(occ_data_vars$bathymetry_20m, probs =0.99, na.rm = TRUE)
# 95% =  9.97939825058 
# 99% =  21.51443481445 
# 1% = -6.653466176987 
# 5% = -3.00764799118 

quantile(occ_data_vars$bathymetry, probs =0.99, na.rm = TRUE)
# 95% =  37.10013580322 
# 99% =  54.61620731354 

ggplot(occ_data_vars, aes(x= bathymetry_20m, y= bathymetry))+
  geom_point()
  
ggplot(occ_data_vars, aes(x= as.factor(Year), y= bathymetry))+
  geom_boxplot()+ theme_bw()

ggplot(occ_data_vars, aes( y= bathymetry_20m))+
  geom_boxplot()+ theme_bw()

occ_data_vars.b<- occ_data_vars%>%
  filter(bathymetry_20m <= Q99_kelp_depth[[1]])

17724 - length(occ_data_vars.b$Year)# 17470 (from 17724)--> 254 records were excluded 


### Extracted cluster values by proximity  =====================================
library(sf)
library(nngeo) 

clusters<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/clusters_MoraSoto2024_3005.csv")

kelp <- st_as_sf(occ_data_vars.b, coords = c("longitude", "latitude"), crs = 3005)
clusters <- st_as_sf(clusters[,c(31:33)], coords = c("lon_3005", "lat_3005"), crs = 3005)  # or your spatial data frame

# Find nearest cluster for each kelp point
nearest_ix <- st_nn(kelp, clusters, k = 1, returnDist = F ) # maxdist = 1000   ---Ensure both are in same CRS (important!)

# Attach cluster data to kelp
kelp$env_cluster <- clusters$Cluster[unlist(nearest_ix)]
head(kelp)



## Plot kelp records by Cluster and check them =================================
mapview::mapview(kelp, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check

as.data.frame(kelp)%>%
  group_by(env_cluster)%>%
  summarize(n_occ= length(Year))
# env_cluster n_occ
# 1           1  8432
# 2           2  6670
# 3           3  1765
# 4           4   572
# 5          NA    31

kelp_na_cluster<- kelp[ which(is.na(kelp$env_cluster)),]
mapview::mapview(kelp_na_cluster, cex=2, zcol= "temperature_summer_mean_blob")#, col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check

### Exclude wrong/suspicious records ===========================================
## Exclude kelp presences in enclosed area near Sooke because there is information about the ENv. CLuster 
kelp<- kelp[-which(is.na(kelp$env_cluster) & kelp$temperature_summer_mean_blob >=12.5),]

## The other records are in Cluster 1 
kelp[ which(is.na(kelp$env_cluster)), "env_cluster"]<- 1

## Remove kelp records without Year information 
kelp_na_year<- kelp[ which(is.na(kelp$Year)),]
mapview::mapview(kelp_na_year, cex=2, add=T)#, col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check
length(kelp_na_year$Class_name) # 1773

mapview::mapview(kelp, cex=2, zcol= "Year")

kelp<- kelp[-which(is.na(kelp$Year)),]


as.data.frame(kelp)%>%
     group_by(env_cluster)%>%
     summarize(n_occ= length(Year))
# env_cluster n_occ
# 1           1  7680
# 2           2  6067
# 3           3  1581
# 4           4   367

mapview::mapview(kelp, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check



## %######################################################%##
#                                                          #
####          Filtering occurrences                     ####
#                                                          #
## %######################################################%##
options(digits=13)


## Filtering kelp presences to keep only one per cell/grid (500 m resolution layers) ====================================
kelp_vect <- vect(kelp)# Convert sf to SpatVector (terra uses this)
coords <- crds(kelp_vect)# Convert SpatVector to coordinates matrix

# Create a raster of cell IDs
bathymetry<- stack_vars[[4]]
grid_id_raster <- setValues(bathymetry, 1:ncell(bathymetry))
names(grid_id_raster) <- "grid_ID"

# Extract grid_ID values at point locations
grid_ids <- extract(grid_id_raster, kelp_vect)

# Combine with point data
points_with_grid <- cbind(kelp_vect, grid_ids["grid_ID"])
plot(kelp_vect)
head(kelp_vect)

length(points_with_grid$grid_ID)
length(unique(points_with_grid$grid_ID))

## Keep only one point per cell
kelp_unique <- points_with_grid[!duplicated(points_with_grid$cell_id), ]
length(kelp_unique$grid_ID) #[1] 1287
length(unique(kelp_unique$grid_ID))# [1] 1287

# set.seed(123)  # for reproducibility
# sampled_points <- points_with_grid %>%
#   group_by(grid_ID) %>%
#   slice_sample(n = 1) %>%   # Select 1 row at random per group
#   ungroup()

# Extract coordinates
coords <- crds(kelp_unique)

# Combine with attribute data
kelp_unique <- cbind(as.data.frame(kelp_unique), coords)


# Initial dataset 
kelp_unique.st<- st_as_sf(kelp_unique, coords = c("x", "y"), crs = 3005)
mapview::mapview(kelp_unique.st, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check



# Calculate the proportion of area occupied by each environmental cluster ====
cl_rel_area<- as.data.frame(clusters) %>%
  group_by(Cluster) %>% 
  summarise(n_segments = length(Cluster)) # same than using n()

# Cluster n_segments
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
# Cluster n_segments prop_area
# 1       1        430      9.98
# 2       2        864     20.1 
# 3       3       1137     26.4 
# 4       4       1753     40.7 
# 5       5        124      2.88

mapview::mapview(clusters, cex=2, zcol= "Cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check

colnames(cl_rel_area)[1]<- "env_cluster"

# Merge kelp data with information of cluster area 
occ_data_cl<- merge(kelp_unique, cl_rel_area, by="env_cluster") 
colnames(occ_data_cl)


### Summarize number of kelp record per Environmental Cluster ==================
occ_data_cl %>%
  group_by(env_cluster, prop_area) %>% 
  summarise(n_occ = n())
# env_cluster prop_area n_occ
# 1           1      9.98   409
# 2           2     20.1    668
# 3           3     26.4    159
# 4           4     40.7     51


### Compare entire kelp records dataset and the filtered one in the map ========
# FInal dataset
occ_data_cl_sf <- occ_data_cl %>%
  st_as_sf(coords = c("x", "y"), crs = 3005)
mapview::mapview(occ_data_cl_sf, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check

# Initial dataset 
mapview::mapview(kelp, cex=2, zcol= "env_cluster", col.regions= c("blue", "green", "red",  "orange", "purple"))# a visual check



### Plot the relationship of number of kelp presences by Cluster
all_pres_plot<- occ_data_cl_sf%>%
  group_by(env_cluster, prop_area) %>% 
  summarise(n_occ = n())%>%
  ggplot(aes(x= prop_area, y= n_occ, color=as.factor(env_cluster)))+
  geom_point(size=3)+
  # ylim(c(0, 700))+
  # geom_abline(intercept = 0, slope = 1, size = 0.5, color="grey") +
  labs(x="Proportion of area (%)", y="Number of presence records", color="Cluster")+
  theme_bw()#+ theme(legend.position="NULL")



write.csv(occ_data_cl, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered.csv")
























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







