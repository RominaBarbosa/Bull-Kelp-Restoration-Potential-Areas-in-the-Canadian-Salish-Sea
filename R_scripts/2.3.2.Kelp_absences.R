##==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
### 2- Absence/speudo-absence data Preparation                  ################
### Author: Romina Barbosa                                      ################
### Date last edition: 25-July-2025                             ################
###==================================================================
# Load libraries
library(terra)
library(sf)

layers_path<- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution"
SDM_path<- "/Volumes/Romina_PSF/PSF/SDM"

plotspath<- "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Plots"
mypath<- "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp"


### Load env layers =============================================================
# source("/Volumes/Romina_PSF/PSF/R_scripts/stack_rasters_path_function.R")
# stack_vars<- stack_rasters_path_funtion(layers_path, terra_class= "Y")
stack_vars<- terra::rast(paste(layers_path,"SalishSeaCast_interp_20m_resolution_FINAL2.tif", sep="/"))
unique(names(stack_vars))

bathymetry_20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathymetry_20m<- crop(bathymetry_20m, stack_vars)


### Determine the area to create absence records ============================
# 1. Mask area by depth < 40 m (shallower than -40)
area<- bathymetry_20m
area[area < -10 | area > 40]<- NA
plot(area)

# Create a raster of cell IDs
grid_id_raster <- setValues(area, 1:ncell(area))
names(grid_id_raster) <- "grid_ID"


# Convert raster to data frame
area.df <- as.data.frame(area, xy = TRUE, cells = TRUE, na.rm = T)
# area.df<- na.exclude(area.df)

# area.df<- flexsdm::sdm_extract(
#   data = area.df,
#   x = "x",
#   y = "y",
#   env_layer = bathymetry_20m, # Raster with environmental variables
#   variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
#   filter_na = F  #NA's   :661, but they should be 315 based on the corrected bathymetry points 
# )

area_sf <- st_as_sf(area.df, coords = c("x", "y"), crs = 3005)
mapview::mapview(area_sf, cex=3, color= "blue")

# Extract grid_ID values at point locations
# Convert sf to SpatVector
area_vect <- vect(area_sf)
grid_ids <- extract(grid_id_raster, area_vect)
area_df <- cbind(area_sf, grid_ids["grid_ID"])

colnames(area.df)<- c("cell_id", "x", "y", "bathymetry_500", "bathymetry_20")

### Load kelp presence points and determine absence records ====================
# 1. Load kelp presence points (assume they have lat/lon)
kelp_pts <- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered.csv")
duplidated(kelp_pts)
length(kelp_pts$cell_id)
length(unique(kelp_pts$cell_id))
(summary(as.factor(kelp_pts$cell_id)))

kelp_sf <- st_as_sf(kelp_pts, coords = c("x", "y"), crs = 3005)

mapview::mapview(area_sf, cex=3, color= "blue")+
 mapview::mapview(kelp_sf[5], cex=3)

# 2. Merge kelp presence points and grids in the bathy raster to exclude kelp sites (assume they have lat/lon)
kelp_pts$kelp<- "1"


length(area.df$cell_id)
length(kelp_pts$cell_id)

kelp_presabs_df<- merge(kelp_pts, area.df, by = c("cell_id"), all=T)
length(kelp_presabs_df$cell_id)
length(unique(kelp_presabs_df$cell_id))


# 3. Make all grids without kelp records of Kelp ABSENCE
kelp_presabs_df[which(is.na(kelp_presabs_df$kelp)), "kelp"]<- "0"

# Remove kelp absences at grids with kelp presence
# df with cell_id and kelp (0 = absence, 1 = presence)
kelp_presabs_df <- kelp_presabs_df %>%
  group_by(cell_id) %>%
  filter(kelp == 1 | all(kelp == 0)) %>% #removes absences from cells that also contain a presence
  ungroup()


length(kelp_presabs_df_2$cell_id)
length(unique(kelp_presabs_df_2$cell_id))


# 4. Save presence and absence dataset - Blob period - 
colnames(kelp_presabs_df)

kelp_presabs_df <- kelp_presabs_df %>%
  mutate(
    x = coalesce(x.x, x.y),
    y = coalesce(y.x, y.y)
  ) %>%
  select(-x.x, -x.y, -y.x, -y.y)  # drop old columns



kelp_presabs_df<- flexsdm::sdm_extract(
  data = kelp_presabs_df[,c(1,3,5:8,28,30:32, 34:36)],
  x = "x",
  y = "y",
  env_layer = stack_vars, # Raster with environmental variables
  variables = NULL, # Vector with the variable names of predictor variables Usage variables. = c("aet", "cwd", "tmin"). If no variable is specified, function will return data for all layers.
  filter_na = F  #NA's   :661, but they should be 315 based on the corrected bathymetry points 
)


# Check and exclude grids without information of env conditions
summary(kelp_presabs_df)
kelp_presabs_df<- kelp_presabs_df%>%
  dplyr::filter(!is.na(northerness_500m_clipped))%>%
  dplyr::filter(!is.na(TPI_500m_clipped))

summary(as.factor(kelp_presabs_df$kelp))
# kelp_presabs_df[which(is.na(kelp_presabs_df$kelp)), "kelp"]<- "0"

# 0      1 
# 14530  1287 


# write.csv(kelp_presabs_df, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_2014_to_2019_filtered_plus_Absences.csv")

head(kelp_presabs_df)




# # 6. Generate 10,000 random points in valid area
# set.seed(123)  # for reproducibility
# random_points <- spatSample(valid_area, size = 10000, method = "random", as.points = TRUE)