


model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres"
terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")


# This chunk was run and its outputs are stored in a RDS file
# # Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_FINAL2.tif")
# names(raster_stack_20m)
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)
terrain_vars<-  rast(tif_files[c(12)]) #,5,12,13,15
# names(terrain_vars)
# summary(terrain_vars)

# Merge rasters of all selected variables including terrain and NEMO
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)

# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask15_40<- bathy20m
bathy20m_mask15_40[!(bathy20m_mask15_40 >= -15 & bathy20m_mask15_40 <= 40)] <- NA # negative values are in land
# plot(bathy20m_mask15_40)

# Mask all layers
glm_mod_s<- readRDS("glm_mod_s.rds")
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)

raster_stack_predict<- raster_stack_20m_all[[vars_selected]] # subset selected variables
# names(raster_stack_predict)
# plot(raster_stack_predict[[1]])
raster_stack_predict<- mask(raster_stack_predict, bathy20m_mask15_40)
# plot(raster_stack_predict[[1]])


# saveRDS(raster_stack_predict, "raster_stack_predict.rds")