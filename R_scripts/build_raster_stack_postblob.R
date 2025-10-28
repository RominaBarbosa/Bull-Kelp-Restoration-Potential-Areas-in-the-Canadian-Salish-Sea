
path_postBlob<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M5"
setwd(path_postBlob)


# =============================================================================
# STEP 1: Load input variables of Post Blob Period of the entire area
# =============================================================================
# Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblobM2.tif")


### Load terrain variables =====================================================
terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)

# Load and stack all tif files
# terrain_vars<- terra::rast(paste(variables_selection_path,"terrain_rasters_selected.tif", sep="/"))
terrain_vars<-  rast(tif_files[c(13)])
names(terrain_vars)<- "slope_5x5"

### Merge rasters of all selected variables including terrain and NEMO =========
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars)

# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask15_40<- bathy20m
bathy20m_mask15_40[!(bathy20m_mask15_40 >= -15 & bathy20m_mask15_40 <= 40)] <- NA # negative values are in land

### Mask all layers by bathymetry  ==============================================
raster_stack_20m_all<- mask(raster_stack_20m_all, bathy20m_mask15_40)
names(raster_stack_20m_all)


### Select final variables  ==============================================
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres"
glm_mod_s<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))

# Get variables used in the model
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)

raster_stack_predict_postBlob<- raster_stack_20m_all[[names(raster_stack_20m_all)%in%vars_selected]]

# saveRDS(raster_stack_predict_postBlob, "raster_stack_predict_postBlob.rds")

