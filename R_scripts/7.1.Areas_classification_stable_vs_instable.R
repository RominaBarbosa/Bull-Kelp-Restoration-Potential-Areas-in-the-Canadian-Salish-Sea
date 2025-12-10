###====================================================================================
###     Species Distribution models    SDMs                                      ######
###                                                                              ######
### 7.1- Assessing the Classified Areas for Restoration                          ######
### Content:                                                                     ######
###   Convert continuous predictions into binary and analyse temporal stability  ######
###   Combine stability and suitability to classify Restoration Categories       ######
### Author: Romina Barbosa                                                       ######
### Date last edition: 24-Nov-2025                                               ######
###====================================================================================

# Load packages
library(terra)
library(dplyr)
library(dismo)
library(ggplot2)
library(raster)
library(tidyr)


model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
setwd(model_results_path)
postblob_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob"
out_dir <- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7"

tab_ensemble_thresholds<- read.csv("Sep2025_M7_weightedPres/tab_ensemble_thresholds_M7_FINAL.csv")
tab_ensemble_thresholds<- tab_ensemble_thresholds[,-1]
# Model   Criterion Threshold Sensitivity Specificity       TSS       AUC
# 1 Ensemble_ave     Max TSS 0.4803167   0.9179523   0.9268092 0.8447615 0.9302505
# 2 Ensemble_ave No omission      -Inf   1.0000000   0.0000000 0.0000000 0.9302505
# 3 Ensemble_ave 1% omission 0.2414677   0.9894811   0.7360197 0.7255008 0.9302505

# FINAL
# Model   Criterion Threshold Sensitivity Specificity       TSS       AUC
# 1 Ensemble_ave     Max TSS 0.4733616   0.9256662   0.9226974 0.8483636 0.9793087
# 2 Ensemble_ave No omission      -Inf   1.0000000   0.0000000 0.0000000 0.9793087
# 3 Ensemble_ave 1% omission 0.2414677   0.9894811   0.7360197 0.7255008 0.9793087


# output dir for any files
# out_dir <- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7"
# dir.create(out_dir, showWarnings = FALSE)


# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
# bathy20m_mask10_40<- bathy20m
# bathy20m_mask10_40[!(bathy20m_mask10_40 >= -10 & bathy20m_mask10_40 <= 40)] <- NA # negative values are in land

bathy20m_mask10_15<- bathy20m
bathy20m_mask10_15[!(bathy20m_mask10_15 >= -10 & bathy20m_mask10_15 <= 15)] <- NA # negative values are in land
plot(bathy20m_mask10_15)

### Load substrate layer =====================================================
substrate_aligned<- rast(paste(model_results_path, "substrate_SOG_aligned_.tif", sep="/"))
# substrate_aligned[substrate_aligned == 2]<- NA
substrate_aligned_masked <- terra::mask(substrate_aligned, bathy20m_mask10_15)
plot(substrate_aligned_masked)
# writeRaster(substrate_aligned_masked, paste(model_results_path, "substrate_SOG_aligned_masked.tif", sep="/"), overwrite=T)


### Load model predictions =====================================================
ens_average<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble.tif", sep="/"))
ens_average_t2<- rast(paste(postblob_path, "M7/tifs/ensemble_postblob.tif", sep="/"))
# ens_average<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif", sep="/"))
# ens_average_t2<- rast(paste(postblob_path, "M7/ens_ave_postblob_M7.tif", sep="/"))
names(ens_average)<- "ens_average"
names(ens_average_t2)<- "ens_average_postblob"

# ------------------------------------------------------------------------------
#  Convert continuous predictions into categoric and analyse overlap
# ------------------------------------------------------------------------------
# thresholds
thresh1 = tab_ensemble_thresholds[which(tab_ensemble_thresholds$Criterion == "1% omission"),"Threshold"] #lower threshold (unsuitable vs moderate)
thresh2 = tab_ensemble_thresholds[which(tab_ensemble_thresholds$Criterion == "Max TSS"),"Threshold"]     #upper threshold (moderate vs optimal)

# -----------------------------------------------------
# 1) Classify continuous predictions into numeric categories (1,2,3)
# unsuitable: 1 
# adequate: 2 
# optimal: 3 
# -----------------------------------------------------
t1_cat <- classify(ens_average, rcl = matrix(c(-Inf, thresh1, 1,
                                      thresh1, thresh2, 2,
                                      thresh2, Inf, 3),
                                    ncol=3, byrow=TRUE))

t2_cat <- classify(ens_average_t2, rcl = matrix(c(-Inf, thresh1, 1,
                                      thresh1, thresh2, 2,
                                      thresh2, Inf, 3),
                                    ncol=3, byrow=TRUE))

# -----------------------------------------------------
# Identify any pixel unsuitable in ANY period
unsuitable_mask <- (t1_cat == 1) | (t2_cat == 1)

# Create an empty restoration class raster
restoration_class <- rast(t1_cat)
restoration_class[] <- NA

# 2) Map categories to numeric scores
# unsuitable: 1 -> -2;
# adequate: 2 -> 2;
# optimal: 3 -> 4
# -----------------------------------------------------
# Categorical habitat suitability maps from both periods were overlaid to identify 
# areas experiencing changes in habitat conditions, to assess where restoration 
# interventions are more likely to succeed. To do that, each habitat suitability 
# category was assigned a numeric score of (unsuitable: -2, marginal: 2, and optimal: 4,) 
# for both periods, and the scores were summed. The resulting index was a 
# combination of habitat quality and stability over time (-4 − 8), which was 
# reclassified to define four types of restoration potential areas: 
# 1) no restoration potential (score: -4 − 2), 2) low restoration potential (score: 4), 
# 3) medium restoration potential (score: 6), and 4) highly restoration potential (score: 8) (Table 1). 

t1_score <- classify(t1_cat, rcl = matrix(c(1, -2,
                                            2, 2,
                                            3, 4), ncol=2, byrow=TRUE))

t2_score <- classify(t2_cat, rcl = matrix(c(1, -2,
                                            2, 2,
                                            3, 4), ncol=2, byrow=TRUE))

# Save rasters
# writeRaster(t1_score, filename = file.path(out_dir, "period1_habitatCategories.tif"), overwrite = TRUE)
# writeRaster(t2_score, filename = file.path(out_dir, "period2_habitatCategories.tif"), overwrite = TRUE)


# --------------------------------------------------
# 5) Classify Restoration Potential classes
# --------------------------------------------------
# Ensure t1_cat and t2_cat are integer rasters with values 1,2,3

# Create output raster
restoration_class <- rast(t1_cat)
restoration_class[] <- NA

# 1) NON-VIABLE: unsuitable in ANY period
restoration_class[(t1_cat == 1) | (t2_cat == 1)] <- 1

# 2) LOW: adequate or marginal (2) in BOTH periods
restoration_class[(t1_cat == 2) & (t2_cat == 2)] <- 2

# 3) MEDIUM: adequate + optimal (2 + 3 in either order)
restoration_class[
  ((t1_cat == 2) & (t2_cat == 3)) |
    ((t1_cat == 3) & (t2_cat == 2))
] <- 3

# 4) HIGH: optimal in BOTH periods
restoration_class[(t1_cat == 3) & (t2_cat == 3)] <- 4


cols <- c("red","orange","yellow","green")  
plot(restoration_class, col = cols, axes = FALSE, legend = FALSE, main="Restoration Potential Map")
legend("topright", legend=c("Non-viable","Low","Medium","High"),
       fill=cols, bty="n")


# writeRaster(restoration_class,
#             filename = file.path(out_dir, "restoration_potential.tif"),
#             overwrite = TRUE)


# ---------------------------
# 3) Sum scores for numeric stability
# ---------------------------
total_score <- t1_score + t2_score  # possible values: -4,0,2,4,6,8
# NA pixels remain NA

# ---------------------------
# 4) Numeric classification into Restoration Suitability classes
# 1= Not suitable (Instable-Unsuitable or Adequate),
# 2= Moderately suitable (Stable-Adequate),
# 3= Suitable (Stable-Adequate to optimal),
# 4= Highly Suitable (Stable-optimal)
# ---------------------------
stability_numeric <- classify(total_score, rcl = matrix(c(-Inf, 2, 1,
                                                          2, 4, 2,
                                                          4, 6, 3,
                                                          6, Inf, 4), 
                                                        ncol=3, byrow=TRUE))

# 1) Save total_score
# writeRaster(total_score,
#             filename = file.path(out_dir, "total_score.tif"),
#             overwrite = TRUE)

# 2) Save stability_numeric
# writeRaster(stability_numeric,
#             filename = file.path(out_dir, "stability_layer.tif"),
#             overwrite = TRUE)



# writeRaster(stability_cat,
#             filename = file.path(out_dir, "stability_layer_categoric.tif"),
#             overwrite = TRUE)

# -------------------------------------------------
# 6) Optional: visualize missing cells
# -------------------------------------------------
na_map <- total_score
na_map[!is.na(na_map)] <- 0
na_map[is.na(na_map)] <- 1
plot(na_map, col=c("transparent","black"), main="Missing cells (NA)")




# Plot
cols <- c("blue","yellow","green", "red")  # match the order of categories

plot(stability_cat, col = cols, axes = FALSE, legend = FALSE, main="Stability (categoric)")
legend("topright", legend=cats$category, fill=cols, bty="n")





###-----------------------------------------------------------------------------
### Identify areas of change from unsuitable to suitable (adequate and optimal together)
###-----------------------------------------------------------------------------
# Inputs: t1 and t2 are continuous rasters (0–1)
# Threshold for suitability
# thresh1 <- threshold_drop

# 0 = unsuitable, 1 = suitable
t1_bin <- classify(ens_average, rcl = matrix(c(-Inf, thresh1, 0,
                                      thresh1, Inf, 1),
                                    ncol=3, byrow=TRUE))
t2_bin <- classify(ens_average_t2, rcl = matrix(c(-Inf, thresh1, 0,
                                      thresh1, Inf, 1),
                                    ncol=3, byrow=TRUE))
# Step 2: Detect changes in status

# Change code = t1_bin * 10 + t2_bin
# Codes:
# 00 = unsuitable → unsuitable
# 01 = unsuitable → suitable  (gain)
# 10 = suitable → unsuitable  (loss)
# 11 = suitable → suitable
change_code <- t1_bin * 10 + t2_bin


# Step 3: Assign numeric codes for change map (optional)
# 0 = no change (unsuitable → unsuitable)
# 1 = gain (unsuitable → suitable)
# 2 = loss (suitable → unsuitable)
# 3 = stable suitable (suitable → suitable)

change_map <- classify(change_code, rcl = matrix(c(
  0, 0, 0,   # 00 → no change
  1, 1, 1,   # 01 → gain
  10,10, 2,  # 10 → loss
  11,11, 3   # 11 → stable suitable
), ncol=3, byrow=TRUE))

# writeRaster(change_map, filename = file.path(out_dir, "change_suitability_M7.tif"), overwrite = TRUE)

# Step 4: Plot changes
cols <- c("gray", "yellow", "maroon2","purple")  # no change, gain, loss, stable suitable
labels <- c("Unsuitable→Unsuitable", "Unsuitable→Suitable", 
            "Suitable→Unsuitable", "Suitable→Suitable")

plot(change_map, col = cols, axes = FALSE, legend = FALSE, main="Suitability change")
legend("bottomright", legend = labels, fill = cols, bty="n")




## Mask layers bu depth and substrate and plot  =================================

t1_cat_maskDepth<- terra::mask(t1_cat, bathy20m_mask10_15)
t1_cat_maskDepthSubstrate<- terra::mask(t1_cat_maskDepth, substrate_aligned)
# writeRaster(t1_cat_maskDepth, paste(model_results_path, "Sep2025_M7_weightedPres/tifs/categoric_masked/ensemble_masked_depth.tif", sep="/"))
# writeRaster(t1_cat_maskDepthSubstrate, paste(model_results_path, "Sep2025_M7_weightedPres/tifs/categoric_masked/ensemble_masked_depthSubstrate.tif", sep="/"))

t2_cat_maskDepth<- terra::mask(t2_cat, bathy20m_mask10_15)
t2_cat_maskDepthSubstrate<- terra::mask(t2_cat_maskDepth, substrate_aligned)
plot(t1_cat_maskDepth)
# writeRaster(t2_cat_maskDepth, paste(postblob_path, "M7/tifs/categoric_masked/ensemble_t2_masked_depth.tif", sep="/"))
# writeRaster(t2_cat_maskDepthSubstrate, paste(postblob_path, "M7/tifs/categoric_masked/ensemble_t2_masked_depthSubstrate.tif", sep="/"))

change_map_maskDepth<- terra::mask(change_map, bathy20m_mask10_15)
change_map_maskDepthSubstrate<- terra::mask(change_map_maskDepth, substrate_aligned)
plot(change_map_maskDepthSubstrate)
# writeRaster(change_map_maskDepth, paste(out_dir, "change_map_masked_depth.tif", sep="/"))
# writeRaster(change_map_maskDepthSubstrate, paste(out_dir, "change_map_maskDepthSubstrate.tif", sep="/"))

restoration_class_maskDepth<- terra::mask(restoration_class, bathy20m_mask10_15)
restoration_class_maskDepthSubstrate<- terra::mask(restoration_class_maskDepth, substrate_aligned)
plot(restoration_class_maskDepthSubstrate)
# writeRaster(restoration_class_maskDepth, paste(out_dir, "restoratio_potential_masked_depth.tif", sep="/"), overwrite=T)
# writeRaster(restoration_class_maskDepthSubstrate, paste(out_dir, "restoratio_potential_maskDepthSubstrate.tif", sep="/"), overwrite=T)


