









model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
setwd(model_results_path)
tab_ensemble_thresholds<- read.csv("Sep2025_M7_weightedPres/tab_ensemble_thresholds_M7.csv")
tab_ensemble_thresholds<- tab_ensemble_thresholds[,-1]


# output dir for any files
out_dir <- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7"
dir.create(out_dir, showWarnings = FALSE)

### Load substrate layer =====================================================
substrate_aligned<- rast(paste(model_results_path, "substrate_SOG_aligned.tif", sep="/"))
substrate_aligned[substrate_aligned == 2]<- NA
substrate_aligned <- terra::mask(substrate_aligned, ens_average)

# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask10_30<- bathy20m
bathy20m_mask10_30[!(bathy20m_mask10_30 >= -10 & bathy20m_mask10_30 <= 15)] <- NA # negative values are in land
plot(bathy20m_mask10_30)

substrate_aligned_masked <- terra::mask(substrate_aligned, bathy20m_mask10_30)
plot(substrate_aligned_masked)
writeRaster(substrate_aligned_masked, paste(model_results_path, "substrate_SOG_aligned.tif", sep="/"), overwrite=T)


### Load model predictions =====================================================
ens_average<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif", sep="/"))
ens_average_t2<- rast(paste(postblob_path, "M7/ens_ave_postblob_M7.tif", sep="/"))
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
# 2) Map categories to numeric scores
# unsuitable: 1 -> -2;
# adequate: 2 -> 2;
# optimal: 3 -> 4
# -----------------------------------------------------
t1_score <- classify(t1_cat, rcl = matrix(c(1, -2,
                                            2, 2,
                                            3, 4), ncol=2, byrow=TRUE))

t2_score <- classify(t2_cat, rcl = matrix(c(1, -2,
                                            2, 2,
                                            3, 4), ncol=2, byrow=TRUE))

# Save rasters
writeRaster(t1_score, filename = file.path(out_dir, "period1_habitatCategories.tif"), overwrite = TRUE)
writeRaster(t2_score, filename = file.path(out_dir, "period2_habitatCategories.tif"), overwrite = TRUE)

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
writeRaster(total_score,
            filename = file.path(out_dir, "total_score.tif"),
            overwrite = TRUE)

# 2) Save stability_numeric
writeRaster(stability_numeric,
            filename = file.path(out_dir, "stability_layer.tif"),
            overwrite = TRUE)

# --------------------------------------------------
# 5) Plot numeric Restoration Suitability map
# --------------------------------------------------
cols <- c("red","orange","yellow","green")  # Instable → Stable-Ideal
plot(stability_numeric, col = cols, axes = FALSE, legend = FALSE, main="Stability Map")
legend("topright", legend=c("Not suitable","Moderately suitable","Suitable","Highly suitable"),
       fill=cols, bty="n")


# -------------------------------------------------
# 6) Optional: visualize missing cells
# -------------------------------------------------
na_map <- total_score
na_map[!is.na(na_map)] <- 0
na_map[is.na(na_map)] <- 1

plot(na_map, col=c("transparent","black"), main="Missing cells (NA)")


# Assume stability_numeric is your numeric raster (1=Instable, 2=Stable-Moderate, 3=Stable-Recommended, 4=Stable-Ideal)
# Convert to factor
stability_cat <- as.factor(stability_numeric)

# Create levels table
cats <- data.frame(
  ID = 1:4,
  category = c("Instable", "Stable-Moderate", "Stable-Recommended", "Stable-Ideal")
)

# Assign levels
levels(stability_cat) <- cats

writeRaster(stability_cat,
            filename = file.path(out_dir, "stability_layer_categoric.tif"),
            overwrite = TRUE)



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

writeRaster(change_map, filename = file.path(out_dir, "change_suitability_M7.tif"), overwrite = TRUE)

# Step 4: Plot changes
cols <- c("gray", "green", "red", "blue")  # no change, gain, loss, stable suitable
labels <- c("Unsuitable→Unsuitable", "Unsuitable→Suitable", 
            "Suitable→Unsuitable", "Suitable→Suitable")

plot(change_map, col = cols, axes = FALSE, legend = FALSE, main="Suitability change")
legend("topright", legend = labels, fill = cols, bty="n")


