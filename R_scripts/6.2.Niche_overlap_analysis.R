


library(terra)
library(ecospat)

# Load environmental predictors (the same used to build models)
env_stack_t1 <- readRDS("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres/raster_stack_predict.rds")
set.seed(123)
env_all_t1 <- spatSample(env_stack_t1, size = 12000, method = "regular", na.rm = TRUE) #
env_all_t1 <- as.data.frame(env_all_t1)# Convert to data.frame

env_stack_t2 <- readRDS("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M5/raster_stack_predict_postBlob.rds")
set.seed(123)
env_all_t2 <- spatSample(env_stack_t2, size = 12000, method = "regular", na.rm = TRUE) #
env_all_t2 <- as.data.frame(env_all_t2)# Convert to data.frame


# Load ensemble suitability rasters for two periods
suit_t1 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres/tifs/ensemble_Average_suitability_blob_M5.tif")
suit_t2 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M5/ensemble_ave_postBlob_M5.tif")


# Define thresholds for moderate and optimal
thr.sui <- 0.27   # moderate
thr.opt <- 0.56   # optimal

# Create masks for moderate and optimal habitat
mask_t1_mod <- (suit_t1 > thr.sui & suit_t1 <= thr.opt)
mask_t1_opt <- (suit_t1 > thr.opt)

mask_t2_mod <- (suit_t2 > thr.sui & suit_t2 <= thr.opt)
mask_t2_opt <- (suit_t2 > thr.opt)

# Sample from moderate and optimal areas
n_samp <- 20000

env_t1_mod <- spatSample(env_stack, size = n_samp, mask = mask_t1_mod, method = "random", na.rm = TRUE)
env_t1_opt <- spatSample(env_stack, size = n_samp, mask = mask_t1_opt, method = "random", na.rm = TRUE)

env_t2_mod <- spatSample(env_stack, size = n_samp, mask = mask_t2_mod, method = "random", na.rm = TRUE)
env_t2_opt <- spatSample(env_stack, size = n_samp, mask = mask_t2_opt, method = "random", na.rm = TRUE)

# Convert to data.frame
env_t1_mod_df <- as.data.frame(env_t1_mod)
env_t1_opt_df <- as.data.frame(env_t1_opt)

env_t2_mod_df <- as.data.frame(env_t2_mod)
env_t2_opt_df <- as.data.frame(env_t2_opt)

### Run PCA-env
env_all_merged <- rbind(env_all_t1, env_all_t2)
pca <- prcomp(env_all_merged, center=TRUE, scale.=TRUE)


# Project moderate and optimal samples into PCA space
scores_t1_mod <- predict(pca, newdata = env_t1_mod_df)[, 1:2]
scores_t1_opt <- predict(pca, newdata = env_t1_opt_df)[, 1:2]

scores_t2_mod <- predict(pca, newdata = env_t2_mod_df)[, 1:2]
scores_t2_opt <- predict(pca, newdata = env_t2_opt_df)[, 1:2]

# PCA background points
scores_bg <- predict(pca, newdata = env_all_merged)[, 1:2]

z_mod_t1 <- ecospat.grid.clim.dyn(glob = scores_bg,
                                  glob1 = scores_bg,
                                  sp = scores_t1_mod, R = 100)

z_mod_t2 <- ecospat.grid.clim.dyn(glob = scores_bg,
                                  glob1 = scores_bg,
                                  sp = scores_t2_mod, R = 100)

# Same for optimal habitat:
z_opt_t1 <- ecospat.grid.clim.dyn(glob = scores_bg,
                                    glob1 = scores_bg,
                                    sp = scores_t1_opt, R = 100)

z_opt_t2 <- ecospat.grid.clim.dyn(glob = scores_bg,
                                  glob1 = scores_bg,
                                  sp = scores_t2_opt, R = 100)

# ------------------------------
# Plot Period 1
# ------------------------------
plot(scores_bg, col = rgb(0.8,0.8,0.8,0.3), pch = 16, cex = 0.5,
     xlab = "PC1", ylab = "PC2", main = "Period 1: Moderate & Optimal Habitat")
points(scores_t1_mod, col = rgb(0,1,0,0.5), pch = 16, cex = 0.6)  # moderate
points(scores_t1_opt, col = rgb(0,0.6,0,0.7), pch = 16, cex = 0.6) # optimal
legend("topright", legend = c("Available","Moderate","Optimal"),
       col = c(rgb(0.8,0.8,0.8,0.3), rgb(0,1,0,0.5), rgb(0,0.6,0,0.7)), pch = 16)

# ------------------------------
# Plot Period 2
# ------------------------------
plot(scores_bg, col = rgb(0.8,0.8,0.8,0.3), pch = 16, cex = 0.5,
     xlab = "PC1", ylab = "PC2", main = "Period 2: Moderate & Optimal Habitat")
points(scores_t2_mod, col = rgb(1,0.8,0,0.5), pch = 16, cex = 0.6)  # moderate
points(scores_t2_opt, col = rgb(1,0,0,0.6), pch = 16, cex = 0.6)    # optimal
legend("topright", legend = c("Available","Moderate","Optimal"),
       col = c(rgb(0.8,0.8,0.8,0.3), rgb(1,0.8,0,0.5), rgb(1,0,0,0.6)), pch = 16)












# Extract values above threshold
mask_t1 <- suit_t1 > thr.sui
mask_t1[] <- as.numeric(mask_t1[])
mask_t2 <- suit_t2 > thr.sui
mask_t2[] <- as.numeric(mask_t2[])

# Sample from suitable areas
n_samp= 20000
env_t1_sample <- spatSample(env_stack, size = n_samp, mask = mask_t1, method = "regular", na.rm = TRUE)
env_t1_df <- as.data.frame(env_t1_sample)

# Period 2
env_t2_sample <- spatSample(env_stack, size = n_samp, mask = mask_t2, method = "regular", na.rm = TRUE)
env_t2_df <- as.data.frame(env_t2_sample)


### Run PCA-env
library(ecospat)
pca <- prcomp(env_all, center=TRUE, scale.=TRUE)

# Project suitable habitats into PCA space
scores_t1 <- predict(pca, newdata=env_t1_df)[, 1:2]
scores_t2 <- predict(pca, newdata=env_t2_df)[, 1:2]
scores_bg <- predict(pca, newdata=env_all)[, 1:2]


z1 <- ecospat.grid.clim.dyn(glob = scores_bg,
                            glob1 = scores_bg,
                            sp = scores_t1, R = 100)

z2 <- ecospat.grid.clim.dyn(glob = scores_bg,
                            glob1 = scores_bg,
                            sp = scores_t2, R = 100)

### Test niche overlap & shifts
# Quantify overlap with Schoener’s D, I, niche equivalency, and similarity tests.
# Schoener's D
ecospat.niche.overlap(z1, z2, cor = TRUE)

# Niche equivalency test
ecospat.niche.equivalency.test(z1, z2, rep = 100)

# Niche similarity test
ecospat.niche.similarity.test(z1, z2, rep = 10)

### Plot the PCA-env
par(mfrow = c(1,2))
ecospat.plot.niche(z1, title = "Period 1")
ecospat.plot.niche(z2, title = "Period 2")

# Centroids of suitable niche (period 1 and 2)
centroid_t1 <- c(mean(z1$z[!is.na(z1$z)]))  # weighted or simple mean of PCA scores
centroid_t1 <- c(colMeans(z1$sp))  # ecospat stores points in $sp

centroid_t2 <- colMeans(z2$sp)

par(mfrow = c(1,2))
ecospat.plot.niche(z1, title = "PCA-env: Niche Shift")
ecospat.plot.niche(z2, title = "PCA-env: Niche Shift")




# PLOTlibrary(terra)
library(grDevices)
# --- Thresholds ---
thr_moderate <- 0.27
thr_optimal  <- 0.56

# --- Example for Period 1 ---
r1 <- z1$z               # raster with predicted suitability for period 1
xy <- as.data.frame(crds(r1))
vals1 <- values(r1)

# Available environment
xy_available <- xy[!is.na(vals1), ]

# Moderate habitat
xy_mod <- xy[vals1 >= thr_moderate & vals1 < thr_optimal, ]

# Optimal habitat
xy_opt <- xy[vals1 >= thr_optimal, ]

# --- Plot ---
plot(1, type="n",
     xlim=range(xy$PC1), ylim=range(xy$PC2),
     xlab="PC1", ylab="PC2",
     main="Period 1: Moderate and Optimal Habitat")

# Full environment in light grey
points(xy_available$PC1, xy_available$PC2, col=rgb(0.8,0.8,0.8,0.3), pch=16, cex=0.5)

# Moderate habitat
points(xy_mod$PC1, xy_mod$PC2, col=rgb(0,1,0,0.5), pch=16, cex=0.6)

# Optimal habitat
points(xy_opt$PC1, xy_opt$PC2, col=rgb(0,0.6,0,0.7), pch=16, cex=0.6)

legend("topright",
       legend=c("Available", "Moderate", "Optimal"),
       col=c(rgb(0.8,0.8,0.8,0.3), rgb(0,1,0,0.5), rgb(0,0.6,0,0.7)),

