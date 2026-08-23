# #===============================================================================
# ### Create Substrate mask            ===========================================
substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/SOG_substrate_20m.tif")
substrate_west<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/WCVI_substrate_20m.tif")
substrate_north<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Substrate/QCS_substrate_20m.tif")

substrate<- merge(substrate, substrate_north, substrate_west)
plot(substrate)
# The predicted raster files are classified as follows:
# 1) Rock,
# 2) Mixed,
# 3) Sand,
# 4) Mud


# Mask substrate to study area
substrate<- crop(substrate, ens_average)

# Align substrate layer to model prediction layer
substrate_aligned <- terra::rast(ens_average)
substrate_aligned <- terra::resample(substrate, substrate_aligned, method = "near")

# Convert categories of substrate in 1 (hard) and 2 (soft) substrate
# substrate_aligned[substrate_aligned == 1]<- 1 # rocky is already 1
substrate_aligned[substrate_aligned == 2]<- 1 # 2 was mixed substrate
substrate_aligned[substrate_aligned == 3]<- 2 # 3 and 4 were mud and sand (I guess)
substrate_aligned[substrate_aligned == 4]<- 2 # 3 and 4 were mud and sand (I guess)

# Save substrate, 1= hard substrate; 2= soft substrate
# writeRaster(substrate_aligned, "substrate_SOG_aligned.tif", overwrite=T)