




library(terra)
library(pbapply)
library(dplyr)

# ---------------------------
# Load or prepare your data
# ---------------------------
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres"

# Models 
glm_mod_s<- readRDS(paste(model_results_path,"glm_mod_s.rds", sep="/"))
gam_mod_s<- readRDS(paste(model_results_path,"gam_mod_s.rds", sep="/"))
rf_mod_s<-  readRDS(paste(model_results_path,"rf_mod_s.rds", sep="/"))
brt_mod_s<- readRDS(paste(model_results_path,"brt_mod_s.rds", sep="/"))

mod_list <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)

# Get variables used in the model
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- gsub("I\\((.+)\\)", "\\1", vars)   # unwrap I()
vars_clean <- gsub("\\^.*", "", vars_clean)      # drop powers (^2, ^3, etc.)
vars_selected <- unique(vars_clean)


# Load environemtnal conditions of bothe periods 
# Ensemble predictions for both periods (blob=env_warm and postblob=env_cold)
# env_warm <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif")   # SpatRaster with predictors as layers
# env_cold <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/ens_ave_postblob_M7.tif")
env_period_1 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/raster_stack_predict.tif")   # SpatRaster with predictors as layers
env_period_2 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/raster_stack_predict_postBlob.tif")


env_period_1<- env_period_1[[names(env_period_1)%in% vars_selected]]
env_period_2<- env_period_2[[names(env_period_2)%in% vars_selected]]
names(env_period_2)

pred_period_1 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/ensemble.tif")
pred_period_2 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/ensemble_postblob.tif")

# $ low suitability areas
pred_period_1_df<- as.data.frame(pred_period_1, xy = TRUE)
pred_period_1_df<- na.exclude(pred_period_1_df)

pred_period_1_df_lowSuitab<- pred_period_1_df%>%
  filter(ensemble <= 0.5)

summary(pred_period_1_df_lowSuitab$ensemble)

##=============================================================================
# --- Stack variables and suitability for each period
stack_period_1 <- c(env_period_1, pred_period_1)
stack_period_2 <- c(env_period_2, pred_period_2)

names(stack_period_1)[dim(stack_period_1)[3]] <- "HS_period_1"
names(stack_period_2)[dim(stack_period_2)[3]] <- "HS_period_2"

bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")

stack_period_1 <- c(stack_period_1, bathy20m)
stack_period_2 <- c(stack_period_2, bathy20m)

stack_period_1_df<- as.data.frame(stack_period_1[[c(2,3,4,7,8)]], xy = TRUE)


vals <- terra::extract( stack_period_1, pred_period_1_df_lowSuitab[, c("x", "y")])
vals<- cbind(pred_period_1_df_lowSuitab[, c("x", "y")], vals)

depper_vals<- vals%>% filter(coastwide_20m >= 10)%>%
  filter(y<= 400000)


ggplot(depper_vals, aes(x= x, y= y, color= HS_period_1))+
  geom_point()

ggplot(depper_vals, aes(x= slope_5x5, y= HS_period_1, color= coastwide_20m))+
  geom_point()+ geom_density(stat="identity")

plot(depper_vals$coastwide_20m, depper_vals$slope_5x5)

df<- depper_vals

# define bins for slope
# get bin labels actually present
bins <- levels(curve_slope$slope_bin)

# extract numeric limits (handles (), [], Inf)
lims <- do.call(
  rbind,
  strsplit(gsub("\\[|\\]|\\(|\\)", "", bins), ",")
)

lims <- apply(lims, 2, as.numeric)

# compute midpoints
curve_slope$slope_mid <- rowMeans(lims)

# remove bins with NA HS
curve_slope <- curve_slope[!is.na(curve_slope$HS_period_1), ]


# indices of non-empty bins
idx <- which(!is.na(curve_slope$HS_period_1))

curve_slope$slope_mid <- (brks[idx] + brks[idx + 1]) / 2

length(brks) - 1 != nrow(curve_slope)


plot(
  curve_slope$slope_mid,
  curve_slope$HS_period_1,
  type = "l",
  xlab = "Slope",
  ylab = "Mean HS",
  main = "HS response to slope"
)
