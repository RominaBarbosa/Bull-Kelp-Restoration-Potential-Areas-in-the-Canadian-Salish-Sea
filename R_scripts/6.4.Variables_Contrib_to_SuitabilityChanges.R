

# sensitivity / gradient mapping approach with an ensemble (average) of models (GLM, GAM, RF, BRT). It:
# Computes the ensemble prediction for the warm baseline (once).
# For each predictor, creates a warm stack where only that predictor is replaced by its cold-period layer.
# Predicts with each model (streaming to disk to avoid memory overload), averages model predictions to produce the ensemble prediction for that “single-variable-changed” scenario.
# Computes per-variable contribution rasters = (ensemble_with_only_var_changed − ensemble_warm).
# Writes contribution rasters to disk and computes a dominant-driver raster (which variable has the largest contribution at each pixel).

# Memory-efficient sensitivity mapping with terra + ensemble averaging
# Requirements: terra, pbapply (for progress), optionally foreach/doParallel if you want parallel
# Assumptions:
#  - env_warm and env_cold are SpatRasters (same extent/resolution/CRS and same layer order)
#  - names(env_warm) give predictor names
#  - mod_list is a list of fitted models. Each element can be:
#       * a model object (glm, gam, randomForest, gbm, etc.) OR
#       * a list with elements: $model (the model object) and $predict_args (a named list passed to predict)
#  - terra::predict works for your model objects (it calls base predict(model, newdata)). If one model needs special args
#    (e.g. gbm: n.trees), put them in predict_args.
#
# Example mod_list item:
#   mod_list <- list(glm_model, gam_model, list(model = rf_model, predict_args = list(type="response")),
#                    list(model = brt_model, predict_args = list(n.trees=1500)))
#
# Edit these paths and names:
# out_dir <- "variables_contrib_change_outputs"   # directory to store temporary rasters and results
# dir.create(out_dir, showWarnings = FALSE)
# terraOptions(tempdir = out_dir)     # make terra use this directory for temp files

library(terra)
library(pbapply)

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

pred_period_1 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif")
pred_period_2 <- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/ens_ave_postblob_M7.tif")

##=============================================================================
# --- Stack variables and suitability for each period
stack_period_1 <- c(env_period_1, pred_period_1)
stack_period_2 <- c(env_period_2, pred_period_2)

names(stack_period_1)[dim(stack_period_1)[3]] <- "HS_period_1"
names(stack_period_2)[dim(stack_period_1)[3]] <- "HS_period_2"

# --- Create a regular grid over raster extent to perform the analysis
spacing= 500
r <- env_period_2[[2]]
ext_r <- ext(r)
crs(r) <- "EPSG:3005"
x_coords <- seq(ext_r[1] + spacing/2, ext_r[2], by = spacing)
y_coords <- seq(ext_r[3] + spacing/2, ext_r[4], by = spacing)

grid <- expand.grid(x = x_coords, y = y_coords)
grid_points <- vect(grid, geom = c("x","y"), crs = crs(r))

values_at_points <- extract(r[[1]], grid_points)
grid_valid <- grid_points[!is.na(values_at_points[,2]), ]  # keep only non-NA cells

plot(r[[1]], main = "0.5 km spaced points over raster")
points(grid_valid, col = "red", pch = 19)

# Extract values
vals_period_1<- terra::extract(stack_period_1, grid_valid)
vals_period_2 <- terra::extract(stack_period_2, grid_valid)
length(vals_period_2$ammonium_spring_mean)

# Get coordinates of points
coords <- crds(grid_valid)   # returns matrix of x and y

# Combine coordinates with raster values
vals_period_1_df <- cbind(x = coords[,1], y = coords[,2], vals_period_1[,-1])  # drop ID
vals_period_2_df <- cbind(x = coords[,1], y = coords[,2], vals_period_2[,-1])


# --- Keep overlapping cells
df_join <- inner_join(vals_period_1_df, vals_period_2_df, by = c("x", "y"),
                      suffix = c("_period_1", "_period_2"))

# --- Compute Δ (differences)
vars <- names(env_period_1)

for (v in vars) {
  df_join[[paste0("d_", v)]] <- df_join[[paste0(v, "_period_2")]] - df_join[[paste0(v, "_period_1")]]
}

# ΔHS = difference in predicted suitability
df_join$delta_HS <- df_join$HS_period_2 - df_join$HS_period_1

# --- Keep only relevant columns
delta_cols <- c("x", "y", "delta_HS", paste0("d_", vars))
delta_df <- df_join[, delta_cols]
delta_df <- na.omit(delta_df)

# ---  Fit regression model
# install.packages("QuantPsyc")  # only once
library(QuantPsyc)

lm_delta <- lm(delta_HS ~ ., data = delta_df[, -c(1,2)])  # exclude coordinates
summary(lm_delta)

# ---  Compute standardized coefficients (comparable importance)
lm_beta <- lm.beta(lm_delta)
summary(lm_beta)

# Convert to dataframe
coef_df <- data.frame(
  variable = names(lm_beta),
  beta = as.numeric(lm_beta)
)


library(broom)

# Standardized model
lm_beta <- lm.beta(lm_delta)

# Extract coefficients + p-values
coef_df <- broom::tidy(lm_delta)
coef_df <- coef_df[coef_df$term != "(Intercept)", ]  # remove intercept

# Add standardized beta
coef_df$beta <- as.numeric(lm_beta)  # use lm_beta directly
coef_df$variable <- names(lm_beta)   # variable names

# Define nice labels (same as for previous plots)
nice_labels <- c(
  d_ammonium_spring_mean     = "Ammonium (mmol m⁻³)",
  d_currentSpeed_summer_mean = "Current Speed (m s⁻¹)",
  d_PAR_summer_mean          = "PAR (µE m⁻² s⁻¹)",
  d_temperature_summer_mean  = "Temperature (°C)",
  d_turbidity_summer_mean    = "Turbidity (NTU)",
  d_salinity_summer_mean     = "Salinity (g kg⁻¹)",
  d_slope_5x5                = "Slope (°)"
)

# Add a column with the nice labels
coef_df <- coef_df %>%
  mutate(variable_label = nice_labels[variable])

# Plot with updated labels
ggplot(coef_df, aes(x = reorder(variable_label, beta), y = beta, fill = beta > 0)) +
  geom_col() + 
  geom_text(aes(label = signif), 
                         hjust = ifelse(coef_df$beta > 0, -0.2, 1.2), # place outside the bar
                         size = 3) +
  ylim(c(-0.23, 0.25))+
  coord_flip() +
  scale_fill_manual(
    values = c("red", "blue"),
    name = "Effect direction",
    labels = c("Negative", "Positive")
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.position = c(0.05, 0.95),      # moves legend inside plot
    legend.justification = c("left", "top")
  ) +
  labs(
    title = "Standardized regression coefficients (ΔHS model)",
    x = "Variable",
    y = "Standardized β (relative importance)"
  )

# ggsave("regression_effect_of_vars_inDeltaSuitability.png", width = 12, height = 8, units="cm",  dpi = 350)




# Compute actual contribution of each variable at each cell
delta_df
delta_vars_df_contrib <- delta_df
for (v in vars) {
  beta_v <- coef_df$beta[coef_df$variable == v]
  delta_vars_df_contrib[[v]] <- delta_df[[v]] * beta_v
}

# Sum contributions for total ΔHS
delta_vars_df_contrib$total <- rowSums(delta_vars_df_contrib[vars])

# For plotting: color by positive/negative contribution
delta_vars_df_contrib$dominant_sign <- ifelse(delta_vars_df_contrib$total >= 0, "Positive", "Negative")





# Vector of variable names (remove "x", "y", "delta_HS")
var_cols <- grep("^d_", colnames(delta_df), value = TRUE)

# Create a contribution dataframe
contrib_df <- delta_df[, c("x", "y", "delta_HS")]
coef_df$variable <- as.character(coef_df$variable)

for (v in vars) {
  delta_col <- paste0("d_", v)        # column in delta_df
  coef_name <- paste0("d_", v)        # name in coef_df$variable
  
  if (!(delta_col %in% colnames(delta_df))) stop("Column not found: ", delta_col)
  
  # extract beta as numeric
  beta_v <- as.numeric(coef_df$beta[coef_df$variable == coef_name])
  if (length(beta_v) == 0) stop("Beta not found for variable: ", v)
  
  # contribution = ΔX * beta
  contrib_df[[v]] <- delta_df[[delta_col]] * beta_v
}

# total predicted ΔHS per location
contrib_df$predicted_deltaHS <- rowSums(contrib_df[, vars])



# ---------------------------
# 1. Reshape to long format for plotting
# ---------------------------
contrib_long <- contrib_df %>%
  pivot_longer(cols = all_of(vars), 
               names_to = "variable", 
               values_to = "contribution")

# Add sign information
contrib_long <- contrib_long %>%
  mutate(sign = ifelse(contribution >= 0, "Positive", "Negative"))

# Optional: scale magnitude by total predicted deltaHS
contrib_long <- contrib_long %>%
  mutate(prop_contrib = contribution / predicted_deltaHS)



# ---------------------------
# 1. Get top 3 contributors per row
# ---------------------------
contrib_only <- contrib_df[, vars]

top_n <- 3

# Initialize list to store results
top_vars_list <- lapply(seq_len(nrow(contrib_only)), function(i) {
  row <- as.numeric(contrib_only[i, ])
  names(row) <- vars
  ord <- order(abs(row), decreasing = TRUE)[1:top_n]
  data.frame(
    rank = 1:top_n,
    variable = names(row)[ord],
    contribution = row[ord],
    x = contrib_df$x[i],
    y = contrib_df$y[i]
  )
})

# Combine into a single dataframe
top_vars_df <- do.call(rbind, top_vars_list)

# Add sign
top_vars_df$sign <- ifelse(top_vars_df$contribution >= 0, "Positive", "Negative")
top_vars_df$rank<- as.numeric(top_vars_df$rank)

# ---------------------------
# Plot top 3 contributors by rank
# ---------------------------
plot_dominats<- top_vars_df%>%
  filter(rank == 1 | rank == 2 | rank ==3)%>%
  ggplot( aes(x = x, y = y, color = contribution)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                        name = "ΔHS contribution") +
  facet_wrap(~ rank, ncol = 1, labeller = labeller(rank = c("1" = "PAR", "2" = "Salinity", "3" = "Turbidity"))) +
  theme_minimal() +
  labs(title = "Top 3 dominant variables per location",
       x = "Longitude", y = "Latitude")



# Find dominant variable per row
dominant_var <- apply(abs(contrib_only), 1, function(x) vars[which.max(x)])

# Extract dominant contribution values
dominant_contrib <- mapply(function(r, v) contrib_only[r, v], 
                           r = seq_len(nrow(contrib_only)), v = dominant_var)

# Determine sign
sign_contrib <- ifelse(dominant_contrib >= 0, "Positive", "Negative")

# Combine into a dataframe
dominant_df <- data.frame(
  x = contrib_df$x,
  y = contrib_df$y,
  dominant_var = dominant_var,
  dominant_contrib = dominant_contrib,
  sign = sign_contrib
)


# ---------------------------
# 2. Compute residuals (ΔHS not explained by variables)
# ---------------------------
residual_df <- contrib_df %>%
  mutate(residual = delta_HS - predicted_deltaHS) %>%
  dplyr::select(x, y, residual) %>%
  mutate(sign = ifelse(residual >= 0, "Positive", "Negative"))

# ---------------------------
# 4. Plot residual contribution
# ---------------------------
p_residual <- ggplot(residual_df, aes(x = x, y = y, color = residual, shape = sign)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                        name = "Residual ΔHS") +
  theme_minimal() +
  labs(title = "Residual habitat suitability change (not explained by dominant variable)",
       x = "Longitude", y = "Latitude")

# ---------------------------
# 5. Print plots
# ---------------------------
plot_dominats +
p_residual


# Prepare a long-format dataframe for plotting
delta_long <- delta_df %>%
  pivot_longer(
    cols = starts_with("d_"), 
    names_to = "variable", 
    values_to = "delta"
  )%>%
  group_by(variable) %>%
  mutate(delta_rescaled = scales::rescale(delta, to = c(-1, 1))) %>%
  ungroup()


# Clean variable names for nicer facet titles
delta_long$variable <- sub("^d_", "", delta_long$variable)

# Plot each variable with its own color scale
ggplot(delta_long, aes(x = x, y = y, color = delta_rescaled)) +
  geom_point(size = 1, alpha = 0.8) +
  facet_wrap(~ variable) +
  scale_color_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Change (Δ) in environmental variables between periods",
    color = "Δ value\n(Period2 - Period1)",
    x = "Longitude",
    y = "Latitude"
  )



## Plot each variable with its own scale of color of change
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)

# Get map 
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Optional: crop to your study area
bbox <- st_bbox(delta_sf_deg)
world_crop <- st_crop(world, bbox)


# Convert coordinates to sf object with EPSG:3005
delta_sf <- st_as_sf(delta_df, coords = c("x", "y"), crs = 3005)

# Transform to geographic coordinates (degrees)
delta_sf_deg <- st_transform(delta_sf, 4326)

# Extract new lon/lat columns
delta_df_deg <- delta_sf_deg %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# List of Δ-variable names
delta_vars <- grep("^d_", names(delta_df_deg), value = TRUE)

# Create a named vector to map variable codes to pretty labels
nice_labels <- c(
  d_ammonium_spring_mean     = "Ammonium\n(mmol m⁻³)",
  d_currentSpeed_summer_mean = "Current Speed\n(m s⁻¹)",
  d_PAR_summer_mean          = "PAR\n(µE m⁻² s⁻¹)",
  d_temperature_summer_mean  = "Temperature\n(°C)",
  d_turbidity_summer_mean    = "Turbidity\n(NTU)",
  d_salinity_summer_mean     = "Salinity\n(g kg⁻¹)",
  d_slope_5x5                = "Slope\n(°)"
)

# Create an empty list to store plots
plot_list <- list()

for (v in delta_vars) {
  # var_label <- sub("^d_", "", v)
  var_label <- nice_labels[[v]] 
  
  p <- ggplot() +
    # 1. Add land polygons
    geom_sf(data = world_crop, fill = "gray90", color = "gray60") +
    
    # 2. Add your sampled points (lon/lat columns)
    geom_point(
      data = delta_df_deg,
      aes(x = lon, y = lat, color = .data[[v]]),
      size = 0.3,
      alpha = 0.8
    ) +
    
    # 3. Color scale for change
    scale_color_gradient2(
      low = "red", mid = "yellow", high = "blue",
      midpoint = 0,
      name = paste0("Δ ", var_label)
    ) +
    
    # 4. Coordinate system and theme
    coord_sf() +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 7),
      legend.title = element_text(size = 5),
      legend.text  = element_text(size = 4),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width  = unit(0.1, "cm"),
      legend.position = c(0.05, 0.05),      # moves legend inside plot
      legend.justification = c("left", "bottom")
    ) +
    labs(
      # title = paste("Change in", var_label, "(Period 2 - Period 1)"),
      x = "Longitude (°)",
      y = "Latitude (°)"
    )
  
  plot_list[[v]] <- p
}

# Combine all plots together
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3, align= "hv")
combined_plot

# # Combine all plots into one grid
# combined_plot <- wrap_plots(plot_list, ncol = 2) +
#   plot_annotation(title = "")

# Save
setwd("/Volumes/Romina_PSF/PSF/SDM/Restoration_areas")
# ggsave("delta_all_variables.pdf", combined_plot, width = 18, height = 20, units="cm",  dpi = 350)





# ---------------------------
# Plot top 3 contributors by rank
# ---------------------------
# Convert coordinates to sf object with EPSG:3005
top_vars_df <- st_as_sf(top_vars_df, coords = c("x", "y"), crs = 3005)

# Transform to geographic coordinates (degrees)
top_vars_sf_deg <- st_transform(top_vars_df, 4326)

# Extract new lon/lat columns
top_vars_sf_deg <- top_vars_sf_deg %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# Filter top 3 contributors
top3 <- top_vars_sf_deg %>% filter(rank %in% 1:3)

# Convert to sf again in degrees for geom_sf
top3_sf <- st_as_sf(top3, coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = world_crop, fill = "gray90", color = "gray60") +
  geom_sf(data = top3_sf, aes(color = contribution), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                        name = "ΔHS contribution") + 
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(labels = function(y) sprintf("%.1f", y)) +
  facet_wrap(~ rank, ncol = 3, labeller = labeller(rank = c("1" = "PAR", 
                                                            "2" = "Salinity", 
                                                            "3" = "Turbidity"))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 5),
    legend.text  = element_text(size = 4),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width  = unit(0.1, "cm"),
    legend.position = c(0.05, 0.05),      # moves legend inside plot
    legend.justification = c("left", "bottom")
  ) +
  labs(title = "Top 3 dominant variables per location",
       x = "Longitude", y = "Latitude")

ggsave("top3_vars_contrib_to_HSuitability_change.png", width = 17, height = 8, units="cm",  dpi = 350)





# Assume your model is lm_delta
# lm_delta <- lm(delta_HS ~ d_ammonium + d_temperature + ..., data = delta_df)

# 1. Residual diagnostics
par(mfrow = c(2,2))
plot(lm_delta)  
# This shows:
# - Residuals vs Fitted: checks linearity & homoscedasticity
# - Normal Q-Q: checks residual normality
# - Scale-Location: checks homoscedasticity
# - Residuals vs Leverage: identifies influential points

# 2. Check multicollinearity
library(car)
vif(lm_delta)
# VIF > 5 or 10 indicates high collinearity

# 3. Check spatial autocorrelation of residuals
library(ape)      # or spdep for more complex analyses
# Create coordinates matrix
coords <- as.matrix(delta_df[, c("x","y")])
# Compute Moran's I
moran_res <- Moran.I(residuals(lm_delta), dnearneigh(coords, 0, 10000))  # example: 0-10 km neighbors
moran_res

# 4. Optional: histogram of residuals
hist(residuals(lm_delta), breaks = 30, main="Residuals distribution", xlab="Residuals")
# Notes:
# plot(lm_delta) gives a quick visual check of linearity, homoscedasticity, and normality.
# vif() identifies predictors that are highly correlated.
# Moran’s I tests whether residuals are spatially autocorrelated. Significant autocorrelation suggests you might need a spatial regression approach.




# Weigthed contribution by the importance of variable to the ensemble model
vars <- c("ammonium_spring_mean", "currentSpeed_summer_mean", "temperature_summer_mean")

contrib_df <- delta_df[, c("x", "y", "delta_HS")]

for (v in vars) {
  delta_col <- paste0("d_", v)
  beta_v <- coef_df$beta[coef_df$variable == paste0("d_", v)]
  importance_v <- importance_df$importance[importance_df$variable == v]
  
  # Weighted contribution
  contrib_df[[v]] <- delta_df[[delta_col]] * beta_v * (importance_v / 100)
}

# Total weighted ΔHS
contrib_df$weighted_deltaHS <- rowSums(contrib_df[, vars])
# Standardized β: shows relative sensitivity of ΔHS to variable changes.
# ΔX: actual change of the variable between periods.
# Importance weighting: respects the SDM’s learned relationships; small variable changes with high SDM importance can be highlighted.
# Spatial mapping: each grid cell can now have contributions from each variable, total weighted ΔHS, and residuals.

















































# Optional: model weights (defaults to equal weighting)
# Provide a numeric vector of length nmodels if you want performance-weighted ensemble
if(!exists("model_weights")){
  model_weights <- rep(1, nmodels)
} else {
  stopifnot(length(model_weights) == nmodels)
}
model_weights <- model_weights / sum(model_weights)

# ---------------------------
# Step 0: warm-period ensemble (baseline)
# ---------------------------
# We compute ensemble prediction for env_warm ONCE and store to disk.
warm_ensemble_file <- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/ensemble_ave_blob_M7.tif"
# cold_ensemble_file <- "/Volumes/Romina_PSF/PSF/SDM//Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/ens_ave_postblob_M7.tif"


# if(!file.exists(warm_ensemble_file)){
#   message("Computing ensemble prediction for the warm baseline (streaming to disk)...")
#   # For each model, predict and save a temporary raster
#   tmp_files <- pbapply::pblapply(seq_len(nmodels), function(i){
#     tf <- file.path(out_dir, paste0("pred_warm_model", i, ".tif"))
#     safe_predict_model(env_warm, mod_list[[i]], filename = tf)
#     return(tf)
#   })
#   # Read predicted rasters and compute weighted average (streamed)
#   preds_rasters <- lapply(tmp_files, rast)
#   # Weighted sum: multiply each raster by weight and sum; write result to disk
#   # We'll compute sum(weights * raster) using terra::app and stack
#   preds_stack <- rast(preds_rasters)
#   # multiply each layer by its weight
#   for(i in seq_len(nlayers(preds_stack))){
#     preds_stack[[i]] <- preds_stack[[i]] * model_weights[i]
#   }
#   # Sum layers into ensemble_warm
#   ensemble_warm <- app(preds_stack, sum, filename = warm_ensemble_file, overwrite=TRUE)
#   # optionally remove per-model warm predictions to save disk (keep only ensemble)
#   file.remove(unlist(tmp_files))
# } else {
#   message("Found existing warm ensemble raster; loading.")
#   ensemble_warm <- rast(warm_ensemble_file)
# }

# ---------------------------
# Step 0: Define sampling locations
# ---------------------------
# You can provide:
#   - coordinates (lon/lat or x/y)
#   - or randomly sample N pixels from env_warm

# Example 1: predefined locations (data.frame or matrix)
# locations <- data.frame(x = c(-127.5, -127.3), y = c(50.6, 50.8))
# Raster in EPSG:3005
r <- env_warm  # example raster

spacing <- 1000  # 1 km in meters

# Create a regular grid over raster extent
ext_r <- ext(r)
x_coords <- seq(ext_r[1] + spacing/2, ext_r[2], by = spacing)
y_coords <- seq(ext_r[3] + spacing/2, ext_r[4], by = spacing)

grid <- expand.grid(x = x_coords, y = y_coords)
grid_points <- vect(grid, geom = c("x","y"), crs = crs(r))

values_at_points <- extract(r[[1]], grid_points)
grid_valid <- grid_points[!is.na(values_at_points[,2]), ]  # keep only non-NA cells

plot(r[[1]], main = "1 km spaced points over raster")
points(grid_valid, col = "red", pch = 19)



# # Example 2: random sample of pixels (for demonstration)
# set.seed(42)
# locations <- spatSample(env_warm, size = 10000, method = "regular", na.rm = TRUE, as.points = TRUE)
# # Plot the raster
# plot(env_warm[[1]], main="Sampled locations on env_warm raster")
# points(locations, col="red", pch=19)
# 
# # Extract the cell indices of those locations
# sample_cells <- cellFromXY(env_warm, crds(locations))

# ---------------------------
# Step 1: Loop variables (but only for selected cells)
# ---------------------------
# Extract values from env_warm/env_cold
# vals_warm <- terra::extract(env_warm, sample_cells)[,, drop=FALSE]
# vals_cold <- terra::extract(env_cold, sample_cells)[,, drop=FALSE]

# Extract values
vals_warm <- terra::extract(env_warm, grid_valid)
vals_cold <- terra::extract(env_cold, grid_valid)
length(vals_warm$ammonium_spring_mean)

# Get coordinates of points
coords <- crds(grid_valid)   # returns matrix of x and y

# Combine coordinates with raster values
vals_warm_df <- cbind(x = coords[,1], y = coords[,2], vals_warm[,-1])  # drop ID
vals_cold_df <- cbind(x = coords[,1], y = coords[,2], vals_cold[,-1])

# Loop through variables, replace one at a time, predict, and compute contributions
source("~/Documents/GitHub/SDMs_Broughton/FINAL_scripts/require_scaling_parameters.R")
source("~/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/predict_model_sdm_function.R")

contrib_list <- vector("list", length = nvars)
model_algos <- sapply(mod_list, function(m) {
  cls <- tolower(class(m)[1])
  
  if (grepl("glm", cls)) return("glm")
  if (grepl("gam", cls)) return("gam")
  if (grepl("randomforest", cls)) return("rf")
  if (grepl("gbm|brt", cls)) return("brt")
  
  return(cls)
})

model_algos
# [1] "glm" "gam" "rf"  "brt"

for (j in seq_len(nvars)) {
  varname <- pred_names[j]
  message(sprintf("[%d/%d] Processing variable: %s", j, nvars, varname))
  
  # replace values of variable j by conditions in period 2
  vals_temp <- vals_warm
  vals_temp[[varname]] <- vals_cold[[varname]] # replace values of variable j by conditions in period 2
  
  # --- Predict for each model with correct scaling ---
  preds_j <- sapply(seq_along(mod_list), function(i) {
    predict_model_sdm(
      model = mod_list[[i]],
      model_name = model_algos[i],
      newdata = vals_temp,
      scale_params = scale_params,
      quad_vars = vars_selected
    )
  })  
  
  preds_warm <- sapply(seq_along(mod_list), function(i) {
    predict_model_sdm(
      model = mod_list[[i]],
      model_name = model_algos[i],
      newdata = vals_warm,
      scale_params = scale_params,
      quad_vars = vars_selected
    )
  })
  
  # --- Build weighted ensemble predictions ---
  ensemble_var  <- preds_j    %*% model_weights
  ensemble_warm <- preds_warm %*% model_weights
  
  # --- Compute contribution (difference) ---
  contrib_vec <- ensemble_var - ensemble_warm
  
  contrib_list[[j]] <- as.numeric(contrib_vec)
  
  # rm(vals_temp, preds_j, preds_warm, ensemble_var, ensemble_warm)
  # gc()
}

# Combine into a data.frame
contrib_df <- as.data.frame(contrib_list)
names(contrib_df) <- vars_selected


# ---------------------------
# Step 2: Combine results into a data frame
# ---------------------------
contrib_df <- as.data.frame(do.call(cbind, contrib_list))
colnames(contrib_df) <- pred_names
contrib_df$x <- vals_cold_df[, "x"]
contrib_df$y <- vals_cold_df[, "y"]

# Reorder columns for clarity
contrib_df <- contrib_df[, c("x", "y", pred_names)]

# ---------------------------
# Step 3: Optional summaries
# ---------------------------
# Extract contribution columns only
contrib_only <- contrib_df[, pred_names]

# 1. Dominant driver index (largest absolute contribution)
dominant_abs_idx <- apply(abs(contrib_only), 1, which.max)
dominant_abs_var <- pred_names[dominant_abs_idx]

# 2. Sign of the dominant contribution
dominant_sign <- mapply(function(row, idx) {
  sign(row[idx])
}, split(contrib_only, seq(nrow(contrib_only))), dominant_abs_idx)

# 3. Magnitude of the dominant contribution (absolute value)
dominant_magnitude <- mapply(function(row, idx) {
  row[idx]
}, split(contrib_only, seq(nrow(contrib_only))), dominant_abs_idx)

# Combine into the dataframe
contrib_df$dominant_var <- dominant_abs_var
contrib_df$dominant_sign <- dominant_sign   # +1 = positive, -1 = negative
contrib_df$dominant_magnitude <- dominant_magnitude


contrib_df$dominant_factor <- factor(
  paste0(ifelse(contrib_df$dominant_sign > 0, "pos_", "neg_"),
         contrib_df$dominant_var),
  levels = c(paste0("pos_", pred_names), paste0("neg_", pred_names))
)

# 
# ggplot(contrib_df, aes(x = x, y = y, color = dominant_factor, size = abs(dominant_magnitude))) +
#   geom_point(alpha = 0.8) +
#   scale_color_brewer(palette = "Paired", name = "Dominant driver") +
#   scale_size_continuous(name = "Contribution magnitude") +
#   theme_minimal() +
#   labs(title = "Dominant driver of suitability change",
#        subtitle = "Magnitude and positive/negative impact",
#        x = "Longitude", y = "Latitude")
# 

### Compute residual contribution
# Extract contributions matrix
contrib_only <- contrib_df[, pred_names]

# Dominant variable index per row
dominant_idx <- apply(abs(contrib_only), 1, which.max)

# Compute residual change = total change minus contribution of dominant variable
residual_change <- mapply(function(row, idx) {
  # total change = sum of all contributions
  sum(row) - row[idx]
}, split(contrib_only, seq(nrow(contrib_only))), dominant_idx)

# Add to dataframe
contrib_df$residual_change <- residual_change

contrib_df$dominant_magnitude <- as.numeric(as.character(contrib_df$dominant_magnitude))
contrib_df$residual_change <- as.numeric(as.character(contrib_df$residual_change))
# Absolute values for plotting
contrib_df$dominant_magnitude_abs <- abs(contrib_df$dominant_magnitude)
# Compute symmetric limits for color scale
max_val <- max(abs(c(contrib_df$dominant_magnitude, contrib_df$residual_change)), na.rm = TRUE)
color_limits <- c(-max_val, max_val)

library(ggplot2)
library(patchwork)

p1 <- ggplot(contrib_df, aes(x = x, y = y,
                             shape = dominant_factor,
                             color = dominant_magnitude)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue",
                        midpoint = 0,
                        limits = color_limits,
                        name = "Contribution") +
  theme_minimal() +
  labs(title = "Dominant driver of suitability change")

p2 <- ggplot(contrib_df, aes(x = x, y = y, color = residual_change)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue",
                        midpoint = 0,
                        limits = color_limits,
                        name = "Residual Contribution") +
  theme_minimal() +
  labs(title = "Residual change not explained by dominant driver")

# Combine plots side by side
p1 + p2


p3 <- contrib_df%>%
  filter(dominant_magnitude >= 0)%>%
  ggplot( aes(x = x, y = y,
              color = dominant_factor,
              size = dominant_magnitude)) +
  geom_point(alpha = 0.8) +
  # scale_color_gradient2(low = "red", mid = "white", high = "blue",
  #                       midpoint = 0,
  #                       limits = color_limits,
  #                       name = "Contribution") +
  theme_bw() +
  facet_wrap(~dominant_factor)+
  labs(title = "Dominant driver of suitability change")

p4 <- contrib_df%>%
  filter(dominant_magnitude < 0)%>%
  ggplot( aes(x = x, y = y,
              color = dominant_factor,
              size = dominant_magnitude)) +
  geom_point(alpha = 0.8) +
  # scale_color_gradient2(low = "red", mid = "white", high = "blue",
  #                       midpoint = 0,
  #                       limits = color_limits,
  #                       name = "Contribution") +
  theme_minimal() +
  labs(title = "Dominant driver of suitability change")


p5 <- ggplot(contrib_df, aes(x = x, y = y, color=dominant_factor,
                             size = residual_change)) +
  geom_point( alpha = 0.8) +
  # scale_color_gradient2(low = "red", mid = "white", high = "blue",
  #                       midpoint = 0,
  #                       limits = color_limits,
  #                       name = "Residual Contribution") +
  theme_minimal() +
  facet_wrap(~dominant_factor)+
  labs(title = "Residual change not explained by dominant driver")

# Combine plots side by side
p1 + p2








# Identify dominant variable per cell
dominant_idx <- apply(abs(contrib_df[, pred_names]), 1, which.max)
dominant_var <- pred_names[dominant_idx]

# Extract the contribution and add to dataframe
dominant_contrib <- mapply(function(row, idx) row[idx],
                           split(contrib_df[, pred_names], seq(nrow(contrib_df))),
                           dominant_idx)

# Compute the corresponding change in the environmental variable
env_change <- mapply(function(varname, idx) {
  vals_cold[idx, varname] - vals_warm[idx, varname]
}, dominant_var, seq_along(dominant_var))

# Add to dataframe
contrib_df$dominant_var <- dominant_var
contrib_df$dominant_contrib <- dominant_contrib
contrib_df$dominant_env_change <- env_change
contrib_df$dominant_contrib <- as.numeric(contrib_df$dominant_contrib)
contrib_df$dominant_env_change <- as.numeric(contrib_df$dominant_env_change)

# dominant_contrib → contribution to suitability change
# dominant_env_change → actual change in the variable driving it

# Summarize magnitude across landscape
summary_df <- contrib_df %>%
  group_by(dominant_var) %>%
  summarise(mean_contribution = mean(dominant_contrib, na.rm = TRUE),
            median_contribution = median(dominant_contrib, na.rm = TRUE),
            mean_env_change = mean(dominant_env_change, na.rm = TRUE),
            median_env_change = median(dominant_env_change, na.rm = TRUE),
            n_cells = n()) %>%
  arrange(desc(abs(mean_contribution)))



# Contribution vs environmental change
plot_df <- contrib_df[!is.na(contrib_df$dominant_contrib) & !is.na(contrib_df$dominant_env_change), ]
ggplot(plot_df, aes(x = dominant_env_change, y = dominant_contrib,
                    color = dominant_var)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_brewer(palette = "Paired", name = "Dominant variable") +
  theme_minimal(base_size = 14) +
  facet_wrap(~dominant_var)+
  labs(
    x = "Change in environmental variable",
    y = "Contribution to suitability change",
    title = "Relationship between environmental change and suitability contribution",
    subtitle = "Points colored by the dominant driver"
  )





# Highlights areas where a single variable explains most of the change vs areas with distributed drivers.
ggplot(contrib_df, aes(x=abs(dominant_contrib), y=abs(residual_change), color=dominant_var)) +
  geom_point(alpha=0.6) +
  geom_abline(slope=1, linetype="dashed") +
  theme_minimal() +
  labs(title="Residual vs dominant contribution",
       x="Dominant contribution magnitude",
       y="Residual contribution magnitude")









####### See ranked variables contribution =====================================
# 1. Convert to long format
contrib_long <- contrib_df %>%
  pivot_longer(cols = all_of(vars_selected),
               names_to = "variable",
               values_to = "contrib") %>%
  group_by(x, y) %>%
  mutate(
    total_change = sum(contrib, na.rm = TRUE),  # ensure we have total change
    abs_contrib = abs(contrib),
    rank = rank(-abs_contrib, ties.method = "first"),      # 1 = most important
    contrib_prop = abs_contrib / sum(abs_contrib, na.rm = TRUE),  # proportion of total magnitude
    change_type = ifelse(total_change >= 0, "Positive", "Negative")
  ) %>%
  ungroup()

contrib_long <- contrib_long %>%
  group_by(x, y) %>%
  mutate(
    total_change = sum(contrib, na.rm = TRUE),
    abs_contrib = abs(contrib),
    abs_total_change = sum(abs(contrib), na.rm = TRUE),
    contrib_prop = abs_contrib / abs_total_change,
    weighted_contrib = contrib_prop * abs(total_change),
    rank = rank(-abs_contrib, ties.method = "first"),
    change_type = ifelse(total_change >= 0, "Positive", "Negative")
  ) %>%
  ungroup()

# Identify top two contributors based on absolute contribution
dominant_df <- contrib_long %>%
  filter(rank == 1) %>%
  select(x, y, variable, weighted_contrib, contrib_prop, change_type) %>%
  rename(first_var = variable,
         first_weighted = weighted_contrib,
         first_prop = contrib_prop)

second_df <- contrib_long %>%
  filter(rank == 2) %>%
  select(x, y, variable, weighted_contrib, contrib_prop, change_type) %>%
  rename(second_var = variable,
         second_weighted = weighted_contrib,
         second_prop = contrib_prop)

## Merge labels for plotting
dominant_df$rank_type <- "First"
second_df$rank_type   <- "Second"

plot_df <- bind_rows(
  dominant_df %>%
    rename(variable = first_var,
           contrib = first_contrib,
           prop = first_prop),
  second_df %>%
    rename(variable = second_var,
           contrib = second_contrib,
           prop = second_prop)
)


ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(aes(color = variable, size = prop), alpha = 0.8) +
  scale_size_continuous(name = "Proportion of total change", range = c(1, 6)) +
  scale_color_brewer(palette = "Set2", name = "Dominant variable") +
  facet_grid(change_type ~ rank_type) +
  theme_minimal() +
  labs(title = "Environmental drivers of suitability change",
       subtitle = "Separated by direction of change and dominance rank",
       x = "Longitude", y = "Latitude")


### Summary by variable =============
summary_by_var <- plot_df %>%
  group_by(variable, change_type, rank_type) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE),
            n_cells = n()) %>%
  arrange(desc(mean_prop))

ggplot(summary_by_var, aes(x = reorder(variable, -mean_prop),
                           y = mean_prop, fill = rank_type)) +
  geom_col(position = "dodge") +
  facet_wrap(~ change_type) +
  coord_flip() +
  labs(x = "Variable", y = "Mean proportion of total change",
       title = "Average relative contribution of variables by change direction") +
  theme_minimal()


