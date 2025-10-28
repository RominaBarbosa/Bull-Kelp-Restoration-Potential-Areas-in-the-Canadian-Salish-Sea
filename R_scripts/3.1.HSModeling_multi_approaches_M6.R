###==================================================================
### Species Distribution models        SDMs          ################
###                                                  ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Last version: 19-Sep-2025                        ################
###==================================================================
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf)

source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")
variables_selection_path<- "/Volumes/Romina_PSF/PSF/SDM/Variables_selection/M6"

# Load kelp records ============================================================
# kelp_presabs_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_2014_to_2019_filtered_plus_Absences.csv")
# names_col<- colnames(kelp_presabs_df)
# names_col<- names_col[-1]
# kelp_presabs_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_2014_to_2019_filtered_plus_AbsencesFINAL.csv")
# kelp_presabs_df<- kelp_presabs_df[,-1]
# colnames(kelp_presabs_df)<- names_col
# write.csv(kelp_presabs_df, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_2014_to_2019_filtered_plus_AbsencesFINAL.csv")


# Load selected environmental variables ========================================
# Load raster stack of variables resampled at 20 m resoltuion
raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_blob.tif")
names(raster_stack_20m)

### Load terrain variables =====================================================
terrain_path<- ("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables")
tif_files <- list.files(terrain_path, pattern = "\\.tif$", full.names = TRUE)
terrain_vars<-  rast(tif_files[c(12,15)])

### Merge rasters of all selected variables including terrain and NEMO, plus substrate type =========
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


# Mask model predictions by substrate 
substrate<- crop(substrate, terrain_vars)
substrate<- as.factor(substrate)

# Mask model predictions by substrate 
substrate_aligned <- terra::rast(terrain_vars)
substrate_aligned <- terra::resample(substrate, substrate_aligned, method = "near")
unique(values(substrate_aligned))

# Force rounding / integer conversion (as a safeguard):
substrate_aligned <- round(substrate_aligned)
substrate_aligned <- as.factor(substrate_aligned)
names(substrate_aligned)<- "substrate_aligned"

# substrate_aligned[substrate_aligned == 1]<- 1
substrate_aligned[substrate_aligned == 2]<- 1
substrate_aligned[substrate_aligned == 3]<- 2
substrate_aligned[substrate_aligned == 4]<- 2

# substrate<- rast("/Volumes/Romina_PSF/PSF/SDM/SDM_results/substrate_SOG_aligned.tif")
# substrate<- mask(substrate, slope)


### Load bathymetry layer ===========
bathy<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")

### Merge layers ====================
raster_stack_20m_all<- c(raster_stack_20m, terrain_vars, substrate_aligned, bathy)
names(raster_stack_20m_all)
plot(raster_stack_20m_all[[10:11]])


### Extract raster values at kelp locations ====================================
kelp_presabs_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&rocky.csv")
# kelp_presabs_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL2.csv")
colnames(kelp_presabs_df)
kelp_presabs_df<- kelp_presabs_df%>%
  select(c(x,y,kelp, Cluster))#[,c(2,3,11:14)]

head(kelp_presabs_df)


# Convert kelp coords to terra SpatVector points
kelp_points <- vect(kelp_presabs_df, geom = c("x", "y"))
crs(kelp_points) <- "EPSG:3005"

extracted_values <- terra::extract(raster_stack_20m_all, kelp_points)

# Combine extracted values with kelp coordinates 
kelp_data_with_variables <- cbind(kelp_presabs_df, extracted_values[, -1])  # remove ID column from extract
colnames(kelp_data_with_variables)

# Explore number of points with NAs (Presences and absences)
kelp_data_with_variables$substrate_aligned<- as.factor(kelp_data_with_variables$substrate_aligned)
summary(kelp_data_with_variables)


# Select the records on the rocky substrate
kelp_data_with_variables%>%
  group_by(substrate_aligned, kelp)%>%
  count()

# substrate_aligned  kelp     n
# 1 1                     0  2604
# 2 1                     1  1190


# View the result
head(kelp_data_with_variables)
summary(kelp_data_with_variables)

kelp_data_with_variables%>%
  group_by(kelp)%>%
  summarize(n= length(kelp), depth_max= max(coastwide_20m, na.rm = T), 
            depth_min= min(coastwide_20m, na.rm = T))
# After filtering by substrate (it's already filter by depth):
#    kelp     n depth_max depth_min
# 1      0  2604      40.0     -14.9
# 2     1  1190      20.7     -14.3


kelp_data_with_variables<- kelp_data_with_variables%>%
  filter(coastwide_20m >= -10)%>%
  filter(coastwide_20m <= 40)

# After filtering by substrate:
# kelp     n depth_max depth_min
# 1     0  2510      40.0     -9.92 # absences go until 40 m depth
# 2     1  1183      20.7     -9.69 # presences go until 20.8 m depth, 10 presences were above 10m (in the coast)


summary(kelp_data_with_variables)
# write.csv(kelp_data_with_variables, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M6.csv")

# kelp_data_with_variables<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M6.csv")
# kelp_data_with_variables<- kelp_data_with_variables[,-1]

#### PERFORM MODELS ############################################################
# Load necessary libraries
library(dplyr)
library(caret)
library(pROC)
library(mgcv)      # GAM
library(randomForest) # RF
library(gbm)       # BRT
# packages for autocorrelation analysis
library(blockCV)
library(terra)  
library(gstat)
library(sp)

# ---- Example data ----
# predictors: RasterStack or SpatRaster of environmental variables
# coords: matrix/data.frame of coordinates (longitude, latitude)
# species: dataframe with presence/absence, column "y"

# ---- Load packages for autocorrelation analysis ----
library(blockCV)
library(raster)   # or terra if you use SpatRaster
library(gstat)
library(sp)

# ---- Example data ----
# predictors: RasterStack or SpatRaster of environmental variables
# coords: matrix/data.frame of coordinates (longitude, latitude)
# species: dataframe with presence/absence, column "y"

# ---- Option 2: Compute variogram for one predictor manually ----
# Example using the first raster layer
coords<- kelp_data_with_variables[,1:2]
  
predictor<- raster(raster_stack_20m_all[[4]])
var_data <- data.frame(
  value = raster::extract(predictor, coords),  # values at species locations
  x = kelp_data_with_variables[,1],
  y = kelp_data_with_variables[,2]
)

coordinates(var_data) <- ~x+y

# Compute empirical variogram
vgm_model <- variogram(value ~ 1, var_data)

# Fit a theoretical model (e.g., spherical)
fit <- fit.variogram(vgm_model, model = vgm("Sph"))

# Plot the variogram
plot(vgm_model, fit)

# ---- Identify block size ----
# The 'range' of the fitted model is the distance at which autocorrelation
# essentially disappears. Use this as 'theRange' in spatialBlock:
block_size <- fit$range[2]   # or read from the plot manually
print(paste("Suggested block size (m):", block_size))
# "Suggested block size (m): 117053.475444824" - temperature mean summer
# "Suggested block size (m): 21452.212637738"  - current speed
# "Suggested block size (m): 95968.126138523"  - nitrate_summer_minimum


# ================================
# STEP 1: Balance dataset
# ================================
# Keep all presences
df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M6.csv")
df<- df[,-1]
colnames(df)[16]<- "bathymetry_20"

# There are much more absences than presence records
df%>%
  group_by(kelp)%>%
  summarize(n= length(kelp))
# M2:
# kelp     n
# 1     0  5217
# 2     1  926

# M3:
# 1     0  5000
# 2     1   911

# M4:
#    kelp     n
# 1     0  2510
# 2     1  1183

# M6:
# kelp     n
# 1     0  2510
# 2     1  1183

# Downsample absences to match presence count
presences <- df %>% filter(kelp == 1)
absences <- df %>% filter(kelp == 0) %>% sample_n(nrow(presences))

summary(presences$bathymetry_20)       
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -9.6917 -0.9426  1.1702  2.1877  4.2456 20.6680 

# M3        # M4
quantile(presences$bathymetry_20, 0.9)  #26.51765  # 8.848001  # 8.145168
quantile(presences$bathymetry_20, 0.95) #32.6647   # 12.02184  # 11.24392
quantile(presences$bathymetry_20, 0.93) #30.14825  # 10.56405  # 9.554493
quantile(presences$bathymetry_20, 0.99) #37.68245  # 17.72347  # 17.48803

df$kelp<- as.factor(df$kelp)
df.p<- df
levels(df.p$kelp)<- c("Absence", "Presence")

ggplot(df.p, aes(x=as.factor(kelp), y=bathymetry_20, fill=as.factor(kelp)))+
  geom_violin()+
  geom_boxplot(width=0.2)+ 
  scale_fill_manual(values=c("blue", "green"))+
  labs(x="", y= "Depth (m)", fill="Kelp Record")+
  theme_bw()


# Combine balanced data
df_balanced <- bind_rows(presences, absences)
df_balanced%>%
  group_by(kelp)%>%
  summarize(n= length(kelp))

#  M2
# kelp     n
# 1     0  926
# 2     1  926

#  M3
# kelp     n
# 1     0   911
# 2     1   911

# M4
# kelp     n
# 1     0  1183
# 2     1  1183

ggplot(df_balanced, aes(x=as.factor(kelp), y=bathymetry_20, fill=as.factor(kelp)))+
  geom_violin()+
  geom_boxplot(width=0.2)+ 
  scale_fill_manual(values=c("blue", "green"))+
  labs(x="", y= "Depth (m)", fill="Kelp Record")+
  theme_bw()



# ================================
# STEP 2: Train-test split (70/30)
# ================================
library(caret)
set.seed(1234) 
train_index <- createDataPartition(df_balanced$kelp, p = 0.7, list = FALSE)
train <- df_balanced[train_index, ]
test  <- df_balanced[-train_index, ]

train_df<- train
train_df$set<- "train"
train_df$cell<- train_index
test_df<- test
test_df$set<- "test"
test_df$cell<- NA
  
train_test_dataset<- rbind(train_df, test_df)

train_test_dataset %>%
  group_by(set)%>%
  summarize(n= length(set))
# set       n
# 1 test    554
# 2 train  1298

# M3
# 1 test    546
# 2 train  1276

# M4
# 1 test    708
# 2 train  1658
# write.csv(train_test_dataset, "/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_M6.csv")


# train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_M6.csv")
train_test_dataset<- (train_test_dataset[,-1])

colnames(train_test_dataset)
predictors <- setdiff(
  names(train_test_dataset),
  c("kelp", "bathymetry_20","x", "y", "set", "cell","Cluster", "substrate_aligned") 
  )

# Subset the data to only predictors
predictors_data <- train_test_dataset %>% select(all_of(predictors))
colnames(predictors_data)

# Calculate correlation matrix (use pairwise complete obs to handle NAs)
cor_mat <- cor(predictors_data, use = "pairwise.complete.obs")

# Plot 
corrplot::corrplot(cor_mat, method = "color", type = "upper", 
                   tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, diag = FALSE)


# Select predictor variables again based on correlations (>0.8 pearson cor)
library(tidyverse)
predictors <- setdiff(
  predictors,
  c("ammonium_winter_mean",
    "ammonium_summer_minimum",
    "nitrate_summer_minimum"))

# Subset the data to only predictors
predictors_data <- train_test_dataset %>% select(all_of(predictors))
colnames(predictors_data)
cor_mat <- cor(predictors_data, use = "pairwise.complete.obs")

# Original long names
original_names <- colnames(cor_mat)

# Define shorter names in the same order
short_names <- c(  "Mean spring ammonium", 
                   "Mean summer current speed",
                   "Mean summer PAR", 
                   "Mean summer SST",    
                   "Mean summer turbidity",
                   "Mean summer SSS",
                    "Slope" ,  
                   "TPI"
)  # adjust as needed

# Assign shorter names to correlation matrix dimnames
colnames(cor_mat) <- short_names
rownames(cor_mat) <- short_names


pdf("/Volumes/Romina_PSF/PSF/SDM/Plots/correlation_plot_for_SelectedVariablesM6.pdf", width = 10, height = 8)
corrplot::corrplot(cor_mat, method = "color", type = "upper", 
                   tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, diag = FALSE)
dev.off()


### Scale predictors for GLM and GAM
# Save scaling parameters for using to scale variables in predictions 
train<- train_test_dataset %>% 
  filter(set== "train")

test<- train_test_dataset %>%
  filter(set== "test")


# Function to apply same scaling
scale_with_params <- function(df, params) {
  df %>%
    mutate(across(where(is.numeric),
                  ~ (.-params[[paste0(cur_column(), "_mean")]]) /
                    params[[paste0(cur_column(), "_sd")]]))
}

scaling_params <- train %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd), na.rm = TRUE))


train <- train %>% select(kelp, all_of(predictors))
test <- test %>% select(kelp, all_of(predictors))
train$kelp<- as.factor(train$kelp)
test$kelp <- as.factor(test$kelp)


train_scaled <- scale_with_params(train, scaling_params)
test_scaled  <- scale_with_params(test,  scaling_params)

train %>%
  pivot_longer(
    cols = -kelp,               # all columns except kelp
    names_to = "variable", 
    values_to = "value"
  ) %>%
  ggplot(aes(x = factor(kelp), y = value)) +
  geom_violin()+
  geom_boxplot(width=0.4) +
  facet_wrap(~ variable, scales = "free_y", ncol=3) +
  labs(x = "Kelp presence (0 = absent, 1 = present)",
       y = "Value") +
  theme_bw()


# ==============================================================================
# STEP 3: Fit models
# ==============================================================================
library(dplyr)
library(cluster)
library(mgcv)
library(randomForest)
library(dismo)
library(ggplot2)
library(cluster)
library(factoextra)

explore_k_clusters <- function(data, k_max = 15, nstart = 25, seed = 42) {
  set.seed(seed)
  
  # ----- 1. Elbow Method -----
  wss <- numeric(k_max)
  for (k in 1:k_max) {
    km <- kmeans(data, centers = k, nstart = nstart)
    wss[k] <- km$tot.withinss
  }
  
  plot(1:k_max, wss, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters K",
       ylab = "Total within-clusters sum of squares",
       main = "Elbow Method")
  
  # ----- 2. Silhouette Method -----
  sil_width <- numeric(k_max)
  for (k in 2:k_max) {
    km <- kmeans(data, centers = k, nstart = nstart)
    ss <- silhouette(km$cluster, dist(data))
    sil_width[k] <- mean(ss[, 3])
  }
  
  plot(2:k_max, sil_width[2:k_max], type = "b", pch = 19,
       xlab = "Number of clusters K",
       ylab = "Average silhouette width",
       main = "Silhouette Method")
  
  # ----- 3. Gap Statistic -----
  gap_stat <- clusGap(data, FUN = kmeans, nstart = nstart, K.max = k_max, B = 50)
  fviz_gap_stat(gap_stat)
  
  return(list(wss = wss, sil_width = sil_width, gap_stat = gap_stat))
}

results <- explore_k_clusters(train[,2:ncol(train)], k_max = 7, nstart = 25, seed = 42)


K <-  5  # <-- set number of clusters
env_pres <- train_scaled 

set.seed(123)
cl <- kmeans(env_pres, centers = K, nstart = 25)

# assign cluster IDs back to full dataset (presences + absences)
train$cluster <- NA
train$cluster<- cl$cluster


# Visualize in biplot:
# Select only numeric columns for PCA
num_vars <- names(train)[sapply(train, is.numeric)]
pca_data <- scale(train[, num_vars[1:8]])

# Run PCA
pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_res$x)  # PCA coordinates
pca_scores$cluster <- as.factor(train$cluster)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA Biplot of Clusters",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Adding variables as arrows
library(factoextra)
fviz_pca_biplot(pca_res,
                label = "var",       # show variable names
                habillage = train$cluster, # color by cluster
                addEllipses = TRUE,  # optional: confidence ellipse
                ellipse.level = 0.95,
                palette = "Set1")


train.xy<- train_test_dataset %>% 
  filter(set== "train")

train.xy$Cluster<- train$cluster
str(train.xy)
ggplot(train.xy, aes(y= y, x= x, color=as.factor(Cluster)))+
  geom_point()

# -----------------------------
# 4. Compute cluster weights  
# -----------------------------
freq <- table(train[which(train$kelp==1), "cluster"])
inv_freq <- 1 / as.numeric(freq)
names(inv_freq) <- names(freq)

# Now lookup by cluster number as character
pres_weights <- inv_freq[as.character(train$cluster[train$kelp==1])]
# pres_weights <- inv_freq[train[which(train$kelp==1), "cluster"]]

# normalize to mean = 1
pres_weights <- pres_weights * (length(pres_weights) / sum(pres_weights))

# add weights to dataset
train_weight<- train
train_weight$weight <- 1
train_weight$weight[train_weight$kelp == 1] <- pres_weights

train_scaled_weight<- train_scaled
train_scaled_weight$weight <- 1
train_scaled_weight$weight[train_scaled_weight$kelp == 1] <- pres_weights

# -----------------------------
# 5. Fit models with weights
# -----------------------------
## GLM
# ==============================================================================
library(car)        # vif
library(broom)      # tidy model output
library(visreg)     # partial effect plots


# Specify the predictors for quadratic terms 
quad_vars <- colnames(train_scaled_weight)[2:9]

# Build the formula with quadratic terms using I(x^2)
#    Use I(x^2) rather than poly() so prediction on new data is straightforward.
terms_quad <- paste0(quad_vars, " + I(", quad_vars, "^2)")
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad), collapse = " + ")))

# Fit the GLM (binomial)
glm_mod_s <- glm(glm_formula,
                 data = train_scaled_weight,
                 family = binomial(link = "logit"),
                 weights = weight)


library(patchwork)
vars <- attr(terms(glm_mod_s), "term.labels")
vars_clean <- vars[!grepl("^I\\(", vars)]

plots <- lapply(vars_clean, function(v) {
  visreg(glm_mod_s, v, scale = "response", gg = TRUE) + 
    ggtitle(paste("Effect of", v))
})

# Combine into one grid
wrap_plots(plots)



# Basic summary & multicollinearity check
summary(glm_mod_s)
vif_vals <- vif(glm_mod_s)     # check for high VIFs (>5 or >10)
print(vif_vals)

# print(vif_vals)
# nitrate_winter_minimum   I(nitrate_winter_minimum^2)       ammonium_winter_minimum  I(ammonium_winter_minimum^2)          ammonium_spring_mean 
# 4.660565                      3.528941                     15.226988                      4.649808                      5.214622 
# I(ammonium_spring_mean^2)       temperature_summer_mean  I(temperature_summer_mean^2)               PAR_summer_mean          I(PAR_summer_mean^2) 
# 3.050223                     17.545143                      3.467860                      5.620510                      1.829542 
# salinity_summer_SD       I(salinity_summer_SD^2)      currentSpeed_summer_mean I(currentSpeed_summer_mean^2)         turbidity_summer_mean 
# 11.353267                      6.549615                      4.488431                      2.774229                     12.578146 
# I(turbidity_summer_mean^2)            ammonium_spring_SD                easterness_3x3               northerness_3x3                     slope_5x5 
# 2.796814                      6.934039                      1.131495                      1.233387                      1.166509 
# TPI_3x3 
# 1.046070 


# ammonium_spring_mean     I(ammonium_spring_mean^2)      currentSpeed_summer_mean I(currentSpeed_summer_mean^2)        nitrate_summer_minimum 
# 3.332036                      1.885950                      4.519733                      2.755218                     25.735134 
# I(nitrate_summer_minimum^2)           nitrate_winter_mean      I(nitrate_winter_mean^2)               PAR_summer_mean          I(PAR_summer_mean^2) 
# 5.608135                      5.026038                      3.274557                      3.562834                      1.847794 
# temperature_summer_mean  I(temperature_summer_mean^2)         turbidity_summer_mean    I(turbidity_summer_mean^2)                easterness_3x3 
# 19.174671                      5.549064                      4.801422                      2.192409                      1.333448 
# I(easterness_3x3^2)               northerness_3x3          I(northerness_3x3^2)                     slope_5x5                I(slope_5x5^2) 
# 2.867631                      1.225061                      2.810393                      2.823536                      2.656227 
# TPI_3x3                  I(TPI_3x3^2) 
# 1.106174                      1.403577 


# M6: 
# ammonium_spring_mean     I(ammonium_spring_mean^2)      currentSpeed_summer_mean I(currentSpeed_summer_mean^2)               PAR_summer_mean 
# 2.597414                      1.859930                      3.825198                      2.313262                      4.396807 
# I(PAR_summer_mean^2)       temperature_summer_mean  I(temperature_summer_mean^2)         turbidity_summer_mean    I(turbidity_summer_mean^2) 
# 1.690142                     10.026678                      2.317880                      6.412240                      2.937120 
# salinity_summer_mean     I(salinity_summer_mean^2)                     slope_5x5                I(slope_5x5^2)                       TPI_3x3 
# 19.017519                      8.635751                      2.830933                      2.435696                      1.416419 
# I(TPI_3x3^2) 
# 1.385584 


## GAM
# ==============================================================================
gam_formula <- as.formula(paste("kelp ~", paste0("s(", predictors, ")", collapse = " + ")))
gam_mod_s <- gam(gam_formula,
                                  data = train_scaled_weight,
                                  family = binomial,
                                  weights = weight)
                 
                 
visreg(gam_mod_s, "salinity_summer_mean", scale = "response")

## Random Forest
# ==============================================================================
train_weight$kelp <- as.factor(train_weight$kelp)
rf_mod <-randomForest(x = train_weight[, predictors], y = train_weight$kelp, ntree = 500)
rf_mod_s <- randomForest(as.formula(paste("kelp ~", paste(predictors, collapse = " + "))),
                         data = train_weight,
                         ntree = 500,
                         importance = TRUE,
                         sampsize = nrow(train_weight),
                         case.weights = train_weight$weight)

# Call:
#   randomForest(formula = as.formula(paste("kelp ~", paste(predictors,      collapse = " + "))), data = train, ntree = 500, importance = TRUE,      sampsize = nrow(train), case.weights = train$weight) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 13.87%
# Confusion matrix:
#   0   1 class.error
# 0 683 146   0.1761158
# 1  84 745   0.1013269

## BRT
# ==============================================================================
train_brt <- train_weight
colnames(train_brt)
train_brt$kelp <- as.numeric(as.character(train_brt$kelp))
brt_mod_s <- dismo::gbm.step(data = train_brt,
                             gbm.x = which(names(train_brt) %in% predictors),
                             gbm.y = which(names(train_brt) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5,
                             site.weights = train_brt$weight)


# M3:
# mean total deviance = 1.386 
# mean residual deviance = 0.43 
# 
# estimated cv deviance = 0.773 ; se = 0.024 
# 
# training data correlation = 0.883 
# cv correlation =  0.732 ; se = 0.012 
# 
# training data AUC score = 0.981 
# cv AUC score = 0.906 ; se = 0.005 
# 
# elapsed time -  0.33 minutes 


# M4
# mean total deviance = 1.386 
# mean residual deviance = 0.432 
# 
# estimated cv deviance = 0.713 ; se = 0.029 
# 
# training data correlation = 0.877 
# cv correlation =  0.747 ; se = 0.013 
# 
# training data AUC score = 0.979 
# cv AUC score = 0.922 ; se = 0.007 
# 
# elapsed time -  0.45 minutes 


# M5 weighted presences
# mean total deviance = 1.386 
# mean residual deviance = 0.475 
# 
# estimated cv deviance = 0.805 ; se = 0.046 
# 
# training data correlation = 0.853 
# cv correlation =  0.717 ; se = 0.014 
# 
# training data AUC score = 0.972 
# cv AUC score = 0.904 ; se = 0.007 
# 
# elapsed time -  0.46 minutes 


# M5 weighted presences and salinity instead of winter nitrate
# fitting final gbm model with a fixed number of 2250 trees for kelp
# 
# mean total deviance = 1.386 
# mean residual deviance = 0.544 
# 
# estimated cv deviance = 0.849 ; se = 0.039 
# 
# training data correlation = 0.828 
# cv correlation =  0.703 ; se = 0.012 
# 
# training data AUC score = 0.963 
# cv AUC score = 0.9 ; se = 0.007 
# 
# elapsed time -  0.44 minutes 



# ================================
# STEP 4: Predictions & Evaluation
# ================================
#-----------------------------
# Plot all response curves in facets
source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/plot_response_curves.R")

library(pdp)
library(visreg)
library(gratia)  # for GAMs if needed

all_curves <- do.call(rbind, lapply(predictors, function(v) {
  grid <- create_var_grid(v, train)   # one-variable-at-a-time grid
  grid$kelp <- NULL  # remove response variable
  rbind(
    get_curve(model=glm_mod_s, v, grid, "GLM", scaling_params),
    get_curve(gam_mod_s,      v, grid, "GAM", scaling_params),
    get_curve(rf_mod_s,       v, grid, "randomForest"),
    get_curve(brt_mod_s,      v, grid, "BRT")
  )
}))



ggplot(all_curves, aes(x = x, y = fit, color = model)) +
  geom_line(size = 1.2) +
  facet_wrap(~var, scales = "free_x", ncol=3) +
  theme_bw() +
  labs(
    y = "Predicted response",
    x = NULL,
    color = "Model",
    title = "Predicted response curves for all variables across models"
  )

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/ResponseCurves_allVariables_FINAL_M4.png", width = 20, height = 17, dpi= 300, units="cm")


#-----------------------------
# Collect results
results <- bind_rows(
  get_metrics_optimized(glm_mod_s, train_scaled, "glm"),
  get_metrics_optimized(gam_mod_s, train_scaled, "gam"),
  get_metrics_optimized(rf_mod_s, train, "rf"),
  get_metrics_optimized(brt_mod_s, train, "brt")
)

# Model Threshold   AUC Sensitivity Specificity   TSS
# 1 glm       0.632 0.886       0.812       0.828 0.639
# 2 gam       0.599 0.911       0.831       0.864 0.694
# 3 rf        0.495 1           1           1     1    
# 4 brt       0.488 0.981       0.955       0.918 0.873

# Model Threshold   AUC Sensitivity Specificity   TSS
# 1 glm       0.601 0.901       0.840       0.836 0.676
# 2 gam       0.585 0.920       0.852       0.855 0.707
# 3 rf        0.494 1           1           1     1    
# 4 brt       0.556 0.979       0.928       0.930 0.858

# M6
# Model Threshold   AUC Sensitivity Specificity   TSS
# <chr>     <dbl> <dbl>       <dbl>       <dbl> <dbl>
# 1 glm       0.550 0.863       0.818       0.755 0.573
# 2 gam       0.530 0.893       0.855       0.800 0.655
# 3 rf        0.498 1           1           1     1    
# 4 brt       0.460 0.963       0.937       0.864 0.801



# Pivot the results for plotting
results_long <- results %>%
  pivot_longer(cols = c("AUC", "Sensitivity", "Specificity"),
               names_to = "Metric", values_to = "Value")
head(results_long)


# Plot
ggplot(results_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set2") +
  ylim(c(0,1))+
  labs(title = "Model Performance Comparison (testing data)",
       y = "Score",
       x = "Model") +
  theme_bw()




# Collect ROC data for each model
roc_df <- bind_rows(
  get_roc_data(glm_mod_s, test_scaled, "GLM"),
  get_roc_data(gam_mod_s, test_scaled, "GAM"),
  get_roc_data(rf_mod_s,  test, "RF"),
  get_roc_data(brt_mod_s,  test, "BRT")
)

# Plot ROC curves
ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "grey") +
  labs(title = "ROC Curves for Model Comparison",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()




# =============================================================================
# STEP 5: Variable Selection based on Variables Importance
# =============================================================================
# each model family measures "importance" differently, so we’ll need to standardize how we extract and plot them:
# GLM & GAM → We can use caret::varImp() which ranks predictors based on their absolute t-statistics (or z-values).
# RF → Use randomForest::importance() or ranger::importance().
# BRT → Use summary() from gbm which returns relative influence of predictors.
# Then we can combine them into one table for side-by-side comparison.


# Use functions to Get variable importance and plot the results
models <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")

p_vars<- plot_var_importance_v2(models, types,  top_n = 16, stacked = TRUE)
p_vars2<- plot_var_importance_v2(models, types,  top_n = 16)


# Explanation for the ranking of variables strategy:
# Standard deviation of the standardized importance (Importance_std) across models for each variable.
# Importance_std = Importance / max(Importance) * 100 rescales each model’s variable importance to a 0–100% scale.
# Then you summarize min, mean, max importance across models.
# Importance_SD → variability of importance across models.
# Weighted_Score → high if a variable is both important and consistent.
# Keep → safeguards against variables that are consistently low.
# Sorting by Weighted_Score gives the final ranking.

# Ranking variables and exclude less important ones using a threshold of 10 weighted score
res <- rank_variables(models, types, avg_threshold = NULL)
res$summary_table 
plot_variable_ranking_v2(res$summary_table, threshold = 20)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_for_SelecBasedOnImp_M4.pdf", width = 15, height = 13, dpi= 300, units="cm")


plot_importance_profiles_v2(res$summary_table, threshold=18)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/VariablesImportance_profile_VarsSelection_FINAL6.png", width = 18, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Plots/VariablesImportance_profile_VarsSelection_FINAL6.pdf", width = 18, height = 12, dpi= 300, units="cm")





# ==============================================================================
# STEP 7: SDMs with selected variables 
# ==============================================================================
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M6_weightedPres")
thresh= 18

vars_selected <- rank_variables(models, types, avg_threshold = 18)
vars_selected <- vars_selected$summary_table$Variable

# vars_selected
# [1] "temperature_summer_mean"  "slope_5x5"                "turbidity_summer_mean"    "ammonium_spring_SD"      
# [5] "PAR_summer_mean"          "currentSpeed_summer_mean" "salinity_summer_SD"      

# 4 Sep 2025
# [1] "slope_5x5"                "nitrate_summer_minimum"   "turbidity_summer_mean"    "nitrate_winter_mean"      "temperature_summer_mean" 
# [6] "ammonium_spring_mean"     "currentSpeed_summer_mean" "PAR_summer_mean"         

#  M6
# [1] "temperature_summer_mean"  "turbidity_summer_mean"    "slope_5x5"                "salinity_summer_mean"     "ammonium_spring_mean"    
# [6] "PAR_summer_mean"          "currentSpeed_summer_mean"


# Refit models using only selected vars
train_sel <- train_weight %>% select(all_of(c("kelp", vars_selected)))
test_sel <- test %>% select(all_of(c("kelp", vars_selected)))

train_sel$kelp<- as.factor(train_sel$kelp)

scaling_params_2 <- train_sel %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, sd = sd), na.rm = TRUE))

train_sel_scaled <- scale_with_params(train_sel, scaling_params_2)
test_sel_scaled  <- scale_with_params(test_sel,  scaling_params_2)


K <-  5  # <-- set number of clusters
env_pres <- train_sel_scaled 
set.seed(42)
cl <- kmeans(env_pres, centers = K, nstart = 25)

# assign cluster IDs back to full dataset (presences + absences)
train_sel_scaled$cluster <- NA
train_sel_scaled$cluster<- cl$cluster


# Visualize in biplot:
# Select only numeric columns for PCA
num_vars <- names(train_sel_scaled)[sapply(train_sel_scaled, is.numeric)]
pca_data <- scale(train_sel_scaled[, num_vars[1:7]])

# Run PCA
pca_res <- prcomp(pca_data, center = TRUE, scale. = F)
pca_scores <- as.data.frame(pca_res$x)  # PCA coordinates
pca_scores$cluster <- as.factor(train_sel_scaled$cluster)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA Biplot of Clusters",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Adding variables as arrows
factoextra::fviz_pca_biplot(pca_res,
                label = "var",       # show variable names
                habillage = train_sel_scaled$cluster, # color by cluster
                addEllipses = TRUE,  # optional: confidence ellipse
                ellipse.level = 0.95,
                palette = "Set1")


train.xy<- train_test_dataset %>% 
  filter(set== "train")

train.xy$Cluster<- train_sel_scaled$cluster
str(train.xy)
ggplot(train.xy, aes(y= y, x= x, color=as.factor(Cluster)))+
  geom_point()

# Compute cluster weights  
# raw weights for presences
# Frequency table already has cluster IDs as names
freq <- table(train_sel_scaled$cluster[train_sel_scaled$kelp == 1])

# Inverse frequency, preserving names
inv_freq <- 1 / freq

# Extract weights for presence rows
pres_weights <- inv_freq[as.character(train_sel_scaled$cluster[train_sel_scaled$kelp == 1])]

# normalize to mean = 1
pres_weights <- pres_weights * (length(pres_weights) / sum(pres_weights))

# add weights to dataset
train_sel_scaled$weight <- 1
train_sel_scaled$weight[train_sel_scaled$kelp == 1] <- pres_weights

train_sel$weight <- 1
train_sel$weight[train_sel$kelp == 1] <- pres_weights


# ==============================================================================
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M6_weightedPres")
# write.csv(train_sel, "train_selected_table_FINALMODELS.csv")
# write.csv(train_sel_scaled, "train_selected_scaled_table_FINALMODELS.csv")
# write.csv(test_sel, "test_selected_table_FINALMODELS.csv")
# write.csv(test_sel_scaled, "test_selected_scaled_table_FINALMODELS.csv")


## GLM
quad_vars_sel <- quad_vars[quad_vars %in% vars_selected]
terms_quad_sel <- paste0(quad_vars_sel, " + I(", quad_vars_sel, "^2)")
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel), collapse = " + ")))
glm_mod_se <- glm(glm_formula, 
                 data = train_sel_scaled, 
                 family = binomial,
                 weights = weight)

saveRDS(glm_mod_se, "glm_mod_s.rds")

## GAM
gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))
gam_mod_se <- gam(gam_formula, 
                 data = train_sel_scaled, 
                 family = binomial,
                 weights = weight)

saveRDS(gam_mod_se, "gam_mod_s.rds")

## Random Forest
train_sel$kelp <- as.factor(train_sel$kelp)
rf_mod_se <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                         data = train_sel,
                         ntree = 500,
                         importance = TRUE,
                         sampsize = nrow(train_sel),
                         case.weights = train_sel$weight)

saveRDS(rf_mod_se, "rf_mod_s.rds")

## BRT
train_brt_sel <- train_sel# %>% select(-env_cluster)
colnames(train_brt_sel)
train_brt_sel$kelp <- as.numeric(as.character(train_brt_sel$kelp))
brt_mod_se <- dismo::gbm.step(data = train_brt_sel, 
                             gbm.x = which(names(train_brt_sel) %in% vars_selected),
                             gbm.y = which(names(train_brt_sel) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5,
                             site.weights = train_brt_sel$weight)

saveRDS(brt_mod_se, "brt_mod_s.rds")

# fitting final gbm model with a fixed number of 2750 trees for kelp
# 
# mean total deviance = 1.386 
# mean residual deviance = 0.211 
# 
# estimated cv deviance = 0.397 ; se = 0.033 
# 
# training data correlation = 0.947 
# cv correlation =  0.875 ; se = 0.012 
# 
# training data AUC score = 0.995 
# cv AUC score = 0.975 ; se = 0.004 
# 
# elapsed time -  0.43 minutes 


# Model without SD_salinity_summer
# fitting final gbm model with a fixed number of 1550 trees for kelp
# 
# mean total deviance = 1.386 
# mean residual deviance = 0.562 
# 
# estimated cv deviance = 0.782 ; se = 0.027 
# 
# training data correlation = 0.826 
# cv correlation =  0.716 ; se = 0.015 
# 
# training data AUC score = 0.959 
# cv AUC score = 0.904 ; se = 0.007 
# 
# elapsed time -  0.2 minutes 

# M4 (only hard substrate)
# mean total deviance = 1.386 
# mean residual deviance = 0.43 
# 
# estimated cv deviance = 0.718 ; se = 0.03 
# 
# training data correlation = 0.878 
# cv correlation =  0.748 ; se = 0.014 
# 
# training data AUC score = 0.979 
# cv AUC score = 0.921 ; se = 0.006 
# 
# elapsed time -  0.44 minutes 


# M5 (only hard substrate + weighted presences)
# mean total deviance = 1.386 
# mean residual deviance = 0.504 
# 
# estimated cv deviance = 0.854 ; se = 0.046 
# 
# training data correlation = 0.837 
# cv correlation =  0.695 ; se = 0.014 
# 
# training data AUC score = 0.967 
# cv AUC score = 0.894 ; se = 0.008 
# 
# elapsed time -  0.39 minutes 

#  M6
# fitting final gbm model with a fixed number of 1600 trees for kelp
# 
# mean total deviance = 1.386 
# mean residual deviance = 0.675 
# 
# estimated cv deviance = 0.971 ; se = 0.03 
# 
# training data correlation = 0.779 
# cv correlation =  0.664 ; se = 0.014 
# 
# training data AUC score = 0.94 
# cv AUC score = 0.876 ; se = 0.009 
# 
# elapsed time -  0.33 minutes 




# save.image("3.1.HSModeling_multi_approaches_M6.RData")


