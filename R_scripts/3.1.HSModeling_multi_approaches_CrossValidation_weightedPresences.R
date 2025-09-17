###==================================================================
### Species Distribution models        SDMs          ################
###                                                  ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Last version: 11-Aug-2025                        ################
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
variables_selection_path<- "/Volumes/Romina_PSF/PSF/SDM/Variables_selection"

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
raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_FINAL2.tif")
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

# Combine extracted values with kelp coordinates (optional)
kelp_data_with_variables <- cbind(kelp_presabs_df, extracted_values[, -1])  # remove ID column from extract
colnames(kelp_data_with_variables)

# Explore number of points with NAs (Presences and absences)
kelp_data_with_variables$substrate_aligned<- as.factor(kelp_data_with_variables$substrate_aligned)
summary(kelp_data_with_variables)
substrate_NAs<- kelp_data_with_variables[which(is.na(kelp_data_with_variables$substrate_aligned)),]
  
substrate_NAs_points <- vect(substrate_NAs, geom = c("x", "y"))
crs(substrate_NAs_points) <- "EPSG:3005"
plot(substrate_NAs_points)

substrate_NAs_points$kelp<- as.factor(substrate_NAs_points$kelp)
as.data.frame(substrate_NAs_points)%>%group_by(kelp)%>%summarize(n= length(kelp))


# Select the records on the rocky substrate
kelp_data_with_variables%>%
  group_by(substrate_aligned, kelp)%>%
  count()
# substrate_aligned     n
# 1 1                  5921
# 2 2                     5  # only 5 records in soft substrate were excluded
# 3 NA                  217  # 217 records with NAs of substrate were excluded

# substrate_aligned  kelp     n
# 1 1                     0  2604
# 2 1                     1  1190


kelp_data_with_variables<- kelp_data_with_variables%>%
  filter(substrate_aligned == 1) # keep only records in hard substrate area



# View the result
head(kelp_data_with_variables)
summary(kelp_data_with_variables)

kelp_data_with_variables%>%
  group_by(kelp)%>%
  summarize(n= length(kelp), depth_max= max(coastwide_20m, na.rm = T), 
            depth_min= min(coastwide_20m, na.rm = T))
#    kelp     n    depth_max depth_min
# 1     0  8761      105.     -104. 
# 2     1  1287      83.0     -45.9

# After filtering by substrate (it's already filter by depth):
#    kelp     n depth_max depth_min
# 1      0  2604      40.0     -14.9
# 2     1  1190      20.7     -14.3


kelp_data_with_variables<- kelp_data_with_variables%>%
  filter(coastwide_20m >= -10)%>%
  filter(coastwide_20m <= 40)

# kelp     n depth_max depth_min
# 1     0  5217      40.0     -9.99
# 2     1   926      39.6     -9.81

# After filtering by substrate:
# kelp     n depth_max depth_min
# 1     0  2510      40.0     -9.92 # absences go until 40 m depth
# 2     1  1183      20.7     -9.69 # presences go until 20.8 m depth, 10 presences were above 10m (in the coast)


summary(kelp_data_with_variables)
# write.csv(kelp_data_with_variables, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M4.csv")


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

# ---- Option 1: Interactive exploration using blockCV ----
# Opens a Shiny app to explore spatial autocorrelation range
rangeExplorer(rasterLayer = predictors,    # your environmental rasters
              speciesData = species,       # your species records
              species = "y")               # response variable column

# In the app:
# 1. Select a predictor
# 2. Look at the variogram
# 3. Identify the distance where semivariance flattens (range)
# 4. Use that distance as 'theRange' in spatialBlock

# ---- Option 2: Compute variogram for one predictor manually ----
# Example using the first raster layer
var_data <- data.frame(
  value = raster::extract(predictors[[1]], coords),  # values at species locations
  x = coords[,1],
  y = coords[,2]
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



# ---- Identify block size ----
# The 'range' of the fitted model is the distance at which autocorrelation
# essentially disappears. Use this as 'theRange' in spatialBlock:
block_size <- fit$range[2]   # or read from the plot manually
print(paste("Suggested block size (m):", block_size))


# ---- Load packages ----
library(blockCV)         
library(caret)           
library(randomForest)    
library(mgcv)            
library(gbm)             
library(pROC)            
library(PresenceAbsence) 
library(dplyr)
library(tidyr)

# ---- Example data (replace with your own) ----
# data: dataframe with predictors (X1, X2, ...) and response (binary 0/1 in "y")
# coords: matrix/data.frame of coordinates (longitude, latitude)
# coords <- data[,c("lon","lat")]

set.seed(123)

# ---- Create spatial blocks ----
sb <- spatialBlock(speciesData = data,
                   species = "y",
                   theRange = 50000,   # adjust distance (m)
                   k = 5,              # number of spatial folds
                   selection = "systematic",
                   iteration = 100)

folds <- sb$foldID

# ---- Function to calculate metrics ----
get_metrics <- function(obs, pred){
  thresh <- optimal.thresholds(data.frame(ID=1:length(obs), obs=obs, pred=pred))$pred
  cm <- cmx(data.frame(ID=1:length(obs), obs=obs, pred=pred), threshold=thresh)
  
  sens <- cm[2,2] / sum(cm[2,])         
  spec <- cm[1,1] / sum(cm[1,])         
  tss  <- sens + spec - 1
  roc  <- auc(obs, pred)                
  
  return(c(TSS=tss, ROC=roc, Sensitivity=sens, Specificity=spec))
}

# ---- Models ----
models <- c("glm","gam","rf","brt")
all_results <- list()

for (m in models){
  metrics <- matrix(NA, nrow=max(folds), ncol=4)
  colnames(metrics) <- c("TSS","ROC","Sensitivity","Specificity")
  
  for (i in 1:max(folds)){
    train_data <- data[folds != i, ]
    test_data  <- data[folds == i, ]
    
    if(m == "glm"){
      mod <- glm(y ~ ., data=train_data, family=binomial)
      pred <- predict(mod, newdata=test_data, type="response")
    }
    if(m == "gam"){
      mod <- gam(y ~ s(X1) + s(X2) + s(X3), data=train_data, family=binomial)
      pred <- predict(mod, newdata=test_data, type="response")
    }
    if(m == "rf"){
      mod <- randomForest(as.factor(y) ~ ., data=train_data, ntree=500)
      pred <- predict(mod, newdata=test_data, type="prob")[,2]
    }
    if(m == "brt"){
      mod <- gbm(y ~ ., data=train_data,
                 distribution="bernoulli",
                 n.trees=2000,
                 interaction.depth=3,
                 shrinkage=0.01,
                 n.minobsinnode=10,
                 verbose=FALSE)
      pred <- predict(mod, newdata=test_data, type="response", n.trees=2000)
    }
    
    metrics[i,] <- get_metrics(test_data$y, pred)
  }
  
  # Collect per-model results
  all_results[[m]] <- data.frame(
    Model = m,
    Metric = colnames(metrics),
    Mean = colMeans(metrics, na.rm=TRUE),
    SD = apply(metrics, 2, sd, na.rm=TRUE)
  )
}

# ---- Combine and pivot to wide table ----
results_table <- bind_rows(all_results) %>%
  pivot_wider(
    names_from = Metric,
    values_from = c(Mean, SD),
    names_glue = "{Metric}_{.value}"
  )

print(results_table)
























set.seed(123) # for reproducibility

# ================================
# STEP 1: Balance dataset
# ================================
# Keep all presences

df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M4.csv")
df<- df[,-1]
colnames(df)[17]<- "bathymetry_20"

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


# Downsample absences to match presence count
presences <- df %>% filter(kelp == 1)
absences <- df %>% filter(kelp == 0) %>% sample_n(nrow(presences))

summary(presences$bathymetry_20)       
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


# ================================
# STEP 2: Train-test split (70/30)
# ================================
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
# write.csv(train_test_dataset, "/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_M4.csv")


train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_M4.csv")
train_test_dataset<- (train_test_dataset[,-1])

colnames(train_test_dataset)
predictors <- setdiff(
  names(train_test_dataset),
  c("kelp", "bathymetry_20","x", "y", "set", "cell","Cluster", "substrate_aligned"#"cell_id", "env_cluster", "easterness_3x3",           "northerness_3x3" 
    ) # "ammonium_spring_mean.1", "bathymetry_20", "summer_hours_above_threshold_18", "summer_cumulated_degrees_18",
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

# predictors <- setdiff(
#   names(train_test_dataset),
#   c("kelp", "cell_id", "env_cluster", "x", "y", 
#      "bathymetry_20",# "ammonium_spring_mean.1","summer_hours_above_threshold_18", "summer_cumulated_degrees_18",
#     "nitrate_summer_minimum",  # correlated with temperature summer mean
#     # "salinity_summer_mean",#"ammonium_winter_mean", # correlated with salinity summer mean
#     "ammonium_summer_SD", #"ammonium_summer_minimum", # correlated with ammonium_summer_minimum and temperature summer mean
#     "PAR_summer_maximum",
#     # "currentDirection_spring_min", #"currentDirection_summer_modal", 
#     "set", "cell")   # correlated with PAR summer mean
# )

# predictors <- setdiff(
#   names(train_test_dataset),
#   c("kelp", "cell_id", "env_cluster", "x", "y", 
#     "bathymetry_20",# "ammonium_spring_mean.1","summer_hours_above_threshold_18", "summer_cumulated_degrees_18",
#     "ammonium_summer_minimum",  # correlated with temperature summer mean
#     # "salinity_summer_mean",#"ammonium_winter_mean", # correlated with salinity summer mean
#     "ammonium_winter_mean", # correlated with nitrate winter mean and temperature summer mean
#     # "currentDirection_spring_min", #"currentDirection_summer_modal", 
#     "set", "cell", "easterness_3x3",   "northerness_3x3"  )
# )


predictors <- setdiff(
  predictors,
  c("ammonium_winter_mean",
    "ammonium_summer_minimum"))

# Subset the data to only predictors
predictors_data <- train_test_dataset %>% select(all_of(predictors))
colnames(predictors_data)
cor_mat <- cor(predictors_data, use = "pairwise.complete.obs")

# Original long names
original_names <- colnames(cor_mat)

# Define shorter names in the same order
# short_names <- c( "nitrate_winter_minimum" ,  "currentSpeed_summer_mean", "salinity_summer_SD",       "ammonium_winter_minimum" , "ammonium_spring_SD",      
#                   "temperature_summer_mean",  "PAR_summer_mean" ,         "ammonium_spring_mean",     "turbidity_summer_mean" ,   "easterness_3x3",          
#                   "northerness_3x3",          "slope_5x5" ,               "TPI_3x3"               
# )  # adjust as needed

# Assign shorter names to correlation matrix dimnames
# colnames(cor_mat) <- short_names
# rownames(cor_mat) <- short_names


# pdf("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/correlation_plot_for_SelecBasedOnImp_FINAL2.pdf", width = 10, height = 8)
corrplot::corrplot(cor_mat, method = "color", type = "upper", 
                   tl.col = "black", tl.srt = 45, addCoef.col = "black", number.cex = 0.7, diag = FALSE)
# dev.off()


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

# -----------------------------
# 1. Select variables
# -----------------------------
# vars_selected <- c("turbidity_summer_mean", "temperature_summer_mean", "nitrate_winter_mean",
#                    "slope_5x5", "currentSpeed_summer_mean",
#                    "ammonium_spring_mean", "PAR_summer_mean")
# 
# train_sel <- train %>% select(all_of(c("kelp", vars_selected)))
# test_sel  <- test  %>% select(all_of(c("kelp", vars_selected)))
# 
# train_sel$kelp <- as.factor(train_sel$kelp)
# test_sel$kelp  <- as.factor(test_sel$kelp)

# -----------------------------
# 2. Scale predictors
# -----------------------------
# scaling_params_2 <- train_sel %>%
#   summarise(across(where(is.numeric),
#                    list(mean = mean, sd = sd), na.rm = TRUE))
# 
# train_sel_scaled <- scale_with_params(train_sel, scaling_params_2)
# test_sel_scaled  <- scale_with_params(test_sel,  scaling_params_2)

# -----------------------------
# 3. Build clusters on presences only
# -----------------------------
# Explore best number of cluster
# wss <- numeric(15)  # up to 15 clusters
# for (k in 1:15) {
#   set.seed(42)
#   km <- kmeans(env_pres, centers = k, nstart = 25)
#   wss[k] <- km$tot.withinss
# }
# 
# plot(1:15, wss, type="b", pch=19, frame=FALSE,
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")
# 
# sil_width <- numeric(15)
# for (k in 2:15) {  # silhouette needs at least 2 clusters
#   set.seed(42)
#   km <- kmeans(env_pres, centers = k, nstart = 25)
#   ss <- silhouette(km$cluster, dist(env_pres))
#   sil_width[k] <- mean(ss[, 3])
# }
# 
# plot(2:15, sil_width[2:15], type="b", pch=19,
#      xlab="Number of clusters K",
#      ylab="Average silhouette width")
# 
# library(factoextra)
# 
# set.seed(42)
# gap_stat <- clusGap(env_pres, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
# fviz_gap_stat(gap_stat)
# 
# # Combine methods
# library(Rtsne)
# 
# tsne_res <- Rtsne(env_pres, perplexity = 30)
# df <- data.frame(tsne_res$Y, cluster = factor(cl$cluster))
# ggplot(df, aes(x = X1, y = X2, color = cluster)) + geom_point()

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

set.seed(42)
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

# raw weights for presences
pres_weights <- inv_freq[train[which(train$kelp==1), "cluster"]]

# normalize to mean = 1
pres_weights <- pres_weights * (length(pres_weights) / sum(pres_weights))

# add weights to dataset
train$weight <- 1
train$weight[train$kelp == 1] <- pres_weights

train_scaled$weight <- 1
train_scaled$weight[train_scaled$kelp == 1] <- pres_weights

# -----------------------------
# 5. Fit models with weights
# -----------------------------
## GLM
# ==============================================================================
library(car)        # vif
library(broom)      # tidy model output
library(visreg)     # partial effect plots


# Specify the predictors for quadratic terms 
quad_vars <- colnames(train_scaled)[2:9]

# If your actual column names differ (e.g. salinity_summer_SD), update quad_vars accordingly.

# Make a training dataframe with only needed columns + response
train_glm_df <- train_scaled %>%
  select(kelp, all_of(quad_vars)) %>%
  # remove rows with NA in any chosen vars (or consider imputation)
  filter(if_all(all_of(quad_vars), ~ !is.na(.)))

# linear_terms<- train_scaled %>%
#   select(kelp, !all_of(quad_vars)) 
# linear_terms<- colnames(linear_terms)
# linear_terms<- linear_terms[-1]

# Build the formula with quadratic terms using I(x^2)
#    Use I(x^2) rather than poly() so prediction on new data is straightforward.
terms_quad <- paste0(quad_vars, " + I(", quad_vars, "^2)")

# glm_formula <- as.formula(paste("kelp ~", paste(predictors, collapse = " + ")))
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad), collapse = " + ")))

# Fit the GLM (binomial)
glm_mod_s <- glm(glm_formula,
                 data = train_scaled,
                 family = binomial(link = "logit"),
                 weights = weight)

# glm_mod_quad <- glm(glm_formula, data = train_scaled, family = binomial(link = "logit"))
# Call:
#   glm(formula = glm_formula, family = binomial(link = "logit"), 
#       data = train_scaled)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    2.34998    0.41224   5.700 1.19e-08 ***
#   nitrate_winter_minimum        -0.47476    0.18991  -2.500 0.012420 *  
#   I(nitrate_winter_minimum^2)   -0.10936    0.08473  -1.291 0.196801    
# ammonium_winter_minimum        0.26249    0.36707   0.715 0.474546    
# I(ammonium_winter_minimum^2)  -0.10236    0.12337  -0.830 0.406718    
# ammonium_spring_mean          -0.30832    0.19227  -1.604 0.108799    
# I(ammonium_spring_mean^2)     -0.12726    0.08901  -1.430 0.152793    
# temperature_summer_mean       -2.41238    0.37761  -6.388 1.68e-10 ***
#   I(temperature_summer_mean^2)  -0.84984    0.16284  -5.219 1.80e-07 ***
#   PAR_summer_mean                0.35366    0.18375   1.925 0.054267 .  
# I(PAR_summer_mean^2)          -0.04635    0.08347  -0.555 0.578691    
# salinity_summer_SD             0.57366    0.26087   2.199 0.027877 *  
#   I(salinity_summer_SD^2)       -0.01221    0.11851  -0.103 0.917962    
# currentSpeed_summer_mean      -0.47663    0.16228  -2.937 0.003313 ** 
#   I(currentSpeed_summer_mean^2)  0.06082    0.05324   1.142 0.253259    
# turbidity_summer_mean          1.04235    0.54997   1.895 0.058053 .  
# I(turbidity_summer_mean^2)    -0.30534    0.14662  -2.082 0.037300 *  
#   northerness_3x3                0.14678    0.08410   1.745 0.080948 .  
# I(northerness_3x3^2)          -0.45839    0.16931  -2.707 0.006781 ** 
#   easterness_3x3                 0.03001    0.08862   0.339 0.734901    
# I(easterness_3x3^2)           -0.54282    0.16171  -3.357 0.000789 ***
#   ammonium_spring_SD             0.12323    0.23770   0.518 0.604162    
# I(ammonium_spring_SD^2)       -0.09609    0.09654  -0.995 0.319580    
# slope_5x5                      0.99895    0.13656   7.315 2.57e-13 ***
#   I(slope_5x5^2)                -0.19364    0.04315  -4.488 7.19e-06 ***
#   TPI_3x3                       -0.06371    0.08169  -0.780 0.435474    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 1799.4  on 1297  degrees of freedom
# Residual deviance: 1092.6  on 1272  degrees of freedom
# AIC: 1144.6


# Model M2
# Call:
#   glm(formula = glm_formula, family = binomial(link = "logit"), 
#       data = train_scaled)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    1.156e+00  2.097e-01   5.512 3.54e-08 ***
#   ammonium_spring_mean          -2.737e-01  1.522e-01  -1.798 0.072148 .  
# I(ammonium_spring_mean^2)     -2.196e-01  7.041e-02  -3.119 0.001817 ** 
#   currentSpeed_summer_mean      -2.312e-01  1.560e-01  -1.482 0.138300    
# I(currentSpeed_summer_mean^2) -5.343e-06  5.206e-02   0.000 0.999918    
# nitrate_summer_minimum        -7.993e-01  3.805e-01  -2.101 0.035647 *  
#   I(nitrate_summer_minimum^2)   -1.853e-01  1.370e-01  -1.353 0.176128    
# nitrate_winter_mean           -7.412e-01  2.019e-01  -3.672 0.000241 ***
#   I(nitrate_winter_mean^2)      -1.015e-01  7.505e-02  -1.353 0.176167    
# PAR_summer_mean                3.343e-01  1.414e-01   2.364 0.018089 *  
#   I(PAR_summer_mean^2)          -9.342e-02  8.174e-02  -1.143 0.253102    
# temperature_summer_mean       -2.645e+00  3.636e-01  -7.274 3.48e-13 ***
#   I(temperature_summer_mean^2)  -3.488e-01  2.012e-01  -1.733 0.083034 .  
# turbidity_summer_mean          1.788e+00  3.057e-01   5.849 4.95e-09 ***
#   I(turbidity_summer_mean^2)    -4.524e-01  1.110e-01  -4.076 4.58e-05 ***
#   slope_5x5                      9.423e-01  1.312e-01   7.181 6.90e-13 ***
#   I(slope_5x5^2)                -1.753e-01  4.416e-02  -3.969 7.22e-05 ***
#   TPI_3x3                       -1.202e-01  8.106e-02  -1.483 0.138035    
# I(TPI_3x3^2)                   2.362e-03  1.942e-02   0.122 0.903187    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 1799.4  on 1297  degrees of freedom
# Residual deviance: 1108.4  on 1279  degrees of freedom
# AIC: 1146.4
# 
# Number of Fisher Scoring iterations: 8


# Model M3
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    0.90061    0.20642   4.363 1.28e-05 ***
#   ammonium_spring_mean          -0.40874    0.14954  -2.733 0.006270 ** 
#   I(ammonium_spring_mean^2)     -0.14151    0.06771  -2.090 0.036633 *  
#   currentSpeed_summer_mean       0.09529    0.15486   0.615 0.538355    
# I(currentSpeed_summer_mean^2) -0.07934    0.04955  -1.601 0.109343    
# nitrate_summer_minimum        -0.90436    0.40922  -2.210 0.027106 *  
#   I(nitrate_summer_minimum^2)   -0.11743    0.14085  -0.834 0.404420    
# nitrate_winter_mean           -0.74201    0.18224  -4.072 4.67e-05 ***
#   I(nitrate_winter_mean^2)      -0.07539    0.05986  -1.260 0.207821    
# PAR_summer_mean                0.59438    0.14146   4.202 2.65e-05 ***
#   I(PAR_summer_mean^2)          -0.04125    0.08879  -0.465 0.642243    
# temperature_summer_mean       -2.45847    0.38907  -6.319 2.64e-10 ***
#   I(temperature_summer_mean^2)  -0.31540    0.21170  -1.490 0.136261    
# turbidity_summer_mean          1.44732    0.28258   5.122 3.03e-07 ***
#   I(turbidity_summer_mean^2)    -0.36466    0.10123  -3.602 0.000315 ***
#   slope_5x5                      0.73639    0.12380   5.948 2.71e-09 ***
#   I(slope_5x5^2)                -0.15267    0.04157  -3.673 0.000240 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 1768.9  on 1275  degrees of freedom
# Residual deviance: 1084.4  on 1259  degrees of freedom
# AIC: 1118.4
# 
# Number of Fisher Scoring iterations: 8


# WEIGHTED PRESENCES
# Call:
#   glm(formula = glm_formula, family = binomial(link = "logit"), 
#       data = train_scaled, weights = weight)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    1.64668    0.18530   8.887  < 2e-16 ***
#   ammonium_spring_mean          -0.32486    0.10467  -3.104 0.001911 ** 
#   I(ammonium_spring_mean^2)     -0.04517    0.04459  -1.013 0.311051    
# currentSpeed_summer_mean      -0.19891    0.11661  -1.706 0.088041 .  
# I(currentSpeed_summer_mean^2) -0.01791    0.03648  -0.491 0.623529    
# nitrate_summer_minimum         0.06987    0.30775   0.227 0.820398    
# I(nitrate_summer_minimum^2)   -0.37891    0.10567  -3.586 0.000336 ***
#   nitrate_winter_mean           -0.80758    0.13934  -5.796 6.80e-09 ***
#   I(nitrate_winter_mean^2)      -0.03732    0.04746  -0.786 0.431632    
# PAR_summer_mean               -0.36830    0.11641  -3.164 0.001558 ** 
#   I(PAR_summer_mean^2)          -0.41050    0.07014  -5.852 4.85e-09 ***
#   temperature_summer_mean       -1.80576    0.29175  -6.189 6.04e-10 ***
#   I(temperature_summer_mean^2)  -0.46703    0.14632  -3.192 0.001414 ** 
#   turbidity_summer_mean          1.84873    0.22802   8.108 5.16e-16 ***
#   I(turbidity_summer_mean^2)    -0.52357    0.12709  -4.120 3.79e-05 ***
#   slope_5x5                      0.65468    0.09806   6.677 2.45e-11 ***
#   I(slope_5x5^2)                -0.10301    0.03310  -3.112 0.001856 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 2298.5  on 1657  degrees of freedom
# Residual deviance: 1613.8  on 1641  degrees of freedom
# AIC: 1452.5
# 
# Number of Fisher Scoring iterations: 8


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


# M3: 


## GAM
# ==============================================================================
gam_formula <- as.formula(paste("kelp ~", paste0("s(", predictors, ")", collapse = " + ")))
# gam_mod <- gam(gam_formula, data = train_scaled, family = binomial)
gam_mod_s <- gam(gam_formula,
                                  data = train_scaled,
                                  family = binomial,
                                  weights = weight)
                 
                 
visreg(gam_mod_s, "ammonium_spring_mean", scale = "response")

## Random Forest
# ==============================================================================
train$kelp <- as.factor(train$kelp)
rf_mod <-randomForest(x = train[, predictors], y = train$kelp, ntree = 500)


rf_mod_s <- randomForest(as.formula(paste("kelp ~", paste(predictors, collapse = " + "))),
                         data = train,
                         ntree = 500,
                         importance = TRUE,
                         sampsize = nrow(train),
                         case.weights = train$weight)

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
# brt_mod <- gbm(glm_formula,
#                data = train,
#                distribution = "bernoulli",
#                n.trees = 2000,
#                interaction.depth = 3,
#                shrinkage = 0.01,
#                bag.fraction = 0.5,
#                train.fraction = 1.0,
#                cv.folds = 5,
#                verbose = FALSE)
# ## BRT (gbm.step)
# train_brt_sel <- train
# train_brt_sel$kelp <- as.numeric(as.character(train_brt_sel$kelp))

train_brt <- train
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



# brt_mod <- dismo::gbm.step(data = train_brt, 
#                     gbm.x = which(names(train_brt) != "kelp"),
#                     gbm.y = which(names(train_brt) == "kelp"),
#                     family = "bernoulli",
#                     tree.complexity = 3,
#                     learning.rate = 0.01,
#                     bag.fraction = 0.5)

# fitting final gbm model with a fixed number of 2700 trees for kelp
# 
# mean total deviance = 1.386 
# mean residual deviance = 0.411 
# 
# estimated cv deviance = 0.783 ; se = 0.033 
# 
# training data correlation = 0.891 
# cv correlation =  0.716 ; se = 0.013 
# 
# training data AUC score = 0.985 
# cv AUC score = 0.908 ; se = 0.008 
# 
# elapsed time -  0.41 minutes 

# M2. Sep 4
# mean total deviance = 1.386 
# mean residual deviance = 0.504 
# 
# estimated cv deviance = 0.774 ; se = 0.031 
# 
# training data correlation = 0.851 
# cv correlation =  0.72 ; se = 0.015 
# 
# training data AUC score = 0.971 
# cv AUC score = 0.908 ; se = 0.007 
# 
# elapsed time -  0.26 minutes 


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

# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_for_SelecBasedOnImp_v4.pdf", width = 12, height = 16, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_for_SelecBasedOnImp_v4_stacked.pdf", width = 12, height = 16, dpi= 300, units="cm")


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
plot_variable_ranking_v2(res$summary_table, threshold = 15)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_for_SelecBasedOnImp_M4.pdf", width = 15, height = 13, dpi= 300, units="cm")


plot_importance_profiles_v2(res$summary_table, threshold=15)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/VariablesImportance_profile_step_1_FINAL4.png", width = 18, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/VariablesImportance_profile_step_1_FINAL3.pdf", width = 18, height = 12, dpi= 300, units="cm")





# ==============================================================================
# STEP 7: SDMs with selected variables 
# ==============================================================================
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M4")
thresh= 15

vars_selected <- rank_variables(models, types, avg_threshold = 15)
vars_selected <- vars_selected$summary_table$Variable

# vars_selected
# [1] "temperature_summer_mean"  "slope_5x5"                "turbidity_summer_mean"    "ammonium_spring_SD"      
# [5] "PAR_summer_mean"          "currentSpeed_summer_mean" "salinity_summer_SD"      

# 4 Sep 2025
# [1] "slope_5x5"                "nitrate_summer_minimum"   "turbidity_summer_mean"    "nitrate_winter_mean"      "temperature_summer_mean" 
# [6] "ammonium_spring_mean"     "currentSpeed_summer_mean" "PAR_summer_mean"         


# Excluded nitrate_summer_minimum because highly correlated with summer mean temperature 
vars_selected<- c("turbidity_summer_mean",    "temperature_summer_mean",  "nitrate_winter_mean", #"nitrate_summer_minimum",
                  "slope_5x5",     "currentSpeed_summer_mean",
                  "ammonium_spring_mean",     "PAR_summer_mean" )


# Refit models using only selected vars
train_sel <- train %>% select(all_of(c("kelp", vars_selected)))
test_sel <- test %>% select(all_of(c("kelp", vars_selected)))

train_sel$kelp<- as.factor(train_sel$kelp)
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
freq <- table(train_sel_scaled[which(train_sel_scaled$kelp==1), "cluster"])
inv_freq <- 1 / as.numeric(freq)

# raw weights for presences
pres_weights <- inv_freq[train_sel_scaled[which(train_sel_scaled$kelp==1), "cluster"]]

# normalize to mean = 1
pres_weights <- pres_weights * (length(pres_weights) / sum(pres_weights))

# add weights to dataset
train_sel_scaled$weight <- 1
train_sel_scaled$weight[train_sel_scaled$kelp == 1] <- pres_weights

train_sel$weight <- 1
train_sel$weight[train_sel$kelp == 1] <- pres_weights


# ==============================================================================
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M5_weightedPres")

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



# ==============================================================================
# STEP 8: Models Evaluation
# ==============================================================================
# ---- Load packages ----
library(blockCV)         
library(caret)           
library(randomForest)    
library(mgcv)            
library(gbm)             
library(pROC)            
library(PresenceAbsence) 
library(dplyr)
library(tidyr)
library(raster)   # or terra if you use SpatRaster
library(gstat)
library(sp)

# ---- Example data (replace with your own) ----
# data: dataframe with predictors (X1, X2, ...) and response (binary 0/1 in "y")
# coords: matrix/data.frame of coordinates (longitude, latitude)
# coords <- data[,c("lon","lat")]

set.seed(123)

# ---- Example data ----
# predictors: RasterStack or SpatRaster of environmental variables
# coords: matrix/data.frame of coordinates (longitude, latitude)
# species: dataframe with presence/absence, column "y"

# ---- Option 1: Interactive exploration using blockCV ----
# Opens a Shiny app to explore spatial autocorrelation range
predictors_rast<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_resolution_FINAL2.tif")
  
# rangeExplorer(rasterLayer = predictors_rast,    # your environmental rasters
#               speciesData = train_sel,       # your species records
#               species = "kelp")               # response variable column
# 
# In the app:
# 1. Select a predictor
# 2. Look at the variogram
# 3. Identify the distance where semivariance flattens (range)
# 4. Use that distance as 'theRange' in spatialBlock

# ---- Option 2: Compute variogram for one predictor manually ----
# Example using the first raster layer
var_data <- train.xy

coordinates(var_data) <- ~x+y

# Compute empirical variogram
vgm_model <- variogram(temperature_summer_mean ~ 1, var_data)

# Fit a theoretical model (e.g., spherical)
fit <- fit.variogram(vgm_model, model = vgm(psill=0.5, model="Sph", nugget=1, range=500))

# Plot the variogram
plot(vgm_model, fit)

# ---- Identify block size ----
# The 'range' of the fitted model is the distance at which autocorrelation
# essentially disappears. Use this as 'theRange' in spatialBlock:
block_size <- fit$range[2]   # or read from the plot manually
print(paste("Suggested block size (m):", block_size))


# ---- Function to calculate metrics ----
get_metrics <- function(obs, pred){
  thresh <- optimal.thresholds(data.frame(ID=1:length(obs), obs=obs, pred=pred))$pred
  cm <- cmx(data.frame(ID=1:length(obs), obs=obs, pred=pred), threshold=thresh)
  
  sens <- cm[2,2] / sum(cm[2,])         
  spec <- cm[1,1] / sum(cm[1,])         
  tss  <- sens + spec - 1
  roc  <- auc(obs, pred)                
  
  return(c(TSS=tss, ROC=roc, Sensitivity=sens, Specificity=spec))
}


# ---- Models ----
models <- c("glm","gam","rf","brt")
all_results <- list()


for (m in models){
  metrics <- matrix(NA, nrow=max(folds), ncol=4)
  colnames(metrics) <- c("TSS","ROC","Sensitivity","Specificity")
  
  for (i in 1:max(folds)){
    train_data <- data[folds != i, ]
    test_data  <- data[folds == i, ]
    
    if(m == "glm"){
      mod <- glm(y ~ ., data=train_data, family=binomial)
      pred <- predict(mod, newdata=test_data, type="response")
    }
    if(m == "gam"){
      mod <- gam(y ~ s(X1) + s(X2) + s(X3), data=train_data, family=binomial)
      pred <- predict(mod, newdata=test_data, type="response")
    }
    if(m == "rf"){
      mod <- randomForest(as.factor(y) ~ ., data=train_data, ntree=500)
      pred <- predict(mod, newdata=test_data, type="prob")[,2]
    }
    if(m == "brt"){
      mod <- gbm(y ~ ., data=train_data,
                 distribution="bernoulli",
                 n.trees=2000,
                 interaction.depth=3,
                 shrinkage=0.01,
                 n.minobsinnode=10,
                 verbose=FALSE)
      pred <- predict(mod, newdata=test_data, type="response", n.trees=2000)
    }
    
    metrics[i,] <- get_metrics(test_data$y, pred)
  }
  
  # Collect per-model results
  all_results[[m]] <- data.frame(
    Model = m,
    Metric = colnames(metrics),
    Mean = colMeans(metrics, na.rm=TRUE),
    SD = apply(metrics, 2, sd, na.rm=TRUE)
  )
}


# ---- Combine and pivot to wide table ----
results_table <- bind_rows(all_results) %>%
  pivot_wider(
    names_from = Metric,
    values_from = c(Mean, SD),
    names_glue = "{Metric}_{.value}"
  )

print(results_table)





















# Get metrics of model performance, AUC, Sensitivity and specificity
results_selected <- bind_rows(
  get_metrics_optimized2(model= glm_mod_se, test_data=train_sel, scale_params= scaling_params_2,
                         model_name = "glm", threshold_type = "youden"),
  get_metrics_optimized2(gam_mod_se, train_sel, scale_params= scaling_params_2, "gam",
                         threshold_type = "youden"),#10pct_omission
  get_metrics_optimized2(rf_mod_se, train_sel, "rf", threshold_type = "youden"),
  get_metrics_optimized2(brt_mod_se, train_sel, "brt", threshold_type = "youden")
)

# write.csv(results_selected, "models_performance_train_Table3_thresholds_Sep2025_M3.csv")

# Pivot the results for plotting
results_long_selected <- results_selected %>%
  pivot_longer(cols = c("AUC", "Sensitivity", "Specificity"),
               names_to = "Metric", values_to = "Value")

head(results_long_selected, 19)


# Plot
ggplot(results_long_selected, aes(fill = Model, y = Value, x = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Model Performance Comparison",
       y = "Score",
       x = "Model") +
  ylim(c(0,1))+
  theme_bw()



# roc_df_sel <- bind_rows(
#   get_roc_data(glm_mod_s, test_sel, "GLM"),
#   get_roc_data(gam_mod_s, test_sel, "GAM"),
#   get_roc_data(rf_mod_s,  test, "RF"),
#   get_roc_data(brt_mod_s,  test, "BRT")
# )


# # Plot ROC curves
# ggplot(roc_df_sel, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
#   geom_line(size = 1) +
#   geom_abline(linetype = "dashed", color = "grey") +
#   labs(title = "ROC Curves for Model Comparison",
#        x = "False Positive Rate (1 - Specificity)",
#        y = "True Positive Rate (Sensitivity)") +
#   theme_minimal()



# Get predictions
pred_df_sel <- bind_rows(
  data.frame(Probability = predict(glm_mod_s, newdata = test_sel_scaled, type = "response"), Model = "GLM", Truth = test_sel$kelp),
  data.frame(Probability = predict(gam_mod_s, newdata = test_sel_scaled, type = "response"), Model = "GAM", Truth = test_sel$kelp),
  data.frame(Probability = predict(rf_mod_s, newdata = test_sel, type = "prob")[, 2], Model = "RF", Truth = test_sel$kelp),
  data.frame(Probability = predict(brt_mod_s, newdata = test_sel, n.trees = brt_mod_s$n.trees, type = "response"), Model = "BRT", Truth = test_sel$kelp)
)

# Plot distributions
ggplot(pred_df_sel, aes(x = Probability, fill = as.factor(Truth))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = c("0" = "red", "1" = "blue"), name = "Observed") +
  labs(title = "Predicted Probability Distributions by Model",
       x = "Predicted Probability", y = "Density") +
  theme_minimal()


models_selected <- list(glm_mod_s, gam_mod_s, rf_mod_s, brt_mod_s)
types <- c("glm", "gam", "rf", "brt")

plot_var_importance_v2(models_selected, types,  top_n = 10, stacked = F)

variablesImportance_models_selected <- rank_variables(models_selected, types, avg_threshold = NULL)
plot<- plot_variable_ranking(variablesImportance_models_selected$summary_table, threshold = NULL)
plot[[2]]


summary_table<- plot[[1]]
summary_table <- summary_table %>%
  arrange(desc(weighted_score)) %>%
  mutate(Variable = factor(Variable, levels = Variable))  # preserve this order

# 2️⃣ Define expression labels
# var_labels_expr <- c(
#   "temperature_summer_mean"   = expression("Mean Summer Temperature ("*degree*C*")"),
#   "slope_5x5"                 = expression("Slope ("*degree*"; 100 x 100 m)"),
#   "turbidity_summer_mean"     = expression("Mean Summer Turbidity"),
#   "ammonium_spring_SD"        = expression("SD of Spring Ammonium"),
#   "PAR_summer_mean"           = expression("Mean Summer PAR"),
#   "currentSpeed_summer_mean"  = expression("Mean Summer Current Speed (m s^-1)"),
#   "salinity_summer_SD"        = expression("SD of Summer Salinity")
# )

var_labels_expr <- c(
  # "nitrate_summer_minimum"    = expression("Min of Summer Nitrate (uM.L^-1)"),
  "temperature_summer_mean"   = expression("Mean Summer Temperature ("*degree*C*")"),
  "slope_5x5"                 = expression("Slope ("*degree*"; 100 x 100 m)"),
  "turbidity_summer_mean"     = expression("Mean Summer Turbidity"),
  "ammonium_spring_mean"        = expression("Mean of Spring Ammonium"),
  "PAR_summer_mean"           = expression("Mean Summer PAR"),
  "currentSpeed_summer_mean"  = expression("Mean Summer Current Speed (m s^-1)"),
  "nitrate_winter_mean"        = expression("Mean of Winter Nitrate")
)

# 3️⃣ Extract the labels in the same order as the factor levels
labels_to_use <- var_labels_expr[levels(summary_table$Variable)]

# 4️⃣ Plot
ggplot(summary_table, aes(x = avg_importance, y = Variable, fill = ColorFlag)) +
  geom_col() +
  scale_fill_manual(values = c("Above" = "grey70", "Below" = "grey80"), guide = "none") +
  scale_y_discrete(labels = labels_to_use) +  # now labels match the reordered factor
  labs(x = "Weighted Score", y = NULL, title = "") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 9, face = "bold")
  )

# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_SelectedModels.pdf", width = 11, height = 6, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_SelectedModels_Sep2025.png", width = 11, height = 6, dpi= 300, units="cm")


plot_importance_profiles(variablesImportance_models_selected$summary_table)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_profile_selectedModels_Sep2025_M4_v2.png", width = 18, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_profile_selectedModels.pdf", width = 18, height = 12, dpi= 300, units="cm")



### PLOT RESPONSE CURVES ======================================================
# Combine response curves across models
create_full_grid <- function(var="PAR_summer_mean", train_sel, model_type = "glm", 
                             quad_vars = NULL, scale_params = NULL, n = 50) {
  
  predictors <- setdiff(names(train_sel), "kelp")
  
  # Sequence for the variable of interest
  x_unscaled <- seq(min(train_sel[[var]], na.rm = TRUE),
                    max(train_sel[[var]], na.rm = TRUE),
                    length.out = n)
  
  # Create a typical row: mean for numeric, mode for factor
  typical_row <- data.frame(matrix(nrow = 1, ncol = 0))
  for (col in predictors) {
    if (col == var) next
    if (is.numeric(train_sel[[col]])) {
      typical_row[[col]] <- mean(train_sel[[col]], na.rm = TRUE)
    } else if (is.factor(train_sel[[col]])) {
      typical_row[[col]] <- names(sort(table(train_sel[[col]]), decreasing = TRUE))[1]
    }
  }
  
  # Replicate the typical row n times and vary the variable of interest
  grid <- typical_row[rep(1, times = n), , drop = FALSE]
  
  grid[[var]] <- x_unscaled
  
  # Handle scaling and quadratic terms
  if (model_type == "glm") {
    if (is.null(scale_params)) stop("GLM requires scale_params.")
    
    scale_grid <- function(grid, params) {
      grid %>%
        mutate(across(where(is.numeric),
                      ~ (.-params[[paste0(cur_column(), "_mean")]]) /
                        params[[paste0(cur_column(), "_sd")]]))
    }
    
    grid_scaled <- scale_grid(grid, scale_params)
    
    # Add quadratic terms for all variables in quad_vars
    if (!is.null(quad_vars)) {
      for (qv in quad_vars) {
        if (qv %in% names(grid_scaled)) {
          if (qv == var) {
            # variable being varied
            grid_scaled[[paste0("I(", qv, "^2)")]] <- grid_scaled[[qv]]^2
          } else {
            # constant variable → square of its mean
            grid_scaled[[paste0("I(", qv, "^2)")]] <- (grid_scaled[[qv]][1])^2
          }
        }
      }
    }
    
    grid_to_predict <- grid_scaled
    
  } else if (model_type == "gam") {
    if (is.null(scale_params)) stop("GAM requires scale_params.")
    grid_to_predict <- as.data.frame(scale_grid(grid, scale_params))
    
  } else if (model_type %in% c("randomForest", "gbm")) {
    # no scaling, no quadratic terms
    grid_to_predict <- grid
    
  } else {
    stop("Unsupported model type")
  }
  
  return(list(grid = grid_to_predict, x_unscaled = x_unscaled))
}



get_curve2 <- function(model, var, grid_info, model_name) {
  
  x_unscaled <- grid_info$x_unscaled
  grid <- grid_info$grid
  cls <- class(model)[1]
  
  if (cls %in% c("glm", "gam")) {
    pred <- predict(model, newdata = grid, type = "response")
    df <- data.frame(
      x = x_unscaled,
      fit = pred,
      model = model_name,
      var = var
    )
    
  } else if (cls %in% c("randomForest", "randomForest.formula")) {
    pd <- pdp::partial(object = model,
                       pred.var = var,
                       train = grid,
                       prob = TRUE)
    df <- data.frame(
      x = pd[[var]],
      fit = 1 - pd$yhat,  # adjust depending on your positive class
      model = model_name,
      var = var
    )
    
  } else if (cls == "gbm") {
    pd <- pdp::partial(object = model,
                       pred.var = var,
                       train = grid,
                       pred.fun = function(object, newdata) {
                         predict(object, newdata, n.trees = model$n.trees, type = "response")
                       },
                       recursive = FALSE)
    df <- data.frame(
      x = pd[[var]],
      fit = pd$yhat,
      model = model_name,
      var = var
    )
    
  } else {
    stop("Unsupported model class")
  }
  
  return(df)
}



# Variables to plot
vars_to_plot <- names(train_sel)[-1]  # all variables except response

# Models list
models_list <- list(
  list(model = glm_mod_s, type = "glm", name = "GLM", quad_vars = quad_vars_sel),
  list(model = gam_mod_s, type = "gam", name = "GAM", quad_vars = NULL),
  list(model = rf_mod_s, type = "randomForest", name = "RF", quad_vars = NULL),
  list(model = brt_mod_s, type = "gbm", name = "GBM", quad_vars = NULL)
)

# Loop over variables and models
all_curves <- lapply(vars_to_plot, function(v) {
  do.call(rbind, lapply(models_list, function(m) {
    # Create the grid
    grid_info <- create_full_grid(
      var = v,
      train_sel = train_sel,
      model_type = m$type,
      quad_vars = m$quad_vars,
      scale_params = scaling_params_2
    )
    # Predict
    get_curve2(
      model = m$model,
      var = v,
      grid_info = grid_info,
      model_name = m$name
    )
  }))
})

# Combine all variables into a single data frame
all_curves_df <- do.call(rbind, all_curves)

# var_labels <- c(
#   "temperature_summer_mean"   = "Mean~Summer~Temperature~(degree*C)", #\n
#   "slope_5x5"                 = "Slope~('°')",
#   "turbidity_summer_mean"     = "Mean~Summer~Turbidity",
#   "ammonium_spring_SD"        = "SD~of~Spring~Ammonium",
#   "PAR_summer_mean"           = "Mean~Summer~PAR",
#   "currentSpeed_summer_mean"  = "Mean~Summer~\nCurrent~Speed~(m.s^{-1})",
#   "salinity_summer_SD"        = "SD~of~Summer~Salinity"
# )

var_labels <- c(
  # "nitrate_summer_minimum"    = "Min of Summer Nitrate (uM.L-1)",
  "temperature_summer_mean"   = "Mean Summer Temperature ('°'C)",
  "slope_5x5"                 = "Slope~('°')~; 100 x 100 m)",
  "turbidity_summer_mean"     = "Mean Summer Turbidity",
  "ammonium_spring_mean"        = "Mean of Spring Ammonium",
  "PAR_summer_mean"           = "Mean Summer PAR",
  "currentSpeed_summer_mean"  = "Mean Summer Current Speed (m s-1)",
  "nitrate_winter_mean"        = "Mean of Winter Nitrate"
)


# Faceted plot for all variables
ggplot(all_curves_df, aes(x = x, y = fit, color = model)) +
  geom_line(size = 1) +
  facet_wrap(~var, scales = "free_x", ncol=4,
             labeller = labeller(var = as_labeller(var_labels_expr, label_parsed))) +
  labs(
    x = "Predictor value",
    y = "Predicted response",
    color = "Model",
    title = "Response Curves for All Variables and Models"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 7),   # facet label size
    axis.title = element_text(size = 9),                 # axis title size
    axis.text = element_text(size = 7),                  # axis tick labels
    plot.title = element_text(size = 9, face = "bold")   # plot title
  )


# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_Sep2025_M4_v2.png", width = 19, height = 9, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_Sep2025.pdf", width = 19, height = 9, dpi= 300, units="cm")


# Plot summary average curves among models: 
summary_curves <- all_curves_df %>%
  group_by(var, x) %>%
  summarise(
    mean_fit = mean(fit),
    sd_fit   = sd(fit),
    ymin = min(fit),   # minimum across models
    ymax = max(fit),   # maximum across models
    .groups = "drop"
  )


# Plots
ggplot(summary_curves, aes(x = x, y = mean_fit)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax),
              fill = "lightblue", alpha = 0.3) +  # shaded area = full range
  geom_line(color = "blue", size = 1) +       # average curve
  facet_wrap(~var, scales = "free_x", ncol=4) +
             # labeller = labeller(var = as_labeller(var_labels, label_parsed))) +
  labs(
    x = "Predictor value",
    y = "Predicted response",
    title = "Average Response Curve with Full Range Across Models"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 7),   # facet label size
    axis.title = element_text(size = 9),                 # axis title size
    axis.text = element_text(size = 7),                  # axis tick labels
    plot.title = element_text(size = 9, face = "bold")   # plot title
  )


# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_average_Sep2025_M4.png", width = 19, height = 9, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_average_Sep2025.pdf", width = 19, height = 9, dpi= 300, units="cm")



