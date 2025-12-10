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
#rast(tif_files[c(2,5,12,15)])
slope<-  rast(tif_files[c(12,15)])

# Load and stack all tif files
# terrain_vars<- terra::rast(paste(variables_selection_path,"terrain_rasters_selected.tif", sep="/"))
# names(terrain_vars)<- c("easterness_3x3",  "northerness_3x3", "slope_5x5" , "TPI_3x3")

plot(terrain_vars)

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
substrate<- crop(substrate, slope)
substrate<- as.factor(substrate)

# Mask model predictions by substrate 
substrate_aligned <- terra::rast(slope)
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
raster_stack_20m_all<- c(raster_stack_20m, slope, substrate_aligned, bathy)
names(raster_stack_20m_all)
plot(raster_stack_20m_all[[10:11]])


### Extract raster values at kelp locations ====================================
kelp_presabs_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL2.csv")
colnames(kelp_presabs_df)
kelp_presabs_df<- kelp_presabs_df%>%
  select(c(x,y,kelp))#[,c(2,3,11:14)]

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
  group_by(substrate_aligned)%>%
  count()
# substrate_aligned     n
# 1 1                  5921
# 2 2                     5  # only 5 records in soft substrate were excluded
# 3 NA                  217  # 217 records with NAs of substrate were excluded


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

# After filtering by substrate (it already filter by depth):
#    kelp     n depth_max depth_min
# 1     0  5000      40.0     -9.99
# 2     1   921      20.8    -47.6 


kelp_data_with_variables<- kelp_data_with_variables%>%
  filter(coastwide_20m >= -10)%>%
  filter(coastwide_20m <= 40)

# kelp     n depth_max depth_min
# 1     0  5217      40.0     -9.99
# 2     1   926      39.6     -9.81

# After filtering by substrate:
#    kelp     n depth_max depth_min
# 1     0  5000      40.0     -9.99 # absences go until 40 m depth
# 2     1   911      20.8     -9.35 # presences go until 20.8 m depth, 10 presences were above 10m (in the coast)


summary(kelp_data_with_variables)
# write.csv(kelp_data_with_variables, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M3.csv")


#### PERFORM MODELS ############################################################
# Load necessary libraries
library(dplyr)
library(caret)
library(pROC)
library(mgcv)      # GAM
library(randomForest) # RF
library(gbm)       # BRT

set.seed(123) # for reproducibility

# ================================
# STEP 1: Balance dataset
# ================================
# Keep all presences

df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/Presence_absences_kelp_variablesSelectedFINAL_M3.csv")
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

# Downsample absences to match presence count
presences <- df %>% filter(kelp == 1)
absences <- df %>% filter(kelp == 0) %>% sample_n(nrow(presences))

summary(presences$bathymetry_20)
quantile(presences$bathymetry_20, 0.9) #26.51765   # 8.848001
quantile(presences$bathymetry_20, 0.95) #32.6647   # 12.02184
quantile(presences$bathymetry_20, 0.93) #30.14825  # 10.56405
quantile(presences$bathymetry_20, 0.99) #37.68245  # 17.72347 

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

# 1 test    546
# 2 train  1276

# write.csv(train_test_dataset, "/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_FINAL3.csv")


train_test_dataset<- read.csv("/Volumes/Romina_PSF/PSF/SDM/SDM_results/training_testing_datasets_blob_FINAL3.csv")
train_test_dataset<- (train_test_dataset[,-1])

colnames(train_test_dataset)
predictors <- setdiff(
  names(train_test_dataset),
  c("kelp", "bathymetry_20","x", "y", "set", "cell","substrate_aligned"#"cell_id", "env_cluster", "easterness_3x3",           "northerness_3x3" 
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

linear_terms<- train_scaled %>%
  select(kelp, !all_of(quad_vars)) 

linear_terms<- colnames(linear_terms)
linear_terms<- linear_terms[-1]

# Build the formula with quadratic terms using I(x^2)
#    Use I(x^2) rather than poly() so prediction on new data is straightforward.
terms_quad <- paste0(quad_vars, " + I(", quad_vars, "^2)")

# glm_formula <- as.formula(paste("kelp ~", paste(predictors, collapse = " + ")))
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad, linear_terms), collapse = " + ")))
# glm_formula_reduced <- kelp ~ 
#   ammonium_spring_mean + I(ammonium_spring_mean^2) +
#   temperature_summer_mean + I(temperature_summer_mean^2) +
#   PAR_summer_mean + I(PAR_summer_mean^2) +
#   salinity_summer_mean + I(salinity_summer_mean^2) +
#   currentSpeed_summer_mean + I(currentSpeed_summer_mean^2) +
#   turbidity_summer_mean + I(turbidity_summer_mean^2) +
#   slope_7x7 + TPI_3x3

# M2 and M3
# kelp ~ ammonium_spring_mean + I(ammonium_spring_mean^2) + currentSpeed_summer_mean + 
#   I(currentSpeed_summer_mean^2) + nitrate_summer_minimum + 
#   I(nitrate_summer_minimum^2) + nitrate_winter_mean + I(nitrate_winter_mean^2) + 
#   PAR_summer_mean + I(PAR_summer_mean^2) + temperature_summer_mean + 
#   I(temperature_summer_mean^2) + turbidity_summer_mean + I(turbidity_summer_mean^2) + 
#   slope_5x5 + I(slope_5x5^2)
summary(train_scaled)

# Fit the GLM (binomial)
glm_mod_quad <- glm(glm_formula, data = train_scaled, family = binomial(link = "logit"))
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


library(patchwork)
vars <- attr(terms(glm_mod_quad), "term.labels")
vars_clean <- vars[!grepl("^I\\(", vars)]

plots <- lapply(vars_clean, function(v) {
  visreg(glm_mod_quad, v, scale = "response", gg = TRUE) + 
    ggtitle(paste("Effect of", v))
})

# Combine into one grid
wrap_plots(plots)



# Basic summary & multicollinearity check
summary(glm_mod_quad)
vif_vals <- vif(glm_mod_quad)     # check for high VIFs (>5 or >10)
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
gam_mod <- gam(gam_formula, data = train_scaled, family = binomial)

visreg(gam_mod, "ammonium_spring_mean", scale = "response")

## Random Forest
# ==============================================================================
train$kelp <- as.factor(train$kelp)
rf_mod <-randomForest(x = train[, predictors], y = train$kelp, ntree = 500)


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


train_brt <- train
colnames(train_brt)
train_brt$kelp <- as.numeric(as.character(train_brt$kelp))
brt_mod <- dismo::gbm.step(data = train_brt, 
                    gbm.x = which(names(train_brt) != "kelp"),
                    gbm.y = which(names(train_brt) == "kelp"),
                    family = "bernoulli",
                    tree.complexity = 3,
                    learning.rate = 0.01,
                    bag.fraction = 0.5)

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




### Save models

saveRDS(glm_mod_s,)

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
    get_curve(model=glm_mod_quad, v, grid, "GLM", scaling_params),
    get_curve(gam_mod,      v, grid, "GAM", scaling_params),
    get_curve(rf_mod,       v, grid, "randomForest"),
    get_curve(brt_mod,      v, grid, "BRT")
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

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/ResponseCurves_allVariables_FINAL_M3.png", width = 20, height = 17, dpi= 300, units="cm")


#-----------------------------
# Collect results
results <- bind_rows(
  get_metrics_optimized(glm_mod_quad, train_scaled, "glm"),
  get_metrics_optimized(gam_mod, train_scaled, "gam"),
  get_metrics_optimized(rf_mod, train, "rf"),
  get_metrics_optimized(brt_mod, train, "brt")
)

# Model Threshold   AUC Sensitivity Specificity   TSS
# 1 glm       0.632 0.886       0.812       0.828 0.639
# 2 gam       0.599 0.911       0.831       0.864 0.694
# 3 rf        0.495 1           1           1     1    
# 4 brt       0.488 0.981       0.955       0.918 0.873


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
  get_roc_data(glm_mod_quad, test_scaled, "GLM"),
  get_roc_data(gam_mod, test_scaled, "GAM"),
  get_roc_data(rf_mod,  test, "RF"),
  get_roc_data(brt_mod,  test, "BRT")
)

# Plot ROC curves
ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "grey") +
  labs(title = "ROC Curves for Model Comparison",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()



# Get predictions
pred_df <- bind_rows(
  data.frame(Probability = predict(glm_mod_quad, newdata = test_scaled, type = "response"), Model = "GLM", Truth = test_scaled$kelp),
  data.frame(Probability = predict(gam_mod, newdata = test_scaled, type = "response"), Model = "GAM", Truth = test_scaled$kelp),
  data.frame(Probability = predict(rf_mod, newdata = test, type = "prob")[, 2], Model = "RF", Truth = test$kelp),
  data.frame(Probability = predict(brt_mod, newdata = test, n.trees = brt_mod$n.trees, type = "response"), Model = "BRT", Truth = test$kelp)
)

# Plot distributions
ggplot(pred_df, aes(x = Probability, fill = as.factor(Truth))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Model) +
  scale_fill_manual(values = c("0" = "red", "1" = "blue"), name = "Observed") +
  labs(title = "Predicted Probability Distributions by Model",
       x = "Predicted Probability", y = "Density") +
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
models <- list(glm_mod_quad, gam_mod, rf_mod, brt_mod)
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
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_for_SelecBasedOnImp_FINAL2.pdf", width = 15, height = 13, dpi= 300, units="cm")


plot_importance_profiles_v2(res$summary_table, threshold=15)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/VariablesImportance_profile_step_1_FINAL3.png", width = 18, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/VariablesImportance_profile_step_1_FINAL3.pdf", width = 18, height = 12, dpi= 300, units="cm")


# =============================================================================
# STEP 6: SDMs with different variables selection based on thresholds
# =============================================================================
# Variables passing the threshold
vars_selected <- rank_variables(models, types, avg_threshold = 15)
vars_selected <- vars_selected$summary_table$Variable
# [1] "ammonium_spring_SD"       "ammonium_summer_mean"     "slope_7x7"                "temperature_summer_mean" 
# [5] "salinity_summer_SD"       "turbidity_summer_mean"    "currentSpeed_summer_mean" "slope_3x3"               
# [9] "ammonium_spring_mean"     "PAR_summer_mean"  

# After running models GLM and GAM with scaled variables
# [1] "slope_7x7"               "ammonium_spring_SD"      "ammonium_winter_mean"    "PAR_summer_mean"         "temperature_summer_mean"
# [6] "turbidity_summer_mean"   "salinity_summer_SD"      "nitrate_summer_mean"     "ammonium_summer_mean"    "ammonium_spring_mean"   

# After removing presence/absences from non-mapped areas 
# [1] "slope_7x7"                "ammonium_spring_SD"       "PAR_summer_mean"          "ammonium_summer_mean"     "slope_3x3"               
# [6] "temperature_summer_mean"  "turbidity_summer_mean"    "salinity_summer_SD"       "salinity_summer_mean"     "currentSpeed_summer_mean"
# [11] "ammonium_winter_mean"  

# After removing presence/absences from non-mapped areas and adding quadratic terms to the GLM
# [1] "slope_7x7"                "temperature_summer_mean"  "turbidity_summer_mean"    "slope_3x3"               
# [5] "ammonium_winter_minimum"  "PAR_summer_mean"          "currentSpeed_summer_mean"

# After correcting slope layer (previously was same as bathymetry)
# [1] "temperature_summer_mean"  "slope_5x5"                "turbidity_summer_mean"    "ammonium_spring_SD"       "PAR_summer_mean"         
# [6] "currentSpeed_summer_mean" "salinity_summer_SD"      


# After removing SD metrics and easterness and northerness --> M2 and M3
# "slope_5x5"                "nitrate_summer_minimum"   "turbidity_summer_mean"    "nitrate_winter_mean"      "temperature_summer_mean" 
# [6] "ammonium_spring_mean"     "currentSpeed_summer_mean" "PAR_summer_mean"         






# # Refit models using only selected vars
# train_sel <- train %>% select(all_of(c("kelp", vars_selected)))
# test_sel <- test %>% select(all_of(c("kelp", vars_selected)))
# 
# scaling_params_2 <- train_sel %>%
#   summarise(across(where(is.numeric),
#                    list(mean = mean, sd = sd), na.rm = TRUE))
# 
# train_sel_scaled <- scale_with_params(train_sel, scaling_params_2)
# test_sel_scaled  <- scale_with_params(test_sel,  scaling_params_2)
# 
# train_sel_scaled$kelp<- as.factor(train_sel_scaled$kelp)
# test_sel_scaled$kelp<- as.factor(test_sel_scaled$kelp)
# 
# 


# Define thresholds to try
thresholds <- c(20, 17, 15)

# Initialize results list
cv_results <- list()

# Run for all thresholds
results2 <- lapply(thresholds, run_models_with_threshold_curves,
                  train_data = train, rank_fun = rank_variables,
                  scaling_params = scaling_params, quad_vars = quad_vars)
# saveRDS(results2, "/Volumes/Romina_PSF/PSF/SDM/SDM_results/comparison_variables_selection_models_FINAL.rds")
  
cv_df <- bind_rows(
  lapply(results2, function(res) {
    do.call(rbind, lapply(c("GLM","GAM","RF","BRT"), function(mod) {
      data.frame(
        Threshold = res$threshold,
        Model = mod,
        Fold = seq_along(res$aucs[[mod]]$aucs),
        AUC = res$aucs[[mod]]$aucs
      )
    }))
  })
)

# 1. Compute Friedman + posthoc letters
summary(cv_df)


library(PMCMRplus)
library(dplyr)
library(PMCMRplus)
library(multcompView)

# Initialize list to store results
model_comp_list <- list()

for(mod in unique(cv_df$Model)) {
  
  # 1. Subset data for this model
  df_mod <- cv_df %>% filter(Model == mod)
  
  thresholds <- sort(unique(df_mod$Threshold))
  
  # 2. Default letters
  letters_full <- setNames(rep("a", length(thresholds)), thresholds)
  
  friedman_p <- NA  # default
  
  if(length(thresholds) > 1) {
    # 3. Friedman test
    ft <- friedman.test(AUC ~ Threshold | Fold, data = df_mod)
    friedman_p <- ft$p.value
    
    # 4. Post-hoc Nemenyi
    ph <- frdAllPairsNemenyiTest(AUC ~ Threshold | Fold, data = df_mod)
    letters <- multcompView::multcompLetters(ph$p.value)$Letters
    letters_full[names(letters)] <- letters
  }
  
  # 5. Summarize mean and SE per threshold
  summary_df <- df_mod %>%
    group_by(Threshold) %>%
    summarise(
      MeanAUC = mean(AUC, na.rm = TRUE),
      SE_AUC = sd(AUC, na.rm = TRUE)/sqrt(n()),
      .groups = "drop"
    )
  
  # 6. Combine all info
  summary_df <- summary_df %>%
    mutate(
      Model = mod,
      Letters = letters_full[as.character(Threshold)],
      Friedman_p = friedman_p
    )
  
  model_comp_list[[mod]] <- summary_df
}

# Combine all models into one dataframe
model_comparison_df <- bind_rows(model_comp_list)

cv_df <- cv_df %>%
  mutate(Threshold = factor(Threshold, levels = sort(unique(Threshold))))

model_comparison_df <- model_comparison_df %>%
  mutate(Threshold = factor(Threshold, levels = sort(unique(Threshold))))

# Plot full CV distribution with letters
p <- ggplot(cv_df, aes(x = Threshold, y = AUC)) +
  geom_boxplot() +
  # Overlay letters from model_comparison_df
  geom_text(data = model_comparison_df, 
            aes(x = Threshold, y = MeanAUC + 0.07, label = Letters), 
            inherit.aes = FALSE, size = 5, vjust = 0) +
  facet_wrap(~Model) +
  ylim(c(0.50,1))+
  theme_bw(base_size = 11) +
  labs(
    x = "Variable Selection Threshold",
    y = "Cross-validated AUC",
    fill = "Threshold"
  ) +
  theme(legend.position = "none")

p









# ==============================================================================
# STEP 7: SDMs with selected variables 
# ==============================================================================
setwd("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M3")
thresh= 15

vars_selected <- rank_variables(models, types, avg_threshold = 15)
vars_selected <- vars_selected$summary_table$Variable

# vars_selected
# [1] "temperature_summer_mean"  "slope_5x5"                "turbidity_summer_mean"    "ammonium_spring_SD"      
# [5] "PAR_summer_mean"          "currentSpeed_summer_mean" "salinity_summer_SD"      

# 4 Sep 2025
# [1] "slope_5x5"                "nitrate_summer_minimum"   "turbidity_summer_mean"    "nitrate_winter_mean"      "temperature_summer_mean" 
# [6] "ammonium_spring_mean"     "currentSpeed_summer_mean" "PAR_summer_mean"         


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



## GLM
quad_vars_sel <- quad_vars[quad_vars %in% vars_selected]
linear_terms_sel <- setdiff(vars_selected, quad_vars_sel)

terms_quad_sel <- paste0(quad_vars_sel, " + I(", quad_vars_sel, "^2)")
glm_formula <- as.formula(paste("kelp ~", paste(c(terms_quad_sel, linear_terms_sel), collapse = " + ")))
glm_mod_s <- glm(glm_formula, data = train_sel_scaled, family = binomial)
# saveRDS(glm_mod_s, "glm_mod_s.rds")

## GAM
gam_formula <- as.formula(paste("kelp ~", paste0("s(", vars_selected, ")", collapse = " + ")))
gam_mod_s <- gam(gam_formula, data = train_sel_scaled, family = binomial)
# saveRDS(gam_mod_s, "gam_mod_s.rds")

## Random Forest
train_sel$kelp <- as.factor(train_sel$kelp)
rf_mod_s <- randomForest(as.formula(paste("kelp ~", paste(vars_selected, collapse = " + "))),
                         data = train_sel, ntree = 500, importance = TRUE)
# saveRDS(rf_mod_s, "rf_mod_s.rds")

## BRT
train_brt_sel <- train_sel# %>% select(-env_cluster)
colnames(train_brt_sel)
train_brt_sel$kelp <- as.numeric(as.character(train_brt_sel$kelp))
brt_mod_s <- dismo::gbm.step(data = train_brt_sel, 
                             gbm.x = which(names(train_brt_sel) != "kelp"),
                             gbm.y = which(names(train_brt_sel) == "kelp"),
                             family = "bernoulli",
                             tree.complexity = 3,
                             learning.rate = 0.01,
                             bag.fraction = 0.5)

# saveRDS(brt_mod_s, "brt_mod_s.rds")

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

dismo::gbm.plot(brt_mod_s, n.trees = brt_mod_s$gbm.call$best.trees)

# save.image("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M3/results_SDM_M3.RData")




# ==============================================================================
# STEP 8: Predictions & Evaluation
# ==============================================================================
# Get metrics of model performance, AUC, Sensitivity and specificity
results_selected <- bind_rows(
  get_metrics_optimized2(model= glm_mod_s, test_data=train_sel, scale_params= scaling_params_2,
                         model_name = "glm", threshold_type = "youden"),
  get_metrics_optimized2(gam_mod_s, train_sel, scale_params= scaling_params_2, "gam",
                         threshold_type = "youden"),#10pct_omission
  get_metrics_optimized2(rf_mod_s, train_sel, "rf", threshold_type = "youden"),
  get_metrics_optimized2(brt_mod_s, train_sel, "brt", threshold_type = "youden")
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
  "nitrate_summer_minimum"    = expression("Min of Summer Nitrate (uM.L^-1)"),
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
ggplot(summary_table, aes(x = weighted_score, y = Variable, fill = ColorFlag)) +
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
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/VariablesImportance_profile_selectedModels_Sep2025.png", width = 18, height = 12, dpi= 300, units="cm")
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
  "nitrate_summer_minimum"    = "Min of Summer Nitrate (uM.L-1)",
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


# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_Sep2025.png", width = 19, height = 9, dpi= 300, units="cm")
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


# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_average_Sep2025.png", width = 19, height = 9, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/plots/ResponseCurves_SelectedVariables_average_Sep2025.pdf", width = 19, height = 9, dpi= 300, units="cm")




# ### TABLES FOR REPORTING RESULTS ========================================================
# # --- 1. Model performance table ---
# 
# results_selected
# 
# model_perf_table <- results_long_selected %>%
#   select(Model, AUC_lower, AUC_upper, Metric, Value) %>%
#   mutate(across(where(is.numeric), ~round(., 3))) 
# 
# # Convert to flextable
# ft_model <- flextable(model_perf_table) %>%
#   autofit() %>%
#   set_header_labels(
#     Model = "Model",
#     AUC = "AUC",
#     Accuracy = "Accuracy",
#     Sensitivity = "Sensitivity",
#     Specificity = "Specificity",
#     Kappa = "Kappa"
#   )
# 
# # --- 2. Variable importance table ---
# var_imp_table <- res$summary_table %>%
#   select(Variable, n_models, min_importance, avg_importance, max_importance,
#          weighted_score, Rank, Category) %>%
#   mutate(across(c(min_importance, avg_importance, max_importance, weighted_score), ~round(., 2))) %>%
#   arrange(desc(weighted_score))
# 
# ft_var <- flextable(var_imp_table) %>%
#   autofit() %>%
#   set_header_labels(
#     Variable = "Variable",
#     n_models = "Models",
#     min_importance = "Min Imp.",
#     avg_importance = "Avg Imp.",
#     max_importance = "Max Imp.",
#     weighted_score = "Weighted Score",
#     Rank = "Rank",
#     Category = "Category"
#   )
# 
# # --- Export to Word ---
# doc <- read_docx() %>%
#   body_add_par("Model Performance", style = "heading 1") %>%
#   body_add_flextable(ft_model) %>%
#   body_add_par("Variable Importance", style = "heading 1") %>%
#   body_add_flextable(ft_var)
# 
# print(doc, target = "Model_Report_Tables.docx")
































# Conver table to long table for plotting and summarizing data =================

kelp_long <- kelp_data_with_rasters[,c(3,7:ncol(kelp_data_with_rasters))] %>%
  pivot_longer(cols = -kelp, names_to = "variable", values_to = "value")

str(kelp_long)
kelp_long$variable_code<- as.factor(kelp_long$variable)
vars<- levels(kelp_long$variable_code)

# levels(kelp_long$variable_code)<- stringr::word(vars, c(1,2), sep = "_")

# Extract first two parts (or just the first if only one exists)
# short_names <- sapply(strsplit(vars, "_"), function(x) {
#   paste(x[1:min(2, length(x))], collapse = "_")
# })

# Combine into a data frame
var_names_df <- data.frame(
  original_name = vars,
  short_name = short_names,
  stringsAsFactors = FALSE
)

# View result
print(var_names_df)



ggplot(kelp_long, aes(x = as.factor(kelp), y = value)) +
  geom_violin() +
  facet_wrap(~ variable_code, nrow= 3) +
  labs(x = "Kelp (0 = absence, 1 = presence)", y = "Value") +
  theme_minimal()












library(ggpmisc)
combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=depth_with_frame , y= bathymetry_brougton_masked20m, color= as.factor(occ), label = site))+
  geom_point(size= 2)+ geom_smooth(method = lm)+
  scale_color_manual(values= c("brown2", "darkgreen"))+
  labs(x= "Measured Depth (m)", y= "Bathymetry depth (m)", color= "presence (green) and absence (red)")+
  stat_poly_line() +
  ylim(c(-5, 45))+
  stat_poly_eq(use_label(c("eq", "adj.R2", "f", "p", "n"))) +
  # geom_text(size=4)+
  theme_bw()


combinePointVarValue$occ<- as.factor(combinePointVarValue$occ)
# write.csv(combinePointVarValue, "C:/Users/romi_/OneDrive - University of Victoria/Kelp_postdoc/SDMs/csvs_to_SDMs/combined_focalsites_VarValues.csv")


b<- combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=occ , y= bathymetry, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="Depth (m)", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none") 
# theme(legend.position = c(0.2, 0.9), legend.direction="horizontal")


s<- combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=occ , y= slope_20m, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="Slope (degrees)", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none")


a<- combinePointVarValue%>%
  filter(bathymetry<= 50)%>%
  ggplot(aes(x=occ , y= aspect_trig, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="Aspect (degrees)", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none")

r<- combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=occ , y= broughton_roughness_20m, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="Roughness", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none")

curve<- combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=occ , y= broughton_profilecurv_20m, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="Prof. curvature", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none")

tri<- combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=occ , y= broughton_TRI_20m, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="TRI", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none")

tpi<- combinePointVarValue%>%
  filter(bathymetry_brougton_masked20m<= 50)%>%
  ggplot(aes(x=occ , y= broughton_TPI_20m, fill= as.factor(occ)))+
  geom_violin()+
  scale_fill_manual(values= c("brown2", "darkgreen" ))+
  geom_boxplot(width= 0.1)+
  theme_bw()+
  labs(y="TPI", x= "", fill= "Kelp Presence")+
  stat_summary(fun=mean, geom='point', shape=20, size = 3)+
  theme(legend.position = "none")

cowplot::plot_grid(b, s, a, r, curve, tri, tpi, ncol = 3)





# points<- as.data.frame(rasterToPoints(btmax))

# btmax_proj<- terra::project(btmax, "EPSG:3005", method = "near")
# btmax<- projectRaster(btmax, crs= ("+proj=longlat +ellps=WGS84"))

# coordinates(coordinates)= ~ lon + lat

# # Extract raster value by points
# lstValue= raster::extract(LST_layers, coordinates)
# 
# # Combine raster values with point and save as a CSV file.
# combinePointlstValue= cbind(coordinates, lstValue)

# write.csv(combinePointlstValue, "combinePoint_LSTValue_atmoorings.csv")