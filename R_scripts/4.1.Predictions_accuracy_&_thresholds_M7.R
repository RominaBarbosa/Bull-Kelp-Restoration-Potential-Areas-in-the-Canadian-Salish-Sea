###==================================================================
### Species Distribution models        SDMs          ################
### Model Predictions - thresholds and accuracy test ################
###
### C0ntent: 
### 1- Plot of predicted habitat suitability         ################
###  and thresholds delimiting the Unsuitable, Adequate,    #########
###  and Optimal habitat                                    #########
### 2- tab_ensemble_thresholds_M7.csv                       #########
### 3- Model val. with testing and independent datasets     #########
### 4- Model validation AFTER masking by substrate
### Author: Romina Barbosa                           ################
### Date: 24-Sep-2025                                ################
### Last edition: 19-Nov-2025                        ################
###==================================================================
library(tidyr)
library(stringr)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(sf) 

model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
postblob_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob"
source("/Users/romina/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/functions_3.1_HSModeling_analyses_&_plots.R")

### Load model predictions =====================================================
ens_average<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble.tif", sep="/"))
ens_average_t2<- rast(paste(postblob_path, "M7/tifs/ensemble_postblob.tif", sep="/"))
names(ens_average)<- "ens_average"
names(ens_average_t2)<- "ens_average_postblob"




# ==============================================================================
#  Identify thresholds based on training dataset and ensemble model (blob)
# ==============================================================================
# Load training and testing dataset 
traintest<- read.csv(paste(model_results_path,"training_testing_datasets_blob_M7.csv", sep="/"))
# train<- read.csv(paste(model_results_path,"Sep2025_M6_weightedPres/train_selected_table_FINALMODELS.csv", sep="/"))
traintest<- traintest[,-1]
colnames(traintest)

summary(traintest)
traintest%>%
  group_by(substrate, set, period, kelp)%>%
  summarize(n_records= length(set))
#   substrate set   period     kelp n_records
# 1         1 test  2014_2019     0       611
# 2         1 test  2014_2019     1       331
# 3         1 test  2020_2022     1       280
# 4         1 train 2014_2019     0      1426
# 5         1 train 2014_2019     1       852
# 6         1 train 2020_2022     1       574


# cols_selected<- c("kelp", "x", "y", "period", "set", "year")
# traintest<- traintest[,colnames(traintest)%in%cols_selected]
# traintest<- #traintest[,c(1,2,3,17,15,4)]

train<- traintest%>%
  filter(set=="train")

train_points<- vect(train, c("x","y"))
  
train$ensemble_t1pred <- terra::extract(ens_average, train_points)[,2]
train$ensemble_t2pred <- terra::extract(ens_average_t2, train_points)[,2]


### Get thresholds  based on training dataset =====
roc_train <- pROC::roc(response = train$kelp,
                         predictor = train$ensemble_t1pred)
plot(roc_train, col = "green", main = "ROC curve: Ensemble model")


tab_ensemble_thresholds <- get_thresh_table(
  roc_obj = roc_train,
  pred = train[which(!is.na(train$ensemble_t1pred)),"ensemble_t1pred"],
  truth = train[which(!is.na(train$ensemble_t1pred)),"kelp"],
  model_name = "Ensemble_ave"
)

#   Model        Criterion   Threshold Sensitivity Specificity   TSS   AUC
# 1 Ensemble_ave Max TSS         0.473       0.926       0.923 0.848 0.979
# 2 Ensemble_ave No omission  -Inf           1           0     0     0.979
# 3 Ensemble_ave 1% omission     0.241       0.989       0.736 0.726 0.979

# threshold_drop <- quantile(pres_pred, 0.01, na.rm = TRUE)  # 1% presence omission
# M7: 1% ==> 0.2414677 
threshold_drop<- tab_ensemble_thresholds[which(tab_ensemble_thresholds$Criterion=="1% omission"), "Threshold"]
threshold_drop<- threshold_drop$Threshold  

thr_maxTSS<- tab_ensemble_thresholds[which(tab_ensemble_thresholds$Criterion=="Max TSS"), "Threshold"]
thr_maxTSS<- thr_maxTSS$Threshold  
# write.csv(tab_ensemble_thresholds, "Sep2025_M7_weightedPres/tab_ensemble_thresholds_M7_FINAL.csv")


### Plot thresholds  - FIGURE SX             ====
# Predicted suitability for presences in training data
train_ave_ens<- train

# Density estimates
pres_pred <- train_ave_ens$ensemble_t1pred[train_ave_ens$kelp == 1]
abs_pred  <- train_ave_ens$ensemble_t1pred[train_ave_ens$kelp == 0]

# Compute density estimates
dens_pres <- density(pres_pred, na.rm = TRUE)
dens_abs  <- density(abs_pred,  na.rm = TRUE)

# Build data frame for ggplot
df_plot <- data.frame(
  x = c(dens_abs$x, dens_pres$x),
  y = c(dens_abs$y, dens_pres$y),
  group = rep(c("Absences", "Presences"), 
              times = c(length(dens_abs$x), length(dens_pres$x)))
)


# Define label positions (x = middle of each region, y = top 90% of max density)
y_max <- 3.7#max(df_plot$y, na.rm = TRUE)
labels <- data.frame(
  label = c("Unsuitable", "Adequate", "Optimal"),
  x = c(threshold_drop / 2,
        (threshold_drop + thr_maxTSS) / 2,
        thr_maxTSS + (max(df_plot$x) - thr_maxTSS)/2),
  y = rep(0.9 * y_max, 3)
)
# labels[1,2]<- 

ggplot(df_plot, aes(x = x, y = y, color = group, fill = group)) +
  # Shaded suitability areas
  # geom_rect(aes(xmin = -Inf, xmax = threshold_drop, ymin = -Inf, ymax = Inf),
  #           fill = "blue", alpha = 0.01, inherit.aes = FALSE) +  # unsuitable
  # geom_rect(aes(xmin = threshold_drop, xmax = thr_maxTSS, ymin = -Inf, ymax = Inf),
  #           fill = "orange", alpha = 0.01, inherit.aes = FALSE) +      # moderate
  # geom_rect(aes(xmin = thr_maxTSS, xmax = Inf, ymin = -Inf, ymax = Inf),
  #           fill = "red", alpha = 0.01, inherit.aes = FALSE) +       # optimal
  geom_line(size = 1.2) +
  # geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.2, color = NA) + # adds semi-transparent fill
  geom_vline(xintercept = threshold_drop, color = "orange", linetype = "dashed", size = 1) +
  geom_vline(xintercept = thr_maxTSS, color = "red", linetype = "dashed", size = 1) +
  # geom_vline(xintercept = thr_maxF1, color = "darkorange", linetype = "dashed", size = 1) +
  # Labels for shaded regions
  lims(y=c(0, 3.5), x=c(0, 1))+
  geom_text(data = labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, color = "black", fontface = "bold") +
  labs(
    x = "Predicted Kelp Habitat Suitability",
    y = "Density",
    color = "",
    fill = ""
  ) +
  scale_color_manual(values = c("black", "grey")) +  # absences = blue, presences = red
  scale_fill_manual(values = c("black", "grey")) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.8, 0.7),
    legend.background = element_rect(fill = alpha("white", 0.1), color = NA),
    legend.key = element_blank()
  ) +
  annotate("text", x = threshold_drop, y = max(df_plot$y)*0.85,
           label = "0.1% Omission", color = "orange", angle = 90, vjust = -0.5) +
  annotate("text", x = thr_maxTSS, y = max(df_plot$y)*0.85,
           label = "Max TSS", color = "red", angle = 90, vjust = -0.5)

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Threshold_selection_plot_M7.pdf", width = 12, height = 8, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/Threshold_selection_plot_M7.png", width = 12, height = 8, dpi= 300, units="cm")



# ==============================================================================
#      Convert predictions in categories based on thresholds and 
#     evaluate its accuracy in both periods: 2014-2019 and 2020-2022
# ==============================================================================
## Testing points from Period 1 =============================
test<- traintest%>%
  filter(set=="test")

# Testing period 1
test_t1_a<- test%>%filter(period == "2014_2019")

test_t1_a%>%
  group_by(kelp, substrate)%>%
  summarize(n= length(kelp))
#    kelp substrate     n
# 1     0         1   611
# 2     1         1   331


test_t1_hardSubstrate<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&rocky.csv")
test_t1_hardSubstrate%>%
  group_by(kelp, substrate)%>%
  summarize(n= length(kelp))
# kelp substrate     n
# 1     0         1  2604
# 2     1         1  1190

colnames(train)
colnames(test_t1_hardSubstrate)

test_t1_hardSubstrate <- test_t1_hardSubstrate %>%
  left_join(train %>% select(x, y, kelp, depth,set),
            by = c("x", "y", "kelp", "depth")) %>%
  mutate(set = ifelse(is.na(set), "test", set))

test_t1_hardSubstrate%>%
  group_by(kelp, set)%>%
  summarize(n= length(kelp))
# kelp set       n
# 1     0 test   1178
# 2     0 train  1426
# 3     1 test    338
# 4     1 train   852

test_t1_hardSubstrate<- test_t1_hardSubstrate%>%
  filter(set== "test")
# kelp set       n
# 1     0 test   1178
# 2     1 test    338


test_t1_softSubstrate<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_data_2014_to_2019_filtered_depth&softSubstrate.csv")
test_t1_softSubstrate%>%
  group_by(kelp, substrate)%>%
  summarize(n= length(kelp))
# kelp substrate     n
# 1     0         2  2497
# 2     1         2    79


# # Downsample absences to balance count (same n of presences and absences)
# presences <- test_t1_softSubstrate %>% filter(kelp == 1)
# absences <- test_t1_softSubstrate %>% filter(kelp == 0) %>% sample_n(nrow(presences)*2)
# 
# test_t1_balanced <- bind_rows(presences, absences)
# test_t1_balanced%>%
#   group_by(kelp, substrate)%>%
#   summarize(n= length(kelp))
# # kelp substrate     n
# # 1     0         2   158
# # 2     1         2    79


# Merge all records of testing period 1
colnames(test_t1_softSubstrate)
test_t1_hardSubstrate$substrate<- "1"
colnames(test_t1_hardSubstrate)

test_t1<- rbind(test_t1_softSubstrate[,which(colnames(test_t1_softSubstrate) %in% colnames(test_t1_hardSubstrate))],
                test_t1_hardSubstrate[,which(colnames(test_t1_hardSubstrate) %in% colnames(test_t1_softSubstrate))])

test_t1%>%  # total 4092 records
  group_by(kelp, substrate)%>%
  summarize(n_records= length(kelp))
# kelp substrate n_records
# 1     0 1              1178
# 2     0 2              2497
# 3     1 1               338
# 4     1 2                79


## Testing points from Period 2 =============================
# Kelp presence records from period 2 were filtered by hard substrate and mixed 
# with period 1 and split in 70/30, whereas all absences were excluded from calibration

test<- traintest%>%
  filter(set=="test")

test%>%group_by(set, period, substrate, kelp)%>%summarize(n_records= length(kelp))
#   set   period    substrate  kelp n_records
# 1 test  2014_2019 1             0       611
# 2 test  2014_2019 1             1       331
# 3 test  2020_2022 1             1       280

test_period_2a<- traintest%>%
  filter(set=="test")%>%
  filter(period=="2020_2022")

traintest%>%
  filter(period=="2020_2022")%>%
  group_by(set, kelp)%>%
  summarize(n_records= length(kelp))


280/(280+574) # 0.32786 % percentage of presences for testing

# thus, Dataset from Period 2: presences in soft substrate and all absences (soft and hard substrate)
# will be used to test the model
test_period_2b<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/MoraSoto_postblob/Presence_absences_kelp_2020_2022_filtered.csv")
colnames(test_period_2b)

test_period_2b%>%
  group_by(kelp, substrate)%>%
  summarize(n_records= length(kelp))
# kelp substrate n_records
# 1     0         1      1795
# 2     0         2      1900
# 3     1         1       863
# 4     1         2        46


# Remove presence records in hard substrate because were used to mix for calibration
test_period_2b <- test_period_2b %>%
  filter(kelp == "0" | (kelp == "1" & substrate == "2"))

test_period_2b%>%
  group_by(kelp, substrate)%>%
  summarize(n_records= length(kelp))
#    kelp substrate n_records
# 1     0         1      1795
# 2     0         2      1900
# 3     1         2        46

test_period_2b$period<- "2020_2022"
test_period_2b$set<- "test"

# Merge all testing records from period 2
test_period_2a$substrate<- "1"

test_period_2b<- test_period_2b[,which(colnames(test_period_2b) %in% colnames(test_period_2a))]
colnames(test_period_2b)
colnames(test_period_2a)

# Merge all records of testing period 2
test_t2<- rbind(test_period_2a[,which(colnames(test_period_2a) %in% colnames(test_period_2b))], test_period_2b)
test_t2%>%  # 4021 
  group_by(kelp, substrate)%>%
  summarize(n_records= length(kelp))
#    kelp substrate n_records
# 1     0 1              1795
# 2     0 2              1900
# 3     1 1               280
# 4     1 2                46

summary(as.factor(test_t2$set)) # test 4021 records in total

# ==============================================================================
# Merge and save table of Testing records from both Periods (from satellite)
# ==============================================================================
colnames(test_t2)
test_t1$gridID<- NA
colnames(test_t1)
test_t1$period<- "2015_2019"
test_t1$set<- "test"

testing_dataset_satellite<- rbind(test_t1[,which(colnames(test_t1)%in%colnames(test_t2))], test_t2)
# write.csv(testing_dataset_satellite, "/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/MoraSoto_postblob/Presence_absences_testing_dataset_satellite_FINAL.csv")

# ==============================================================================
# Testing both Periods - Binary results & Evaluation with presence/absence data from satellite
# ==============================================================================
library(terra)
library(dismo)   # for evaluate()
library(pROC)    # for AUC
library(biomod2) # for Boyce Index if needed
library(dplyr)

### Calculate Prevalence of the kelp in the area training + testing data ==========
traintest%>%
  # filter(set=="train")%>%
  filter(kelp==1)%>%
  summarize(bathy= quantile(bathymetry_20, 0.97)) #15.62306

records_total_t<- traintest%>%
  filter(set=="train")%>%
  group_by(set, kelp)%>%
  summarize(n_records= length(kelp))
# set    kelp n_records
# 1 train     0      1426
# 2 train     1      1426

records_total_test<- testing_dataset_satellite%>%
  filter(set=="test")%>%
  group_by(set, kelp)%>%
  summarize(n_records= length(kelp))
# set    kelp n_records
# 1 test      0      7370
# 2 test      1       743

records_total<- rbind(records_total_t, records_total_test)
records_total<-  testing_dataset_satellite%>%
  group_by(kelp)%>%
  summarize(n_records= length(kelp))

prev_PeriodsCombined<- records_total[which(records_total$kelp==1),2]/sum(records_total[,2])
# 0.09158141

# separated prevalence from period 1
records_total_t1<-  testing_dataset_satellite%>%
  filter(period=="2015_2019")%>%
  group_by(kelp)%>%
  summarize(n_records= length(kelp))

prev_Period1<- records_total_t1[which(records_total_t1$kelp==1),2]/sum(records_total_t1[,2])
# 0.1019062


# separated prevalence from period 2
records_total_t2<-  testing_dataset_satellite%>%
  filter(period=="2020_2022")%>%
  group_by(kelp)%>%
  summarize(n_records= length(kelp))

prev_Period2<- records_total_t2[which(records_total_t2$kelp==1),2]/sum(records_total_t2[,2])
# 0.08107436


### Calibrate model predictions and thresholds  (save tif files) =================================
# calibration function
# calibrate_prob <- function(pred_vals, prev_train, prev_target) {
#   # clip to avoid 0/1
#   p <- pmin(pmax(pred_vals, 1e-8), 1 - 1e-8)
#   odds <- p / (1 - p)
#   odds_adj <- odds * ((prev_target / (1 - prev_target)) / (prev_train / (1 - prev_train)))
#   p_cal <- odds_adj / (1 + odds_adj)
#   return(p_cal)
# }

# ---- inputs ----
prev_train <- 0.5                     # training prevalence used to train models
threshold_raw <- threshold_drop           # 1% drop threshold = 0.2414677, derived from raw training predictions
ensemble_raster_p1<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble.tif", sep="/"))
ensemble_raster_p2<- rast(paste(postblob_path, "M7/tifs/ensemble_postBlob.tif", sep="/"))


# prevalences to test
prev_vals <- seq(0.07, 0.11, 0.01)  # add prev_all if available

# function to compute calibrated threshold from raw threshold
calc_calibrated_threshold <- function(t_raw, prev_train, prev_target) {
  calibrate_prob(t_raw, prev_train, prev_target)
}

# Function to compute area above threshold (km² and % of valid area)
area_above_threshold <- function(rast_obj, thresh) {
  # ensure raster is in projected coordinates (meters)
  cell_area_m2 <- abs(prod(res(rast_obj)))
  
  # identify valid (non-NA) cells
  valid_cells <- !is.na(values(rast_obj))
  n_valid <- sum(valid_cells)
  
  # identify cells above threshold among valid ones
  pred_vals <- values(rast_obj)[valid_cells]
  n_above <- sum(pred_vals >= thresh)
  
  # compute total and proportional areas
  total_area_km2 <- (n_valid * cell_area_m2) / 1e6
  area_above_km2 <- (n_above * cell_area_m2) / 1e6
  perc_area <- (n_above / n_valid) * 100
  
  # return both absolute and relative results
  return(data.frame(
    Threshold = thresh,
    Area_km2 = area_above_km2,
    Total_km2 = total_area_km2,
    Perc_valid = perc_area
  ))
}


# load rasters
r1 <- (ensemble_raster_p1)
r2 <- (ensemble_raster_p2)

# iterate prevalences
results <- expand.grid(Period = c("Period1","Period2"), Prevalence = prev_vals, stringsAsFactors = FALSE)
results$Threshold_cal <- NA
results$Area_km2 <- NA

t_cal <- calc_calibrated_threshold(threshold_raw, prev_train, prev_PeriodsCombined)$n_records
# 0.03109478

for (i in seq_len(nrow(results))) {
  pv <- results$Prevalence[i]
  results$Threshold_cal[i] <- t_cal

  if (results$Period[i] == "Period1") {
    r_cal <- app(r1, function(x) calibrate_prob(x, prev_train, pv))
    # Binary raster
    r_bin <- app(r_cal, function(x) as.numeric(x >= t_cal))
    # writeRaster(r_bin, paste("/Volumes/Romina_PSF/PSF/SDM/SDM_results/calibrate_predictions/ensemble_t1", "calibrated", pv, ".tif", sep="_"),overwrite=TRUE)
    results$Area_km2[i] <- area_above_threshold(r_cal, t_cal)$Area_km2
    # optionally remove r_cal from memory after use
    rm(r_cal); gc()
  } else {
    r_cal <- app(r2, function(x) calibrate_prob(x, prev_train, pv))
    # Binary raster
    r_bin <- app(r_cal, function(x) as.numeric(x >= t_cal))
    # writeRaster(r_bin, paste("/Volumes/Romina_PSF/PSF/SDM/SDM_results/calibrate_predictions/ensemble_t2", "calibrated", pv, ".tif", sep="_"),overwrite=TRUE)
    results$Area_km2[i] <- area_above_threshold(r_cal, t_cal)$Area_km2
    rm(r_cal); gc()
  }
}




# plot comparison
# results_0.07_0.12<- results
results$Period <- as.factor(results$Period)
results$Prevalence <- as.numeric(results$Prevalence)  # ensure numeric
results$Area_km2 <- as.numeric(unlist(results$Area_km2))


ggplot(results, aes(x = Prevalence, y = Area_km2, color = Period, group = Period)) +
  # geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(
    x = "Prevalence",
    y = expression("Area above threshold (km"^2*")"),
    title = "Sensitivity to prevalence calibration",
    color = "Period"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/calibrate_predictions/sensitivity_test_prevalence.pdf",
#        width = 10, height = 9, units="cm", dpi = 300)




### Sensitive analysis: predictions calibration/scaling with diff. prevalence (save tif files and SD among results) ====
# range of prevalence values from 0.07 to 0.11, every 0.01




# Load process_period and calc_metrics functions
source("~/Documents/GitHub/PSF_kelp-HSM/R_scripts/functions/models_validation_metrics_function.R")

# -----------------------------
#  Parameters
# -----------------------------
prev_train <- 0.5
prev_PeriodsCombined<- prev_PeriodsCombined$n_records
prevalence=prev_PeriodsCombined
thresholds <- list(
  drop1 = 0.2414677,   # threshold_drop (1% drop)
  maxTSS = 0.4733616   # thr_maxTSS
)

# -----------------------------
#  Process both periods
# -----------------------------
# Period 1
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p1 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/renamed"
metrics_p1 <- process_period(raster_folder_p1, test_t1, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p1$Period <- "Period1"

metrics_p1_v2 <- process_period_v2(raster_folder_p1, test_t1, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p1_points<- metrics_p1_v2$Classification
metrics_p1_points$Period <- "Period1"
metrics_p1_points$dataset<- "30pct_satellite"
  
# Period 2
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed"
metrics_p2 <- process_period(raster_folder_p2, test_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p2$Period <- "Period2"


metrics_p2_v2 <- process_period_v2(raster_folder_p2, test_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p2_points<- metrics_p2_v2$Classification
metrics_p2_points$Period <- "Period2"
metrics_p2_points$dataset<- "30pct_satellite"


metrics_points_satellite<- rbind(metrics_p2_points, metrics_p1_points)
# write.csv(metrics_points_satellite, "/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_points_satellite.csv")

# -----------------------------
#  Combine results
# -----------------------------
metrics_all <- bind_rows(metrics_p1, metrics_p2) %>%
  dplyr::select(Period, Model, Threshold, everything())

print(metrics_all)
metrics_all$dataset<- "satellite_testing"

# -------------------------------------------------------------------------------
#  Save Table of results Model validation without masking by substrate
# -------------------------------------------------------------------------------
# write.csv(metrics_all, "/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_satellite_evaluation_results_without_masking.csv")

# Boyce_point evaluates model discrimination at the testing points, using presence vs absence scores.
# Boyce_raster evaluates the ecological and spatial realism of the model, testing whether presences fall in the high-suitability regions of the projected map.
# These two indices often differ because point-based evaluation can show good discrimination 
# even when the spatial prediction surface is flat, overextended, or poorly structured. 
# Boyce_raster is generally considered the more ecologically meaningful measure for SDM projection maps, but testing points should represent the entire area

metrics_summary_thresInd <- metrics_all %>%
  dplyr::select(Period, Model, AUC, Boyce_point, Boyce_raster) %>%
  distinct() %>%  
  arrange(Period, Model)%>%
  mutate(
    AUC = round(AUC, 2),
    Boyce_point = round(Boyce_point, 2)
  )

# write.csv(metrics_summary_thresInd, "/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_threshold_independent.csv", row.names = FALSE)
# write.csv(metrics_summary_thresInd,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_threshold_independent_satellite_Periods1&2.csv")

metrics_summary <- metrics_all %>%
  dplyr::select( Period,Model, Threshold, Accuracy, TSS, Sensitivity, Specificity) %>%
  arrange(Period, Model, Threshold)%>%
  mutate(
    Accuracy = round(Accuracy, 2),
    TSS = round(TSS, 2),
    Sensitivity = round(Sensitivity, 2),
    Specificity = round(Specificity, 2)
  )

# write.csv(metrics_summary,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_threshold_dependent_satellite_Periods1&2.csv")




# -------------------------------------------------------------------------------
#  Independent datasets for validation of Period 2
# -------------------------------------------------------------------------------
# Mapillary dataset -------------------------------------------------------------
# Capillary dataset has only presences, thus, AUC, Boyce and TSS, as well as the Specificity is not possible to be calculated

# Load raster
ensemble_raster_p1<- rast(paste(model_results_path, "Sep2025_M7_weightedPres/tifs/ensemble.tif", sep="/"))
ensemble_raster_p2<- rast(paste(postblob_path, "M7/tifs/ensemble_postBlob.tif", sep="/"))
capillary_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/fieldwork_2025/KelpData_BoatTransect_CourtneyNorth_GitaN_xy3005.csv")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed"

# Load buffer polygons
# buf <- vect("/Volumes/Romina_PSF/PSF/SDM/fieldwork_2025/Kelp_Surveys_Shapefile/mapillary_kelp_presences_bufferPolygons.shp")

capillary_df$kelp<- 1
colnames(capillary_df)[11:12]<- c("x", "y") 

# Filter records
min_dist <- 180 * 2  # meters
pts_sf <- st_as_sf(capillary_df, coords = c("x","y"), crs = 3005)
keep <- rep(TRUE, nrow(pts_sf))  # start by keeping all

for (i in seq_len(nrow(pts_sf))) {
  if (!keep[i]) next
  dists <- st_distance(pts_sf[i,], pts_sf[(i+1):nrow(pts_sf),])
  too_close <- which(as.numeric(dists) < min_dist)
  if (length(too_close) > 0) keep[i + too_close] <- FALSE
}

capillary_df_t2_filtered <-capillary_df[ which(capillary_df$Presence.b %in% pts_sf[keep, ]$Presence.b),]
length(capillary_df_t2_filtered$Presence.b) # 35
length(pts_sf$Presence.b) # 147 --> 112 were removed


model_names <- c("brt", "ensemble", "gam", "glm", "rf")
buffer_capillary_pred_df <- build_buffer_pred_df(
  model_names=model_names,
  xy_df = capillary_df_t2_filtered,     # your table of ID, kelp, x_3005, y_3005
  raster_folder = raster_folder_p2,
  buffer_radius = 180,
  threshold = thresholds$drop1[1]        # example threshold
)


model_names <- c("brt", "ensemble", "gam", "glm", "rf")
eval_capillary <- process_buffer_predictions(
  buffer_pred_df = buffer_capillary_pred_df,   # your averaged predictions per ID
  thresholds = thresholds,                    # named list e.g. list(opt = 0.3, maxSens = 0.25)
  model_names = c("brt", "ensemble", "gam", "glm", "rf")
)

eval_capillary_metrics<- eval_capillary$Metrics
eval_capillary_metrics<- eval_capillary_metrics[,c(1,2,6:9)]
eval_capillary_metrics$period <- "Period2"
eval_capillary_metrics$dataset<- "Capillary_2024"

metrics_capillary_p2_points<- eval_capillary$Classification
metrics_capillary_p2_points



# Aerial and shore images dataset -------------------------------------------------------------
test_images_t2<- read.csv("/Volumes/Romina_PSF/PSF/SDM/fieldwork_2025/KelpData_fromshore_&_aerial_2025_NicoleC_correct_xyKelp.csv", sep=";")
colnames(test_images_t2)[12:13]<- c("x", "y")
colnames(test_images_t2)[3]<- "kelp"
summary(as.factor(test_images_t2$kelp))
# 0  1 (absence presence)
# 3 34 

# Filter records
min_dist <- 60  # meters
pts_sf <- st_as_sf(test_images_t2, coords = c("x","y"), crs = 3005)
keep <- rep(TRUE, nrow(pts_sf))  # start by keeping all

for (i in seq_len(nrow(pts_sf))) {
  if (!keep[i]) next
  dists <- st_distance(pts_sf[i,], pts_sf[(i+1):nrow(pts_sf),])
  too_close <- which(as.numeric(dists) < min_dist)
  if (length(too_close) > 0) keep[i + too_close] <- FALSE
}

test_images_t2_filtered <-test_images_t2[ which(test_images_t2$picture %in% pts_sf[keep, ]$picture),]
length(test_images_t2_filtered$picture)
length(test_images_t2$picture) # 3 were removed

# create buffer area 
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed"

buffer_images_pred_df <- build_buffer_pred_df(
  model_names= model_names,
  xy_df = test_images_t2_filtered,     # your table of ID, kelp, x_3005, y_3005
  raster_folder = raster_folder_p2,
  buffer_radius = 30,
  threshold = thresholds$drop1[1]        # example threshold
)

# calculate metrics of evaluation
eval_images <- process_buffer_predictions(
  buffer_pred_df = buffer_images_pred_df,   # your averaged predictions per ID
  thresholds = thresholds,                    # named list e.g. list(opt = 0.3, maxSens = 0.25)
  model_names = c("brt", "ensemble", "gam", "glm", "rf")
)

metrics_images_p2<- eval_images$Metrics[,c(1,2,6:9)]
metrics_images_p2$period <- "Period2"
metrics_images_p2$dataset<- "Images_2025"
head(metrics_images_p2)

metrics_images_p2_points<- eval_images$Classification


# Kayak surveys dataset -------------------------------------------------------------
test_kayak_t2<- read.csv("/Volumes/Romina_PSF/PSF/SDM/fieldwork_2025/Kelp_Surveys_3_centroidxy3005.csv")
colnames(test_kayak_t2)[16:17]<- c("x", "y")
test_kayak_t2$kelp<- test_kayak_t2$Bull_Kelp +test_kayak_t2$Giant_Kelp # presence of any kelp (bull or giant kelp)

summary(as.factor(test_kayak_t2$Bull_Kelp))
# 0  1 (absence presence)
# 88 10 
summary(as.factor(test_kayak_t2$Giant_Kelp))
# 0  1 (absence presence)
# 96  2 

summary(as.factor(test_kayak_t2$kelp))
# 0  1 (absence presence)
# 86 12

model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed"
eval_kayak_p2 <- process_period_v2(raster_folder_p2, test_kayak_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_kayak_p2<- eval_kayak_p2$Metrics
metrics_kayak_p2$period <- "Period2"
metrics_kayak_p2$dataset<- "kayak_2025"


metrics_kayak_p2_points<- eval_kayak_p2$Classification
metrics_kayak_p2_points$Classification<- as.factor(metrics_kayak_p2_points$Classification)


# Kayak plus images evaluation =====
colnames(test_images_t2)
test_images_t2_filtered$dataset<- "Images_2025"
test_kayak_t2$dataset<- "Kayak_2025"
colnames(test_kayak_t2)

dataset_2025<- rbind(test_kayak_t2[, which(colnames(test_kayak_t2)%in% colnames(test_images_t2_filtered))], test_images_t2_filtered[, which(colnames(test_images_t2_filtered)%in% colnames(test_kayak_t2))])
summary(as.factor(dataset_2025$kelp))
# 0  1 (absence presence)
# 89 43

model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed"
buffer_pred_df <- build_buffer_pred_df(
  model_names= model_names,
  xy_df = dataset_2025,     # your table of ID, kelp, x_3005, y_3005
  raster_folder = raster_folder_p2,
  buffer_radius = 30,
  threshold = thresholds$drop1[1]        # example threshold
)

eval_dataset_2025 <- process_buffer_predictions(
  buffer_pred_df = buffer_pred_df,   # your averaged predictions per ID
  thresholds = thresholds,                    # named list e.g. list(opt = 0.3, maxSens = 0.25)
  model_names = c("brt", "ensemble", "gam", "glm", "rf")
)

metrics_dataset_2025_p2<- eval_dataset_2025$Metrics[,c(1,2,6:9)]
metrics_dataset_2025_p2$period <- "Period2"
metrics_dataset_2025_p2$dataset<- "dataset_2025"

metrics_dataset_2025_p2_points<- eval_dataset_2025$Classification
metrics_dataset_2025_p2_points$Classification<- as.factor(metrics_dataset_2025_p2_points$Classification)


# Maine Island Conservancy Society surveys dataset -------------------------------------------------------------
setwd("/Volumes/Romina_PSF/PSF/SDM/fieldwork_2025/Maine_Island_Conservancy/Shapefiles_Mayne")
dir()
MICS2019<- read.csv("points_2019_MayneC.csv")
MICS2019$year<- 2019

MICS2020<- read.csv("points_2020_MayneC.csv")
MICS2020$year<- 2020
sf_MICS2020 <- st_as_sf(MICS2020, coords = c("x", "y"), crs = 26910)
sf_MICS2020_3005 <- st_transform(sf_MICS2020, 3005)# Transform to EPSG:3005
coords <- st_coordinates(sf_MICS2020_3005)
sf_MICS2020_3005$x <- coords[,1]
sf_MICS2020_3005$y <- coords[,2]
MICS2020<- as.data.frame(sf_MICS2020_3005)
MICS2020<- MICS2020[,-7]

MICS2021<- read.csv("points_2021_MayneC.csv")
MICS2021$year<- 2021

MICS2022<- read.csv("points_2022_MayneC.csv")
MICS2022$year<- 2022
sf_MICS2022 <- st_as_sf(MICS2022, coords = c("longitude", "latitude"), crs = 26910)
sf_MICS2022_3005 <- st_transform(sf_MICS2022, 3005)# Transform to EPSG:3005
coords <- st_coordinates(sf_MICS2022_3005)
sf_MICS2022_3005$x <- coords[,1]
sf_MICS2022_3005$y <- coords[,2]
MICS2022<- as.data.frame(sf_MICS2022_3005)
MICS2022<- MICS2022[,-7]


MICS_dataset<- rbind(MICS2019, MICS2020, MICS2021, MICS2022)
MICS_dataset$dataset<- "MICS"
MICS_dataset$kelp<- 1
MICS_dataset<- na.exclude(MICS_dataset)

# write.csv(MICS_dataset, "MICS_dataset_EPSG3005.csv")
MICS_dataset<- read.csv("MICS_dataset_EPSG3005.csv")

sf_MICS_dataset<- st_as_sf(MICS_dataset, coords = c("x", "y"), crs = 3005)
  
# Extract cell ID for each point
sf_MICS_dataset$cell_id <- raster::cellFromXY(ensemble_raster_p2, sf::st_coordinates(sf_MICS_dataset))

# Keep one point per raster cell
sf_MICS_dataset_unique <- sf_MICS_dataset %>%
  group_by(cell_id) %>%
  slice(1) %>% 
  ungroup()

summary(as.factor(sf_MICS_dataset_unique$kelp))# 1243 presences 


coords <- sf::st_coordinates(sf_MICS_dataset_unique)

test_MICS_t2 <- cbind(
  sf_MICS_dataset_unique |> st_drop_geometry(),
  x = coords[,1],
  y = coords[,2]
)

head(test_MICS_t2)

model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed"
eval_MICS_p2 <- process_period_v2(raster_folder_p2, test_MICS_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_MICS_p2<- eval_MICS_p2$Metrics
metrics_MICS_p2$period <- "Period2"
metrics_MICS_p2$dataset<- "MICSIndependent_dataset"


metrics_MICS_p2_points<- eval_MICS_p2$Classification
metrics_MICS_p2_points$Classification<- as.factor(metrics_MICS_p2_points$Classification)


# Merge all classification metrics
colnames(metrics_all)[1]<- "period"
metrics_evaluation_all_datasets<- rbind(metrics_all[which(colnames(metrics_all)%in% colnames(metrics_kayak_p2))],
                                        metrics_kayak_p2, 
                                        metrics_images_p2, 
                                        eval_capillary_metrics, 
                                        metrics_dataset_2025_p2,
                                        metrics_MICS_p2)

metrics_evaluation_all_datasets_period2<- rbind(
                                        eval_capillary_metrics, 
                                        metrics_dataset_2025_p2,
                                        metrics_MICS_p2)


# write.csv(metrics_evaluation_all_datasets_period2,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_evaluation_all_datasets_period2.csv")
# write.csv(metrics_evaluation_all_datasets,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_evaluation_all_datasets.csv")

metrics_dataset_2025_p2_points$dataset<- "2025_dataset"
metrics_MICS_p2_points$dataset<- "Independent_MICSperiod2"
metrics_capillary_p2_points$dataset<- "2024_dataset"

evaluation_points_all_datasets_period2<- rbind(metrics_dataset_2025_p2_points,
      metrics_MICS_p2_points, 
      metrics_capillary_p2_points)
    
# write.csv(evaluation_points_all_datasets_period2,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/points_evaluation_matrix_Independent_datasets_period2.csv")






# -------------------------------------------------------------------------------
#  4- Model validation AFTER masking by substrate
# -------------------------------------------------------------------------------
###  Mask ensemble model predictions by substrate  =======
# source("create_substrate_mask.R")

# ens_averaget1_masked <- terra::mask(ens_average, substrate_aligned)
# names(ens_averaget1_masked)<- "ens_average_t1"
# ens_averaget2_masked <- terra::mask(ens_average_t2, substrate_aligned)
# names(ens_averaget2_masked)<- "ens_average_t2"
# writeRaster(ens_averaget1_masked, "ens_average_masked_M7.tif", overwrite=T)
# writeRaster(ens_averaget2_masked, "ens_average_masked_M7.tif", overwrite=T)

#---------------------------------------------------------
# INPUTS
#---------------------------------------------------------
pred_path   <- raster_folder_p1          # folder with prediction rasters
out_path    <- paste(raster_folder_p1, "tifs_masked", sep="/")      # folder to save masked versions
pred_path_t2   <- raster_folder_p2          # folder with prediction rasters
out_path_t2    <- paste(raster_folder_p2, "tifs_masked", sep="/")      # folder to save masked versions
substrate_r <- paste(model_results_path, "substrate_SOG_aligned_.tif", sep="/")    # 1 = suitable, 2/other = unsuitable
depth_r     <- "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif"
depth_threshold <- 15             # example

# dir.create(out_path_t2, showWarnings = FALSE)

#---------------------------------------------------------
# 1. LOAD RASTERS FOR CONDITIONS
#---------------------------------------------------------
r_sub <- rast(substrate_r)          # substrate raster
r_dep <- rast(depth_r)              # depth raster

#---------------------------------------------------------
# 2. BUILD BINARY MASK (1=suitable, 0=unsuitable)
#---------------------------------------------------------

# substrate suitability
m_sub <- r_sub == 1            

# depth suitability
m_dep <- r_dep < depth_threshold   

# combine
mask_binary <- m_sub & m_dep

# convert from logical (TRUE/FALSE) to numeric (1/0)
mask_binary <- as.numeric(mask_binary)

# SAVE mask for reference
# writeRaster(mask_binary, 
#             file.path(out_path, "suitability_mask.tif"), 
#             overwrite = TRUE)

#---------------------------------------------------------
# 3. MASK EACH MODEL RASTER
#---------------------------------------------------------
pred_files <- list.files(pred_path_t2, pattern = "\\.tif$", full.names = TRUE)

for (f in pred_files) {
  
  cat("Masking:", basename(f), "\n")
  
  r <- rast(f)
  
  # apply the mask (multiply instead of terra::mask to retain spatial coverage)
  r_masked <- r * mask_binary        # suitable areas keep values, unsuitable become 0
  
  out_file <- file.path(out_path_t2, paste0("masked_", basename(f)))
  
  writeRaster(r_masked, out_file, overwrite = TRUE)
}

cat("\n✔️ DONE: All masked rasters saved to:", out_path, "\n")



#---------------------------------------------------------
# Test predictions with substrate mask 
#---------------------------------------------------------
# Using satellite dataset (3~0% split plus absences period-2)
# Period 1
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p1 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_results/Sep2025_M7_weightedPres/tifs/renamed/tifs_masked"
metrics_p1_maskedpreds <- process_period(raster_folder_p1, test_t1, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p1_maskedpreds$Period <- "Period1"

metrics_p1_maskedpreds_v2 <- process_period_v2(raster_folder_p1, test_t1, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p1_maskedpreds_points<- metrics_p1_maskedpreds_v2$Classification
metrics_p1_maskedpreds_points$Period <- "Period1"
metrics_p1_maskedpreds_points$dataset<- "30pct_satellite"

# Period 2
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed/tifs_masked"
metrics_p2_maskedpreds <- process_period(raster_folder_p2, test_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p2_maskedpreds$Period <- "Period2"
metrics_p2_maskedpreds

metrics_p2_maskedpreds_v2 <- process_period_v2(raster_folder_p2, test_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_p2_maskedpreds_points<- metrics_p2_maskedpreds_v2$Classification
metrics_p2_maskedpreds_points$Period <- "Period2"
metrics_p2_maskedpreds_points$dataset<- "30pct_satellite"


## Save results validation with making predictions by substrate and depth
metrics_maskedpreds_satellite<- rbind(metrics_p1_maskedpreds, metrics_p2_maskedpreds)
# write.csv(metrics_maskedpreds_satellite, "/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_satellite_evaluation_results_maskedpreds.csv")

metrics_points_maskedpreds_satellite<- rbind(metrics_p2_maskedpreds_points, metrics_p1_maskedpreds_points)
# write.csv(metrics_points_maskedpreds_satellite, "/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/points_satellite_maskedpreds.csv")

metrics_summary_thresInd_maskedpreds <- metrics_maskedpreds_satellite %>%
  dplyr::select(Period, Model, AUC, Boyce_point, Boyce_raster) %>%
  distinct() %>%  
  arrange(Period, Model)%>%
  mutate(
    AUC = round(AUC, 2),
    Boyce_point = round(Boyce_point, 2)
  )

# write.csv(metrics_summary_thresInd_maskedpreds, "/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_threshold_independent_satellite_maskedpreds.csv", row.names = FALSE)

metrics_summary_maskedpreds <- metrics_maskedpreds_satellite %>%
  dplyr::select( Period,Model, Threshold, Accuracy, TSS, Sensitivity, Specificity) %>%
  arrange(Period, Model, Threshold)%>%
  mutate(
    Accuracy = round(Accuracy, 2),
    TSS = round(TSS, 2),
    Sensitivity = round(Sensitivity, 2),
    Specificity = round(Specificity, 2)
  )

# write.csv(metrics_summary_maskedpreds,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_threshold_dependent_satellite_maskedpreds.csv")



# Dataset 2024: Capillary data ======
raster_folder_p2 <- out_path_t2

model_names <- c("brt", "ensemble", "gam", "glm", "rf")
buffer_capillary_maskedpred_df <- build_buffer_pred_df(
  model_names=model_names,
  xy_df = capillary_df_t2_filtered,     # your table of ID, kelp, x_3005, y_3005
  raster_folder = raster_folder_p2,
  buffer_radius = 180,
  threshold = thresholds$drop1[1]        # example threshold
)


model_names <- c("brt", "ensemble", "gam", "glm", "rf")
eval_capillary_maskedpreds <- process_buffer_predictions(
  buffer_pred_df = buffer_capillary_maskedpred_df,   # your averaged predictions per ID
  thresholds = thresholds,                    # named list e.g. list(opt = 0.3, maxSens = 0.25)
  model_names = c("brt", "ensemble", "gam", "glm", "rf")
)

metrics_capillary_maskedpreds<- eval_capillary_maskedpreds$Metrics
metrics_capillary_maskedpreds<- metrics_capillary_maskedpreds[,c(1,2,6:9)]
metrics_capillary_maskedpreds$period <- "Period2"
metrics_capillary_maskedpreds$dataset<- "Capillary_2024"

metrics_capillary_maskedpreds_points<- eval_capillary_maskedpreds$Classification
metrics_capillary_maskedpreds_points


## Dataset 2025:  Images   ====
# create buffer area 
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- out_path_t2

buffer_images_pred_masked_df <- build_buffer_pred_df(
  model_names= model_names,
  xy_df = test_images_t2_filtered,     # your table of ID, kelp, x_3005, y_3005
  raster_folder = raster_folder_p2,
  buffer_radius = 30,
  threshold = thresholds$drop1[1]        # example threshold
)
# calculate metrics of evaluation
eval_images_maskedpreds <- process_buffer_predictions(
  buffer_pred_df = buffer_images_pred_masked_df,# averaged predictions per ID
  thresholds = thresholds,                      
  model_names = c("brt", "ensemble", "gam", "glm", "rf")
)

metrics_images_maskedpreds_p2<- eval_images_maskedpreds$Metrics[,c(1,2,6:9)]
metrics_images_maskedpreds_p2$period <- "Period2"
metrics_images_maskedpreds_p2$dataset<- "Images_2025"
head(metrics_images_maskedpreds_p2)

metrics_images_maskedpreds_p2_points<- eval_images_maskedpreds$Classification


## Dataset 2025:  kayak surveys data   ====
raster_folder_p2 <- out_path_t2
colnames(test_kayak_t2)[16:17]<- c("x", "y")

eval_kayak_maskedpreds_p2 <- process_period_v2(raster_folder_p2, test_kayak_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_kayak_maskedpreds_p2<- eval_kayak_maskedpreds_p2$Metrics
metrics_kayak_p2$period <- "Period2"
metrics_kayak_p2$dataset<- "kayak_2025"


metrics_kayak_p2_points<- eval_kayak_p2$Classification
metrics_kayak_p2_points$Classification<- as.factor(metrics_kayak_p2_points$Classification)


# Kayak plus images evaluation =====
dataset_2025<- rbind(test_kayak_t2[, which(colnames(test_kayak_t2)%in% colnames(test_images_t2_filtered))], test_images_t2_filtered[, which(colnames(test_images_t2_filtered)%in% colnames(test_kayak_t2))])
summary(as.factor(dataset_2025$kelp))
# 0  1 (absence presence)
# 89 43

model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed/tifs_masked"
buffer_maskedpred_df <- build_buffer_pred_df(
  model_names= model_names,
  xy_df = dataset_2025,     # your table of ID, kelp, x_3005, y_3005
  raster_folder = raster_folder_p2,
  buffer_radius = 30,
  threshold = thresholds$drop1[1]        # example threshold
)

eval_dataset_2025_maskedpreds <- process_buffer_predictions(
  buffer_pred_df = buffer_maskedpred_df,   # your averaged predictions per ID
  thresholds = thresholds,                    # named list e.g. list(opt = 0.3, maxSens = 0.25)
  model_names = c("brt", "ensemble", "gam", "glm", "rf")
)

metrics_dataset_2025_maskedpreds_p2<- eval_dataset_2025_maskedpreds$Metrics[,c(1,2,6:9)]
metrics_dataset_2025_maskedpreds_p2$period <- "Period2"
metrics_dataset_2025_maskedpreds_p2$dataset<- "dataset_2025"

metrics_dataset_2025_maskedpreds_p2_points<- eval_dataset_2025_maskedpreds$Classification
metrics_dataset_2025_maskedpreds_p2_points$Classification<- as.factor(metrics_dataset_2025_maskedpreds_p2_points$Classification)


# Dataset Mayne Island Conservancy Society ====
model_names <- c("brt", "ensemble", "gam", "glm", "rf")
raster_folder_p2 <- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob/M7/tifs/renamed/tifs_masked"
eval_MICS_maskedpreds_p2 <- process_period_v2(raster_folder_p2, test_MICS_t2, prev_train, thresholds, model_names, prevalence= prev_PeriodsCombined)
metrics_MICS_maskedpreds_p2<- eval_MICS_maskedpreds_p2$Metrics
metrics_MICS_maskedpreds_p2$period <- "Period2"
metrics_MICS_maskedpreds_p2$dataset<- "MICSIndependent_dataset"


metrics_MICS_maskedpreds_p2_points<- eval_MICS_maskedpreds_p2$Classification
metrics_MICS_maskedpreds_p2_points$Classification<- as.factor(metrics_MICS_maskedpreds_p2_points$Classification)


# Merge all classification metrics
# colnames(metrics_all)[1]<- "period"
# metrics_evaluation_all_datasets_maskedpreds<- rbind(#metrics_all[which(colnames(metrics_all)%in% colnames(metrics_kayak_p2))],
#                                         metrics_kayak_maskedpreds_p2, 
#                                         metrics_images_maskedpreds_p2, 
#                                         metrics_capillary_maskedpreds, 
#                                         metrics_dataset_2025_maskedpreds_p2,
#                                         metrics_MICS_maskedpreds_p2)

metrics_evaluation_all_datasets_period2_maskedpreds<- rbind(
  metrics_capillary_maskedpreds, 
  metrics_dataset_2025_maskedpreds_p2,
  metrics_MICS_maskedpreds_p2)


# write.csv(metrics_evaluation_all_datasets_period2_maskedpreds,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_evaluation_all_datasets_period2_maskedpreds.csv")
# write.csv(metrics_evaluation_all_datasets,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_evaluation_all_datasets.csv")

metrics_dataset_2025_maskedpreds_p2_points$dataset<- "2025_dataset"
metrics_MICS_maskedpreds_p2_points$dataset<- "Independent_MICSperiod2"
metrics_capillary_maskedpreds_points$dataset<- "2024_dataset"

evaluation_points_all_datasets_period2_maskedpreds<- rbind(metrics_dataset_2025_maskedpreds_p2_points,
                                               metrics_MICS_maskedpreds_p2_points, 
                                               metrics_capillary_maskedpreds_points)

# write.csv(evaluation_points_all_datasets_period2_maskedpreds,"/Volumes/Romina_PSF/PSF/SDM/Models_validation/points_evaluation_matrix_Independent_datasets_period2_maskedpreds.csv")


