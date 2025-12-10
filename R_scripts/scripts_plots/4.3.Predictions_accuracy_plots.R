###====================================================================================
###     Species Distribution models    SDMs                                      ######
###                                                                              ######
### 4.3- Evaluation of Predictions Accuracy (models performance)                 ######
### Content:                                                                     ######
###   Maps with presence and absence records from each dataset                   ######
###   Maps with presence and absence classified as TP, TN, FP and FN X period    ######
### Author: Romina Barbosa                                                       ######
### Date last edition: 24-Nov-2025                                               ######
###====================================================================================


figures_output_path<- "/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures"


# -------------------------------------------------------------------------------
#  Plot MAP with testing points (satellite data)
# -------------------------------------------------------------------------------
library(rnaturalearth)
library(rnaturalearthhires)  # high-resolution coastlines
library(rnaturalearthdata)

data_val_satellite<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/MoraSoto_postblob/Presence_absences_testing_dataset_satellite_FINAL.csv")


# Testing period 1
test_t1<- data_val_satellite %>%
  filter(Period == "2015_2019")

test_t1_sf <- st_as_sf(
  test_t1,
  coords = c("x", "y"),
  crs = 3005
)


# Testing period 2
test_t2<- data_val_satellite %>%
  filter(Period == "2020_2022")


test_t2_sf <- st_as_sf(
  test_t2_balanced,
  coords = c("x", "y"),
  crs = 3005
)


# Coast shapefile 
coast <- ne_coastline(scale = "large", returnclass = "sf")
coast_3005 <- st_transform(coast, 3005)
coast_vect <- vect(coast_3005)
coast_crop<- crop(coast_vect, ens_average)
coast_crop_sf <- st_as_sf(coast_crop)

# Get high-res countries polygons
land <- ne_countries(scale = "large", returnclass = "sf")
land_3005 <- st_transform(land, 3005)

# Crop to your raster area
# Get raster extent and convert to sf bbox
ext_bbox <- st_bbox(ext(ens_average))

# Crop land polygons
land_crop <- st_crop(land_3005, ext_bbox)

# Filter data by kelp
test_t1_kelp1 <- test_t1_sf %>% filter(kelp == 1)
test_t1_kelp0 <- test_t1_sf %>% filter(kelp == 0)
test_t2_kelp1 <- test_t2_sf %>% filter(kelp == 1)
test_t2_kelp0 <- test_t2_sf %>% filter(kelp == 0)

# Function to create a map with specific colors for kelp
create_map <- function(data, title, pints_size= .4, legend=F) {
  # p<- 
    ggplot() +
    geom_sf(data = land_crop, fill = "grey85", color = "black", linewidth = 0.3) +
    geom_sf(data = data, aes(color = factor(kelp)), shape = 16, size = pints_size) +
    scale_color_manual(values = c("0" = "blue", "1" = "green")) +
    coord_sf(crs = st_crs(3005)) +
    theme(panel.background = element_rect(fill = "aliceblue")) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "aliceblue"),
      axis.text.x = element_text(size = 7),   # smaller x-axis labels
      axis.text.y = element_text(size = 7)    # smaller y-axis labels
    ) +
    labs(title = title, color = "Kelp")
  
  # if( legend==T){
  #   p + theme(
  #     legend.position = c(0.8, 0.8))
  # }
  # return(p)
}

# Create the 4 maps
map_t1_kelp1 <- create_map(test_t1_kelp1, "Period 1: Kelp = Presence", pints_size= .4)
map_t1_kelp0 <- create_map(test_t1_kelp0, "Period 1: Kelp = Absence", pints_size= .4)
map_t2_kelp1 <- create_map(test_t2_kelp1, "Period 2: Kelp = Presence", pints_size= .4)
map_t2_kelp0 <- create_map(test_t2_kelp0, "Period 2: Kelp = Absence", pints_size= .4)

# Combine into a 2x2 grid
cowplot::plot_grid(map_t1_kelp1, map_t1_kelp0,
                   map_t2_kelp1, map_t2_kelp0,
                   ncol = 2)

# ggsave(paste(model_results_path, "Figure_testing_records_maps.pdf", sep="/"), width = 15, height = 17, units="cm", dpi = 300)



# -- Kayak dataset map -----------------------------------------------------------------------------
test_kayak_t2<- read.csv("/Volumes/Romina_PSF/PSF/SDM/fieldwork_2025/Kelp_Surveys_3_centroidxy3005.csv")
colnames(test_kayak_t2)[16:17]<- c("x", "y")
test_kayak_t2$kelp<- test_kayak_t2$Bull_Kelp +test_kayak_t2$Giant_Kelp # presence of any kelp (bull or giant kelp)

test_kayak_t2_sf <- st_as_sf(
  test_kayak_t2,
  coords = c("x", "y"),
  crs = 3005
)

# Plot
create_map(test_kayak_t2_sf, "Period 2: Kayak observations", pints_size= .8, legend=T) + theme(
  legend.position = c(0.8, 0.8))


## Images datset map 
test_images_t2_sf <- st_as_sf(
  test_images_t2,
  coords = c("x", "y"),
  crs = 3005
)

# Plot
create_map(test_images_t2_sf, "Period 2: Kayak observations", pints_size= .8, legend=T) + theme(
  legend.position = c(0.8, 0.8))

# -------------------------------------------------------------------------------
### Classification results evaluation ==========================================
# -------------------------------------------------------------------------------
metrics_p2_evaluation<- read.csv("")


create_map_confusion <- function(data, title, pints_size= .4, legend=F) {
  p<- data%>%
    # filter(Classification=="TP")%>%
    filter(Model=="ensemble")%>%
    filter(ThresholdName=="drop1")%>%
    ggplot() +
    geom_sf(data = land_crop, fill = "grey85", color = "black", linewidth = 0.3) +
    geom_sf( aes(color = factor(Classification)), shape = 16, size = pints_size) +
    facet_wrap(~Classification)+
    scale_color_manual(values = c("TN" = "blue", "FN" = "cyan", "TP"= "darkgreen", "FP"= "orange")) +
    coord_sf(crs = st_crs(3005)) +
    # theme(panel.background = element_rect(fill = "aliceblue")) +
    theme_bw() +
    theme(
      # legend.position = c(0.8, 0.8),
      panel.background = element_rect(fill = "aliceblue"),
      axis.text.x = element_text(size = 7),   # smaller x-axis labels
      axis.text.y = element_text(size = 7)    # smaller y-axis labels
    ) +
    labs(title = "", color = "Kelp")
  
  return(p)
}



# Plot satellite testing data classified ==========
metrics_points_satellite<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Models_validation/metrics_points_satellite.csv")

# Period 1
metrics_p1_points<- metrics_points_satellite%>%
  filter(Period=="Period1")

metrics_p1_points_sf <- st_as_sf(
  metrics_p1_points,
  coords = c("x", "y"),
  crs = 3005
)

create_map_confusion(metrics_p1_points_sf, "Period 1: Testing Dataset", 
                     pints_size= .5, legend=T) + theme(
                       legend.position = "right")
# ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_SatellitePeriod2_maps.jpg", sep="/"), width = 18, height = 17, units="cm", dpi = 300)


# Period 2
metrics_p2_points<- metrics_points_satellite%>%
  filter(Period=="Period2")

metrics_p2_points_sf <- st_as_sf(
  metrics_p2_points,
  coords = c("x", "y"),
  crs = 3005
)

create_map_confusion(metrics_p2_points_sf, "Period 2: Testing Dataset", 
                     pints_size= .5, legend=T) + theme(
                       legend.position = "right")
# ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_SatellitePeriod1_maps.pdf", sep="/"), width = 18, height = 17, units="cm", dpi = 300)


## Plot Independent period-2 dataset (Mayne Island Conservancy data) ========
metrics_MICS_p2_points<- metrics_p2_evaluation%>%
  filter(dataset=="MICS")

classif_MICSdatset_t2_sf <- st_as_sf(
  metrics_MICS_p2_points,
  coords = c("x", "y"),
  crs = 3005
)

classif_MICSdatset_t2_sf$Classification<- as.factor(classif_MICSdatset_t2_sf$Classification)
classif_MICSdatset_t2_sf

# Plot
# Get bounding box of your points
bb <- st_bbox(sf_MICS_dataset_unique)
expand_dist <- 20000   # e.g. 2 km buffer
  
create_map_confusion(classif_MICSdatset_t2_sf, "Period 2: Independent Period-2 Dataset", 
                     pints_size= .8, legend=T, facet=F) + theme(
  legend.position = "right")

# create_map_confusion(classif_MICSdatset_t2_sf, "Period 2: Independent Period-2 Dataset", pints_size= .8, legend=T) + theme(
#   legend.position = "right") + coord_sf(
#     xlim = c(bb$xmin - expand_dist, bb$xmax + expand_dist),
#     ylim = c(bb$ymin - expand_dist, bb$ymax + expand_dist),
#     expand = FALSE
#   )
# ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_IndepDataset_maps.pdf", sep="/"), width = 9, height = 8, units="cm", dpi = 300)

## Plot 2024 Dataset (Capillary dataset)  ========
metrics_capillary_p2_points<- metrics_p2_evaluation%>%
  filter(dataset=="Capillary_2024")

classif_capillary_t2_sf <- st_as_sf(
  metrics_capillary_p2_points,
  coords = c("x", "y"),
  crs = 3005
)

# Plot
create_map_confusion(classif_capillary_t2_sf, "Period 2: Capillary observations", pints_size= .8, legend=T) + theme(
  legend.position = c(0.9, 0.8))
# ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_2024Dataset_maps.jpg", sep="/"), width = 15, height = 8, units="cm", dpi = 300)



## Plot dataset 2025 ========
metrics_dataset_2025_p2_points<- metrics_p2_evaluation%>%
  filter(dataset=="Kayak_2025")

classif_datset2025_t2_sf <- st_as_sf(
  metrics_dataset_2025_p2_points,
  coords = c("x", "y"),
  crs = 3005
)

classif_datset2025_t2_sf$Classification<- as.factor(classif_datset2025_t2_sf$Classification)
classif_datset2025_t2_sf

# Plot
create_map_confusion(classif_datset2025_t2_sf, "Period 2: Dataset 2025", 
                           pints_size= .8, legend=T)  + theme(legend.position = "right")

# ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_datset2025_maps.pdf", sep="/"), width = 18, height = 17, units="cm", dpi = 300)


## Plot aerial images dataset  ========
metrics_images_p2_points<- metrics_p2_evaluation%>%
  filter(dataset=="Images_2025")

classif_images_t2_sf <- st_as_sf(
  metrics_images_p2_points,
  coords = c("x", "y"),
  crs = 3005
)


create_map_confusion(classif_images_t2_sf, "Period 2: Aerial images observations", pints_size= .8, legend=T) + theme(
  legend.position = "right")
ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_images_maps.pdf", sep="/"), width = 18, height = 17, units="cm", dpi = 300)


## Plot kayak dataset  ========
metrics_kayak_p2_points<- metrics_p2_evaluation%>%
  filter(dataset=="Kayak_2025")

classif_kayak_t2_sf <- st_as_sf(
  metrics_kayak_p2_points,
  coords = c("x", "y"),
  crs = 3005
)

classif_kayak_t2_sf$Classification<- as.factor(classif_kayak_t2_sf$Classification)
classif_kayak_t2_sf

# Plot
create_map_confusion(classif_kayak_t2_sf, "Period 2: Kayak observations", pints_size= .8, legend=T) + theme(
  legend.position = "right")

# ggsave(paste(figures_output_path, "Figure_classification_error_ensemble_kayak_maps.pdf", sep="/"), width = 18, height = 17, units="cm", dpi = 300)





# -------------------------------------------------------------------------------
#  Plot results of model performance based on testing dataset (satellite data)
# -------------------------------------------------------------------------------
metrics_all<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_satellite_evaluation_results_without_masking.csv")
metrics_all_masked<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_satellite_evaluation_results_maskedpreds.csv")

metrics_all<- metrics_all_masked

# compare AUC and Boyce_point across models and periods
metrics_all$Model <- factor(metrics_all$Model,
                            levels = c("ensemble", "brt", "gam", "glm", "rf"))


plot_Boyce<- ggplot(metrics_all, aes(x = Model, y = Boyce_point, fill = Period)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  # geom_text(aes(label = round(Boyce_point, 2), y = Boyce_point - 0.1),
  #           position = position_dodge(width = 0.8), 
  #           vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("Period2" = "#1f77b4", "Period1" = "#ff7f0e")) +
  theme_bw(base_size = 11) +
  theme(
    # legend.position = "none",
    # panel.background = element_rect(fill = "aliceblue"),
    axis.text.x = element_text(size = 8),   # smaller x-axis labels
    axis.text.y = element_text(size = 8)    # smaller y-axis labels
  ) +
  ylim(c(0,1))+
  labs(#title = "Model performance (Boyce Index) by period and threshold",
    y = "Boyce Index", x = "Model")


plot_AUC<- ggplot(metrics_all, aes(x = Model, y = AUC, fill = Period)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  ylim(c(0,1))+
  # geom_text(aes(label = round(AUC, 2), y = AUC - 0.1),
  #           position = position_dodge(width = 0.8), 
  #           vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("Period2" = "#1f77b4", "Period1" = "#ff7f0e")) +
  theme_bw(base_size = 11) +
  theme(
    legend.position= ,
    axis.text.x = element_text(size = 8),  
    axis.text.y = element_text(size = 8)    
  ) +
  labs(#title = "Model performance (AUC) by period and threshold",
    y = "AUC", x = "Model") 


cowplot::plot_grid(plot_AUC, plot_Boyce, ncol=2, align = "hv", labels = c("A)", "B)"), label_size = 11)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/AUC_Boyce_testing_barplot.jpg", width = 19, height = 6, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/AUC_Boyce_testing_maskedpreds_barplot.jpg", width = 19, height = 6, dpi= 300, units="cm")



## Threshold-dependent metrics 
metrics_all$TSS<- as.numeric(as.character(metrics_all$TSS))
metrics_all$Threshold<- as.factor(metrics_all$Threshold)
metrics_long <- metrics_all %>%
  pivot_longer(
    cols = c(TSS, Sensitivity, Specificity),
    names_to = "Metric",
    values_to = "Value"
  )

metrics_long$Model <- factor(metrics_long$Model,
                             levels = c("ensemble", "brt", "gam", "glm", "rf"))

plot_metrics <- ggplot(metrics_long, 
                       aes(x = Model, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Period + Threshold) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0, 1) +
  theme_bw(base_size = 13) +
  labs(
    # title = "Model performance by Period and Metric",
    y = "Value",
    x = "Model",
    fill = "M etric"
  )

plot_metrics
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/Sens_Spe_TSS_testing_barplot.pdf", width = 17, height = 13, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/Sens_Spe_TSS_testing_maskedpreds_barplot.jpg", width = 17, height = 13, dpi= 300, units="cm")



# Select threshold-dependent metrics
df_bar <- metrics_all %>%
  dplyr::select(Period, Model, Threshold, Accuracy) 

# Make Model a factor to control order
# df_bar$Model <-  factor(df_bar$Model,
#                         levels = c("ensemble", "brt", "gam", "glm", "rf"))


ggplot(df_bar, aes(x = Model, y = Accuracy, fill = Period)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Period1" = "#1f77b4", "Period2" = "#ff7f0e")) +
  facet_wrap(~ Threshold)+
  ylim(0, 1) +
  labs(title = "", y = "Accuracy", x = "Model") +
  theme_bw() +
  theme(
    legend.position= ,
    axis.text.x = element_text(size = 8),  
    axis.text.y = element_text(size = 8)    
  )

# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/Accuracy_testing_barplot.jpg", width = 17, height = 6, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/Accuracy_testing_maskedpreds_barplot.pdf", width = 17, height = 6, dpi= 300, units="cm")




##------------------------------------------------------------------------------
# Plot comparison of metrics with and without masking by substrate and depth
metrics_all<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_satellite_evaluation_results_without_masking.csv")
metrics_all$predictions<- "original"

metrics_all_masked<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Models_validation/validation_satellite/metrics_satellite_evaluation_results_maskedpreds.csv")
metrics_all_masked$predictions<- "masked"

metrics_all<- rbind(metrics_all_masked[,which(colnames(metrics_all_masked)%in% colnames(metrics_all))], 
                    metrics_all[,which(colnames(metrics_all)%in% colnames(metrics_all_masked))])

# compare AUC and Boyce_point across models and periods
metrics_all$Model <- factor(metrics_all$Model,
                            levels = c("ensemble", "brt", "gam", "glm", "rf"))

metrics_all$predictions <- factor(metrics_all$predictions,
                            levels = c("original", "masked"))


plot_Boyce_compare<- ggplot(metrics_all, aes(x = Model, y = Boyce_point, fill = predictions)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Boyce_point, 2), y = Boyce_point - 0.1),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c( "darkgrey", "lightgrey")) +
  theme_bw(base_size = 11) +
  facet_wrap(~Period)+
  theme(
    # legend.position = "none",
    # panel.background = element_rect(fill = "aliceblue"),
    axis.text.x = element_text(size = 8),   # smaller x-axis labels
    axis.text.y = element_text(size = 8)    # smaller y-axis labels
  ) +
  ylim(c(0,1))+
  labs(#title = "Model performance (Boyce Index) by period and threshold",
    y = "Boyce Index", x = "Model")


plot_AUC_compare<- ggplot(metrics_all, aes(x = Model, y = AUC, fill = predictions)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  ylim(c(0,1))+
  geom_text(aes(label = round(AUC, 2), y = AUC - 0.1),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c( "darkgrey", "lightgrey")) +
  theme_bw(base_size = 11) +
  facet_wrap(~Period)+
  theme(
    legend.position= ,
    axis.text.x = element_text(size = 8),  
    axis.text.y = element_text(size = 8)    
  ) +
  labs(#title = "Model performance (AUC) by period and threshold",
    y = "AUC", x = "Model") 


cowplot::plot_grid(plot_AUC_compare, plot_Boyce_compare, ncol=1, align = "hv", labels = c("A)", "B)"), label_size = 11)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/AUC_Boyce_testing_comparing_masked_barplot.jpg", width = 19, height = 12, dpi= 300, units="cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/AUC_Boyce_testing_comparing_masked_barplot.jpg", width = 19, height = 6, dpi= 300, units="cm")


## Threshold-dependent metrics 
metrics_all$TSS<- as.numeric(as.character(metrics_all$TSS))
metrics_all$Threshold<- as.factor(metrics_all$Threshold)
metrics_long <- metrics_all %>%
  pivot_longer(
    cols = c(TSS, Sensitivity, Specificity),
    names_to = "Metric",
    values_to = "Value"
  )

metrics_long$Model <- factor(metrics_long$Model,
                             levels = c("ensemble", "brt", "gam", "glm", "rf"))

plot_metrics_compare <- ggplot(metrics_long, 
                       aes(x = Model, y = Value, fill = Metric, color= predictions)) +
  geom_col(position = position_dodge(width = 0.8)) +
  # geom_text(aes(label = round(Value, 2), y = Value - 0.1),
  #           position = position_dodge(width = 0.8),
  #           vjust = -0.3, size = 2.5) +
  facet_wrap(~ Period + Threshold) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_manual(values= c("black", "grey")) +
  ylim(0, 1) +
  theme_bw(base_size = 13) +
  labs(
    # title = "Model performance by Period and Metric",
    y = "Value",
    x = "Model",
    fill = "M etric"
  )

plot_metrics_compare
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/Sens_Spe_TSS_testing_comparing_predictions_barplot.jpg", width = 17, height = 13, dpi= 300, units="cm")


# Select threshold-dependent metrics
df_bar <- metrics_all %>%
  dplyr::select(Period, Model, Threshold, Accuracy, predictions) 

ggplot(df_bar, aes(x = Model, y = Accuracy, fill = predictions)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Accuracy, 2), y = Accuracy - 0.1),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c( "darkgrey", "lightgrey")) +
  facet_wrap(~ Period + Threshold)+
  ylim(0, 1) +
  labs(title = "", y = "Accuracy", x = "Model") +
  theme_bw() +
  theme(
    legend.position= ,
    axis.text.x = element_text(size = 8),  
    axis.text.y = element_text(size = 8)    
  )

ggsave("/Volumes/Romina_PSF/PSF/SDM/Models_validation/Figures/Accuracy_testing_comparing_predictions_barplot.jpg", width = 17, height = 12, dpi= 300, units="cm")

