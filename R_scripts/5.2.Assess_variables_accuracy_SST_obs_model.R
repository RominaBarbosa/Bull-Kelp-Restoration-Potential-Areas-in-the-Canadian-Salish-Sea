##==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
### 5- Assessing the quality/accuracy of modeled variables      ################
### to be used in the Habitat suitability models                ################
### Author: Romina Barbosa                                      ################
### Date last edition: 31-July-2025                             ################
###==================================================================
# Analysis of correlation of satellite SST from Mora-Soto et al. 2024 and modeled data
# Plots of correlations 
# Plots of SST values in map, modeled vs satellite and in-situ


library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(purrr)

# Inport dataset of nutients
path_data<- "/Volumes/Romina_PSF/PSF"
shp_summer <- st_read(paste(path_data, "Mora-Soto_SalishSea/SST/Summer1984_2022_Complete/Summer1984_2022_Complete.shp", sep="/"))
# shp_spring <- st_read(paste(path_data, "Mora-Soto_SalishSea/SST/Spring1984_2022_Complete/Spring1984_2022_Complete.shp", sep="/"))

st_geometry_type(shp_summer)
shp_summer_3005 <- st_transform(shp_summer, crs = 3005)

centroids <- st_centroid(shp_summer_3005)
coords <- st_coordinates(centroids)
shp_summer_3005$longitude <- coords[, 1]
shp_summer_3005$latitude <- coords[, 2]


# Spring data 
shp_spring <- st_read(paste(path_data, "Mora-Soto_SalishSea/SST/Spring1984_2022_Complete/Spring1984_2022_Complete.shp", sep="/"))

st_geometry_type(shp_spring)
shp_spring_3005 <- st_transform(shp_spring, crs = 3005)

validity <- st_is_valid(shp_spring_3005)
table(validity)
sum(st_is( shp_spring, "MULTIPOLYGON"))

centroids <- st_centroid(shp_spring_3005)
coords <- st_coordinates(centroids)

shp_spring_3005$longitude <- coords[, 1]
shp_spring_3005$latitude <- coords[, 2]

shp_spring_3005df <- shp_spring_3005 %>% st_drop_geometry()

# plot(st_geometry(shp_spring), col = 'lightblue')
plot(st_geometry(centroids), col = 'red',  pch = 20, cex = 1.5)


# Merge by stacking the features
# LST_merged_Blob <- rbind(shp_spring, shp_summer)



### Extract info from CLusters to group by env. region =========================
cluster_points<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/clusters_MoraSoto2024_3005.csv")
coastal_segments_LST <- st_as_sf(cluster_points[,c(13:33)], coords = c("lon_3005", "lat_3005"), crs = 3005)  # or your spatial data frame

colnames(coastal_segments_LST)
plot(coastal_segments_LST[18])
plot(coastal_segments_LST[19])



### Upload metrics and compare with in situ measurements spatial pattern ========
my_path<- ("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics")
files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE)
files<- files[77:92]

# Load file 1 to start building the merged dataset
df_all=read.csv(files[1])
head(df_all)
df_all<- df_all[,c(2:4)]
filename <- strsplit(tools::file_path_sans_ext(files[1]), "/")[[1]][7]
filename <- strsplit(filename, "_")[[1]][c(1,4,3)]
filename <- paste0(filename, collapse ="_")
colnames(df_all)[3]<- filename

# Merge all files by adding a column of the metric to the df_all dataset
for (i in 2:length(files)) {
  
  df= read.csv(files[i])
  df= df[,c(2:4)]

    filename <- strsplit(tools::file_path_sans_ext(files[i]), "/")[[1]][7]
    filename <- strsplit(filename, "_")[[1]][c(1,4,3)]
    filename <- paste0(filename, collapse ="_")
    
  colnames(df)[3]<- filename
  df_all= merge(df_all, df, by=c("latitude" , "longitude"))
  
}


df_all



# Convert to spatial data and plot in map
df_sf <- st_as_sf(df_all, coords = c("longitude", "latitude"), crs = 4326)  # WGS84

# Transform df_sf to match CRS of stations_with_metrics
df_sf_3005 <- st_transform(df_sf, crs = st_crs(coastal_segments_LST))




### Merge data from modeled variable at locations of in situ data ==============
# Find index of nearest modeled point for each station
nearest_indices <- st_nearest_feature(coastal_segments_LST, df_sf_3005)
length(nearest_indices)# 4727 segments/stations
length((coastal_segments_LST$SprTSM_mea)) # 4727, values were extracted at all segments positions 

# Extract only the nearest modeled rows (i.e., exclude values not pairing LST locations)
modeled_nearest <- df_sf_3005[nearest_indices, ]

modeled_nearest<- modeled_nearest%>%
mutate(
  x = st_coordinates(.)[,1],
  y = st_coordinates(.)[,2]
) %>%
st_drop_geometry()


# combine columns to add coordinates and station name
stations_modeled_data_nearest<- bind_cols(coastal_segments_LST, modeled_nearest)
library(stringr)


################################################################################
### Standardize values to compare spatial patterns =============================
library(ggpubr)
library(patchwork)


LST_and_modeledSST<- stations_modeled_data_nearest%>%
  select(Long, Lat, Spring_mea, Summer_mea, temperature_mean_spring, temperature_mean_summer, Cluster,  x, y)

# Standardize value within nutrient-metric-season groups
data_standardized <- LST_and_modeledSST %>%
  mutate(across(3:6, ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE)), .names = "{.col}_std"))

# With Min-Max Scaling, we scale the data values between a range of 0 to 1 only. 
# Due to this, the effect of outliers on the data values suppresses to a certain extent. 
# Moreover, it helps us have a smaller value of the standard deviation of the data scale.
# https://www.digitalocean.com/community/tutorials/normalize-data-in-r


### SAVE std values of combined DATA ========
# write.csv(data_standardized, "/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/LSTclimatology_Vs_modeled_SST_standardized.csv")

# Correlation for summer
cor_summer <- cor(data_standardized$Summer_mea_std, data_standardized$temperature_mean_summer_std, use = "complete.obs")

# Correlation for spring
cor_spring <- cor(data_standardized$Spring_mea_std, data_standardized$temperature_mean_spring_std, use = "complete.obs")

cor_summer
# [1] 0.8963167
cor_spring
# [1] 0.8777167


### PLOT Correlation ===========================================================

# Build plot for SUMMER conditions
test <- cor.test(data_standardized$Summer_mea_std, data_standardized$temperature_mean_summer_std)
r <- round(test$estimate, 3)
p <- signif(test$p.value, 3)
annotation_text1 <- paste0("r = ", r, "\n", "p = ", p)


p_summer<- na.exclude(data_standardized)%>%
ggplot( aes(x = Summer_mea_std, y = temperature_mean_summer_std, color = as.factor(Cluster))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 1:1 red dashed line
  labs(
    y = "Standardized Summer SST (modeled)",
    x = "Standardized Summer SST (satellite)",
    color= "Env. Cluster"
  ) +
  annotate(
    "text",
    x = min(data_standardized$Summer_mea_std, na.rm = TRUE),   # Adjust x, y position as needed
    y = max(data_standardized$temperature_mean_summer_std, na.rm = TRUE),
    label = annotation_text1,
    hjust = 0,
    vjust = 1,
    size = 4,
    color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )



# Build plot fro SPRING standardized conditions
test <- cor.test(data_standardized$Spring_mea_std, data_standardized$temperature_mean_spring_std)
r <- round(test$estimate, 3)
p <- signif(test$p.value, 3)
annotation_text2 <- paste0("r = ", r, "\n", "p = ", p)


p_spring<- na.exclude(data_standardized)%>%
  ggplot( aes(x = Spring_mea_std, y = temperature_mean_spring_std, color = as.factor(Cluster))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 1:1 red dashed line
  labs(
    y = "Standardized Spring SST (modeled)",
    x = "Standardized Spring SST (satellite)",
    color= "Env. Cluster"
  ) +
  annotate(
    "text",
    x = min(data_standardized$Summer_mea_std, na.rm = TRUE),   # Adjust x, y position as needed
    y = max(data_standardized$temperature_mean_summer_std, na.rm = TRUE),
    label = annotation_text2,
    hjust = 0,
    vjust = 1,
    size = 4,
    color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )



# Extract the legend from the first plot
legend <- cowplot::get_legend(p_summer)

# Remove the legend from both plots for alignment
p2 <- p_summer + theme(legend.position = "none")
p1 <- p_spring

# Combine plots side by side with aligned widths
combined <- cowplot::plot_grid(p1, p2, nrow = 1, align = "hv")

# Add the legend next to them (optional)
cowplot::plot_grid(combined, legend, nrow = 1, rel_widths = c(1, 0.2))

# ggsave("/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/correlation_LSTclimatology_SSTmodel_Plots.pdf",
#        width = 26, height = 10, dpi = 300, units = "cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/correlation_LSTclimatology_SSTmodel_Plots.png",
#        width = 26, height = 10, dpi = 300, units = "cm")
# 
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/correlation_LSTclimatology_SSTmodel_Plots_black.pdf",
#        width = 26, height = 10, dpi = 300, units = "cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/correlation_LSTclimatology_SSTmodel_Plots_black.png",
#        width = 26, height = 10, dpi = 300, units = "cm")


# Build plot of SPRING conditions (absolute values)
test <- cor.test(data_standardized$Spring_mea, data_standardized$temperature_mean_spring)
r <- round(test$estimate, 3)
p <- signif(test$p.value, 3)
annotation_text3 <- paste0("r = ", r, "\n", "p = ", p)

p_spring_abs<- na.exclude(data_standardized)%>%
  ggplot( aes(x = Spring_mea, y = temperature_mean_spring, color = as.factor(Cluster))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 1:1 red dashed line
  labs(
    y = "Spring SST (modeled)",
    x = "Spring SST (satellite)",
    color= "Env. Cluster"
  ) +
  annotate(
    "text",
    x = min(data_standardized$Summer_mea, na.rm = TRUE),   # Adjust x, y position as needed
    y = max(data_standardized$temperature_mean_spring, na.rm = TRUE),
    label = annotation_text3,
    hjust = 0,
    vjust = 1,
    size = 4,
    color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )


test <- cor.test(data_standardized$Summer_mea, data_standardized$temperature_mean_summer)
r <- round(test$estimate, 3)
p <- signif(test$p.value, 3)
annotation_text4 <- paste0("r = ", r, "\n", "p = ", p)

p_summer_abs<- na.exclude(data_standardized)%>%
  ggplot( aes(x = Summer_mea, y = temperature_mean_summer, color = as.factor(Cluster))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 1:1 red dashed line
  labs(
    y = "Summer SST (modeled)",
    x = "Summer SST (satellite)",
    color= "Env. Cluster"
  ) +
  annotate(
    "text",
    x = min(data_standardized$Summer_mea, na.rm = TRUE),   # Adjust x, y position as needed
    y = max(data_standardized$temperature_mean_summer, na.rm = TRUE),
    label = annotation_text4,
    hjust = 0,
    vjust = 1,
    size = 4,
    color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

# Extract the legend from the first plot
legend <- cowplot::get_legend(p_summer_abs)
p2 <- p_summer_abs + theme(legend.position = "none")
p1 <- p_spring_abs

combined <- cowplot::plot_grid(p1, p2, nrow = 1, align = "hv")
cowplot::plot_grid(combined, legend, nrow = 1, rel_widths = c(1, 0.2))

# ggsave("/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/correlation_LSTclimatology_SSTmodel_Plots_absolute.pdf",
#        width = 26, height = 10, dpi = 300, units = "cm")
# ggsave("/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel/correlation_LSTclimatology_SSTmodel_Plots_absolute.png",
#        width = 26, height = 10, dpi = 300, units = "cm")



#### Plot standardized values of modeled and satellite SST in a map ====
library(sf)
library(purrr)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

### Function to plot map with nutrients values =================================

# Load coastline and land basemap (if not already loaded)
coastline <- ne_coastline(scale = "large", returnclass = "sf")
coastline <- st_transform(coastline, st_crs(combined_data))
land_highres <- ne_download(scale = "large", type = "land", category = "physical", returnclass = "sf")
land_highres <- st_transform(land_highres, st_crs(combined_data))

# Generate maps with basemap
ylim <- c(354740.5, 629504.7)  # Example East-West bounds
xlim <- c(970000.7, 1258935.2)    # Example North-South bounds

# bb_3005 <- st_bbox(c(xmin = 970000.7, xmax = 1258935.2, ymin = 354740.5, ymax = 629504.7),
#                    crs = st_crs(3005))


# Plotting

map_SST_plotfunction<- function(data= data_standardized, plot_standardized= T, 
                                save_format= "png", output_path= "/Volumes/Romina_PSF/PSF/PSF_nutrients_data",
                                output_name= "map_Satellite_&_Modeled_SST"){
  
  # Select columns based on the flag
  columns_to_pivot <- names(data) %>%
    keep(~ {
      is_std <- str_ends(.x, "_std")
      is_satellite <- str_starts(.x, "Spring_mea") | str_starts(.x, "Summer_mea")
      is_modeled <- str_starts(.x, "temperature")
      
      if (plot_standardized) {
        (is_std && (is_satellite | is_modeled))
      } else {
        (!is_std && (is_satellite | is_modeled))
      }
    })
  
  # Reshape the data
  data_long <- data_standardized %>%
    pivot_longer(
      cols = all_of(columns_to_pivot),
      names_to = "variable",
      values_to = "SST"
    ) %>%
    mutate(
      data_type = case_when(
        str_starts(variable, "temperature") ~ "Modeled",
        TRUE ~ "Satellite"
      ),
      season = case_when(
        str_detect(variable, "spring|Spring") ~ "Spring",
        str_detect(variable, "summer|Summer") ~ "Summer",
        TRUE ~ "Unknown"
      ),
      standardized = str_ends(variable, "_std")
    )
  
  if (plot_standardized) {
    color_scale <- scale_color_distiller(palette = "RdYlBu", direction = -1) #
  } else {
    color_scale <- scale_color_distiller(palette = "Spectral", direction = -1)
  }
  
  p<- ggplot(data_long) +
    geom_sf(data = land_highres, fill = "grey95", color = NA) +
    geom_sf(data = coastline, color = "grey40", size = 0.3) +
    geom_sf(aes(color = SST), size = 0.5) +
    color_scale +# scale_color_distiller(palette = "Spectral", direction = -1)
    facet_grid(data_type ~ season) +
    theme_bw() +
    scale_x_continuous(breaks = seq(-126, -122, by = 1)) +
    scale_y_continuous(breaks = seq(47, 51, by = 1))+
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +  # ← ZOOM HERE
    labs(
      title = "Sea Surface Temperature (SST)",
      subtitle = ifelse(plot_standardized, "Standardized Values", "Absolute Values"),
      color = ifelse(plot_standardized, "Standardized SST", "SST (°C)"))+
    theme(panel.border = element_rect(colour = "black", fill = NA),
              axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
              plot.title = element_text(size = 10, face = "bold"),
              # plot.text = element_text(size = 9, face = "bold"),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 7),
              panel.grid = element_blank()
        )
    
  
  if (plot_standardized) {
  output_name= paste(paste(output_name, "std", sep="_"), save_format, sep=".")
  }else{
    output_name= paste(output_name, save_format, sep=".") 
  }
  
  # Save  
  ggsave(paste(output_path, output_name, sep="/"), #map_Standardized_inSitu_Vs_Modeled_Nitrate.pdf",
         width = 17, height = 20, dpi = 300, units = "cm")
  
  return(p)
}



### True values MAPs ===========================================================
# Define common color scale
plot_std<- map_SST_plotfunction(data= data_standardized, plot_standardized= T, 
                                save_format= "png", output_path= "/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel",
                                output_name= "map_Satellite_&_Modeled_SST")

plot_std<- map_SST_plotfunction(data= data_standardized, plot_standardized= T, 
                                save_format= "pdf", output_path= "/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel",
                                output_name= "map_Satellite_&_Modeled_SST")

plot_abs<- map_SST_plotfunction(data= data_standardized, plot_standardized= F, 
                                save_format= "png", output_path= "/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel",
                                output_name= "map_Satellite_&_Modeled_SST")

plot_abs<- map_SST_plotfunction(data= data_standardized, plot_standardized= F, 
                                save_format= "pdf", output_path= "/Volumes/Romina_PSF/PSF/SDM/Correlation_LSTclimatology_SSTmodel",
                                output_name= "map_Satellite_&_Modeled_SST")



