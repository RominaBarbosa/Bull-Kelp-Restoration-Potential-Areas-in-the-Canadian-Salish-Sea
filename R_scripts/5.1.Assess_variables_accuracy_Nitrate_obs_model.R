##==================================================================
###     Species Distribution models    SDMs                     ################
###                                                             ################
### 5- Assessing the quality/accuracy of modeled variables      ################
### to be used in the Habitat suitability models                ################
### Author: Romina Barbosa                                      ################
### Date last edition: 31-July-2025                             ################
###==================================================================
# Analysis of correlation against in-situ data of Nitrate from PSF surveys (available online)
# Plots of correlations (standardized values from selected stations processed in
# "exploring_insitu_nutrients_PSF.R")
# Plots of nitrate values in map, modeled vs in-situ


library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(purrr)

# Inport dataset of nutients
library(readxl)
# Nutrients_2015_2023 <- read_excel("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/original_data_nutrients/2015-2023_Nutrients_240604.xlsx", 
#                                   col_types = c("text", "text", "numeric", 
#                                                 "numeric", "numeric", "numeric", 
#                                                 "numeric", "numeric", "numeric"))
# colnames(Nutrients_2015_2023)<- c("Station_name", "Date", "Time_PST", "lat", "lon", "Depth_m","Nitrate_um", "Phosphate_um", "Silicon_um")

Nutrients_2015_2019 <-read.csv("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/filtered_data_with_status3YEARSfiltering.csv")

Nutrients_2015_2019<- Nutrients_2015_2019 %>%
  filter(status == "Kept")

# Ensure date column is in Date format
Nutrients_2015_2019$Date <- ymd(Nutrients_2015_2019$Date)

# Keep only one coordinate per site to plot 
stations_df <- Nutrients_2015_2019 %>%
  distinct(Station_name, .keep_all = TRUE) %>%
  select(Station_name, Latitude_dd = lat, Longitude_dd = lon)

# Plot stations
stations_sf <- stations_df %>%
  st_as_sf(coords = c("Longitude_dd", "Latitude_dd"), crs = 4326)

mapview::mapview(stations_sf)

### Extract info from CLusters to group by env. region =========================
# Convert stations to sf
stations_sf <- stations_df %>%
  st_as_sf(coords = c("Longitude_dd", "Latitude_dd"), crs = 4326) %>%
  st_transform(crs = 3005)

cluster_points<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Variables_selection/clusters_MoraSoto2024_3005.csv")
cluster_sf <- st_as_sf(cluster_points[,c(31:33)], coords = c("lon_3005", "lat_3005"), crs = 3005)  # or your spatial data frame

# Get index of nearest cluster point for each station
nearest_index <- st_nearest_feature(stations_sf, cluster_sf)

# Bind cluster info to stations
stations_with_clusters <- stations_sf %>%
  bind_cols(cluster_id = cluster_sf$Cluster[nearest_index])


# Count most frequent cluster_id per Station_name
station_cluster_mode <- stations_with_clusters %>%
  st_drop_geometry() %>%
  count(Station_name, cluster_id, sort = TRUE) %>%
  group_by(Station_name) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%  # get only the top cluster_id
  ungroup()

length(station_cluster_mode$Station_name) #36
# station_cluster_mode<- station_cluster_mode[-duplicated(station_cluster_mode$Station_name)]

stations_with_coords <- stations_sf %>%
  mutate(
    Longitude = st_coordinates(.)[,1],
    Latitude = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

# merge with station coordinates
stations_with_cluster_mode <- stations_with_coords %>%
  left_join(station_cluster_mode, by = "Station_name")

stations_with_cluster_mode%>%
  group_by(Station_name, Longitude, Latitude)%>%
  summarize(n= length(cluster_id))

# Add cluster id to the entire dataset 
nutrients_with_clusters <- Nutrients_2015_2019 %>%
  left_join(stations_with_cluster_mode, by = "Station_name", relationship = "many-to-many")

summary(as.factor(stations_with_cluster_mode$cluster_id))
length(Nutrients_2015_2019$X)
length(nutrients_with_clusters$X)
length(unique(nutrients_with_clusters$Station_name))



# Pivot longer to make metrics
long_df <- nutrients_with_clusters #%>%
colnames(long_df)[9]<- "Concentration"  
colnames(long_df)[12]<-"season" 
long_df$Station_name<- as.factor(long_df$Station_name)


# Prepare dataframe with season and filter for surface depth
long_df <- long_df %>%
  mutate(
    Month = lubridate::month(Date)
  )


# Summarise statistics: mean, sd, min, max
cluster_season_stats <- long_df %>%
  group_by(cluster_id, Year, season, Nutrient) %>%
  summarise(
    mean = mean(Concentration, na.rm = TRUE),
    SD = sd(Concentration, na.rm = TRUE),
    minimum = if (all(is.na(Concentration))) NA_real_ else min(Concentration, na.rm = TRUE),
    maximum = if (all(is.na(Concentration))) NA_real_ else max(Concentration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Nutrient,
    values_from = c(mean, SD, minimum, maximum),
    names_glue = "{.value}_{Nutrient}"
  )


# Make dataframe long
cluster_stats_long <- cluster_season_stats %>%
  pivot_longer(
    cols = -c(cluster_id, season, Year),
    names_to = c("metric", "nutrient"),
    names_sep = "_",
    values_to = "Value"
  )


# Join metric values to stations by cluster_id and season:
stations_with_metrics <- stations_with_clusters %>%
  left_join(cluster_stats_long, by = "cluster_id", relationship = "many-to-many") #%>%
  # filter(season %in% c("summer", "spring", "fall", "winter")) # optional season filter



### Upload metrics and compare with in situ measurements spatial pattern ========
my_path<- ("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics")
files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE)
files<- files[33:48]

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
df_sf_3005 <- st_transform(df_sf, crs = st_crs(stations_with_metrics))


# Plot a single variable with mapview
# mapview::mapview(df_sf, zcol = "nitrate_mean_summer", legend = TRUE)# takes a long time because it's all points in the grid!



### Merge data from modeled variable at locations of in situ data ==============
# Find index of nearest modeled point for each station
nearest_indices <- st_nearest_feature(stations_with_clusters, df_sf_3005)
length(nearest_indices)# 36 stations
length(unique(stations_with_metrics$Station_name))

# 2. Extract only the nearest modeled rows
modeled_nearest <- df_sf_3005[nearest_indices, ]

# combine columns to add coordinates and station name
stations_modeled_data_nearest<- bind_cols(stations_with_clusters, modeled_nearest)
names(stations_modeled_data_nearest)[3]<- "geometry"
names(stations_modeled_data_nearest)[20]<- "geometry_nearest"

library(stringr)

# --- Step 1: Reshape modeled data (wide → long) ---
# Assume 'stations_modeled_data_nearest' is your wide-format modeled dataset
# Do NOT drop geometry — pivot inside the sf object
modeled_long <- stations_modeled_data_nearest %>%
  pivot_longer(
    cols = starts_with("nitrate"),  # all metric columns
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    nutrient = str_extract(variable, "^[^_]+"),
    metric = str_extract(variable, "(?<=_)[^_]+"),
    season = str_extract(variable, "[^_]+$"),
    data_type = "modeled"
  ) %>%
  select(Station_name, cluster_id, nutrient, metric, season, value, data_type, geometry, geometry_nearest)


# In situ data average for hte Blob period (2015-2019)
in_situ_data <- stations_with_metrics %>%
  filter(Year<= 2019)%>%
  group_by(Station_name, geometry, cluster_id, season, metric, nutrient)%>%
  summarize(value= mean(Value))%>%
  mutate(data_type = "in_situ") %>%
  # select(Station_name, cluster_id, nutrient, metric, season, value, geometry, data_type) %>%
  st_as_sf(crs = 3005)  # Ensure same CRS

# Make sure metric names are the same in both datasets
levels(as.factor(modeled_long$metric))
in_situ_data$metric<- as.factor(in_situ_data$metric)
levels(in_situ_data$metric)
modeled_long$nutrient<- as.factor(modeled_long$nutrient)
levels(modeled_long$nutrient)<- "Nitrate"
modeled_long$season<- as.factor(modeled_long$season)
levels(modeled_long$season)<- c( "Fall"  , "Spring" ,"Summer" ,"Winter")

# Combine
combined_data <- bind_rows(in_situ_data, modeled_long)
combined_data<- combined_data[,1:8]


unique(trimws(in_situ_data$Station_name)) %in% # Dataset 1
unique(trimws(modeled_long$Station_name))  # Dataset 2


################################################################################
### Standardize values to compare spatial patterns =============================
library(ggpubr)
library(patchwork)

# Standardize value within nutrient-metric-season groups
data_standardized <- combined_data %>%
  group_by(nutrient, metric, season, data_type) %>%
  mutate(value_std = (value-min(value,na.rm=TRUE))/(max(value,na.rm=TRUE)-min(value,na.rm=TRUE))) %>%  # scale returns a matrix
  ungroup()

# With Min-Max Scaling, we scale the data values between a range of 0 to 1 only. 
# Due to this, the effect of outliers on the data values suppresses to a certain extent. 
# Moreover, it helps us have a smaller value of the standard deviation of the data scale.
# https://www.digitalocean.com/community/tutorials/normalize-data-in-r


### SAVE std values of combined DATA ========
# write.csv(data_standardized, "/Volumes/Romina_PSF/PSF/PSF_nutrients_data/obs_modeled_nutrients_data_standardized.csv")




### PLOT Correlation ===========================================================
# Only modeled ammonium and in situ nutrients
data_filtered <- data_standardized %>%
  filter(
    (nutrient == "Nitrate" & data_type == "modeled") |
      (nutrient %in% c("Nitrate") & data_type == "in_situ")
  ) %>%
  # filter(metric== "mean", season == "spring") %>%
  mutate(nutrient_id = paste0(nutrient, "_", data_type)) %>%
  select(Station_name, geometry, cluster_id, metric, season, nutrient_id, value_std)

data_filtered_with_coords <- data_filtered %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()  # remove geometry column, keep coords


# Reshape to wide format to have data type in different columns
data_wide <- data_filtered_with_coords[,-2] %>%
  pivot_wider(names_from = nutrient_id, values_from = value_std)


data_wide<- na.exclude(data_wide)
length(!is.na(data_wide$Nitrate_modeled))# 424 pairing values were extracted 
summary(data_wide)

### Summary of stations number for correlations ======================================
data_wide%>%
  group_by(metric, season)%>%
  summarise(n_stations= length(Station_name))

# metric  season n_stations
# 1 SD      Fall           31
# 2 SD      Spring         36
# 3 SD      Summer         36
# 4 maximum Fall           35
# 5 maximum Spring         36
# 6 maximum Summer         36
# 7 mean    Fall           35
# 8 mean    Spring         36
# 9 mean    Summer         36
# 10 minimum Fall           35
# 11 minimum Spring         36
# 12 minimum Summer         36




# List of in situ nutrient columns
nutrients <- c("Nitrate_in_situ")#, "Phosphate_in_situ", "Silicon_in_situ")

plot_correlation <- function(df, x_var, y_var, season_name, metric_name) {
  # Pearson correlation test
  test <- cor.test(df[[x_var]], df[[y_var]], method = "pearson", use = "complete.obs")
  
  # Format correlation and p-value labels
  cor_label <- paste0("r = ", round(test$estimate, 2))
  p_label <- paste0("p = ", signif(test$p.value, 2))
  annotation_text <- if (test$p.value < 0.05) paste(cor_label, p_label, sep = "\n") else cor_label
  
  # Build plot
  ggplot(df, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(
      title = paste(season_name, "| Metric:", metric_name),
      x = "Standardized Nitrate (modeled)",
      y = paste("Standardized", y_var)
    ) +
    annotate(
      "text",
      x = min(df[[x_var]], na.rm = TRUE),
      y = max(df[[y_var]], na.rm = TRUE),
      label = annotation_text,
      hjust = 0,
      vjust = 1.5,
      size = 2.5,
      color = if (test$p.value < 0.05) "red" else "gray40"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
}




# This assumes `data_wide` and `nutrients` are already in memory
plots <- data_wide %>%
  filter(!is.na(Nitrate_modeled)) %>%
  filter(season != "Winter") %>%
  # filter(metric != "maximum" & metric != "minimum") %>%
  group_by(season, metric) %>%
  group_split() %>%
  lapply(function(group_data) {
    season_name <- unique(group_data$season)
    metric_name <- unique(group_data$metric)
    
    lapply(nutrients, function(nutrient_col) {
      plot_correlation(group_data, "Nitrate_modeled", nutrient_col, season_name, metric_name)
    })
  }) %>%
  purrr::flatten()  # flattens to a simple list of ggplot objects




library(patchwork)
wrap_plots(plots[c(1:12)], ncol = 4)
# ggsave("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/correlation_Nitrate_inSitu_Vs_Modeled_Plots.pdf",
#        width = 22, height = 17, dpi = 300, units = "cm")
# ggsave("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/correlation_Nitrate_inSitu_Vs_Modeled_Plots.png",
#        width = 22, height = 17, dpi = 300, units = "cm")



# Initialize storage
correlation_results <- list()

# Plot and collect correlation stats
data_wide.2<- na.exclude(data_wide)

# plots <- 
  data_wide.2 %>%
  filter(!is.na(Nitrate_modeled)) %>%
  group_by(season, metric) %>%
  group_split() %>%
  lapply(function(group_data) {
    season_name <- unique(group_data$season)
    metric_name <- unique(group_data$metric)
    n_stations  <- group_data %>%
      summarize(n_stations= length(Station_name))
    n_stations<- n_stations[,"n_stations"]
    
    lapply(nutrients, function(nutrient_col) {
      df <- group_data %>%
        filter(!is.na(.data[[nutrient_col]]), !is.na(Nitrate_modeled))
      
      # Compute correlation
      test <- cor.test(df[[nutrient_col]], df$Nitrate_modeled, method = "pearson")
      
      # Store in summary list
      correlation_results[[length(correlation_results) + 1]] <<- tibble(
        season = season_name,
        metric = metric_name,
        nutrient = nutrient_col,
        correlation = test$estimate,
        p_value = test$p.value,
        significance = case_when(
          test$p.value < 0.001 ~ "***",
          test$p.value < 0.01 ~ "**",
          test$p.value < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        n_stations= n_stations$n_stations
      )
      
  #     # Build plot
  #     ggplot(df, aes_string(x = "Nitrate_modeled", y = nutrient_col)) +
  #       geom_point(alpha = 0.6) +
  #       geom_smooth(method = "lm", se = TRUE, color = "blue") +
  #       labs(
  #         title = paste("Season:", season_name, "| Metric:", metric_name),
  #         subtitle = paste("r =", round(test$estimate, 2), ", p =", signif(test$p.value, 2)),
  #         x = "Standardized Nitrate (modeled)",
  #         y = paste("Standardized", nutrient_col)
  #       ) +
  #       theme_minimal() +
  #       theme(panel.border = element_rect(colour = "black", fill = NA),
  #                             axis.text = element_text(size=8),
  #                             plot.title = element_text(size = 10, face = "bold"),
  #                             # axis.text = element_text(size = 9, face = "bold"),
  #             axis.title.y = element_text(size = 9),  # Smaller x-axis text
  #             axis.title.x = element_text(size = 9),   # Smaller y-axis text
  #                             panel.grid = element_blank()
  #       )
    })
  }) #%>%
  # purrr::flatten()


correlation_summary <- bind_rows(correlation_results)
# write.csv(correlation_summary, "/Volumes/Romina_PSF/PSF/PSF_nutrients_data/correlation_inSitu_NItrate_inSitu_Vs_Modeled_Table.csv", row.names = FALSE)

wrap_plots(plots[c(1:12)], ncol = 4)
ggsave()



#### Plot standardized values of selected variables and in situ nitrate in a map ====
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

data_standardized <- data_standardized %>%
  mutate(plot_label = paste(nutrient, metric, season, data_type, sep = "_"))

map_list_Std <- data_standardized %>%
  filter(season != "Winter")%>% # filter(season == "Spring"| season=="Summer")%>%
  filter(metric == "mean")%>%
  group_split(plot_label) %>%
  map(~{
    label <- unique(.x$plot_label)
    ggplot() +
      geom_sf(data = land_highres, fill = "grey95", color = NA) +
      geom_sf(data = coastline, color = "grey40", size = 0.3) +
      geom_sf(data = .x, aes(color = value_std), size = 2, alpha = 0.8) +
      scale_color_viridis_c(option = "plasma", name = "Standardized\nValue") +
      labs(title = label) +
      scale_x_continuous(breaks = seq(-126, -122, by = 1)) +
      scale_y_continuous(breaks = seq(47, 51, by = 1))+
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +  # ← ZOOM HERE
      theme_minimal() +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(size = 10, face = "bold"),
        # plot.text = element_text(size = 9, face = "bold"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank()
      )
  })


library(patchwork)
wrap_plots(map_list_Std, ncol = 2)

# ggsave("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/map_Standardized_inSitu_Vs_Modeled_Nitrate.pdf",
#        width = 17, height = 20, dpi = 300, units = "cm")
# ggsave("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/map_Standardized_inSitu_Vs_Modeled_Nitrate.png",
#        width = 17, height = 20, dpi = 300, units = "cm")



### True values MAPs ===========================================================
# Define common color scale
summary(data_standardized$value)

common_color_scale <- scale_color_viridis_c(
  option = "plasma",
  name = "Nitrate (uM)",
  limits = c(-1, 27)  # <- This makes the color scale fixed for all plots
)

map_list_spring_<- data_standardized %>%
  filter(season != "Winter")%>% # filter(season == "Spring"| season=="Summer")%>%
  filter(metric == "mean")%>%
  # filter(data_type == "in_situ")%>%
  group_split(plot_label) %>%
  map(~{
    label <- unique(.x$plot_label)
    ggplot() +
      geom_sf(data = land_highres, fill = "grey95", color = NA) +
      geom_sf(data = coastline, color = "grey40", size = 0.3) +
      geom_sf(data = .x, aes(color = value), size = 2, alpha = 0.8) +
      # scale_color_viridis_c(option = "plasma", name = "Nitrate (uM.L-1)") +
      common_color_scale+
      labs(title = label) +
      scale_x_continuous(breaks = seq(-126, -122, by = 1)) +
      scale_y_continuous(breaks = seq(47, 51, by = 1))+
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +  # ← ZOOM HERE
      theme_minimal() +
      theme(panel.border = element_rect(colour = "black", fill = NA),
            axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
            plot.title = element_text(size = 10, face = "bold"),
            # plot.text = element_text(size = 9, face = "bold"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            panel.grid = element_blank()
      )
  })


wrap_plots(map_list_spring_, ncol = 2)
# ggsave("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/map_inSitu_vs_modeled_Nitrate.pdf",
#        width = 17, height = 20, dpi = 300, units = "cm")
# ggsave("/Volumes/Romina_PSF/PSF/PSF_nutrients_data/map_inSitu_vs_modeled_Nitrate.png",
#        width = 17, height = 20, dpi = 300, units = "cm")




