###====================================================================================
###     Species Distribution models    SDMs                                      ######
###                                                                              ######
### 7.1- Assessing the Classified Areas for Restoration                          ######
### Content:                                                                     ######
###   Plot maps of Restoration Areas Classification                              ######
###   Barplot of area per Restoration category along the latitudinal gradient    ######
### Author: Romina Barbosa                                                       ######
### Date last edition: 24-Nov-2025                                               ######
###====================================================================================

# Load packages
library(rnaturalearth)
library(rnaturalearthhires)  # high-resolution coastlines
library(rnaturalearthdata)

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
land_crop_4326 <- st_transform(land_crop, crs = 4326)



create_map <- function(land_crop= land_crop_4326, data, title, point_size = 0.4, show_legend = FALSE) {
    ggplot() +
      # Base map
      geom_sf(data = land_crop, fill = "grey85", color = "black", linewidth = 0.3) +
      # Raster / tile layer
      geom_tile(data = data, aes(x = x, y = y, fill = factor(category))) +
      scale_fill_manual(values = c("blue", "yellow",  "orange", "red")) +
      # Coordinate system
      coord_sf(crs = st_crs(3005)) +
      # Theme
      theme_bw() +
      theme(
        legend.position = if (show_legend) "right" else "none",
        panel.background = element_rect(fill = "aliceblue"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7)
      ) +
      # Labels
      labs(title = title, fill = "Class")
  }
  

    
r_df <- as.data.frame(r, xy = TRUE)  # adds x and y columns
colnames(r_df)[3] <- "category"
r_df$category <- as.factor(r_df$category)
summary(r_df$y)

# Create maps
p_mapa <- create_map(land_crop= land_crop, data=r_df, "Period 1: Kelp = Presence")



### Barplot 
# Steps in the function:
# Takes a raster (SpatRaster) in EPSG:3005
# Extracts cells → dataframe
# Removes NA
# Computes area per latitude bin × class
# Converts y-bin into latitude in degrees correctly using sf::st_transform()
# Handles any number of factor levels
# Lets you supply a color palette
# Returns a ggplot barplot

make_latitude_area_plot <- function(
    r,                   # SpatRaster in CRS 3005
    y_bin_size = 2000,   # bin height in meters
    cell_size = 20,      # raster cell size (m)
    palette = NULL,      # vector of colors (optional)
    flip = TRUE ,         # flip axes for final plot
    facet = FALSE,       # TRUE = facet by class
    axis = "lat"        # "lat" or "y"
) {
  
  # -------------------------
  # 1️⃣ Extract raster values
  # -------------------------
  pts <- as.data.frame(cbind(
    terra::xyFromCell(r, 1:ncell(r)),
    value = values(r)
  ))
  colnames(pts) <- c("x", "y", "class")
  pts <- pts %>% filter(!is.na(class))
  
  pts$class <- as.factor(pts$class)
  
  # -------------------------
  # 2️⃣ Bin in y direction (meters)
  # -------------------------
  pts <- pts %>%
    mutate(y_bin = floor(y / y_bin_size) * y_bin_size)
  
  # -------------------------
  # 3️⃣ Area calculation
  # -------------------------
  cell_area_km2 <- (cell_size * cell_size) / 1e6
  
  lat_summary <- pts %>%
    group_by(y_bin, class) %>%
    summarise(area_km2 = n() * cell_area_km2, .groups = "drop")
  
  # -------------------------
  # 4️⃣ If axis == "lat", convert y_bin → latitude
  # -------------------------
  if (axis == "lat") {
    
    # Dummy X column; real Y from y_bin
    lat_coords <- data.frame(x = 0, y = lat_summary$y_bin)
    
    lat_sf <- st_as_sf(lat_coords, coords = c("x", "y"), crs = 3005) %>%
      st_transform(4326)
    
    lat_summary$lat <- st_coordinates(lat_sf)[,2]
  }
  
  # -------------------------
  # 5️⃣ Set default palette if not provided
  # -------------------------
  n_classes <- length(levels(lat_summary$class))
  
  if (is.null(palette)) {
    palette <- scales::hue_pal()(n_classes)
  }
  
  # -------------------------
  # 6️⃣ Base ggplot
  # -------------------------
  x_var <- if (axis == "lat") "lat" else "y_bin"
  
  p <- ggplot(lat_summary,
              aes_string(x = x_var,
                         y = "area_km2",
                         fill = "class")) +
    geom_col(position = "identity", alpha = 0.6) +
    scale_fill_manual(values = palette) +
    labs(
      x = if (axis == "lat") "Latitude (°)" else "Y coordinate (m)",
      y = "Area (km²)",
      fill = "Class"
    ) +
    theme_bw() +
    theme(
      # panel.background = element_rect(fill = "aliceblue"),
      legend.position = "right",
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    )
  
  # -------------------------
  # 7️⃣ Optional facetting
  # -------------------------
  if (facet) {
    p <- p + facet_wrap(~ class, nrow = 1)
  }
  
  # -------------------------
  # 8️⃣ Flip axes if requested
  # -------------------------
  if (flip) {
    p <- p + coord_flip()
  }
  
  return(p)
}


r<- rast(paste(out_dir, "stability_layer_categoric.tif", sep="/"))
p_bar <- make_latitude_area_plot(
  r,
  y_bin_size = 2000,
  cell_size = 20,
  palette = c("blue", "yellow", "orange", "red"),   # optional
  axis = "y",
  facet = T
)

p_bar


category = c("Instable", "Stable-Moderate", "Stable-Recommended", "Stable-Ideal")

figures_path<- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas"
# ggsave(paste(figures_path, "Figure_stability_class_barplot.pdf", sep="/"), width = 15, height = 16, units="cm", dpi = 300)



library(patchwork)
summary(lat_summary$y_bin)

# Combine 
combined_plot <- p_map | p_bar + 
  plot_layout(widths = c(2, 1)) # map wider than bar plot

# Display
combined_plot


# ggsave(paste(model_results_path, "Figure_testing_records_maps.pdf", sep="/"), width = 15, height = 17, units="cm", dpi = 300)

