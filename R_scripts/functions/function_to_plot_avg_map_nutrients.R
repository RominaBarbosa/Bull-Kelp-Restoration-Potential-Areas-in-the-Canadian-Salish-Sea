library(dplyr)
library(ggplot2)
library(sf)
library(ggspatial)
library(stringr)

# FUNCTION TO PLOT FACETED SEASONAL AVERAGES (2014–2019)
plot_seasonal_avg_map <- function(nutrient = "nitrate",
                                  data_sf,
                                  land,
                                  coastline,
                                  bbox) {
  
  # Select the proper column for each nutrient
  nutrient_col <- case_when(
    nutrient == "nitrate" ~ "avg_nitrate",
    nutrient == "phosphate" ~ "avg_phosphate",
    nutrient == "silicate" ~ "avg_silicate",
    TRUE ~ stop("Invalid nutrient name")
  )
  
  # Filter years and calculate seasonal averages
  seasonal_avg <- data_sf %>%
    filter(Year >= 2014, Year <= 2019) %>%
    select(Station_name, cluster_id, geometry,
           starts_with("avg_"), Year) %>%
    pivot_longer(cols = starts_with("avg_"),
                 names_to = c(".value", "season"),
                 names_pattern = "avg_(.*)_(.*)") %>%
    filter(!is.na(!!sym(nutrient))) %>%
    group_by(Station_name, cluster_id, season) %>%
    summarise(nutrient_mean = mean(!!sym(nutrient), na.rm = TRUE),
              geometry = first(geometry),
              .groups = "drop") %>%
    st_as_sf()
  
  # Plot faceted map
  ggplot() +
    geom_sf(data = land, fill = "lightgray", color = "gray") +
    geom_sf(data = coastline, color = "darkgrey", size = 0.2) +
    geom_sf(data = seasonal_avg,
            aes(fill = nutrient_mean, shape = factor(cluster_id)),
            color = "black", size = 2, stroke = 0.8) +
    scale_fill_viridis_c(name = paste("Mean", nutrient, "(µM)"),
                         option = "plasma", direction = -1) +
    scale_shape_manual(values = 0:(length(unique(seasonal_avg$cluster_id)) - 1),
                       name = "Cluster ID") +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3) +
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                           style = north_arrow_fancy_orienteering) +
    labs(title = paste("Average", str_to_title(nutrient), "by Season (2014–2019)"),
         x = "Longitude", y = "Latitude") +
    facet_wrap(~ season) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    )
}
