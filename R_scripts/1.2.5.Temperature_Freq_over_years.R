

library(tidyverse)

# Load the data
df <- read.csv("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_blob_/temperature_blob_summer_mean_surface_.csv")
dfpostblob<- read.csv("/Volumes/Romina_PSF/PSF/modeled_variables_original/climatology_metrics_post_blob/temperature_post_blob_summer_mean_surface_.csv")
head(df)

df_temp<- merge(df, dfpostblob, by = c("gridY", "gridX", "bathymetry", "latitude", "longitude"))

df_long <- df_temp %>%
  pivot_longer(
    cols = starts_with("value_"),
    names_to = "year",
    names_prefix = "value_",
    values_to = "temperature"
  ) %>%
  filter(!is.na(temperature))  # remove NA values


ggplot(df_long, aes(x=temperature)) +
  geom_histogram(binwidth=.25, fill="steelblue", color="black", alpha=0.7) +
  facet_wrap(~year, ncol=1, scales="free_y") +  # one panel per year
  labs(title="Temperature distribution by year",
       x="Temperature (°C)",
       y="Number of cells") +
  theme_minimal()


ggsave("/Volumes/Romina_PSF/PSF/modeled_variables_original/Temperature_summerMean_over_time.png", width = 20, height = 25, dpi= 300, units="cm")


ggplot(df_long, aes(x=temperature, y=..density.., fill=year)) +
  geom_density(alpha=0.5) +
  facet_wrap(~year,ncol=1) +
  labs(title="Temperature density by year",
       x="Temperature (°C)",
       y="Density") +
  theme_minimal() +
  theme(legend.position="none")
