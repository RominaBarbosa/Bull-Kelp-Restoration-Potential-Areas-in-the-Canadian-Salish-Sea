


library(terra)     # for raster handling
library(ggplot2)   # for plotting
library(dplyr)     # for data manipulation
library(tidyr)     # to reshape data

# Adjust the path to your raster files
my_path<-("/Volumes/Romina_PSF/PSF/modeled_variables_original/seasonal_metrics_0.5m_depth")
raster_files <- list.files(path = my_path, pattern = "\\.csv$", full.names = TRUE)
# raster_files<- raster_files[]
# raster_files <- strsplit(tools::file_path_sans_ext(raster_files), "/")[[1]][7]



metric=  "18"
threshold= "hours_above_threshold"
season= "summer"

metrics<- c()
for (i in 2014:2024) {
  file_name_i= paste("temperature", i, season, threshold,  metric, "surface", ".csv", sep="_")
  # file_name_i<- paste(my_path, file_name_i, sep="/")
  metrics<- c(metrics, file_name_i)
}


raster_files<- raster_files[basename(raster_files) %in% metrics]

all_data <- lapply(raster_files, function(f) {
  df <- read.csv(f)
  df<- df[,-1]
  
  # Extract year from filename
  year <- gsub(".*temperature_(\\d+)_.*", "\\1", basename(f))
  
  # Rename metric column if needed
  value_col <- setdiff(names(df), c("lat", "lon", "latitude", "longitude", "year"))[1]
  names(df)[names(df) == value_col] <- "value"
  
  # Add year column
  df$year <- year
  return(df)
}) %>% bind_rows()

summary(as.factor(all_data$year))

all_data$year<- as.factor(all_data$year)

data_20214<- all_data%>%
  filter(year==2014)

plot(data_20214$value)

all_data%>%
  filter(value>=1)%>%
  filter(value<=2000)%>%
ggplot( aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ year, scales = "free_y", ncol=1) +  # scales = "free_y" allows different y-axis scales
  labs(x = "Hours above threshold", y = "Frequency",
       title = "Distribution of Metric by Year") +
  theme_minimal()



