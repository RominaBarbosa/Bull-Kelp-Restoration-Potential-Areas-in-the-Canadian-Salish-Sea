###====================================================================================
###     Species Distribution models    SDMs                                      ######
###                                                                              ######
### 7.2- Assessing the Env conditions at the Areas for Restoration               ######
### Content:                                                                     ######
### Statistical differences among Restoration categories, for each env. variable ######
### Export tables of stats results and summary env conditions                    ######
### Figure 7, boxplot of env conditons ~ Restoration category (Facet_wrap ~ Variable) #
### Author: Romina Barbosa                                                       ######
### Date last edition: 26-Nov-2025                                               ######
###====================================================================================
# Load packages
library(terra)
library(dplyr)
library(dismo)
library(ggplot2)
library(raster)
library(tidyr)


# Paths
model_results_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_results"
setwd(model_results_path)
postblob_path<- "/Volumes/Romina_PSF/PSF/SDM/SDM_predict_postBlob"
out_dir <- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7"

### Load substrate layer =====================================================
substrate_aligned<- rast(paste(model_results_path, "substrate_SOG_aligned.tif", sep="/"))
substrate_aligned[substrate_aligned == 2]<- NA
substrate_aligned <- terra::mask(substrate_aligned, ens_average)

# Create Mask from bathymetry
bathy20m<- rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/Topographic_Variables/topographic_variables_20mres/coastwide_20m.tif")
bathy20m_mask10_30<- bathy20m
bathy20m_mask10_30[!(bathy20m_mask10_30 >= -10 & bathy20m_mask10_30 <= 30)] <- NA # negative values are in land

bathy20m_mask10_15<- bathy20m
bathy20m_mask10_15[!(bathy20m_mask10_15 >= -10 & bathy20m_mask10_15 <= 15)] <- NA # negative values are in land


### Summarize environmental conditions during each periods at each Restoration category ====
# Load raster stack of variables resampled at 20 m resoltuion
r<- rast( paste(out_dir, "restoratio_potential_masked_depth.tif", sep="/"))
r<- terra::mask(r, bathy20m_mask10_15)
r_df <- as.data.frame(r, xy = TRUE)  # adds x and y columns
colnames(r_df)[3] <- "category"
r_df$category <- as.factor(r_df$category)


raster_stack_20m<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_blob.tif")
names(raster_stack_20m)

cummulated_18SST<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated_Blob/summer_cumulated_degrees_18_interpolated_output.tif")
hourAbove_18SST<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interpolated_Blob/summer_hours_above_threshold_18_interpolated_output.tif")
temperature_layers<- c(cummulated_18SST, hourAbove_18SST)
names(temperature_layers)<- c("cummulated_18SST", "hourAbove_18SST")


raster_stack_20m_masked <- terra::mask(raster_stack_20m, bathy20m_mask10_15)
plot(raster_stack_20m_masked[[1]])
# writeRaster(raster_stack_20m_masked, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/raster_stack_blob_masked.tif", overwrite=T)

raster_stack_postblob<- terra::rast("/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/SalishSeaCast_interp_20m_postblob.tif")
names(raster_stack_postblob)

raster_stack_postblob_masked <- terra::mask(raster_stack_postblob, bathy20m_mask10_15)
plot(raster_stack_postblob_masked[[1]])
# writeRaster(raster_stack_postblob_masked, "/Volumes/Romina_PSF/PSF/SDM/environmental_layers/SalishSeaCast_interp_20m_resolution/raster_stack_postblob_masked.tif", overwrite=T)


# Extract env conditions at each cell and merge with Restoration category
# Select the period, period 1 or period 2 (commented lines)
# df_env <- cbind(
#   r_df,
#   terra::extract(raster_stack_20m, r_df[, c("x", "y")])
# )
# 
# df_env <- cbind(
#   df_env,
#   terra::extract(temperature_layers, df_env[, c("x", "y")])
# )

df_env <- cbind(
  r_df,
  terra::extract(raster_stack_postblob_masked, r_df[, c("x", "y")])
)

# Remove NA classes
df_env <- df_env[!is.na(df_env$category), ]
df_env$class <- as.factor(df_env$category)


##===============================================================================##
### Assess differences on environmental conditions among restoration categories ===
# Convert from wide to long table
df_long_raw <- df_env %>% 
  pivot_longer(cols = !c(x, y, category, ID, class), names_to = "variable", values_to = "value")

levels(df_long_raw$category)<- c("Non-viable",  "Low", "Medium", "High" )


df_long_raw$class<- as.factor(df_long_raw$class)
levels(df_long_raw$class)<- c("Non-viable",  "Low", "Medium", "High" )

summary(df_long_raw)

box_stats <- df_long_raw %>%
  group_by(category, class, variable) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    data_min = min(value, na.rm = TRUE),
    data_max = max(value, na.rm = TRUE),
    lower_whisker = max(data_min, q1 - 1.5 * iqr),
    upper_whisker = min(data_max, q3 + 1.5 * iqr),
    .groups = "drop"
  )

# write.csv(box_stats, "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7/Table_boxplot_stats_Figure8_depthMasked.csv")
# write.csv(box_stats, "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas/stability_outputs_M7/Table_boxplot_stats_postblob_Figure8_depthMasked.csv")

box_stats%>% 
  mutate(mean= round(mean, 2), sd= round(sd,2))%>%
  dplyr::select(category, variable, mean, sd)

### Perform statistical test of differences among RESTORATION CATEGORIES =======
# SInce htere are too many rows, I performed sampling before testing for differences
library(rstatix)
library(data.table)
library(PMCMRplus)
setDT(df_long_raw)

# Function to process one variable
# # process_variable <- function(varname= i, dt = df_long_raw, alpha = 0.05) {
# 
#   library(PMCMRplus)
#   library(data.table)
#   library(multcompView)
# 
#   # 1. Subset data
#   df_var <- dt[variable == varname]
# 
#   # 2. Kruskal-Wallis Dunn Test
#   dunn_res <- kwAllPairsDunnTest(x = df_var$value,
#                                  g = as.factor(df_var$category),
#                                  p.adjust.method = "holm")
# 
#   # 3. Convert p-value matrix to long format
#   p_mat <- as.matrix(dunn_res$p.value)
#   df_pvals <- as.data.table(as.table(p_mat))
#   setnames(df_pvals, c("group1","group2","p.value"))
# 
#   # remove self-comparisons
#   df_pvals <- df_pvals[group1 != group2]
# 
#   # 4. Ensure matrix is symmetric for multcompLetters
#   all_groups <- unique(df_var$category)
#   p_mat_full <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
#                        dimnames = list(all_groups, all_groups))
# 
#   for(i in 1:nrow(df_pvals)) {
#     g1 <- df_pvals$group1[i]
#     g2 <- df_pvals$group2[i]
#     p_val <- df_pvals$p.value[i]
#     p_mat_full[g1, g2] <- p_val
#     p_mat_full[g2, g1] <- p_val
#   }
#   
#   p_mat_full[is.na(p_mat_full)] <- 1
#   
#   # 5. Generate compact letter display
# 
#   cld <- multcompLetters(p_mat_full)
# 
#   letters_dt <- data.table(
#     variable = varname,
#     category = names(cld$Letters),
#     letter = as.vector(cld$Letters)
#   )
# 
#   # 6. Return both letters and p-values
#   return(list(
#     letters = letters_dt,
#     p_values = df_pvals
#   ))
# }

process_variable <- function(varname, dt = df_long_raw, alpha = 0.05) {

  library(PMCMRplus)
  library(agricolae)
  library(data.table)
  library(multcompView)

  # 1. Subset and sample data
  df_var <- dt[variable == varname]
  # df_var <- df_var[, .SD[sample(.N, min(.N, max_n))], by = category]

  # 2. Kruskal-Wallis
  pval.matrix <- kwAllPairsDunnTest(x = df_var$value,
                                    g = as.factor(df_var$category), p.adjust.method = "holm")

  #square and diagonalize the p-value matrix
  new.pval.matrix <- rbind(1,pval.matrix$p.value)
  new.pval.matrix <- cbind(new.pval.matrix, 1)
  diag(new.pval.matrix) <- 1
  library(Matrix)

  new.pval.matrix <- as.matrix(forceSymmetric(new.pval.matrix, "L"))

  #Add Instable to the row and column names
  rownames(new.pval.matrix)[dim(pval.matrix$p.value)+1] <-
    rownames(pval.matrix$p.value)[dim(pval.matrix$p.value)[1]]
  colnames(new.pval.matrix)[dim(pval.matrix$p.value)+1] <-
    rownames(pval.matrix$p.value)[dim(pval.matrix$p.value)[1]]


  #  Convert p-value matrix to long format
  df_pvals <- as.data.table(as.table(new.pval.matrix))
  setnames(df_pvals, c("group1","group2","p.value"))

  # remove self-comparisons
  df_pvals <- df_pvals[group1 != group2]


  # 5. Generate compact letter display
  # 1. Get all pairwise combinations
  categories <- c("Non-viable", "Low", "Medium", "High")
  all_combs <- t(combn(categories, 2))
  all_combs <- rbind(all_combs, all_combs[,2:1])  # make symmetric
  all_combs <- as.data.frame(all_combs)
  colnames(all_combs) <- c("group1", "group2")

  # 2. Merge with actual p-values
  df_full <- merge(all_combs, df_pvals, by=c("group1","group2"), all.x=TRUE)

  # 3. Fill missing p-values with a small number (optional)
  # df_full$p.value[is.na(df_full$p.value)] <- 1e-6


  # 4. Convert to matrix
  p_mat <- reshape2::acast(df_full, group1 ~ group2, value.var="p.value")
  p_mat <- as.matrix(forceSymmetric(p_mat, "L"))  # now proper symmetric

  # 5. Generate letters
  cdl <- multcompLetters(p_mat)

  # Convert to a proper data.table
  letters_dt <- data.table(
    variable = varname,
    category = names(cdl$Letters),
    letter = as.vector(cdl$Letters)
  )

  df_full$p.value <- format(df_full$p.value, scientific = TRUE)
  return(list(df_full, letters_dt))
}

# Loop over all variables
variables <- unique(df_long_raw$variable)
letters_list <- vector("list", length(variables))
pairwise_DunnTest_list<-  vector("list", length(variables))

for(i in seq_along(variables)) {
  out <- process_variable(variables[i], df_long_raw)
  letters_list[[i]] <- out[[2]]   # out must be a data.table
  pairwise_DunnTest_list[[i]] <- out[[1]]   # out must be a data.table
  cat("Processed:", variables[i], "\n")
}


# Combine all letters
options(digits = 13)
letters_df <- rbindlist(letters_list)
pairwise_DunnTest_df <- rbindlist(pairwise_DunnTest_list)

box_stats2 <- merge(
  box_stats,
  letters_df,
  by = c("variable", "category"),
  all.x = TRUE
)%>%mutate(DunnTest_Letter=letter)



# Function to compute pairwise Cliff's delta for a variable

library(data.table)
library(rstatix)
library(PMCMRplus)
library(rcompanion)  # for cldList
library(effsize)     # for Cliff’s delta
library(DescTools)  # for Cliff's delta


# Define categories and threshold for effect size


cliffs_delta_fast <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  # Convert to numeric early to avoid integer overflow
  nx <- as.numeric(length(x))
  ny <- as.numeric(length(y))
  
  xy <- c(x, y)
  r <- data.table::frank(xy, ties.method = "average")  # fast ranking
  
  rx <- r[seq_len(nx)]
  ry <- r[(nx + 1):(nx + ny)]
  
  # Mann–Whitney U statistic (use numeric to avoid integer overflow)
  U <- sum(rx) - nx * (nx + 1) / 2
  
  # Convert U → Cliff's delta (safe from overflow)
  delta <- (2 * U) / (nx * ny) - 1
  
  return(delta)
}

compute_effect_sizes <- function(dt_var, 
                                 categories = c("Non-viable", "Low", "Medium", "High")) {
  
  cats <- intersect(categories, unique(dt_var$category))
  all_pairs <- t(combn(cats, 2))
  
  res_list <- vector("list", nrow(all_pairs))
  
  for(i in seq_len(nrow(all_pairs))) {
    g1 <- all_pairs[i,1]
    g2 <- all_pairs[i,2]
    
    x <- dt_var[category == g1, value]
    y <- dt_var[category == g2, value]
    
    delta <- cliffs_delta_fast(x, y)
    
    res_list[[i]] <- data.table(
      group1 = g1,
      group2 = g2,
      cliff_delta = delta
    )
  }
  
  rbindlist(res_list)
}

# Example: compute effect sizes for one variable

all_vars <- unique(df_long_raw$variable)
all_effects <- list()

for(v in all_vars){
  dt_var <- df_long_raw[variable == v]
  es <- compute_effect_sizes(dt_var)
  es[, variable := v]
  all_effects[[v]] <- es
}

effect_sizes_all <- rbindlist(all_effects)


effect_sizes_all[, magnitude := fifelse(abs(cliff_delta) < 0.147, "negligible",
                                        fifelse(abs(cliff_delta) < 0.33, "small",
                                                fifelse(abs(cliff_delta) < 0.474, "medium", "large")))]


library(data.table)


# Example threshold for practical significance
delta_threshold <- 0.147  # Cliff's delta small effect threshold

# Step 1: Mask p-values where effect size is too small
effect_sizes_all[, masked_p := ifelse(abs(cliff_delta) >= delta_threshold, 0.05, 1)]
# 0.05 here is arbitrary to indicate "significant", 1 = "not significant"

# Step 2: Loop over variables to generate letters
categories <- c("Non-viable", "Low", "Medium", "High")
letters_list <- list()

for(v in variables) {
  df_var <- effect_sizes_all[variable == v]
  
  # Initialize p-matrix
  p_mat <- matrix(1, nrow = length(categories), ncol = length(categories),
                  dimnames = list(categories, categories))
  
  # Fill with masked_p
  for(i in 1:nrow(df_var)) {
    g1 <- df_var$group1[i]
    g2 <- df_var$group2[i]
    p_mat[g1, g2] <- df_var$masked_p[i]
    p_mat[g2, g1] <- df_var$masked_p[i]
  }
  
  # Convert to binary for multcompLetters
  p_mat_bin <- ifelse(p_mat <= 0.05, 0, 1)
  
  # Generate letters
  letters <- multcompLetters(p_mat_bin)$Letters
  letters_list[[v]] <- data.table(category = names(letters), Letter = letters, variable = v)
}

letters_df_effect <- rbindlist(letters_list)
colnames(letters_df_effect)[2]<- "letter_effectSize"

# Step 2: Merge by variable and category
effect_sizes_df<- as.data.frame(effect_sizes_all)


letters_combined <- merge(
  box_stats2,
  letters_df_effect,
  by = c("variable", "category"),
  all.x = TRUE
)

# colnames(letters_combined)[15]<- "letter_DunnTest_holmCorrect"


### Save TABLE for Figure 7 and supplementary information (pairwise comparison results)
formatted_table <- letters_combined %>%
  mutate(
    mean_sd = sprintf("%.2f ± %.2f", mean, sd),
    mean    = round(mean, 2),
    sd      = round(sd, 2),
    area_km2= round(n*0.02*0.02, 2) 
    # median  = round(median, 2)
  ) %>%
  select(variable, category, class, n, mean_sd, area_km2)

# write.csv(letters_combined, paste(figures_path, "table_Figure7_RestorationHabitats_EnvPredictors_barplot_v2.csv", sep="/"))
# write.csv(effect_sizes_all, paste(figures_path, "table_Figure7_cliff_delta_depthMasked.csv", sep="/"))
# write.csv(pairwise_DunnTest_df, paste(figures_path, "table_Figure7_pairwise_DunnTest_depthMasked.csv", sep="/"))
# write.csv(formatted_table, paste(figures_path, "table_Figure8_formated_depthMasked.csv", sep="/"))

# write.csv(letters_combined, paste(figures_path, "table_Figure8_RestorationHabitats_EnvPredictors_barplot_postblob_depthMasked.csv", sep="/"))
# write.csv(effect_sizes_all, paste(figures_path, "table_Figure8_cliff_delta_postblob_depthMasked.csv", sep="/"))
# write.csv(pairwise_DunnTest_df, paste(figures_path, "table_Figure8_pairwise_DunnTest_postblob_depthMasked.csv", sep="/"))
# write.csv(formatted_table, paste(figures_path, "table_Figure8_postblob_formated_depthMasked.csv", sep="/"))


figures_path<- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas"
# letters_combined_period1<- read.csv(paste(figures_path, "table_Figure7_RestorationHabitats_EnvPredictors_barplot_v2.csv", sep="/"))
letters_combined_period1<- read.csv(paste(figures_path, "table_Figure7_RestorationHabitats_EnvPredictors_barplot_postblob.csv", sep="/"))


# BOXPLOTS WITH LETTERS
labels <- c(
  "PAR_summer_mean"           = "Mean Summer PAR (W m-2)",
  "ammonium_spring_mean"      = "Mean Spring Ammonium (μM)",
  # "ammonium_summer_minimum"   = "Minimum Summer Ammonium",
  "currentSpeed_summer_mean"  = "Mean Summer Current Speed (m s-1)",
  "nitrate_summer_minimum"    = "Minimum Summer Nitrate (μM)",
  "temperature_summer_mean"   = "Mean Summer SST (°C)",
  "turbidity_summer_mean"     = "Mean Summer Turbidity (NTU)",
  "salinity_summer_mean"      = "Mean Summer Salinity (psu)"
  # "hourAbove_18SST"           = "# of hours above 18 °C (hours)",
  # "cummulated_18SST"          = "Accumulated degree hours above 18 °C (°C)"
)


levels(as.factor(letters_combined_period1$variable))

letters_combined_period1<- letters_combined_period1 %>%
  mutate(category = recode(category,
                           "Instable"   = "Non-viable",
                           "Moderate"   = "Low",
                           "Recommended"= "Medium",
                           "Ideal"      = "High"))


letters_combined_period1 %>%
  mutate(category = factor(category,
                           levels = c("Non-viable", "Low", "Medium", "High")))%>%
  filter(variable != "nitrate_winter_mean")%>%
  # filter(variable != "nitrate_summer_minimum")%>%
  filter(variable != "ammonium_summer_minimum")%>%
  filter(variable != "ammonium_winter_mean")%>%
  mutate(variable = factor(variable,
                           levels = c("turbidity_summer_mean", "temperature_summer_mean",
                                      "PAR_summer_mean", "salinity_summer_mean", "ammonium_spring_mean",
                                      "currentSpeed_summer_mean", "nitrate_summer_minimum",
                                      "hourAbove_18SST", "cummulated_18SST")))%>%
ggplot(aes(x = category, fill = category)) +
  geom_boxplot(
    aes(
      ymin = lower_whisker,
      lower = q1,
      middle = median,
      upper = q3,
      ymax = upper_whisker
    ),
    stat = "identity", width =0.75
  ) +
  geom_text(
    aes(x = category, y = upper_whisker, label = letter_effectSize),
    hjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(values = c(
    "Non-viable"  = "#882E72",    # purple
    "Low"    = "#A6D854",    # light green
    "Medium"      = "#1B9E77",    # colorblind-friendly green-teal
    "High"       = "#1B7837"     # dark green
  )) +
  #scale_fill_manual(values=c("purple", "greenyellow", "springgreen3", "forestgreen")) + #"yellow", "orange"
  geom_point(aes(x = category, y = mean)) +
  labs(title = "", fill = "Restoration Category", y= "Predictor value", x="") +
  facet_wrap(~ variable, scales = "free_y", nrow = 3, labeller = as_labeller(labels))  +
  theme_bw() +
  theme(
    legend.position = "null",#"bottom",
    # legend.direction = "horizontal",
    axis.text.x = element_text(size = 8, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 7, face = "bold", color = "black")
  )

figures_path<- "/Volumes/Romina_PSF/PSF/SDM/Restoration_areas"
# ggsave(paste(figures_path, "Figure_8_RestorationHabitats_EnvPredictors_barplot_postblob_depthMasked.pdf", sep="/"),
#        width = 19, height = 12, units="cm", dpi = 300)

# ggsave(paste(figures_path, "Figure_8_RestorationHabitats_EnvPredictors_barplot_blob_v2.pdf", sep="/"),
# width = 19, height = 12, units="cm", dpi = 300)





letters_combined_period1$period<- "Period1"
letters_combined_period2$period<- "Period2"

colnames(letters_combined_period1)
colnames(letters_combined_period2)[17]<- "letter_DunnTest_holmCorrect"
letters_combined_period2<- letters_combined_period2[,-16]

letters_combined_2periods<- rbind(letters_combined_period1, letters_combined_period2)
str(letters_combined_2periods)

library(ggpattern)

letters_combined_2periods %>%
  filter(variable != "nitrate_winter_mean",
         variable != "nitrate_summer_minimum",
         variable != "ammonium_summer_minimum",
         variable != "ammonium_winter_mean") %>%
  mutate(period = factor(period)) %>%   # ensure period is a factor
  # reorder category here
  mutate(
    period   = factor(period),
    category = factor(category,
                      levels = c("Instable", "Moderate", "Recommended", "Ideal"))
  ) %>%
  ggplot(aes(x = category)) +
  
  geom_boxplot_pattern(
    aes(
      ymin   = lower_whisker,
      lower  = q1,
      middle = median,
      upper  = q3,
      ymax   = upper_whisker,
      fill   = category,   # ← keep your original colors
      pattern = period     # ← texture for each period
    ),
    stat     = "identity",
    width    = 0.7,
    position = position_dodge(width = 0.9),
    pattern_fill = "black",          # outline color of texture
    pattern_colour = "black",
    pattern_density = 0.3,           # adjust
    pattern_spacing = 0.02
  ) +
  
  # Text (Dunn letters)
  geom_text(
    aes(y = upper_whisker, label = letter_effectSize),
    position = position_dodge(width = 0.9),
    hjust = -0.2,
    size = 3
  ) +
  
  # mean points
  geom_point(
    aes(y = mean, group = period),
    position = position_dodge(width = 0.9),
    size = 1.8
  ) +
  
  scale_fill_manual(values = c(
    "Instable"    = "lightblue",
    "Moderate"    = "green",
    "Recommended" = "darkgreen",
    "Ideal"       = "purple"
  )) +
  
  scale_pattern_manual(values = c(
    "Period1" = "stripe",
    "Period2" = "crosshatch"
  )) +
  
  facet_wrap(~ variable, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    strip.text = element_text(size = 8, face = "bold")
  )


