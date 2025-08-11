###==================================================================
### Species Distribution models        SDMs          ################
###                                                  ################
###                                                  ################
### Author: Romina Barbosa                           ################
### Last version: 11-Aug-2025                        ################
###==================================================================
library(tidyr)
library(stringr)

variables_selection_path<- "/Volumes/Romina_PSF/PSF/SDM/Variables_selection"

# Load kelp records ============================================================
kelp_presabs_df<- read.csv("/Volumes/Romina_PSF/PSF/SDM/Presence_absences_kelp/presence_2014_to_2019_filtered_plus_Absences.csv")
colnames(kelp_presabs_df)
kelp_presabs_df<- kelp_presabs_df[,c(2,3,11:14)]
head(kelp_presabs_df)


# Load selected environmental variables ========================================



# Interpolate variables ===







# SSCast_variables<- df_all
# # Transform your spatial object 
# coordinates(SSCast_variables) <- ~longitude + latitude  # Convert to spatial object
# proj4string(SSCast_variables) <- CRS("+proj=longlat +datum=WGS84")  # Set CRS (WGS84)
# 
# # Convert to sf with original long/lat and WGS84 CRS first
# sf_points <- st_as_sf(SSCast_variables, coords = c("longitude", "latitude"), crs = 4326)
# 
# bc_albers <- CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs")
# SSCast_bc <- spTransform(SSCast_variables, bc_albers)





stack_nemo_vars<- list()
for (i in selected_variables){
  
  variable<- read.csv()
  
  stack_nemo_vars[i]<- variable
}










terrain_vars<- terra::rast(paste(variables_selection_path,"terrain_rasters_selected.tif", sep="/"))
names(terrain_vars)
# "easterness_3x3"  "northerness_3x3" "slope_3x3"       "slope_7x7"       "TPI_3x3"   










# Conver table to long table for plotting and summarizing data =================

kelp_long <- kelp_presabs_df %>%
  pivot_longer(cols = -kelp, names_to = "variable", values_to = "value")

str(kelp_long)
kelp_long$variable_code<- as.factor(kelp_long$variable)
vars<- levels(kelp_long$variable_code)


levels(kelp_long$variable_code)<- stringr::word(vars, c(1,2), sep = "_")

# Extract first two parts (or just the first if only one exists)
short_names <- sapply(strsplit(vars, "_"), function(x) {
  paste(x[1:min(2, length(x))], collapse = "_")
})

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