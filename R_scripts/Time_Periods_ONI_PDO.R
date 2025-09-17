





### ANALYSE period of time based on ONI and PDO indices --> from Mora-Soto et al. (2024)
# data downloaded from: PDO: https://psl.noaa.gov/pdo/ and ONI: https://psl.noaa.gov/data/timeseries/month/DS/ONI/ 

# Mora-Soto et al. 2024 - methods: "The global scale was represented by ONI (NOAA, 2023a) 
# and PDO (NOAA, 2023b). In general, we consider that optimal periods for kelp were when 
# the PDO or ONI was negative, while positive PDO or ONI were classified as suboptimal 
# periods (Dayton et al., 1999; Bell et al., 2015; Pfister et al., 2018; Cavanaugh et al., 2019; Gendall, et al., in prep). 
# To implement this classification, ONI and PDO values from spring (May to June) 
# and summer (July to August) 2002 to 2022 were rescaled using z-scores. 
# The median value of z-scored ONI (-0.05) and PDO (-0.26) were the parameters to distinguish optimal versus suboptimal periods. 
# The Kruskal-Wallis test (Kruskal and Wallis, 1952) was used to identify if the 
# difference among the periods was statistically significant, whereas the Dunn test (Dunn, 1961) 
# was used as a post hoc analysis to define the groups that were statistically different."


ONI_PDO<- read.csv("/Volumes/Romina_PSF/PSF/PDO_ONI_data/indices_RVB2025.csv", sep=";")
ONI_PDO<- ONI_PDO[,c(1:5,9,10)]
colnames(ONI_PDO)[6:7]<- c("PDO", "ONI") 
str(ONI_PDO)

spring_indices<- ONI_PDO%>%
  filter(month == 5 | month == 6 )%>%
  mutate(season = "Spring")

summer_indices<- ONI_PDO%>%
  filter(month == 7 | month == 8 )%>%
  mutate(season = "Summer")


ONI_PDO_season<- rbind(spring_indices, summer_indices)
ONI_PDO_season<- ONI_PDO_season%>%
  group_by(year, season)%>%
  summarize(ONI= mean(ONI), PDO= mean(PDO))

baseline <- subset(ONI_PDO_season, year >= 2002 & year <= 2022)
future   <- subset(ONI_PDO_season, year >= 2023 & year <= 2025)


#Compute mean & SD from baseline
# ONI baseline stats
oni_mean <- mean(baseline$ONI, na.rm = TRUE)
oni_sd   <- sd(baseline$ONI, na.rm = TRUE)

# PDO baseline stats
pdo_mean <- mean(baseline$PDO, na.rm = TRUE)
pdo_sd   <- sd(baseline$PDO, na.rm = TRUE)

#Rescale both baseline and future applying the same parameters:
# Baseline z-scores
baseline$ONI_z <- (baseline$ONI - oni_mean) / oni_sd
baseline$PDO_z <- (baseline$PDO - pdo_mean) / pdo_sd

# Future z-scores (using SAME mean & sd)
future$ONI_z <- (future$ONI - oni_mean) / oni_sd
future$PDO_z <- (future$PDO - pdo_mean) / pdo_sd

# Combine datasets
df_all <- rbind(baseline, future) # df_all has z-scored values for 2002–2025, all rescaled consistently using 2002–2022 as the reference period.

### Classification of time periods with thresholds ============
# I apply the same thresholds from Mora-Soto et al. (2024) (−0.05 for ONI, −0.26 for PDO):
df_all$ONI_group <- ifelse(df_all$ONI_z <= -0.05, "Optimal", "Suboptimal")
df_all$PDO_group <- ifelse(df_all$PDO_z <= -0.26, "Optimal", "Suboptimal")

df_all$season<- as.factor(df_all$season)

# PLOT
# Reshape to long format
df_long <- df_all %>%
  select(year, season, ONI_z, PDO_z, ONI_group, PDO_group) %>%
  pivot_longer(cols = c(ONI_z, PDO_z),
               names_to = "Index",
               values_to = "Z_value") %>%
  mutate(Index = factor(Index, levels = c("ONI_z","PDO_z"), labels = c("ONI","PDO")),
         xyear = year + ifelse(season == "Spring", -0.25, 0.25))

# Quick check
df_long <- df_long %>%
  mutate(Period = case_when(
    year >= 2002 & year <= 2006 ~ "Suboptimal 1",
    year >= 2007 & year <= 2013 ~ "Optimal 1",
    year >= 2014 & year <= 2019 ~ "Suboptimal 2",
    year >= 2020 & year <= 2025 ~ "Optimal 2",
    TRUE ~ NA_character_
  ))



# compute period bounds for rect/labels (in year units)
period_bounds <- df_long %>%
  group_by(Period) %>%
  summarize(
    xmin = min(year) - 0.5,
    xmax = max(year) + 0.5,
    .groups = "drop"
  )

period_bounds <- period_bounds %>%
  arrange(xmin) 

period_vlines <- period_bounds %>%
  arrange(xmin) %>%
  slice(-n()) %>%          # remove last period
  pull(xmax)

ggplot(df_long, aes(x = year, y = Index, fill = Z_value)) +
  geom_tile(width = 1, height = 0.9, color = "grey70") +
  geom_rect(data = period_bounds,
            aes(xmin = xmin, xmax = xmax, ymin = 2.25, ymax = 2.5),
            inherit.aes = FALSE, fill = "grey90", color = NA, alpha = 0.6) +
  scale_fill_gradient2(low = "#1f78b4",    # strong blue
                       mid = "#ffffff",    # white
                       high = "#e31a1c",   # strong red, 
                       midpoint = 0) +
  scale_x_continuous(breaks = seq(min(df_long$year), max(df_long$year), by = 1)) +
  # Add vertical dashed lines at period boundaries
  geom_vline(xintercept = period_vlines, linetype = "dashed", color = "black", size = 0.8) +
  geom_text(data = period_bounds,
            aes(x = (xmin + xmax)/2, y = 2.5, label = Period),
            inherit.aes = FALSE, size = 3.5) +
  theme_minimal(base_size = 11.5) +
  facet_wrap(~season, nrow=2)+
  labs(title = "",
       x = "Year", y = "", fill = "Z-score") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# save
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/plot_ONI_PDO.png", width = 17, height = 7, units="cm",dpi = 300)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/plot_ONI_PDO.pdf", width = 17, height = 7, units="cm",dpi = 300)



ggplot(df_long, aes(x = year, y = Z_value, color = Index, group = interaction(Index, season))) +
  geom_line(aes(linetype = season), size = 0.8) +
  geom_point(aes(shape = season), size = 2) +
  scale_color_manual(values = c("ONI" = "#1f78b4", "PDO" = "#e31a1c")) +
  theme_minimal(base_size = 12) +
  geom_text(data = period_bounds,
            aes(x = (xmin + xmax)/2, y = 2.5, label = Period),
            inherit.aes = FALSE, size = 3.5) +
  labs(title = "ONI and PDO z-scores over time",
       x = "Year",
       y = "Z-score",
       color = "Index",
       linetype = "Season",
       shape = "Season") +
  scale_x_continuous(breaks = seq(min(df_long$year), max(df_long$year), by = 1)) +
  geom_vline(xintercept = period_vlines, linetype = "dashed", color = "black", size = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/plot_ONI_PDO_linesInTime.png", width = 17, height = 7, units="cm",dpi = 300)
# ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/plot_ONI_PDO_linesInTime.pdf", width = 17, height = 7, units="cm",dpi = 300)




# install.packages("FSA")   # if not installed
# install.packages("rcompanion")
# install.packages("dunn.test")
library(FSA)
library(dunn.test)
library(rcompanion)

# ONI
kruskal_ONI <- df_long %>%
  filter(Index == "ONI") %>%
  kruskal.test(Z_value ~ Period, data = .)
kruskal_ONI
# Kruskal-Wallis rank sum test
# data:  Z_value by Period
# Kruskal-Wallis chi-squared = 12.225, df = 3, p-value = 0.006651

# PDO
kruskal_PDO <- df_long %>%
  filter(Index == "PDO") %>%
  kruskal.test(Z_value ~ Period, data = .)
kruskal_PDO
# Kruskal-Wallis rank sum test
# data:  Z_value by Period
# Kruskal-Wallis chi-squared = 30.729, df = 3, p-value = 9.696e-07

# ONI
dunn_ONI <- dunn.test(x = df_long$Z_value[df_long$Index=="ONI"],
          g = df_long$Period[df_long$Index=="ONI"],
          method = "bonferroni")
cld_ONI <- cldList(P.adjusted ~ comparisons, data=dunn_ONI, threshold=0.05)

# PDO
dunn_PDO <- dunn.test(x = df_long$Z_value[df_long$Index=="PDO"],
          g = df_long$Period[df_long$Index=="PDO"],
          method = "bonferroni")
cld_PDO <- cldList(P.adjusted ~ comparisons, data=dunn_PDO, threshold=0.05)
cld_PDO

# ONI Plut PDO
dunn_<- dunn.test(x = df_long$Z_value,
                      g = df_long$Period,
                      method = "bonferroni")
cld_ONIPDO <- cldList(P.adjusted ~ comparisons, data=dunn_, threshold=0.05)



# Insert a space before the digit in Group
cld_PDO$Group <- sub("([a-zA-Z]+)(\\d)", "\\1 \\2", cld_PDO$Group)
cld_ONI$Group <- sub("([a-zA-Z]+)(\\d)", "\\1 \\2", cld_ONI$Group)

df_letters <- df_long %>%
  distinct(Period, Index) %>%
  left_join(cld_ONI %>% select(Group, Letter), by=c("Period"="Group")) %>%
  rename(ONI_letter = Letter) %>%
  left_join(cld_PDO %>% select(Group, Letter), by=c("Period"="Group")) %>%
  rename(PDO_letter = Letter)

# reorder periods
period_levels <- c("Suboptimal 1", "Optimal 1", "Suboptimal 2", "Optimal 2")
df_long$Period <- factor(df_long$Period, levels = period_levels)
df_letters$Period <- factor(df_letters$Period, levels = period_levels)
period_bounds$Period <- factor(period_bounds$Period, levels = period_levels)


df_long$period_kelp<- as.factor(df_long$Period)
levels(df_long$period_kelp)<- c("Suboptimal", "Optimal", "Suboptimal", "Optimal")

# PDO boxplot
a<- ggplot(df_long %>% filter(Index=="PDO"), aes(x=Period, y=Z_value)) +
  geom_boxplot(aes(fill=period_kelp), width=0.6) +
  geom_text(data=df_letters, aes(x=Period, y=max(df_long$Z_value)+0.1, label=PDO_letter),
            inherit.aes = FALSE, size=3.5) +
  scale_fill_manual(values=c("#fc8d59","#91bfdb")) +
  theme_minimal(base_size = 11) +
  labs(title="PDO", y="Z-score", x="Period") +
  theme(legend.position="none")


b<- ggplot(df_long %>% filter(Index=="ONI"), aes(x=Period, y=Z_value)) +
  geom_boxplot(aes(fill=period_kelp), width=0.6) +
  geom_text(data=df_letters, aes(x=Period, y=max(df_long$Z_value)+0.1, label=PDO_letter),
            inherit.aes = FALSE, size=3.5) +
  scale_fill_manual(values=c("#fc8d59","#91bfdb")) +
  theme_minimal(base_size = 11) +
  labs(title="ONI", y="Z-score", x="Period") +
  theme(legend.position="none")

cowplot::plot_grid(a, b, labels=c("A)", "B)"))


ggsave("/Volumes/Romina_PSF/PSF/SDM/SDM_results/boxplot_stats_ONI_PDO.png", width = 18, height = 7, units="cm",dpi = 300)



df_long %>%
  group_by(Index, Period) %>%
  summarize(
    mean_z = mean(Z_value),
    sd_z   = sd(Z_value),
    n      = n()
  )

# Index Period         mean_z  sd_z     n
# 1 ONI   Suboptimal 1  0.503 0.596    10
# 2 ONI   Optimal 1    -0.525 0.801    14
# 3 ONI   Suboptimal 2  0.657 1.02     12
# 4 ONI   Optimal 2    -0.142 1.15     11
# 5 PDO   Suboptimal 1  0.585 0.526    10
# 6 PDO   Optimal 1    -0.635 0.832    14
# 7 PDO   Suboptimal 2  0.809 0.589    12
# 8 PDO   Optimal 2    -1.48  0.589    11

