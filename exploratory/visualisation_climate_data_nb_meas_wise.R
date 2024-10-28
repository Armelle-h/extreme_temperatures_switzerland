setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

gc()
rm(list = ls())
library(rnaturalearth) #library for map of switzerland
library(tidyverse)
library(sf)

data<- read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv", header=TRUE)

legend<- read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv", header=TRUE)

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

#count nb of measurement associated with each station
nb_meas_df <- data %>%
  group_by(stn) %>%
  summarise(count = n())

merged_data <- legend %>%
  left_join(nb_meas_df, by = "stn")

points_sf <- st_as_sf(merged_data, coords = c("longitude", "latitude"), crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightgrey", color = "black") +  # Plot Switzerland
  geom_sf(data = points_sf, aes(size = count, color = count), alpha=1.0) +  # Plot points sized by measurements
  scale_size_continuous(range = c(1, 4)) +  # Adjust size range
  scale_color_viridis_c(name = "Number of Measurements", option = "viridis") +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")