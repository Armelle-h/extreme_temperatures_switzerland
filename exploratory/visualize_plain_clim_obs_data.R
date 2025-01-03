library(tidyverse)
library(sf)
library(rnaturalearth)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_data_loc_id.csv")

legend = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")%>%
  select(stn, Nom, altitude, longitude, latitude)

obs_data = obs_data %>%
  left_join(legend, by="stn")

pts_sf <- st_as_sf(obs_data, coords = c("longitude", "latitude"), crs = 4326)

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = pts_sf) + 
  #  ggtitle("plain vs mountain") +
  labs(x = "longitude", y = "latitude") +
  theme_minimal()

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

pts_sf <- st_as_sf(I_plain, coords = c("longitude", "latitude"), crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = pts_sf) + 
  labs(x = "longitude", y = "latitude") +
  theme_minimal()

