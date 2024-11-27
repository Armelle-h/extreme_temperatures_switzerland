
library(tidyverse)
library(sf)
library(rnaturalearth)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>%
  rename(altitude = Altitude.m.)%>%
  select(stn, Nom, altitude, longitude, latitude)

obs_data = obs_data %>%
  left_join(legend, by="stn")

obs_data_plain = obs_data %>%
  filter(altitude<1200)%>%
  filter(!(longitude>=6.8 & latitude<=46.4))%>%
  filter(!(longitude>=9.5 & latitude<=47.1))%>%
  filter(!(longitude>=8.5 & latitude<=46.7))%>%
  filter(!(stn %in% c("GST", "KALTB", "GWA", "ENG", "GTT", "INT", "BEA", "MMBIS", "ELM", "VAE", "ENNESF", "ILZ", "MMLEN", "FRU", "MMGTT", "MER", "MMLIN", "TFD", "ALT", "INNESF")))

#combining some stations that don't have the same name but are at exactly the same place, and record temperatures at non intersecting times


obs_data_plain <- obs_data_plain %>%
  mutate(stn = if_else(stn == "AGAAR", "AAR", stn))%>%
  mutate(stn = if_else(stn == "BZN", "BEZ", stn))%>%
  mutate(stn = if_else(stn == "WSLLAE", "NABLAE", stn))%>%
  mutate(stn = if_else(stn == "TGFRA", "FRF", stn))

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

# Step 3: Convert the dataframe to an sf object, assuming coordinates in WGS84 (EPSG:4326)
pts_sf <- st_as_sf(obs_data_plain, coords = c("longitude", "latitude"), crs = 4326)

# Step 4: Ensure the CRS for both datasets match (WGS84)
switzerland <- st_transform(switzerland, crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = pts_sf) + 
  ggtitle("plain vs mountain") +
  theme_minimal()

legend_filtered <- legend %>%
  filter(stn %in% obs_data_plain$stn)

write.csv(legend_filtered, "Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv", row.names = FALSE)

write.csv(obs_data_plain%>%select(stn, date, maxtp, id), "Data/Observed_data/plain_1971_2022_JJA_obs_data_loc_id.csv", row.names=FALSE)

