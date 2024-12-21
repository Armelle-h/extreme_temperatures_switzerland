gc() 
rm(list = ls())

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(data.table)
library(sf)
library(rnaturalearth)

obs_data = read.csv("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])
  
  clim_data$maxtp = as.integer(clim_data$maxtp)
  
  sing_clim_thresh_values = clim_data %>%
    filter(id %in% unique(obs_data$id))
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_thresh_values = do.call(rbind, clim_thresh_values_list)

date_id = read.csv("Data/id_date_correspondance.csv")

clim_thresh_values = clim_thresh_values %>%
  left_join(date_id, by = "date_id")%>%
  select(-date_id)%>%
  rename(maxtp_clim = maxtp)

obs_clim_data = obs_data %>%
  left_join(clim_thresh_values, by = c("id", "date"))


obs_clim_data = obs_clim_data %>%
  mutate(tp_diff = abs(maxtp-maxtp_clim))


legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>%
  rename(altitude = Altitude.m.)%>%
  select(stn, altitude)

obs_clim_data = obs_clim_data %>%
  left_join(legend_data, by = "stn")

T = obs_clim_data %>% filter(tp_diff>=5) %>% select(stn, tp_diff) %>%count(stn)
T = T  %>% left_join(legend_data%>%select(stn, longitude, latitude, Nom, Altitude.m.)  , by = "stn" )

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

points_sf <- st_as_sf(T, coords = c("longitude", "latitude"), crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = points_sf, aes(size = n), alpha=0.7) + 
  ggtitle("Stations where the climate model mis-estimated obs data \n by more that 8 degrees, size wrt number of mismeasurements") +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")



