#plots the climate data generated for the 1st of august 2020 on a grid of 1km over Switzerland
#plots the climate quantiles for the 90th quantile computed in file prep_data_for_bulk_models


gc() 
rm(list = ls())

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

#plotting first of August 2020

clim_data = read.csv("Data/Climate_data/2018_2022_JJA_climate_data.csv")

legend_data = read.csv("Data/id_lon_lat_correspondance.csv")

date_id_corres =  read.csv("Data/id_date_correspondance.csv")

id_2020 = date_id_corres %>%
  filter(date == "2020-08-01") %>%
  pull(date_id)

clim_data_lon_lat= clim_data %>%
  left_join(legend_data, by="id")

#clim_data_lon_lat$date <- as.Date(clim_data_lon_lat$date)

clim_data_lon_lat = clim_data_lon_lat %>%
  filter(date_id == id_2020)

ggplot(clim_data_lon_lat, aes(x = longitude, y = latitude, color = maxtp)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  theme_minimal()

glob_anom = read.csv("Data/global_tp_anomaly_JJA.csv")

ggplot(glob_anom, aes(x=year, y=JJA))+
  geom_line()+
  labs(x = "year", y = "global anomaly") +
  theme_minimal()