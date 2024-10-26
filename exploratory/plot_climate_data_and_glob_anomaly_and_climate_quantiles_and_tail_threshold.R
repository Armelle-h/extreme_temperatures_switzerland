#plots the climate data generated for the 1st of august 2020 on a grid of 1km over Switzerland
#plots the average over months June, July and August of the global anomaly for each year
#plots the climate quantiles for 4 different quantiles computed in file prep_data_for_bulk_models


gc() 
rm(list = ls())

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

#plotting first of August 2020

clim_data = read.csv("Data/Climate_data/2018_2022_JJA_climate_data.csv")

lon_lat_corres = read.csv("Data/id_lon_lat_correspondance.csv")

date_id_corres =  read.csv("Data/id_date_correspondance.csv")

id_2020 = date_id_corres %>%
  filter(date == "2020-08-01") %>%
  pull(date_id)

clim_data_lon_lat= clim_data %>%
  left_join(lon_lat_corres, by="id")

#clim_data_lon_lat$date <- as.Date(clim_data_lon_lat$date)

clim_data_lon_lat = clim_data_lon_lat %>%
  filter(date_id == id_2020)

ggplot(clim_data_lon_lat, aes(x = longitude, y = latitude, color = maxtp)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  theme_minimal()

# plotting global anomaly

glob_anom = read.csv("Data/global_tp_anomaly_JJA.csv")

ggplot(glob_anom, aes(x=year, y=JJA))+
  geom_line()+
  labs(x = "year", y = "global anomaly") +
  theme_minimal()

#plotting climate quantiles for the 1st of August 2020 (doesn't change anything since they are spatially and not temporally dependent)
num_quantiles = 30

clim_quantile = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

clim_quantile = clim_quantile %>%
  mutate(
    quantile_1 = map_dbl(value, ~ .x[1] %||% NA),
    quantile_10 = map_dbl(value, ~ .x[10] %||% NA),
    quantile_20 = map_dbl(value, ~ .x[20] %||% NA),
    quantile_30 = map_dbl(value, ~ .x[length(.x)] %||% NA)
  )

clim_quantile_loc = clim_quantile%>%
  left_join(lon_lat_corres, by="id")

#1st quantile
ggplot(clim_quantile_loc, aes(x = longitude, y = latitude, color = quantile_1)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  ggtitle("0.1th quantile") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#10th quantile
ggplot(clim_quantile_loc, aes(x = longitude, y = latitude, color = quantile_10)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  ggtitle("30.8th quantile") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#20th quantile
ggplot(clim_quantile_loc, aes(x = longitude, y = latitude, color = quantile_20)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  ggtitle("64.9th quantile") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#30th quantile 
ggplot(clim_quantile_loc, aes(x = longitude, y = latitude, color = quantile_30)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  ggtitle("99th quantile") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#plotting the tail threshold u_o(s)

#-----------probably to delete
obs_threshold = read.csv("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")

obs_threshold = obs_threshold %>%
  select(c("id", "threshold_9"))%>%
  unique()

obs_threshold_loc = obs_threshold %>%
  left_join(lon_lat_corres, by="id")

ggplot(obs_threshold_loc, aes(x = longitude, y = latitude, color = threshold_9)) +
  geom_point(size = 3, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "°C") +
  theme_minimal()
#end of probably to delete ---------------


threshold_model_9 = readRDS("output/threshold_model_9")

print(threshold_model_9$location$coefficients)

#compute the 0.9 quantile of climate data for all locations 

library(data.table)

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])
  
  clim_data$maxtp = as.integer(clim_data$maxtp)
  
  sing_clim_thresh_values = clim_data %>%
    group_by(id) %>%
    summarise(clim_thresh_value_9 = quantile(maxtp, 0.9, na.rm = TRUE)) #ensures missing values are ignored
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_thresh_values = do.call(rbind, clim_thresh_values_list)

beta_0 = threshold_model_9$location$coefficients[1]
beta_1 =  threshold_model_9$location$coefficients[2] 

clim_thresh_values$threshold_obs = beta_0 + beta_1*clim_thresh_values$clim_thresh_value_9

clim_thresh_values_loc = clim_thresh_values %>%
  left_join(lon_lat_corres, by="id")

#using the threshold regression, estimate the observed threshold

ggplot(clim_thresh_values_loc, aes(x = longitude, y = latitude, color = threshold_obs)) +
  geom_point(size = 1, alpha = 0.6) +  # Adjust size and transparency as needed
  scale_color_viridis_c() +  # Optional: For color scale
  labs(x = "Longitude", y = "Latitude", color = "u_o") +
  theme_minimal()


