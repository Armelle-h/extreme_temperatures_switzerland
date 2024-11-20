gc()
rm(list = ls())
library(tidyverse)
library(mgcv)
library(data.table)

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

raw_obs_data<- read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv", header=TRUE)
obs_data <- raw_obs_data %>%
  filter(maxtp != "-")#filtering out rows with missing values
obs_data$maxtp <- as.integer(obs_data$maxtp) #converting last column elements as integers

#adding a column year
obs_data = obs_data %>%
  mutate(year = lubridate::year(date))

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

obs_data = obs_data %>%
  left_join(clim_thresh_values) #associates 0.9 climate quantile to obs data

obs_data_ = obs_data %>%
  group_by(stn) %>%
  mutate(obs_quant_9 = quantile(maxtp, 0.9))


obs_data_$obs_threshold = obs_data_$clim_thresh_value_9/2 + obs_data_$obs_quant_9/2

write.csv(obs_data_, "Data/processed/obs_threshold_1971_2022_JJA_obs_data_bulk_model.csv", row.names = FALSE)