#TL DR; Observational data threshold model -- get threshold at each site

#in this file, I am creating a duplicate of the observed data with a column threshold_9
#corresponding to the 0.9 quantile of the climate date associate with the observed data's site.
#I am adding a column year associated with the year of the observed data
#I am also creating and saving a quantile regression model with the ALD for the 0.9 quantile
#where the regression model has the shape beta_0 + beta_1 climate_quantile

#takes 30 sec to run

gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(data.table)
library(raster) # package for netcdf manipulation

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

#fits a regression model to the 0.9 quantile 
quantile_model <-  maxtp ~ clim_thresh_value_9
quantile_model_fit <- evgam(quantile_model, obs_data, family = "ald", ald.args = list(tau = 0.9))
obs_data$threshold_9 = quantile_model_fit$location$fitted

#saves the estimates 0.9 quantile of obs data, estimates using regression based on 0.9 climate quantile
quantile_model_fit %>% saveRDS("output/threshold_model_9") #keeping the folder name output for now
  
write.csv(obs_data, "Data/processed/1971_2022_JJA_obs_data_bulk_model.csv", row.names = FALSE)