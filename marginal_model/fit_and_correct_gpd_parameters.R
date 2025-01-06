# Description --- This script fits all potential marginal gpd models to each 
# bootstrap and predicts scale and shape parameter from the GPD on the climate 
# model grid

gc()
rm(list=ls()) 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(gridExtra)

#loading the custom functions
source('marginal_model/gpd_models.R') 

num_quantiles = 30

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

#adding the altitude

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

#adding threshold 9 

threshold_9_df = vroom::vroom("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom(paste0("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv") )

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_9_df, by="id")

#doing it only for model 0, only sigma_c (=scale_9) as covariate
covars = obs_data %>%
  dplyr::select(stn, year, scale_9, glob_anom, thresh_exceedance_9, threshold_9) %>%
  unique()

# ----- Fit model
  
#extracting the excess data
extreme_dat_true = obs_data %>%
  mutate(excess = maxtp - threshold_9) %>%
  filter(excess > 0)
  
#fit and save model 0, 1, 2
  
#reminder: scale_9 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_9, c(0.8, -0.05, -0.1))  %>%
  matrix() %>% t() %>% as.data.frame() %>%
  write_csv("output/gpd_model_fits/model_0_true.csv")
  
fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_9,
          extreme_dat_true$glob_anom, c(0.8, -0.05, 0.001, -0.1))  %>%
  matrix() %>% t() %>% as.data.frame() %>%
  write_csv("output/gpd_model_fits/model_1_true.csv")
  
fit_mod_2(extreme_dat_true$excess, extreme_dat_true$scale_9,
          extreme_dat_true$glob_anom, extreme_dat_true$altitude, c(0.03, -0.03,  0.08,  1, -0.08, -0.2))  %>%
  matrix() %>% t() %>% as.data.frame() %>%
  write_csv("output/gpd_model_fits/model_2_true.csv")
  
fit_mod_3(extreme_dat_true$excess, extreme_dat_true$scale_9,
          extreme_dat_true$altitude, c(0.5, 0.01, 0.01, -0.1))  %>%
  matrix() %>% t() %>% as.data.frame() %>%
  write_csv("output/gpd_model_fits/model_3_true.csv")