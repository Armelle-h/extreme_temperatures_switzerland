# Description --- This script fits all potential marginal gpd models to each 
# bootstrap and predicts scale and shape parameter from the GPD on the climate 
# model grid

gc()
rm(list=ls()) 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(gridExtra)

#loading the custom functions
source('gpd_models.R') #for now, might need to do a bit of relocating


fit_uncorrected_models = T
fit_true_models = T

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

#adding threshold 9 

threshold_9_df = vroom::vroom("Data/Observed_data/1971_2023_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom("Data/processed/thresh_exceedance_lambda_num_quantiles_30.csv") 

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
               names_to = "month", 
               values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

obs_data = obs_data %>% 
  mutate(year = year(date), month = month(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = c("year", "month")) %>%
  left_join(threshold_9_df, by="id")

#doing it only for model 0, only sigma_c (=scale_9) as covariate
covars = obs_data %>%
  dplyr::select(stn, year, scale_9, glob_anom, thresh_exceedance_9, threshold_9) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)
  
  #fit and save model 0, 1, 2
  
  #reminder: scale_9 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
  fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_9)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_0_true.csv")
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_9,
            extreme_dat_true$glob_anom)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_1_true.csv")
}

#loading the fitted models -- for later

#model_0_true = readr::read_csv("output/gpd_model_fits/model_0_true.csv") %>%
#  rename(b0 = V1, b1 = V2, xi = V3)

#model_1_true = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
#  rename(b0 = V1, b1 = V2, b2 = V3, xi = V4)
