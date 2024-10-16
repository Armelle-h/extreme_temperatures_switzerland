# Description --- This script fits all potential marginal gpd models to each 
# bootstrap and predicts scale and shape parameter from the GPD on the climate 
# model grid

gc()
rm(list=ls()) 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(gridExtra)

fit_uncorrected_models =T
fit_true_models = T

#loading the custom functions
source('gpd_models.R') #for now, might need to do a bit of relocating

#observed data that has exceeded the threshold
obs_data = vroom::vroom("data/processed/obs_data_gpd_model.csv") %>%
  left_join(vroom::vroom("Data/processed/thresh_exceedance_lambda_num_quantiles_30.csv")) #lambda associated with observed data, vary the quantile model

#doing it only for model 0, only sigma_c as covariate
covars = obs_data %>%
  dplyr::select(Station, year, scale_9, thresh_exceedance_9, threshold_9) %>%
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
}

#loading the fitted models

model_0_true = readr::read_csv("output/gpd_model_fits/model_0_true.csv") %>%
  rename(b0 = V1, b1 = V2, xi = V3)

