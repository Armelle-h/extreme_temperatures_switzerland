# This script fits and saves quantile regression model and lambda estimates

#until spline models for climate data, takes 12 minutes to run

gc()
rm(list = ls())
library(tidyverse)
library(evgam)
library(car)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

# ---- get covariates for prediction
temporal_covariates = obs_data %>%
  dplyr::select(year, altitude) %>% 
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # bool, re-estimate clim quantiles?

if(fit_clim_quants){
  
  file.remove(paste0("Data/processed/complex_altitude_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
  # Loop over each quantile and fit the quantile regression model
  for(q in seq_along(quantiles_to_estimate_bulk)){
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    print(paste0("Fitting  tau = ", zeta))
    
    #in prep data bulk model we had associated to each data a quantile climate. 
    #regression only in terms of the climate quantile covariate
    
    # Set the 'value' column for current quantile
    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    
    # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
    quantile_model_fit <- evgam(maxtp ~ value + altitude + I(altitude^2) + log(altitude), obs_data_for_quant_reg,
                                family = "ald", ald.args = list(tau = zeta))
    
    # Save the fitted parameter estimates for each quantile
    tibble(tau = zeta,
           beta_0 = quantile_model_fit$location$coefficients[1],
           beta_1 = quantile_model_fit$location$coefficients[2],
           beta_2 = quantile_model_fit$location$coefficients[3],
           beta_3 = quantile_model_fit$location$coefficients[4],
           beta_4 = quantile_model_fit$location$coefficients[5]) %>%
      write_csv(paste0("Data/processed/complex_altitude_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("Data/processed/complex_altitude_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3', 'beta_4'))

cat("beta_0", mean(quant_reg_model_pars$beta_0, na.rm=TRUE),  
    "beta_1", mean(quant_reg_model_pars$beta_1, na.rm=TRUE),
    "beta_2", mean(quant_reg_model_pars$beta_2, na.rm=TRUE),
    "beta_3", mean(quant_reg_model_pars$beta_3, na.rm=TRUE),
    "beta_4", mean(quant_reg_model_pars$beta_4, na.rm=TRUE) )     