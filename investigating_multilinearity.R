gc()
library(car)
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


num_quantiles = 30
quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)
#joining obs data with global anomaly

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
              names_to = "month", 
              values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

obs_data = obs_data %>% 
  mutate(month = month(date)) %>%
  left_join(glob_anomaly_reshaped, by = c("year", "month"))
index = 0

for(q in seq_along(quantiles_to_estimate_bulk)){
  
  index = index + 1
  
  if (index>4){break}

  obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
  
  zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
  print(paste0("Fitting  tau = ", zeta))
  
  obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist

  # Fit a linear model with the same covariates
  lm_model_quant_glob_anom <- lm(maxtp ~ value + glob_anom, data = obs_data_for_quant_reg)
  lm_model_quant_log_altitude <- lm(maxtp ~ value + log(altitude), data = obs_data_for_quant_reg)
  lm_model_quant_altitude_complex <- lm(maxtp ~ value + altitude + I(altitude^2) + log(altitude), data = obs_data_for_quant_reg)
  lm_model_quant_altitude_squared_log <- lm(maxtp ~ value + I(altitude^2)+log(altitude), data = obs_data_for_quant_reg)
  lm_model_quant_glob_anom_log_altitude <- lm(maxtp ~ value + glob_anom + log(altitude), data = obs_data_for_quant_reg)

  # Calculate VIF
  vif_model_quant_glob_anom <- vif(lm_model_quant_glob_anom)
  vif_model_quant_log_altitude <- vif(lm_model_quant_log_altitude)
  vif_model_quant_altitude_complex <- vif(lm_model_quant_altitude_complex)
  vif_model_quant_altitude_squared_log <- vif(lm_model_quant_altitude_squared_log)
  vif_model_quant_glob_anom_log_altitude <- vif(lm_model_quant_glob_anom_log_altitude)


  print(vif_model_quant_glob_anom)
  print(vif_model_quant_log_altitude)
  print(vif_model_quant_altitude_complex)
  print(vif_model_quant_altitude_squared_log)
  print(vif_model_quant_glob_anom_log_altitude)
}


