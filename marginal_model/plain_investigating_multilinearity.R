gc()
rm(list = ls())
library(tidyverse)
library(evgam)
library(car)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


num_quantiles = 30
quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
obs_data = readRDS(paste0("Data/processed/plain_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, altitude), by="stn")
#joining obs data with global anomaly

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

df_quant_glob_anom = data.frame(vif_value = numeric(), vif_glob_anom = numeric())

df_quant_log_alt = data.frame(vif_value = numeric(), vif_log_alt = numeric())

df_quant_glob_log_alt = data.frame(vif_value = numeric(), vif_glob_anom = numeric(), vif_log_alt = numeric())

index = 0

for(q in seq_along(quantiles_to_estimate_bulk)){
  
  obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
  
  zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
  print(paste0("Fitting  tau = ", zeta))
  
  obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
  
  # Fit a linear model with the same covariates
  lm_model_quant_glob_anom <- lm(maxtp ~ value + glob_anom, data = obs_data_for_quant_reg)
  lm_model_quant_log_altitude <- lm(maxtp ~ value + log(altitude), data = obs_data_for_quant_reg)
  lm_model_quant_glob_anom_log_altitude <- lm(maxtp ~ value + glob_anom + log(altitude), data = obs_data_for_quant_reg)
  
  
  # Calculate VIF
  vif_model_quant_glob_anom <- vif(lm_model_quant_glob_anom)
  vif_model_quant_log_altitude <- vif(lm_model_quant_log_altitude)
  vif_model_quant_glob_anom_log_altitude <- vif(lm_model_quant_glob_anom_log_altitude)
  
  
  df_quant_glob_anom <- rbind(df_quant_glob_anom, data.frame(value = vif_model_quant_glob_anom["value"], glob_anom = vif_model_quant_glob_anom["glob_anom"]))
  df_quant_log_alt <- rbind(df_quant_log_alt, data.frame(value = vif_model_quant_log_altitude["value"], log_alt = vif_model_quant_log_altitude["log(altitude)"]))
  df_quant_glob_log_alt <- rbind(df_quant_glob_log_alt, data.frame(value = vif_model_quant_glob_anom_log_altitude["value"], glob_anom = vif_model_quant_glob_anom_log_altitude["glob_anom"], log_alt = vif_model_quant_glob_anom_log_altitude["log(altitude)"]))
}
