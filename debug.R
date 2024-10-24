#testing if bulk_models_function_cv works well

gc()
rm(list = ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('bulk_models_function_cv.R')

library(tidyverse)
library(spatialsample)
library(scoringRules)

# ========  ========  Global parameters ========

# ========  ========  Read in data  ========  ========
# Observational data with covariates

#new code
num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
               names_to = "month", 
               values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

obs_data = obs_data %>% 
  mutate(year = year(date), month = month(date), week=week(date)) %>% #week is used for the temporal cross validation
  left_join(glob_anomaly_reshaped, by = c("year", "month"))

# onlykeep full years for cv
obs_data <- obs_data %>%
  group_by(year, stn) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n == 92) %>% 
  left_join(obs_data) %>%
  ungroup()

#temporal_covariates = obs_data %>%
# dplyr::select(year, altitude, glob_anom) %>% 
# unique() %>%
# arrange(year)

obs_data_raw = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

num_quantiles_to_estimate = 50
quantiles_to_estimate = seq(0.001,0.99,length.out = num_quantiles_to_estimate)

#empirical quantiles of obs_data
#obs_data_quantile <- obs_data_raw %>%
# select(-id) %>%
# mutate(year = year(date)) %>%
# select(-date) %>%
# group_by(stn, year) %>%
# summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))

num_quantiles_to_fit = 30
fitting_quantiles = seq(0.001,0.99,length.out = num_quantiles_to_fit)

model_name = "quant_log_alt"

rmse_no_covar = estimate_RMSE(fitting_quantiles, obs_data, model_name, quantiles_to_estimate, obs_data_raw)
