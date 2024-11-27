gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30

obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_", num_quantiles, ".csv"))

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")%>%
  mutate(year = lubridate::year(date))


obs_quantile = obs_data %>%
  group_by(stn, year) %>%
  mutate(obs_q_25 = quantile(maxtp, 0.25))%>%
  mutate(obs_q_50 = quantile(maxtp, 0.5))%>%
  mutate(obs_q_75 = quantile(maxtp, 0.75))%>%
  select(stn, year, obs_q_25, obs_q_50, obs_q_75)



est_quantile <- obs_smoothed_quantiles %>%
  rowwise() %>%
  mutate(
    est_q_25 = tau_to_temp(0.25),  
    est_q_50 = tau_to_temp(0.5),  
    est_q_75 = tau_to_temp(0.75)   
  ) %>%
  ungroup() %>%
  select(stn, year, est_q_25, est_q_50, est_q_75)

quantile = obs_quantile %>%
  left_join(est_quantile, by= c("stn", "year"))%>%
  mutate(q_25 = abs(obs_q_25-est_q_25))%>%
  mutate(q_50 = abs(obs_q_50-est_q_50))%>%
  mutate(q_75 = abs(obs_q_75-est_q_75))


col = quantile$q_25

print(min(col))
print(max(col))
print(mean(col))
print(quantile(col, 0.25))
print(quantile(col, 0.5))
print(quantile(col, 0.75))
print(quantile(col, 0.90))


