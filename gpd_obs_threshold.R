gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30

obs_data = read.csv("Data/processed/obs_threshold_1971_2022_JJA_obs_data_bulk_model.csv")

temporal_covariates_year = sort(unique(obs_data$year))

obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

lambda_thresh_ex = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    # Calculate exceedance probability for threshold_9 
    #(threshold_9 using a regression model is an estimate of the 0.9 obs quantile)
    
    #using the temp_to_tau function computed above to find the tau associated with threshold_9
    thresh_exceedance_9 = obs_smoothed_quantiles%>%
      filter(stn == .x$stn[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$obs_threshold[1], x)) #threshold defined using empirical obs and clim quantile
    
    
    tibble(stn = .x$stn[1],
           year = temporal_covariates_year,
           thresh_exceedance_9 = 1-thresh_exceedance_9)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

#replacing the negative proba by 0

lambda_thresh_ex = lambda_thresh_ex %>%
  mutate(thresh_exceedance_9 = pmax(thresh_exceedance_9, 0))


lambda_thresh_ex %>%
  write_csv(paste0("Data/processed/obs_threshold_glob_anomaly_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))


