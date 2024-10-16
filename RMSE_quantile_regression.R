#takes a long time to read RDS file
#can't change the format as spline functions are stored

#makes sense that the tau_to-temp function don't vary with time as in the regression model, 
#max_tp DOESN'T depend on temporal covariates !!

#if I'm putting all the obs_smoothed quantiles in the same folder, I could do a for loop !

gc()
rm(list = ls())
library(tidyverse)
library(lubridate)
library(purrr)
library(future)
library(future.apply)

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#num_cores <- availableCores() - 1
#plan(multisession, workers = num_cores)

#reading RDS files takes a long time, parallelising it.
#rds_files <- list.files(path = "output", pattern = "*.csv", full.names = TRUE)

#data_list <- future_lapply(rds_files, readRDS)

#plan(sequential)

#ran until here

#no_covariate_obs_smoothed_quantiles = data_list[[1]]
#obs_smoothed_quantiles = data_list[[2]]
#glob_anom_obs_smoothed_quantiles = data_list[[3]]

no_covariate_obs_smoothed_quantiles = readRDS("output/no_covariate_quant_models_num_quantiles_30.csv")
obs_smoothed_quantiles = readRDS("output/quant_models_num_quantiles_30.csv")
glob_anom_obs_smoothed_quantiles = readRDS("output/glob_anomaly_quant_models_num_quantiles_30.csv")

altitude_obs_smoothed_quantiles = readRDS("output/altitude_quant_models_num_quantiles_30.csv")

obs_data = read.csv("Data/Observed_data/1971_2023_JJA_obs_data_loc_id.csv")

num_quantiles = 50 
quantiles_to_estimate = seq(0.001,0.99,length.out = num_quantiles)

#empirical quantiles of obs_data
obs_data_quantile <- obs_data %>%
  select(-id) %>%
  mutate(year = year(date)) %>%
  select(-date) %>%
  group_by(stn, year) %>%
  summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))

extract_quantile <- function(df){
  df_pred_quantiles = df %>%
    mutate(quantiles = map(tau_to_temp, ~ sapply(quantiles_to_estimate, .x))) %>%
    select(year, stn, quantiles)
}

no_covariate_pred_obs_quantiles <- extract_quantile(no_covariate_obs_smoothed_quantiles) 

pred_obs_quantiles <- extract_quantile(obs_smoothed_quantiles)

glob_anom_pred_obs_quantiles <- extract_quantile(glob_anom_obs_smoothed_quantiles)

altitude_pred_obs_quantiles <- extract_quantile(altitude_obs_smoothed_quantiles)

#the code works because all the lists have the same length

# Function to calculate Root Mean Squared Error
calculate_rmse <- function(list1, list2) {
  if (length(list1) == length(list2)) {
    return(sqrt(mean((list1 - list2)^2, na.rm = TRUE)))  # Calculate RMSE
  } else {
    return(NA)  # Return NA if the lengths are different
  }
}

total_rmse <- function(df){
  RMSE_df = df%>%
    full_join(obs_data_quantile, by = c("stn", "year"), suffix = c("_pred", "_empirical"))
  
  # Calculate the RMSE for corresponding quantile lists
  results <- RMSE_df %>%
    mutate(rmse = map2_dbl(quantiles_pred, quantiles_empirical, calculate_rmse))
  
  # Sum the RMSE over all years and stations
  total_rmse_df <- results %>%
    summarise(total_rmse = sqrt(sum(rmse^2, na.rm = TRUE) / sum(!is.na(rmse)))) 
  
  total_rmse_value <- total_rmse_df$total_rmse
}

rmse_no_covariates = total_rmse(no_covariate_pred_obs_quantiles)

rmse_quantiles = total_rmse(pred_obs_quantiles)

rmse_glob_anom = total_rmse(glob_anom_pred_obs_quantiles)

rmse_altitude = total_rmse(altitude_pred_obs_quantiles)

print(rmse_no_covariates)
print(rmse_quantiles)
print(rmse_glob_anom)
print(rmse_altitude)