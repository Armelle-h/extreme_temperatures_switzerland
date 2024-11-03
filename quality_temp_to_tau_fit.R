#By default, the function temp to tau doesn't necessarily have images in (0,1).
#In this file, we investigate techniques to enforce this condition

gc()
rm(list = ls())
library(tidyverse)
library(job)
library(lubridate)

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 40

glob_anom_obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_", num_quantiles, ".csv"))

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

num_quantiles = 75 
quantiles_to_estimate = seq(0.001,0.99,length.out = num_quantiles)

#empirical quantiles of obs_data
obs_data_quantile <- obs_data %>%
  select(-id) %>%
  mutate(year = year(date)) %>%
  select(-date) %>%
  group_by(stn, year) %>%
  summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))

extract_tau <- function(df){
  merged_df <- df %>%
    inner_join(obs_data_quantile, by = c("stn", "year"))
  
  # Apply the function from col3_func to the array in col3_array for each row
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(est_tau = list(temp_to_tau(quantiles))) %>%
    ungroup()

  # Select only the columns you need
  final_df <- merged_df %>%
    select(stn, year, est_tau)

  return(final_df)
}

glob_anom_pred_obs_quantiles <- extract_tau(glob_anom_obs_smoothed_quantiles)

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
  df <- df %>%
    rowwise() %>%
    mutate(rmse_result = calculate_rmse(est_tau, quantiles_to_estimate)) %>%
    ungroup()
  
  mean_RMSE = mean(df$rmse_result) 
  return(mean_RMSE)
}

#function to compute the MAE 

calculate_mae <- function(list1, list2) {
  if (length(list1) == length(list2)) {
    return(mean(abs(list1 - list2), na.rm = TRUE))  # Calculate RMSE
  } else {
    return(NA)  # Return NA if the lengths are different
  }
}

total_mae <- function(df){
  df <- df %>%
    rowwise() %>%
    mutate(mae_result = calculate_mae(est_tau, quantiles_to_estimate)) %>%
    ungroup()
  
  mean_mae = mean(df$mae_result) 
  return(mean_mae)
}


rmse_glob_anom = total_rmse(glob_anom_pred_obs_quantiles)

mae_glob_anom = total_mae(glob_anom_pred_obs_quantiles)

print(rmse_glob_anom)

print(mae_glob_anom)
