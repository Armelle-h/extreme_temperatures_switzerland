gc()
rm(list=ls()) 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(data.table)
library(evgam)

num_quantiles = 30

obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

obs_data = obs_data%>%
  select(-c("threshold_9", "clim_thresh_value_9"))

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])
  
  clim_data$maxtp = as.integer(clim_data$maxtp)
  
  sing_clim_thresh_values = clim_data %>%
    group_by(id) %>%
    summarise(clim_thresh_value_95 = quantile(maxtp, 0.95, na.rm = TRUE)) #ensures missing values are ignored
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_thresh_values = do.call(rbind, clim_thresh_values_list)

obs_data = obs_data %>%
  left_join(clim_thresh_values, by="id") #associates 0.95 climate quantile to obs data

#fits a regression model to the 0.95 quantile 
quantile_model <-  maxtp ~ clim_thresh_value_95
quantile_model_fit <- evgam(quantile_model, obs_data, family = "ald", ald.args = list(tau = 0.95))
obs_data$threshold_95 = quantile_model_fit$location$fitted

#saves the estimates 0.95 quantile of obs data, estimates using regression based on 0.95 climate quantile
quantile_model_fit %>% saveRDS("output/threshold_model_95") #keeping the folder name output for now

obs_data %>%
  saveRDS(paste0("Data/processed/obs_data_95_for_bulk_model_num_quantiles_",num_quantiles,".csv")) #adding the 95 to avoid mistakes
#messing everything up

#obs_data = readRDS(paste0("Data/processed/obs_data_95_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

# ---- get covariates for prediction
temporal_covariates = obs_data %>%  #issue, be wary of, glob anom is defined for June July and August so each year is associated with 3 different values
  dplyr::select(year) %>% 
  unique() %>%
  arrange(year)

#computing the threshold function 

obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

# Calculate exceedance probability (lambda) for a threshold (threshold_95)
lambda_thresh_ex = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    print(.x$stn[1])
    
    # Calculate exceedance probability for threshold_95 
    #(threshold_95 using a regression model is an estimate of the 0.95 obs quantile)
    
    #using the temp_to_tau function computed above to find the tau associated with threshold_95
    thresh_exceedance_95 = obs_smoothed_quantiles%>%
      filter(stn == .x$stn[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_95[1], x))
    
    
    tibble(stn = .x$stn[1],
           year = temporal_covariates$year,
           thresh_exceedance_95 = 1-thresh_exceedance_95)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


lambda_thresh_ex %>%
  write_csv(paste0("Data/processed/glob_anomaly_95_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))

#BEFORE RUNNING BELOW; NEED TO HAVE RUN CLIM_DATA_GPD_MODEL_THRESHOLD_95 !!!!!

#computing the gpd parameters 

source('gpd_models.R') #for now, might need to do a bit of relocating


fit_uncorrected_models = T
fit_true_models = T

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_95_gpd_model.csv") #computed in clim_data_gpd_model_threshold_95

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 
obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

#adding threshold 95 

threshold_95_df = readRDS(paste0("Data/processed/obs_data_95_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  dplyr::select(threshold_95, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom("Data/processed/glob_anomaly_95_thresh_exceedance_lambda_num_quantiles_30.csv") 

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_95_df, by="id")

#loading the covariates
covars = obs_data %>%
  dplyr::select(stn, year, scale_95, glob_anom, thresh_exceedance_95, threshold_95, altitude) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_95) %>%
    filter(excess > 0)
  
  #fit and save model 0, 1, 2  #you fit on the OBSERVED data. So enough to have altitude info only on the observed points :)
  
  #reminder: scale_95 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
  fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_95, c(0.4, 0.1, -0.1))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_0_true_95.csv")
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_95,
            extreme_dat_true$glob_anom, c(0.1, 0.1, 0.5, -0.1))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_1_true_95.csv")
  
  fit_mod_2(extreme_dat_true$excess, extreme_dat_true$scale_95,
            extreme_dat_true$glob_anom, extreme_dat_true$altitude, c(-0.4, 0.1, 0.1, 1.2, -0.1, -0.07))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_2_true_95.csv")
  
  fit_mod_3(extreme_dat_true$excess, extreme_dat_true$scale_95,
            extreme_dat_true$altitude, c(0.7, 0.01, -0.03, -0.1))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_3_true_95.csv")
  
}

#TO FINISH; RUN THE FILE CV_MODEL_01_THRESHOLD_95 TO GET THE PERFORMANCE METRICS