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
    summarise(clim_thresh_value_8 = quantile(maxtp, 0.8, na.rm = TRUE)) #ensures missing values are ignored
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_thresh_values = do.call(rbind, clim_thresh_values_list)

obs_data = obs_data %>%
  left_join(clim_thresh_values, by="id") #associates 0.8 climate quantile to obs data

#fits a regression model to the 0.8 quantile 
quantile_model <-  maxtp ~ clim_thresh_value_8
quantile_model_fit <- evgam(quantile_model, obs_data, family = "ald", ald.args = list(tau = 0.8))
obs_data$threshold_8 = quantile_model_fit$location$fitted

#saves the estimates 0.8 quantile of obs data, estimates using regression based on 0.8 climate quantile
quantile_model_fit %>% saveRDS("output/threshold_model_8") #keeping the folder name output for now

obs_data %>%
  saveRDS(paste0("Data/processed/obs_data_8_for_bulk_model_num_quantiles_",num_quantiles,".csv")) #adding the 8 to avoid mistakes
#messing everything up

obs_data = readRDS(paste0("Data/processed/obs_data_8_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

# ---- get covariates for prediction
temporal_covariates = obs_data %>%  #issue, be wary of, glob anom is defined for June July and August so each year is associated with 3 different values
  dplyr::select(year) %>% 
  unique() %>%
  arrange(year)

#computing the threshold function 

obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

# Calculate exceedance probability (lambda) for a threshold (threshold_8)
lambda_thresh_ex = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    print(.x$stn[1])
    
    # Calculate exceedance probability for threshold_8 
    #(threshold_8 using a regression model is an estimate of the 0.8 obs quantile)
    
    #using the temp_to_tau function computed above to find the tau associated with threshold_8
    thresh_exceedance_8 = obs_smoothed_quantiles%>%
      filter(stn == .x$stn[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_8[1], x))
    
    
    tibble(stn = .x$stn[1],
           year = temporal_covariates$year,
           thresh_exceedance_8 = 1-thresh_exceedance_8)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


lambda_thresh_ex %>%
  write_csv(paste0("Data/processed/glob_anomaly_8_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))

#BEFORE RUNNING BELOW; NEED TO HAVE RUN CLIM_DATA_GPD_MODEL_THRESHOLD_8 !!!!!

#computing the gpd parameters 

source('gpd_models.R') #for now, might need to do a bit of relocating


fit_uncorrected_models = T
fit_true_models = T

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_8_gpd_model.csv") #computed in clim_data_gpd_model_threshold_8

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 
obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

#adding threshold 8 

threshold_8_df = readRDS(paste0("Data/processed/obs_data_8_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  dplyr::select(threshold_8, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom("Data/processed/glob_anomaly_8_thresh_exceedance_lambda_num_quantiles_30.csv") 

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_8_df, by="id")

#loading the covariates
covars = obs_data %>%
  dplyr::select(stn, year, scale_8, glob_anom, thresh_exceedance_8, threshold_8, altitude) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_8) %>%
    filter(excess > 0)
  
  #fit and save model 0, 1, 2  #you fit on the OBSERVED data. So enough to have altitude info only on the observed points :)
  
  #reminder: scale_8 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
  fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_8)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/change_init_param_model_0_true_8.csv")
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_8,
            extreme_dat_true$glob_anom)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_1_true_8.csv")
  
  fit_mod_2(extreme_dat_true$excess, extreme_dat_true$scale_8,
            extreme_dat_true$glob_anom, extreme_dat_true$altitude)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_2_true_8.csv")
  
  fit_mod_3(extreme_dat_true$excess, extreme_dat_true$scale_8,
            extreme_dat_true$altitude)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_3_true_8.csv")
  
}

#TO FINISH; RUN THE FILE CV_MODEL_01_THRESHOLD_8 TO GET THE PERFORMANCE METRICS
