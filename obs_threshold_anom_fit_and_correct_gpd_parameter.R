# Description --- This script fits all potential marginal gpd models to each 
# bootstrap and predicts scale and shape parameter from the GPD on the climate 
# model grid

gc()
rm(list=ls()) 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(gridExtra)

#loading the custom functions
source('gpd_models.R') #for now, might need to do a bit of relocating

num_quantiles = 30

fit_uncorrected_models = T
fit_true_models = T

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 
obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

#adding threshold 9 

threshold_9_df = vroom::vroom("Data/processed/obs_threshold_1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(obs_threshold, stn)%>%
  rename(threshold_9 = obs_threshold)%>%
  unique()

#lambda associated with observed data, vary the quantile model
#lambda_df = vroom::vroom(paste0("Data/processed/obs_threshold_glob_anomaly_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv") )

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  #left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_9_df, by="stn")


obs_data = obs_data %>% 
  group_by(stn)%>% 
  filter(maxtp < quantile(maxtp, 0.999))%>% 
  filter(maxtp > quantile(maxtp, 0.001))%>% 
  ungroup


#loading the covariates
covars = obs_data %>%
  #dplyr::select(stn, year, scale_9, glob_anom, thresh_exceedance_9, threshold_9, altitude) %>%
  dplyr::select(stn, year, scale_9, glob_anom, threshold_9, altitude) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)
  
  #fit and save model 0, 1, 2  #you fit on the OBSERVED data. So enough to have altitude info only on the observed points :)
  
  #reminder: scale_9 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
  
  
  fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_9, initial_pars = c(0.3, 0.7, -0.1) )  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/full_obs_threshold_model_0_true.csv")
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_9,
            extreme_dat_true$glob_anom, initial_pars = c(0.2, 0.5, 0.2, -0.1))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/full_obs_threshold_model_1_true.csv")
  
  
  #not done for full threshold
  fit_mod_2(extreme_dat_true$excess, extreme_dat_true$scale_9,
            extreme_dat_true$glob_anom, extreme_dat_true$altitude, initial_pars = c(0.1, 0.4, 0.01, 0.2, 0.01, -0.1))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/full_obs_threshold_model_2_true.csv")
  
  #not done for full threshold
  fit_mod_3(extreme_dat_true$excess, extreme_dat_true$scale_9,
            extreme_dat_true$altitude, initial_pars = c(0.2,  0.5,  0.02, -0.1))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/full_obs_threshold_model_3_true.csv")
  
}























list_ip = list(c(0.2, 0.5, 0.2, -0.1), c(1, 0.1, 0.5, -0.1), c(0.5, 0.5, 0,5, -0.1)  )

for (ip in list_ip){
  
  ip = c(0.5, 0.15, 0.4, -0.25)
  
  par = fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_9,
                  extreme_dat_true$glob_anom, ip)
  
  sol=ngll_1(par)
  
  print(par)
  print(sol)
}





