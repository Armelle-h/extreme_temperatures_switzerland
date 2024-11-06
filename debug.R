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

  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_95) %>%
    filter(excess > 0)
  
  #fit and save model 0, 1, 2  #you fit on the OBSERVED data. So enough to have altitude info only on the observed points :)
  
  #reminder: scale_95 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
  
  L_0 = list(c(0.5, 0.1, -0.03, -0.1), c(-1, 0.5, 0.3, -0.1), c(-0.1, 0.01, -0.01, -0.1))
  
  for (init_par in L_0){
    p = fit_mod_3(extreme_dat_true$excess, extreme_dat_true$scale_95,
                  extreme_dat_true$altitude, init_par)
    print(p)
    print(ngll_3(p))
  }
  
  fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_95)  
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_95,
            extreme_dat_true$glob_anom)  
  
  fit_mod_2(extreme_dat_true$excess, extreme_dat_true$scale_95,
            extreme_dat_true$glob_anom, extreme_dat_true$altitude)  
  
  fit_mod_3(extreme_dat_true$excess, extreme_dat_true$scale_95,
            extreme_dat_true$altitude)  