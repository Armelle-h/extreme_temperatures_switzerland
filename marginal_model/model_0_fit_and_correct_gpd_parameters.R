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

num_quantiles = 40

fit_uncorrected_models = T
fit_true_models = T

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

#adding threshold 9 

threshold_9_df = vroom::vroom("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom(paste0("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv") )

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
               names_to = "month", 
               values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

obs_data = obs_data %>% 
  mutate(year = year(date), month = month(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = c("year", "month")) %>%
  left_join(threshold_9_df, by="id")

#doing it only for model 0, only sigma_c (=scale_9) as covariate
covars = obs_data %>%
  dplyr::select(stn, year, scale_9, glob_anom, thresh_exceedance_9, threshold_9) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)
  
  #fit and save model 0, 1, 2
  
  #reminder: scale_9 is the scaling parameter of the climate model (computed in clim_data_gpd_model)
  fit_mod_0(extreme_dat_true$excess, extreme_dat_true$scale_9)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_0_true.csv")
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_9,
            extreme_dat_true$glob_anom)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_1_true.csv")
}

#loading the fitted models

model_0_true = readr::read_csv("output/gpd_model_fits/model_0_true.csv") %>%
  rename(b0 = V1, b1 = V2, xi = V3)

model_1_true = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
  rename(b0 = V1, b1 = V2, b2 = V3, xi = V4)

model_2_true = readr::read_csv("output/gpd_model_fits/model_2_true.csv") %>%
  rename(b0 = V1, b1 = V2, b2 = V3, b3 = V4, b4 = V5, xi = V6)

#fitting and saving gpd models to bootstrap data

bts_files = list.files("Data/processed/bootstrap_data/bts_under_gpd_models/")

if(fit_uncorrected_models){
  
  # ----- fit all models
  for(file_name in bts_files){
    
    print(paste0("fitting to bootstrap ", file_name))
    
    #load bootstrap data and merge with covariates
    dat = read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/", file_name))%>%
      setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3")) %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(covars %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station'))
    
    # --- remove observations poorly interpolated by quantile model (i.e. <0)
    dat = dat %>% filter(maxtp_0 > 0,
                         maxtp_1 > 0,
                         maxtp_2 > 0)
    
    # ----- fit Model 0
    obs_data_to_pred = dat %>%
      mutate(excess = maxtp_0 - threshold_9) %>%
      filter(excess > 0)
    
    this_fit_mod_0 = fit_mod_0(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
                               initial_pars = as.numeric(unlist(model_0_true)))
    
    c(file_name, this_fit_mod_0) %>% matrix() %>% t() %>% as.data.frame() %>%
      write_csv("output/gpd_model_fits/model_0_uncorrected.csv", append = T)
    
    
    # ----- Model 1
    obs_data_to_pred = dat %>%
      mutate(excess = maxtp_1 - threshold_9) %>%
      filter(excess > 0)
    this_fit_mod_1 = fit_mod_1(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
                               obs_data_to_pred$loess_temp_anom,
                               initial_pars = as.numeric(unlist(model_1_true)))
    c(file_name, this_fit_mod_1)  %>% matrix() %>% t() %>% as.data.frame() %>%
      write_csv("output/gpd_model_fits/model_1_uncorrected.csv", append = T)
    
    
    # ----- Model 2
    obs_data_to_pred = dat %>%
      mutate(excess = maxtp_2 - threshold_9) %>%
      filter(excess > 0)
    this_fit_mod_2 = fit_mod_2(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
                               obs_data_to_pred$loess_temp_anom,
                               obs_data_to_pred$dist_sea,
                               initial_pars = as.numeric(unlist(model_2_true)))
    c(file_name, this_fit_mod_2)  %>% matrix() %>% t() %>% as.data.frame() %>%
      write_csv("output/gpd_model_fits/model_2_uncorrected.csv", append = T)
  }
}


read_csv("output/gpd_model_fits/model_0_uncorrected.csv")



# ----  calculate correction terms
model_0_xi = model_0_true %>% pull(xi)
model_1_xi = model_1_true %>% pull(xi)
model_2_xi = model_2_true %>% pull(xi)

model_0_correction = model_0_xi - (read_csv("output/gpd_model_fits/model_0_uncorrected.csv",
                                            col_names = c('bts', 'b0', 'b1', 'xi')) %>%
                                     pull(xi) %>% mean)

model_1_correction = model_1_xi - (read_csv("output/gpd_model_fits/model_1_uncorrected.csv",
                                            col_names = c('bts', 'b0', 'b1',  'b2', 'xi')) %>%
                                     pull(xi) %>%
                                     mean)

model_2_correction = model_2_xi - (read_csv("output/gpd_model_fits/model_2_uncorrected.csv",
                                            col_names = c('bts', 'b0', 'b1', 'b2','b3','b4', 'xi')) %>%
                                     pull(xi) %>% mean)



#fit GPD models with corrected parameters
for(file_name in bts_files){
  
  print(paste0("fitting to bootstrap ", file_name))
  dat = readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/", file_name))%>%
    mutate(year = lubridate::year(date)) %>%
    left_join(covars %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station'))
  
  # --- remove observations poorly interpolated by quantile model (i.e. <0)
  dat = dat %>% filter(maxtp_0 > 0,
                       maxtp_1 > 0,
                       maxtp_2 > 0)
  
  #Fit corrected Model 0
  obs_data_to_pred = dat %>%
    mutate(excess = maxtp_0 - threshold_9) %>%
    filter(excess > 0)
  
  uncorrected_0_xi = read_csv("output/gpd_model_fits/model_0_uncorrected.csv",
                              col_names = c('bts', 'b0', 'b1', 'xi'))  %>%
    filter(bts == file_name) %>%
    pull(xi)
  
  this_fit_mod_0_corrected = fit_mod_0_fix_shape(obs_data_to_pred$excess,
                                                 obs_data_to_pred$scale_9,
                                                 this_shape_est = (uncorrected_0_xi + model_0_correction),
                                                 initial_pars = as.numeric(unlist(model_0_true))[-ncol(model_0_true)] + rnorm(1, 0, 0.01))
  c(file_name, this_fit_mod_0_corrected, (uncorrected_0_xi + model_0_correction)) %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_0_pars_corrected.csv", append = T)
  
  
  # # ----- Model 1
  obs_data_to_pred = dat %>%
    mutate(excess = maxtp_1 - threshold_9) %>%
    filter(excess > 0)
  
  uncorrected_1_xi = read_csv("output/gpd_model_fits/model_1_uncorrected.csv",
                              col_names = c('bts', 'b0', 'b1', 'b2', 'xi'))  %>%
    filter(bts == file_name) %>%
    pull(xi)
  
  this_fit_mod_1_corrected = fit_mod_1_fix_shape(obs_data_to_pred$excess,
                                                 obs_data_to_pred$scale_9,
                                                 obs_data_to_pred$loess_temp_anom,
                                                 this_shape_est = (uncorrected_1_xi + model_1_correction),
                                                 initial_pars =as.numeric(unlist(model_1_true))[-ncol(model_1_true)])
  
  c(file_name, this_fit_mod_1_corrected, (uncorrected_1_xi +  model_1_correction)) %>%
    matrix() %>% t() %>% as.data.frame() %>% write_csv("output/gpd_model_fits/model_1_pars_corrected.csv", append = T)
  
  # ----- Model 2
  obs_data_to_pred = dat %>%
    mutate(excess = maxtp_2 - threshold_9) %>%
    filter(excess > 0)
  
  uncorrected_2_xi = read_csv("output/gpd_model_fits/model_2_uncorrected.csv",
                              col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi'))  %>%
    filter(bts == file_name) %>% pull(xi)
  
  this_fit_mod_2_corrected = fit_mod_2_fix_shape(obs_data_to_pred$excess,
                                                 obs_data_to_pred$scale_9,
                                                 obs_data_to_pred$loess_temp_anom,
                                                 obs_data_to_pred$dist_sea,
                                                 this_shape_est = (uncorrected_2_xi + model_2_correction),
                                                 initial_pars =as.numeric(unlist(model_2_true))[-ncol(model_2_true)])
  
  c(file_name, this_fit_mod_2_corrected, (uncorrected_2_xi +  model_2_correction)) %>%
    matrix() %>% t() %>% as.data.frame() %>% write_csv("output/gpd_model_fits/model_2_pars_corrected.csv", append = T)
}


# # # # ----- Visualise results
bias_correction_mod_1 = gridExtra::grid.arrange(read_csv("output/gpd_model_fits/model_2_uncorrected.csv",
                                                         col_names = c('bts', 'β0', 'β1', 'β2', 'β3', 'β4', 'ξ')) %>%
                                                  pivot_longer(-bts) %>%
                                                  ggplot()+
                                                  geom_density(aes(value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_2_true.csv") %>%
                                                               rename(`β0` = V1, `β1` = V2, `β2` = V3, `β3` = V4, `β4` = V5, `ξ` = V6) %>%
                                                               pivot_longer(everything()),
                                                             aes(xintercept = value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  labs(y= "Density", x = "Parameter value")+
                                                  theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1),
                                                                          axis.title.x=element_blank()),
                                                read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                                                         col_names = c('bts', 'β0', 'β1', 'β2', 'β3', 'β4', 'ξ')) %>%
                                                  pivot_longer(-bts) %>%
                                                  ggplot()+
                                                  geom_density(aes(value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_2_true.csv") %>%
                                                               rename(`β0` = V1, `β1` = V2, `β2` = V3, `β3` = V4, `β4` = V5, `ξ` = V6) %>%
                                                               pivot_longer(everything()),
                                                             aes(xintercept = value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  labs(y= "Density", x = "Parameter value")+
                                                  theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))

bias_correction_mod_0 = gridExtra::grid.arrange(read_csv("output/gpd_model_fits/model_0_uncorrected.csv",
                                                         col_names = c('bts', 'β0', 'β1', 'ξ')) %>%
                                                  pivot_longer(-bts) %>%
                                                  ggplot()+
                                                  geom_density(aes(value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_0_true.csv") %>%
                                                               rename(`β0` = V1, `β1` = V2, `ξ` = V3) %>%
                                                               pivot_longer(everything()),
                                                             aes(xintercept = value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  labs(y= "Density", x = "")+
                                                  theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1),
                                                                          axis.title.x=element_blank()),
                                                read_csv("output/gpd_model_fits/model_0_pars_corrected.csv",
                                                         col_names = c('bts', 'β0', 'β1', 'ξ')) %>%
                                                  pivot_longer(-bts) %>%
                                                  ggplot()+
                                                  geom_density(aes(value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_0_true.csv") %>%
                                                               rename(`β0` = V1, `β1` = V2,  `ξ` = V3) %>%
                                                               pivot_longer(everything()),
                                                             aes(xintercept = value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  labs(y= "Density", x = "Parameter value")+
                                                  theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))

bias_correction_mod_0 = gridExtra::grid.arrange(read_csv("output/gpd_model_fits/model_1_uncorrected.csv",
                                                         col_names = c('bts', 'β0', 'β1', 'β2', 'ξ')) %>%
                                                  pivot_longer(-bts) %>%
                                                  ggplot()+
                                                  geom_density(aes(value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
                                                               rename(`β0` = V1, `β1` = V2, `β2` = V3,  `ξ` = V4) %>%
                                                               pivot_longer(everything()),
                                                             aes(xintercept = value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  labs(y= "Density", x = "")+
                                                  theme_minimal(10)+
                                                  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1),
                                                        axis.title.x=element_blank()),
                                                read_csv("output/gpd_model_fits/model_1_pars_corrected.csv",
                                                         col_names = c('bts', 'β0', 'β1', 'β2', 'ξ')) %>%
                                                  pivot_longer(-bts) %>%
                                                  ggplot()+
                                                  geom_density(aes(value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
                                                               rename(`β0` = V1, `β1` = V2, `β2` = V3,  `ξ` = V4) %>%
                                                               pivot_longer(everything()),
                                                             aes(xintercept = value))+
                                                  facet_wrap(~name, scales = 'free', nrow = 1)+
                                                  labs(y= "Density", x = "Parameter value")+
                                                  theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))

ggsave(plot = bias_correction_mod_1, filename = "output/figs/bias_correction_mod_1.png", height = 3.75, width = 8)
ggsave(plot = bias_correction_mod_0, filename = "output/figs/bias_correction_mod_0.png", height = 3.75, width = 5.33333)
ggsave(plot = bias_correction_mod_0, filename = "output/figs/bias_correction_mod_0.png", height = 3.75, width = 4)
