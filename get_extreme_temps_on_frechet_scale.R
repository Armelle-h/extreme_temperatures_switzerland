gc()
rm(list = ls())
library(tidyverse)

num_quantiles = 30

marg_mod = "mod_1"
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')

obs_sites = read_csv("Data/Observed_data/plain_obs_data_gpd_model.csv") %>%
  dplyr::select(stn, date, scale_9, id) %>%
  unique()

threshold_9_df = vroom::vroom("Data/processed/plain_1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, stn)%>%
  unique()

lambda_df = vroom::vroom(paste0("Data/processed/plain_glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv") )

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_grid = obs_sites %>% 
  mutate(year = year(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_9_df, by="stn")

#only doing for 2 models

if(marg_mod == "mod_0"){
  # ----- model 0
  obs_grid$scale = my_predict_0(as.numeric(unlist(read_csv("output/gpd_model_fits/model_0_true.csv"))), obs_grid$scale_9)$scale
  obs_grid$shape = my_predict_0(as.numeric(unlist(read_csv("output/gpd_model_fits/model_0_true.csv"))), obs_grid$scale_9)$shape
}else if(marg_mod == "mod_1"){
  # # ----- model 1
  obs_grid$scale = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/model_1_true.csv"))), obs_grid$scale_9, obs_grid$glob_anom)$scale
  obs_grid$shape = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/model_1_true.csv"))), obs_grid$scale_9, obs_grid$glob_anom)$shape
}

extreme_temps = seq(22.9, 36, by = 0.25)  #will need to be adapted

res = c()
for(tmp in extreme_temps){
  res = cbind(res,(-1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - obs_grid$threshold_9),  scale = obs_grid$scale, shape =  obs_grid$shape[1])))))
}

res = as_tibble(res)
names(res) = as.character(extreme_temps)

res %>%
  mutate(year = obs_grid$year,
         stn = obs_grid$stn) %>%
  pivot_longer(-c(year, stn)) %>%
  rename(temp = name,
         frechet_value = value) %>%
  mutate(temp = as.numeric(temp)) %>%
  write_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv"))


# # ---- Climate data ----------------------

#NOT THE SMAE AS OBS SMOOTHED QUANTILE !!!

num_quantiles = 30

#defined for id %%10 ==0
clim_smoothed_quantiles = readRDS(paste0("output/plain_glob_anomaly_quant_models_clim_num_quantiles_",num_quantiles,".csv"))

#I have id and scale_9, defined for all id
clim_grid = read.csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")%>%
  filter(id %% 10 == 0)

marg_mod = "mod_1"

clim_grid = clim_grid %>%
  left_join(clim_smoothed_quantiles %>%
              dplyr::select(id, glob_anom, tau_to_temp, threshold_9, thresh_exceedance_9, year), by = "id")

this_fit_mod = read_csv("output/gpd_model_fits/model_1_true.csv") %>%
  unlist() %>% as.numeric
clim_grid$shape = this_fit_mod[length(this_fit_mod)]
clim_grid$scale = my_predict_1(this_fit_mod, clim_grid$scale_9, clim_grid$glob_anom)$scale


extreme_temps = seq(28, 38, by = 0.25)  #used to be seq(26, 36, by = 0.25)  | no clue, would have to read the theory on Pareto processes

res = c()
for(tmp in extreme_temps){
  res = cbind(res,(-1/log(1-clim_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - clim_grid$threshold_9),  scale = clim_grid$scale, shape =  clim_grid$shape[1])))))
}

res = as_tibble(res)
names(res) = as.character(extreme_temps)

res %>%
  mutate(year = clim_grid$year,
         id = clim_grid$id) %>%
  pivot_longer(-c(year, id)) %>%
  rename(temp = name,
         frechet_value = value) %>%
  mutate(temp = as.numeric(temp)) %>%
  write_csv(paste0("output/plain_clim_grid_extreme_temps_frechet_scale_",marg_mod,".csv"))