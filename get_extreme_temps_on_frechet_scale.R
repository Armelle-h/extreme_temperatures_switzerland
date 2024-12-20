gc()
rm(list = ls())
library(tidyverse)

num_quantiles = 30

marg_mod = "mod_0"
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')

#all stations are associated with a value
obs_sites = read_csv("Data/Observed_data/plain_obs_data_gpd_model_025.csv") %>%
  dplyr::select(stn, scale_9) %>%
  unique()

#all stations are associated with a value
threshold_9_df = vroom::vroom("Data/processed/plain_1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, stn)%>%
  unique()

#thresh_exceedance proba defined for all stations and all years
lambda_df = vroom::vroom(paste0("Data/processed/plain_glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv") )

#defined for all years
glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_grid = lambda_df %>% 
  left_join(obs_sites, by = "stn") %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_9_df, by="stn")

#only doing for 2 models

if(marg_mod == "mod_0"){
  # ----- model 0
  obs_grid$scale = my_predict_0(as.numeric(unlist(read_csv("output/gpd_model_fits/plain_model_0_true_025.csv"))), obs_grid$scale_9)$scale
  obs_grid$shape = my_predict_0(as.numeric(unlist(read_csv("output/gpd_model_fits/plain_model_0_true_025.csv"))), obs_grid$scale_9)$shape
}else if(marg_mod == "mod_1"){
  # # ----- model 1
  obs_grid$scale = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/plain_model_1_true_025.csv"))), obs_grid$scale_9, obs_grid$glob_anom)$scale
  obs_grid$shape = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/plain_model_1_true_025.csv"))), obs_grid$scale_9, obs_grid$glob_anom)$shape
}

extreme_temps = seq(30, 46, by = 0.25)  #converting values of interest to their frechet and pareto scale depending on location and time

res_frechet = c()
res_pareto = c()

for(tmp in extreme_temps){
  #converts to uniform scale using the distribution function for tail data in the bulk model then converts to frechet
  #or pareto scale
  #this is done on tmp-threshold_9 which is the definition domain
  res_frechet = cbind(res_frechet,(-1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - obs_grid$threshold_9),  scale = obs_grid$scale, shape =  obs_grid$shape[1])))))
  res_pareto = cbind(res_pareto,(1/(obs_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - obs_grid$threshold_9),  scale = obs_grid$scale, shape =  obs_grid$shape[1])))))
}

res_frechet = as_tibble(res_frechet)
res_pareto = as_tibble(res_pareto)
names(res_frechet) = as.character(extreme_temps)
names(res_pareto) = as.character(extreme_temps)


#converting the temperature into the frechet scale or pareto scale depending on the year and the location
res_frechet = res_frechet %>%
  mutate(year = obs_grid$year,
         stn = obs_grid$stn) %>%
  pivot_longer(-c(year, stn)) %>%
  rename(temp = name,
         frechet_value = value) %>%
  mutate(temp = as.numeric(temp))

res_pareto = res_pareto %>%
  mutate(year = obs_grid$year,
         stn = obs_grid$stn) %>%
  pivot_longer(-c(year, stn)) %>%
  rename(temp = name,
         pareto_value = value) %>%
  mutate(temp = as.numeric(temp))

if (identical(res_frechet[, c("year", "stn", "temp")], res_pareto[, c("year", "stn", "temp")])) {
  # Add col4 from df1 to df2
  res_frechet$pareto_value <- res_pareto$pareto_value
} else {
  stop("Dataframes do not align. Check year, stn, temp.")
}

#defined for all years, all stations, all temperatures in "extreme_temps"
res_frechet %>%
  write_csv(paste0("output/plain_obs_sites_extreme_temps_frechet_pareto_scale_",marg_mod,".csv"))


# # ---- Climate data ----------------------

#NOT THE SAME AS OBS SMOOTHED QUANTILE !!!

num_quantiles = 30

#defined for id %%10 ==0
clim_smoothed_quantiles = readRDS(paste0("output/plain_glob_anomaly_quant_models_clim_num_quantiles_",num_quantiles,".csv"))

#I have id and scale_9, defined for all id
clim_grid = read.csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")%>%
  filter(id %% 10 == 0)

clim_grid = clim_grid %>%
  left_join(clim_smoothed_quantiles %>%
              dplyr::select(id, glob_anom, tau_to_temp, threshold_9, thresh_exceedance_9, year), by = "id")

if(marg_mod == "mod_1"){
  
this_fit_mod = read_csv("output/gpd_model_fits/plain_model_1_true_025.csv") %>%
  unlist() %>% as.numeric
clim_grid$shape = this_fit_mod[length(this_fit_mod)]
clim_grid$scale = my_predict_1(this_fit_mod, clim_grid$scale_9, clim_grid$glob_anom)$scale

}
if(marg_mod == "mod_0"){
  
  this_fit_mod = read_csv("output/gpd_model_fits/plain_model_0_true_025.csv") %>%
    unlist() %>% as.numeric
  clim_grid$shape = this_fit_mod[length(this_fit_mod)]
  clim_grid$scale = my_predict_0(this_fit_mod, clim_grid$scale_9)$scale
  
}

extreme_temps = seq(30, 46, by = 0.25)

res_frechet = c()
res_pareto = c()
for(tmp in extreme_temps){
  #converts to uniform scale using the distribution function for tail data in the bulk model then converts to frechet
  #or pareto scale
  #this is done on tmp-threshold_9 which is the definition domain 
  res_frechet = cbind(res_frechet,(-1/log(1-clim_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - clim_grid$threshold_9),  scale = clim_grid$scale, shape =  clim_grid$shape[1])))))
  res_pareto = cbind(res_pareto,(1/(clim_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - clim_grid$threshold_9),  scale = clim_grid$scale, shape =  clim_grid$shape[1])))))
}

res_frechet = as_tibble(res_frechet)
res_pareto = as_tibble(res_pareto)
names(res_frechet) = as.character(extreme_temps)
names(res_pareto) = as.character(extreme_temps)

res_frechet = res_frechet %>%
  mutate(year = clim_grid$year,
         id = clim_grid$id) %>%
  pivot_longer(-c(year, id)) %>%
  rename(temp = name,
         frechet_value = value) %>%
  mutate(temp = as.numeric(temp))

res_pareto = res_pareto %>%
  mutate(year = clim_grid$year,
         id = clim_grid$id) %>%
  pivot_longer(-c(year, id)) %>%
  rename(temp = name,
         pareto_value = value) %>%
  mutate(temp = as.numeric(temp))

if (identical(res_frechet[, c("year", "id", "temp")], res_pareto[, c("year", "id", "temp")])) {
  # Add col4 from df1 to df2
  res_frechet$pareto_value <- res_pareto$pareto_value
} else {
  stop("Dataframes do not align. Check year, stn, temp.")
}

res_frechet %>%
  write_csv(paste0("output/plain_clim_grid_extreme_temps_frechet_pareto_scale_",marg_mod,".csv"))