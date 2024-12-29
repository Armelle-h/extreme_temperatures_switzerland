gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

obs_data_extreme = obs_data %>%
  group_by(stn)%>%
  filter(maxtp>quantile(maxtp, 0.999))%>%
  ungroup


obs_data_extreme = obs_data %>%
  group_by(stn)%>%
  mutate(quantile_0.999 = quantile(maxtp, 0.999))%>%
  ungroup 

obs_data_extreme_filter = obs_data_extreme %>%
  filter(maxtp>quantile_0.999)


n_out = 100 

tau_list = seq(0.5, 1, length.out=n_out)[-n_out]


nb_measurements = c()

for (tau in tau_list){
  
  T = obs_data %>%
    group_by(stn) %>%
    filter(maxtp>quantile(maxtp, tau))%>% 
    ungroup
  
  nb_measurements = c(nb_measurements, nrow(T) ) 
  
}

plot(tau_list, nb_measurements)


n_out = 100 

tau_list = seq(0, 0.05, length.out=n_out)[-1]


nb_measurements = c()

for (tau in tau_list){
  
  T = obs_data %>%
    group_by(stn) %>%
    filter(maxtp<quantile(maxtp, tau))%>% 
    ungroup
  
  nb_measurements = c(nb_measurements, nrow(T) ) 
  
}

plot(tau_list, nb_measurements)