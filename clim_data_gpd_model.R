#estimating the shape and scale parameter of the tail of the climate data
#the tail being modeled as a GPD
#In the paper, we're only interested in the scale parameter so the scale is not saved

gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(evd)
library(mgcv)
library(data.table)
library(raster) # package for netcdf manipulation

#think how to do the for loop 


# ---- Climate model

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_data_extreme_9_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])

  #keeping only temperatures above the associated threshold 0.9 quantile
  sing_clim_data_extreme_9 = clim_data %>%
    group_by(id) %>%
    mutate(threshold = quantile(maxtp, 0.9),
         excess = maxtp - threshold) %>%
    filter(excess > 0) %>%
    ungroup()

  clim_data_extreme_9_list[[i]] = sing_clim_data_extreme_9

  #to free memory
  rm(clim_data)
  gc()
}

clim_data_extreme_9 =  do.call(rbind, clim_data_extreme_9_list)

# Negative log-likelihood function for GPD
ngll = function(par){
  if(par <= 0) return(2^30) # Return large penalty if the scale parameter is non-positive
  if(par > -1/shape_param) return(2^30) # Return penalty if the scale does not satisfy GPD constraint
  if(any((1+shape_param*this.dat/par)< 0)) return(2^30) # Check if GPD condition is violated
  
  # Compute negative log-likelihood using the GPD density function (log=True for log-likelihood)
  -sum(evd::dgpd(x = this.dat, loc=0, scale = par, shape=shape_param, log=T))
}

# Function to estimate the scale parameter with a fixed shape parameter
estimate_scale_fixed_shape = function(x,shape_c){
  this.dat <<- x
  shape_param <<- shape_c
  #optimizing negative log likelihood using Brent
  optim(par = c(0), fn = ngll, method = 'Brent', lower=0, upper = 5)
}

# fit scale parameter to each location with constant shape
num_sites = clim_data_extreme_9$id %>% unique()
scales = c()
loglik_sum = c()

#Computes log-likelihood value for potential shape parameter values
for(potential_shape in potential_shape_values_climate){
  print(potential_shape)
  
  loglik = c()
  
  for(i in (clim_data_extreme_9$id %>% unique())){
    this_clim_extm_irel = clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)
    
    #estimate scale parameter and computes associated log-likelihood
    model_fit = estimate_scale_fixed_shape(this_clim_extm_irel, potential_shape)
    
    scales = c(scales, model_fit$par)
    loglik = c(loglik, model_fit$value)
    
  }
  loglik_sum = c(loglik_sum, sum(loglik))
}

#return the shape parameter minimizing the negative log-likelihood
optimal_shape = potential_shape_values_climate[which.min(loglik_sum)]

scales_9 = c()
for(i in (clim_data_extreme_9$id %>% unique())){
  this_clim_extm_irel_9 = clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)
  #estimate scale parameter with optimal shape parameter
  model_fit_9 = estimate_scale_fixed_shape(this_clim_extm_irel_9, optimal_shape_9)
  scales_9 = c(scales_9, model_fit_9$par)
}

# --- save estimates on climate grid

#keeping only location where the temperature is above the associated threshold 0.9 quantile
clim_data_extreme_9 %>%
  dplyr::select(longitude, latitude, id) %>%
  unique() %>%
  mutate(scale_9 = scales_9) %>%
  write_csv("Data/Climate_data/clim_scale_grid_gpd_model.csv")

obs_sites = obs_sites %>%
  left_join(clim_data_extreme_9 %>%
              dplyr::select(longitude, latitude, id) %>%
              unique() %>%
              mutate(scale_9 = scales_9) %>%
              dplyr::select(id, scale_9), by = 'id')

obs_data %>%
  left_join(obs_sites) %>% 
  write_csv("Data/Observed_data/obs_data_gpd_model.csv")