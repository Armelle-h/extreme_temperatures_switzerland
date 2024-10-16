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
library(parallel)
library(tictoc)
library(raster) # package for netcdf manipulation

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

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

potential_shape_values_climate = seq(-0.23, -0.18, length.out = 50)

sub_clim_data_extreme_9 = clim_data_extreme_9 %>%
  dplyr::select(id, excess)%>%
  filter(id %% 5 == 0)

#Computes log-likelihood value for potential shape parameter values in parallel
num_cores <- detectCores() - 1

# Create a cluster
cl <- makeCluster(num_cores, type = "PSOCK")

clusterEvalQ(cl, library(dplyr))

# Export necessary objects and functions to the workers
clusterExport(cl, c("sub_clim_data_extreme_9", "estimate_scale_fixed_shape", "ngll"))

# Initialize variables for results
scales <- c()
loglik_sum <- c()

# Outer loop over potential shapes
for (potential_shape in potential_shape_values_climate) {
  print(potential_shape)
  
  # Define function to process each 'id' in parallel
  process_id <- function(i, potential_shape) {
    this_clim_extm_irel <- sub_clim_data_extreme_9 %>%
      filter(id == i) %>%
      pull(excess)
    
    model_fit <- estimate_scale_fixed_shape(this_clim_extm_irel, potential_shape)
    
    list(scale = model_fit$par, loglik = model_fit$value)
  }
  
  # Parallelize over unique ids using parLapply
  results <- parLapply(cl, unique(sub_clim_data_extreme_9$id), process_id, potential_shape)
  
  # Combine results from each worker
  scales <- c(scales, sapply(results, function(res) res$scale))
  loglik_sum <- c(loglik_sum, sum(sapply(results, function(res) res$loglik)))
}

# Stop the cluster after computations are done
stopCluster(cl)

#return the shape parameter minimizing the negative log-likelihood
optimal_shape = potential_shape_values_climate[which.min(loglik_sum)]

print(potential_shape_values_climate)
saveRDS(potential_shape_values_climate, file = "shape_candiates.rds")
print(loglik_sum)
saveRDS(loglik_sum, file = "associated_loglikelihood.rds")
print(optimal_shape)
saveRDS(optimal_shape, file = "optimal_shape.rds")
min_loglik_sum = min(loglik_sum)
print(min_loglik_sum)
saveRDS(min_loglik_sum, file = "optimal_loglikelihood.rds")