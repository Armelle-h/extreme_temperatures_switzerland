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
clim_data_extreme_95_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])
  
  #keeping only temperatures above the associated threshold 0.95 quantile
  sing_clim_data_extreme_95 = clim_data %>%
    group_by(id) %>%
    mutate(threshold = quantile(maxtp, 0.95),
           excess = maxtp - threshold) %>%
    filter(excess > 0) %>%
    ungroup()
  
  clim_data_extreme_95_list[[i]] = sing_clim_data_extreme_95
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_data_extreme_95 =  do.call(rbind, clim_data_extreme_95_list)

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
num_sites = clim_data_extreme_95$id %>% unique()

#don't forget to do some tests to change the interval
potential_shape_values_climate = seq(-0.26, -0.22, length.out = 50) #used to be length.out = 50, -0.23, -0.18

sub_clim_data_extreme_95 = clim_data_extreme_95 %>%
  dplyr::select(id, excess)%>%
  filter(id %% 5 == 0) 

#Computes log-likelihood value for potential shape parameter values in parallel (the cluster takes 4 hours to run)
num_cores <- detectCores() - 1

# Create a cluster
cl <- makeCluster(num_cores, type = "PSOCK")

clusterEvalQ(cl, library(dplyr))

# Export necessary objects and functions to the workers
clusterExport(cl, c("sub_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll"))

# Initialize variables for results
scales <- c()
loglik_sum <- c()

# Outer loop over potential shapes
for (potential_shape in potential_shape_values_climate) {
  print(potential_shape)
  
  # Define function to process each 'id' in parallel
  process_id <- function(i, potential_shape) {
    this_clim_extm_irel <- sub_clim_data_extreme_95 %>%
      filter(id == i) %>%
      pull(excess)
    
    model_fit <- estimate_scale_fixed_shape(this_clim_extm_irel, potential_shape)
    
    list(scale = model_fit$par, loglik = model_fit$value)
  }
  
  # Parallelize over unique ids using parLapply
  results <- parLapply(cl, unique(sub_clim_data_extreme_95$id), process_id, potential_shape)
  
  # Combine results from each worker
  scales <- c(scales, sapply(results, function(res) res$scale))
  loglik_sum <- c(loglik_sum, sum(sapply(results, function(res) res$loglik)))
}

# Stop the cluster after computations are done
stopCluster(cl)

#return the shape parameter minimizing the negative log-likelihood
optimal_shape = potential_shape_values_climate[which.min(loglik_sum)]

print(potential_shape_values_climate)
saveRDS(potential_shape_values_climate, file = "shape_candiates_95.rds")
print(loglik_sum)
saveRDS(loglik_sum, file = "associated_loglikelihood_95.rds")
print(optimal_shape)
saveRDS(optimal_shape, file = "optimal_shape_95.rds")
min_loglik_sum = min(loglik_sum)
print(min_loglik_sum)
saveRDS(min_loglik_sum, file = "optimal_loglikelihood_95.rds")

#DO THE SPLINE THING TO GET THE MAX AND THEN RUN THE CODE BELOW

potential_shape_values_climate = readRDS("shape_candiates_95.rds")
loglik_sum = readRDS("associated_loglikelihood_95.rds")

spline_loglik <- splinefun(potential_shape_values_climate, loglik_sum, method = "natural")

# Find the minimum of the spline in the range of potential_shape_values_climate
result <- optimize(spline_loglik, range(potential_shape_values_climate))

# Extract the minimum value and the pre-image (argmin)
estimated_optimal_loglik <- result$objective
optimal_shape <- result$minimum

# Display the results
cat("The minimum loglikelihood is:", estimated_optimal_loglik, "\n") #with this technique the optimal loglik go from 
cat("The optimal shape is:", optimal_shape, "\n")  #with this technique the optimal shape go from -0.2014286 (before) to -0.2019049 (after). 

#computing the associated loglikelihood ----------------------------------------


id_clim_data_extreme_95 = clim_data_extreme_95 %>%
  dplyr::select(id, excess)%>%
  filter(id %% 5 == 0)  #for comparison purpose with the other loglikelihoods computed 

process_id = function(i, optimal_shape){
  this_clim_extm_irel_95 <- id_clim_data_extreme_95 %>% filter(id == i) %>% pull(excess)
  # Estimate scale parameter with the optimal shape parameter
  model_fit_95 <- estimate_scale_fixed_shape(this_clim_extm_irel_95, optimal_shape)
  return(list(par = model_fit_95$par, loglik= model_fit_95$value))
}

job_process_id = function (indices, optimal_shape){
  scales = c()
  loglik_sum = 0
  for (i in indices){
    fun_output = process_id(i, optimal_shape)
    scales = c(scales, fun_output$par)
    loglik_sum = loglik_sum + fun_output$loglik
  }
  
  return (list(scales_=scales, loglik=loglik_sum))
}

indices <- unique(id_clim_data_extreme_95$id)

n <- length(unique(id_clim_data_extreme_95$id))

# Create a sequence of indices


# Split the indices into 5 chunks
chunk_size <- ceiling(n / 5)
chunks <- split(indices, ceiling(seq_along(indices) / chunk_size))

job::job ({
  result_1_loglik = job_process_id(chunks[[1]], optimal_shape)
  saveRDS(result_1_loglik, "result_1_loglik.rds")
  job::export(result_1_loglik)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_2_loglik = job_process_id(chunks[[2]], optimal_shape)
  saveRDS(result_2_loglik, "result_2_loglik.rds")
  job::export(result_2_loglik)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_3_loglik = job_process_id(chunks[[3]], optimal_shape)
  saveRDS(result_3_loglik, "result_3_loglik.rds")
  job::export(result_3_loglik)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_4_loglik = job_process_id(chunks[[4]], optimal_shape)
  saveRDS(result_4_loglik, "result_4_loglik.rds")
  job::export(result_4_loglik)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_5_loglik = job_process_id(chunks[[5]], optimal_shape)
  saveRDS(result_5_loglik, "result_5_loglik.rds")
  job::export(result_5_loglik)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

loglik_optimal_shape = result_1_loglik$loglik + result_2_loglik$loglik + result_3_loglik$loglik + result_4_loglik$loglik + result_5_loglik$loglik

saveRDS(loglik_optimal_shape, "optimal_loglikelihood_5_index_95.rds") #loglikelihood computed on one out of 5 indexes





#end of computing the associated loglikelihood ----------------------------------

#new version , the cluster takes two hours to run --> no use computing the associated loglikelihood to compare it to the others 
#since the other was computed on one id out of 5 and here we consider all ids. It would be useful to have a result it and of 
#itself 

id_clim_data_extreme_95 = clim_data_extreme_95 %>%
  dplyr::select(id, excess)

process_id = function(i, optimal_shape){
  this_clim_extm_irel_95 <- id_clim_data_extreme_95 %>% filter(id == i) %>% pull(excess)
  # Estimate scale parameter with the optimal shape parameter
  model_fit_95 <- estimate_scale_fixed_shape(this_clim_extm_irel_95, optimal_shape)
  return(list(par = model_fit_95$par, loglik= model_fit_95$value))
}

job_process_id = function (indices, optimal_shape){
  scales = c()
  loglik_sum = 0
  for (i in indices){
    fun_output = process_id(i, optimal_shape)
    scales = c(scales, fun_output$par)
    loglik_sum = loglik_sum + fun_output$loglik
  }
  
  return (list(scales_=scales, loglik=loglik_sum))
}

indices <- unique(id_clim_data_extreme_95$id)

n <- length(unique(id_clim_data_extreme_95$id))

# Create a sequence of indices


# Split the indices into 5 chunks
chunk_size <- ceiling(n / 5)
chunks <- split(indices, ceiling(seq_along(indices) / chunk_size))

job::job ({
  result_1 = job_process_id(chunks[[1]], optimal_shape)
  saveRDS(result_1, "result_1.rds")
  job::export(result_1)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_2 = job_process_id(chunks[[2]], optimal_shape)
  saveRDS(result_2, "result_2.rds")
  job::export(result_2)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_3 = job_process_id(chunks[[3]], optimal_shape)
  saveRDS(result_3, "result_3.rds")
  job::export(result_3)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_4 = job_process_id(chunks[[4]], optimal_shape)
  saveRDS(result_4, "result_4.rds")
  job::export(result_4)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

job::job ({
  result_5 = job_process_id(chunks[[5]], optimal_shape)
  saveRDS(result_5, "result_5.rds")
  job::export(result_5)
}, import=c("job_process_id", "chunks", "optimal_shape", "process_id", "id_clim_data_extreme_95", "estimate_scale_fixed_shape", "ngll")
, packages = c("tidyverse"))

loglik_optimal_shape = result_1$loglik + result_2$loglik + result_3$loglik + result_4$loglik + result_5$loglik
# Combine results from all chunks into one data frame
scales_95 <- c(result_1$scales_, result_2$scales_, result_3$scales_, result_4$scales_, result_5$scales_)


#end of new version

saveRDS(loglik_optimal_shape, "loglik_optimal_shape_95.rds") #loglikelihood computed for all indexes 

saveRDS(scales_95, "scales_95.rds")

# --- save estimates on climate grid

#keeping only location where the temperature is above the associated threshold 0.95 quantile
clim_data_extreme_95 %>%
  dplyr::select(id) %>%
  unique() %>%
  mutate(scale_95 = scales_95) %>%
  write_csv("Data/Climate_data/clim_scale_grid_95_gpd_model.csv")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

obs_sites = obs_data %>%
  dplyr::select(stn, id) %>%
  unique()

obs_sites = obs_sites %>%
  left_join(clim_data_extreme_95 %>%
              dplyr::select(id) %>%
              unique() %>%
              mutate(scale_95 = scales_95) %>%
              dplyr::select(id, scale_95), by = 'id')

obs_data %>%
  left_join(obs_sites) %>% 
  write_csv("Data/Observed_data/obs_data_95_gpd_model.csv")