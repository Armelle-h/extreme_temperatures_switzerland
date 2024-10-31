#In this file I'm computing the empirircal climate quantiles for quantile values 
#in quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

#for everything, takes 25 minutes to run :)

#All good, don't forget to readRDS the files

#25 minutes for 40 quantiles, 34 for 50 quantiles

gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(raster) # package for netcdf manipulation
library(plyr)
library(purrr)
library(furrr)
library(future)
library(data.table)
library(tictoc)
# 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

# # i. ========  ========  Global parameters ========  ========  ========
# marginal threshold 

# Start timing
tic.clearlog()

tic("Timing the group processing")

num_quantiles = 50
quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

raw_obs_data<- read.csv("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv", header=TRUE)
obs_data <- raw_obs_data %>%
  filter(maxtp != "-")#filtering out rows with missing values
obs_data$maxtp <- as.integer(obs_data$maxtp) #converting last column elements as integers


files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)

clim_quant_list = list()

for (i in seq_along(files)){

  clim_data = fread(files[[i]])

  #----- up until here okay ! 

  #parallelized version
  available_cores = parallel::detectCores()

  plan(multisession, workers = available_cores -1 )

  clim_quantiles = clim_data %>%
    group_by(id) %>%
    group_split() %>%
    future_map(~{
      res = c()
      for(q in quantiles_to_estimate_bulk){
        res = rbind(res, tibble(id = .x$id[1],
                              date_id = .x$date_id[1],
                              quantile = q,
                              value = as.numeric(quantile(.x$maxtp, q))))
    }
    res
  }) %>%
  bind_rows() %>%
  as_tibble()

  #summarizes clim_quantiles by putting all computed quantiles in a list associated with location id
  sing_clim_quantiles_subset = clim_quantiles %>%
    group_by(id) %>%
    group_split() %>%
    future_map(~{
      tibble(id = .x$id[1], 
           quantile = list(.x$quantile), 
           value = list(.x$value))
    }) %>%
    bind_rows() %>%
    as_tibble()
  
  clim_quant_list[[i]] = sing_clim_quantiles_subset

  # After the parallel work, revert to sequential processing
  plan(sequential)

  #end of parallelized version
  #removing clim_data from memory as the file is quite heavy 
  rm(clim_data)
  gc()
}

clim_quantiles_subset = do.call(rbind, clim_quant_list)

clim_quantiles_subset %>%
  saveRDS(paste0("Data/processed/debug_clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

#all good until here -----------------------------------------------

#need to set a column id corresponding to unique pair of longitude latitude
#in obs data 
obs_data = obs_data %>%
  left_join(clim_quantiles_subset)

obs_data %>%
  saveRDS(paste0("Data/processed/debug_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

toc(log = TRUE)

# Retrieve the logged time
time_log <- tic.log(format = TRUE)

# Print the time
print(time_log)