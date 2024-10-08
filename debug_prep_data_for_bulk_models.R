#testing if works for year 2017
#takes 11~15 minutes to run

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
setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

# # i. ========  ========  Global parameters ========  ========  ========
# marginal threshold 

# Start timing
tic.clearlog()

tic("Timing the group processing")

num_quantiles = 30
quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

raw_obs_data<- read.csv("Data/Observed_data/1971_2023_JJA_obs_data_loc_id.csv", header=TRUE)
obs_data <- raw_obs_data %>%
  filter(maxtp != "-")#filtering out rows with missing values
obs_data$maxtp <- as.integer(obs_data$maxtp) #converting last column elements as integers

obs_data <- obs_data %>%
  filter(
    (substr(date, 1, 4) %in% c("2017"))    # Check for years 2016 or 2017
  )

clim_data = fread("Data/Climate_data/2017_JJA_climate_data.csv")

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

#summarizes clim_quantiles by putting all computed quantiles in a list asociated with location id
clim_quantiles_subset = clim_quantiles %>%
  group_by(id) %>%
  group_split() %>%
  future_map(~{
    tibble(id = .x$id[1], 
           quantile = list(.x$quantile), 
           value = list(.x$value))
  }) %>%
  bind_rows() %>%
  as_tibble()

# After the parallel work, revert to sequential processing
plan(sequential)

#end of parallelized version

clim_quantiles_subset %>%
  saveRDS(paste0("Data/Climate_data/2016_2017_clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

#all good until here -----------------------------------------------

#need to set a column id corresponding to unique pair of longitude latitude
#in obs data 
obs_data = obs_data %>%
  left_join(clim_quantiles_subset)

obs_data %>%
  saveRDS(paste0("Data/Observed_data/2016_2017_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

toc(log = TRUE)

# Retrieve the logged time
time_log <- tic.log(format = TRUE)

# Print the time
print(time_log)