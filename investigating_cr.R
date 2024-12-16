gc()
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(data.table)
library(evd)

#spatial_threshold = readRDS("Data/spatial_threshold.rds")
spatial_threshold = 8

filter_id = 10

num_quantiles = 30

#quantiles are only defined with respect to space. Not great for functin approximation but doing with that 
#could be interesting to do a QQplot to see how good this all is.
clim_grid = read_csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")%>%
  filter(id %% filter_id == 0)


clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  filter(id %% filter_id == 0)


clim_smooth = clim_quantiles_subset %>%
  group_by(id) %>%
  group_map(~{
    tibble( id = .x$id[[1]],
            tau_to_temp = list(splinefun(unlist(.x$quantile),unlist(.x$value),  method = 'monoH.FC')),
           temp_to_tau = list(splinefun(unlist(.x$value),unlist(.x$quantile),  method = 'monoH.FC')))
  }, .keep = T)%>%
  plyr::rbind.fill() %>%
  as_tibble() 


files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])%>%
    filter(id %% filter_id == 0)%>%
    filter(id %in% I_plain$id)
  
  # Calculate the climate threshold values (0.9 quantile of 'maxtp')
  sing_clim_thresh = clim_data %>%
    group_by(id) %>%
    summarise(clim_thresh_value_9 = quantile(maxtp, 0.9, na.rm = TRUE))
  
  clim_thresh_values_list[[i]] = sing_clim_thresh
  
  #to free memory
  rm(clim_data)
  gc()
  
}

#give access to clim_thresh_value_9
clim_thresh = do.call(rbind, clim_thresh_values_list)

# Calculate lambda (exceedance probability) for the threshold
# (difference compared to before is that now we're working with the regression estimated climate quantile)
lambda_thresh_ex = clim_thresh %>%
  group_by(id) %>%
  group_map(~{
    
    thresh_exceedance_9 = clim_smooth%>%
      filter(id == .x$id[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$clim_thresh_value_9[1], x))# Apply the threshold function
    
    # Create a tibble with exceedance probabilities
    tibble(id = .x$id[1],
           thresh_exceedance_9 = 1-thresh_exceedance_9)# Inverse exceedance probability  
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

shape = readRDS("Data/clim_data_gpd_model/plain_optimal_shape.rds")




#division needs to be done time_wise !!!!!

exceedance_indicator = 0

files = list.files(path = "Data/Climate_data/By_year", full.names = TRUE)

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

I_date = read.csv("Data/id_date_correspondance.csv")%>%
  mutate(year = lubridate::year(date))

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])%>%
    filter(id %% filter_id == 0)%>%
    filter(id %in% I_plain$id)
  
  I_date_reduced = I_date %>%
    filter(date_id %in% clim_data$date_id)
    
  clim_data = clim_data %>%
    left_join(I_date, by = "date_id")%>%
    select(-date_id)
  
  lambda_thresh_ex_reduced = lambda_thresh_ex %>%
    filter(id %in% clim_data$id)
  
  clim_smooth_reduced = clim_smooth %>%
    filter(id %in% clim_data$id)
  
  clim_grid_reduced = clim_grid %>%
    filter(id %in% clim_data$id)
  
  clim_thresh_reduced = clim_thresh %>%
    filter(id %in% clim_data$id)
  
  clim_data = clim_data %>%
    left_join(lambda_thresh_ex_reduced, by = "id")%>%
    left_join(clim_smooth_reduced, by = "id")%>%
    left_join(clim_thresh_reduced, by = "id")%>%
    left_join(clim_grid_reduced, by = "id")
  
  
  clim_data_standardised = clim_data %>%  #will have to do a for loop over the clim data files
    group_by(id, year) %>% #the function is just defined with respect to space, not time
    group_map(~{
      data = .x$maxtp
      threshold = .x$clim_thresh_value_9
      res = rep(NA, length(data))
      num_extremes = sum(data>threshold)
      
      if(num_extremes >0){ # if extreme obs in this year at this site
        scle = .x$scale_9
        shpe = shape
        my_lambda = .x$thresh_exceedance_9
        
        res[data > threshold] = 1 - my_lambda[data > threshold]*(1-evd::pgpd((data[data > threshold] - threshold[data > threshold]), loc = 0,scale = scle[data > threshold], shape = shpe[1]))
      }
      
      this_quant_mod = clim_smooth_reduced %>%
        filter(id == .x$id[[1]]) %>%
        pull(temp_to_tau) %>%
        .[[1]]
      
      res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
      res[res<0]=0
      res[res>0.999999]=0.999999
      
      .x$unif = res
      #.x$frechet_marg = -1/log(.x$unif)
      .x$pareto_marg = 1/(1-.x$unif)
      .x
    },.keep = T) %>%
    plyr::rbind.fill()%>%
    as_tibble() 
  

  #would be to heavy to save as a list so doing cr
  
  
  extreme_dates = clim_data_standardised %>%
    dplyr::group_by(date) %>%
    #dplyr::summarise(cost = mean(frechet_marg)) %>%
    dplyr::summarise(cost = median(pareto_marg)) %>%  #dplyr::summarise(cost = mean(pareto_marg)) %>%
    ungroup() %>%
    arrange(desc(cost))
  
  exceedance_indicator =  exceedance_indicator + nrow(extreme_dates %>% filter(cost > spatial_threshold))
  
  #to free memory
  rm(clim_data, lambda_thresh_ex_reduced, clim_smooth_reduced, clim_data_standardised, extreme_dates)
  gc()
  
}


estimate = exceedance_indicator * spatial_threshold / nrow(I_date) #supposed to be between 0 and 1   #idea is then to make spatial threshold vary to see if the estimate cr holds

print(estimate)