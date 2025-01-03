gc() 
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')
  
  num_quantiles = 30
  obs_data = readRDS(paste0("Data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles)) %>% arrange(date)
  
  quantile_models = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))%>%
    select(-c("temp_to_tau"))
  
  obs_data = obs_data %>% left_join(quantile_models) 
  
  rm('quantile_models')
  
  sites_by_date = readRDS("data/processed/bootstrap_data/sites_by_date")
  bootstrap_dates = readRDS("data/processed/bootstrap_data/bootstrapped_dates")
  
  bts_rng = seq(1, 100)
  
  for(b in bts_rng){
    
    if(b<=75){next}
    print(b)
    
    bts_data <- rlang::duplicate(obs_data, shallow = FALSE) # make a deep copy
    sites_by_date$bts.dat = bootstrap_dates[[b]]
    
    dates_and_their_data = bts_data %>%
      group_by(date) %>%
      group_map(~{
        tibble(date = .x$date[1],
               stn = list(.x$stn),
               unif_0 = list(.x$unif_0),
               unif_1 = list(.x$unif_1),
               unif_2 = list(.x$unif_2),
               unif_3 = list(.x$unif_3))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    rm('bts_data')
    
    dates_and_their_data_copy = rlang::duplicate(dates_and_their_data, shallow = FALSE)
    
    for(i in seq(nrow(sites_by_date))){
      data_to_replace = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]
      sampled_data = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$bts.dat,]
      sites_to_replace = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$stn[[1]]
      sites_in_bts_date = sampled_data$stn[[1]]
      ind_of_sites_to_rep = which(sites_in_bts_date %in% sites_to_replace)
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_0 = list(sampled_data$unif_0[[1]][ind_of_sites_to_rep])
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_1 = list(sampled_data$unif_1[[1]][ind_of_sites_to_rep])
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_2 = list(sampled_data$unif_2[[1]][ind_of_sites_to_rep])
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_3 = list(sampled_data$unif_3[[1]][ind_of_sites_to_rep])
    }
    
    rm('dates_and_their_data')
    
    expanded_data = dates_and_their_data_copy %>%
      group_by(date) %>%
      group_map(~{
        tibble(date = .x$date,
               stn = .x$stn[[1]],
               unif_0 = .x$unif_0[[1]],
               unif_1 = .x$unif_1[[1]],
               unif_2 = .x$unif_2[[1]],
               unif_3 = .x$unif_3[[1]])
      }, .keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    rm('dates_and_their_data_copy')
    
    covars = obs_data %>%
      dplyr::select(stn, year, glob_anom,
                    threshold_9, scale_9,thresh_exceedance_9,
                    "tau_to_temp", 
                    "scale_0", 
                    "scale_1",
                    "scale_2",
                    "scale_3",
                    "shape_0",
                    "shape_1",
                    "shape_2",
                    "shape_3") %>% unique()
    
    expanded_data = expanded_data %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(covars)
    
    # ---- back transform
    expanded_data%>%
      group_by(stn) %>%
      group_map(~{
        
        scle_0 = .x$scale_0
        shpe_0 = .x$shape_0
        
        scle_1 = .x$scale_1
        shpe_1 = .x$shape_1
        
        scle_2 = .x$scale_2
        shpe_2 = .x$shape_2
        
        scle_3 = .x$scale_3
        shpe_3 = .x$shape_3
        
        threshold = .x$threshold_9
        my_lambda = .x$thresh_exceedance_9
        
        .x$maxtp_0 = NA
        .x$maxtp_1 = NA
        .x$maxtp_2 = NA
        
        for(i in seq(nrow(.x))){
          
          #if the proba associated with the element is above the exceedance threshold, then compute the tp0 using the fact 
          #we know the data follows a gpd distribution 
          #else use the tau_to_temp function which is defined on the body of the distribution 
          
          # --- model 0
          if(.x$unif_0[i] > (1-my_lambda[i])){
            .x$maxtp_0[i] =  evd::qgpd(p =(1 + (.x$unif_0[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_0[i],shape=shpe_0[i])
          }else{
            .x$maxtp_0[i] = .x$tau_to_temp[i][[1]](.x$unif_0[i])
          }
          
          # --- model 1
          if(.x$unif_1[i] > (1-my_lambda[i])){
            .x$maxtp_1[i] =  evd::qgpd( p =(1 + (.x$unif_1[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_1[i],shape=shpe_1[i])
          }else{
            .x$maxtp_1[i] = .x$tau_to_temp[i][[1]](.x$unif_1[i])
          }
          
          # --- model 2
          if(.x$unif_2[i] > (1-my_lambda[i])){
            .x$maxtp_2[i] =  evd::qgpd( p =(1 + (.x$unif_2[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_2[i],shape=shpe_2[i])
          }else{
            .x$maxtp_2[i] = .x$tau_to_temp[i][[1]](.x$unif_2[i])
          }
          # --- model 3
          if(.x$unif_3[i] > (1-my_lambda[i])){
            .x$maxtp_3[i] =  evd::qgpd( p =(1 + (.x$unif_3[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_3[i],shape=shpe_3[i])
          }else{
            .x$maxtp_3[i] = .x$tau_to_temp[i][[1]](.x$unif_3[i])
          }
        }
        .x
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      dplyr::select(stn, date, scale_9, threshold_9, maxtp_0, maxtp_1, maxtp_2, maxtp_3) %>%
      write.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_",b, ".csv"), row.names=FALSE, col.names = FALSE)
  }