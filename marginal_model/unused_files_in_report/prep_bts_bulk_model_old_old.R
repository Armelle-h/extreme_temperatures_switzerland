#in the second part, 1 file takes 10 minutes

gc() 
rm(list = ls())
library(tidyverse)
library(job)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')


standardise_data = T

#independent of the bootstrap thing
if(standardise_data){
  
  num_quantiles = 30
  
  obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")
  
  legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")
  
  #joining obs data with the altitude 
  obs_data = obs_data %>%
    left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
    rename(altitude = Altitude.m.)
  
  #adding threshold 9 
  
  threshold_9_df = vroom::vroom("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")%>%
    dplyr::select(threshold_9, id) %>%
    unique()
  
  #lambda associated with observed data, vary the quantile model
  lambda_df = vroom::vroom(paste0("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv")) 
  
  glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")
  
  glob_anomaly_reshaped = glob_anomaly %>%
    select(c("year", "JJA"))%>%
    rename(glob_anom = JJA)
  
  obs_data = obs_data %>% 
    mutate(year = year(date)) %>%
    left_join(lambda_df, by=c("stn", "year")) %>%
    left_join(glob_anomaly_reshaped, by = "year") %>%
    left_join(threshold_9_df, by="id")
  
  rm("lambda_df", "glob_anomaly", "glob_anomaly_reshaped", "threshold_9_df", "legend_data") #free as much space as possible
  
  quantile_models = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv")) #the regression model selected
  
  quantile_models = quantile_models %>%
    select(-c(tau_to_temp))
  
  extreme_data = obs_data %>%
    group_by(stn) %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)
  
  this_fit_mod_0 = fit_mod_0(extreme_data$excess, extreme_data$scale_9)
  this_fit_mod_1 = fit_mod_1(extreme_data$excess, extreme_data$scale_9, extreme_data$glob_anom)
  this_fit_mod_2 = fit_mod_2(extreme_data$excess, extreme_data$scale_9, extreme_data$glob_anom, extreme_data$altitude)
  this_fit_mod_3 = fit_mod_3(extreme_data$excess, extreme_data$scale_9, extreme_data$altitude)
  
  #freeing memory
  rm('extreme_data')
  
  pred_0 = my_predict_0(this_fit_mod_0, obs_data$scale_9)
  pred_1 = my_predict_1(this_fit_mod_1, obs_data$scale_9, obs_data$glob_anom)
  pred_2 =my_predict_2(this_fit_mod_2, obs_data$scale_9, obs_data$glob_anom, obs_data$altitude)
  pred_3 =my_predict_3(this_fit_mod_3, obs_data$scale_9, obs_data$altitude)
  
  obs_data = obs_data %>%
    mutate(scale_0 = pred_0$scale,
           scale_1 = pred_1$scale,
           scale_2 = pred_2$scale,
           scale_3 = pred_3$scale,
           shape_0 = pred_0$shape,
           shape_1 = pred_1$shape,
           shape_2 = pred_2$shape,
           shape_3 = pred_3$shape)
  
  rm('this_fit_mod_0', 'this_fit_mod_1',  'this_fit_mod_2',  'this_fit_mod_3','pred_0', 'pred_1', 'pred_2', 'pred_3')
  
  print("Joining bulk models")
  obs_data = obs_data %>%
    left_join(quantile_models)
  
  obs_data = obs_data %>%
    group_by(stn) %>%
    group_map(~{
      print(.x$stn[1])
      data = .x$maxtp
      threshold = .x$threshold_9
      
      scle_0 = .x$scale_0
      shpe_0 = .x$shape_0
      
      scle_1 = .x$scale_1
      shpe_1 = .x$shape_1
      
      scle_2 = .x$scale_2
      shpe_2 = .x$shape_2
      
      scle_3 = .x$scale_3
      shpe_3 = .x$shape_3
      
      yr = .x$year
      my_lambda = .x$thresh_exceedance_9
      
      res_0 = rep(NA, length(data))
      res_1 = rep(NA, length(data))
      res_2 = rep(NA, length(data))
      res_3 = rep(NA, length(data))
      
      for(i in seq(length(data))){
        if(data[i] > threshold[i]){ # transform tail
          
          res_0[i] = 1 - (my_lambda[i]) *(1+shpe_0[i]* ((data[i] - threshold[i])/scle_0[i]))^(-1/(shpe_0[i]))
          res_1[i] = 1 - (my_lambda[i]) *(1+shpe_1[i]* ((data[i] - threshold[i])/scle_1[i]))^(-1/(shpe_1[i]))
          res_2[i] = 1 - (my_lambda[i]) *(1+shpe_2[i]* ((data[i] - threshold[i])/scle_2[i]))^(-1/(shpe_2[i]))
          res_3[i] = 1 - (my_lambda[i]) *(1+shpe_3[i]* ((data[i] - threshold[i])/scle_3[i]))^(-1/(shpe_3[i]))
          
        }else{ 
          res_0[i] = .x$temp_to_tau[i][[1]](data[i])
          res_1[i] = res_0[i]
          res_2[i] = res_0[i]
          res_3[i] = res_0[i]
        }
      }
      .x$unif_0 = res_0
      .x$unif_1 = res_1
      .x$unif_2 = res_2
      .x$unif_3 = res_3
      .x
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble()
  
  #before 305 GB heavy, after 99.1 MB heavy 
  obs_data = obs_data  %>%
    select(-c(temp_to_tau))
  
  obs_data$unif_0[obs_data$unif_0 < 0] = 0 
  obs_data$unif_1[obs_data$unif_1 < 0] = 0
  obs_data$unif_2[obs_data$unif_2 < 0] = 0
  obs_data$unif_3[obs_data$unif_3 < 0] = 0
  
  #not awesome bypass but should do the trick 
  obs_data$unif_0[obs_data$unif_0 == 1] = 0.9999999
  obs_data$unif_1[obs_data$unif_1 == 1] = 0.9999999
  obs_data$unif_2[obs_data$unif_2 == 1] = 0.9999999
  obs_data$unif_3[obs_data$unif_3 == 1] = 0.9999999
  
  obs_data %>%
    saveRDS(paste0("Data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles))
}

#Sanity check if unif_0,..., unif_3 does have values in (0,1) 

invalid_values <- obs_data %>%
  filter(if_any(c(unif_0, unif_1, unif_2, unif_3), ~ . < 0 | . > 1))

print(nrow(invalid_values)) #0 rows of invalid values :)


#END OF FIXED PART

apply_tau_to_temp = function (df, val1, val3, val){
  
  spline_function <- df$tau_to_temp[df$year == val1 & df$stn == val3]
  result <- spline_function[[1]](val)
  return(result)
}

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')

num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles)) %>% arrange(date)

quantile_models = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv")) #only need tau_to_temp

quantile_models = quantile_models %>%
  select(-c("temp_to_tau"))

print("Finished reading quantile_models")

#obs_data = obs_data %>% left_join(quantile_models)  nope, otherwise memory issues. Designing functions to directly access the data I'm interested in

#rm('quantile_models')

sites_by_date = readRDS("Data/processed/bootstrap_data/sites_by_date") 
bootstrap_dates = readRDS("Data/processed/bootstrap_data/bootstrapped_dates") 

print("Finished reading sites and bootstrap_dates")

covars = obs_data %>%
  dplyr::select(stn, year, glob_anom, altitude,
                threshold_9, scale_9,thresh_exceedance_9,
                # "tau_to_temp", 
                "scale_0", 
                "scale_1",
                "scale_2",,
                "scale_3",
                "shape_0",
                "shape_1",
                "shape_2",
                "shape_3") %>% unique()

unique_values <- unique(obs_data$stn)

chunk_size <- ceiling(length(unique_values) / 4)
group_chunks <- split(unique_values, ceiling(seq_along(unique_values) / chunk_size))


df_1 = obs_data %>% filter(stn %in% group_chunks[[1]])
df_2 = obs_data %>% filter(stn %in% group_chunks[[2]])
df_3 = obs_data %>% filter(stn %in% group_chunks[[3]])
df_4 = obs_data %>% filter(stn %in% group_chunks[[4]])

cv_1 = covars %>% filter(stn %in% group_chunks[[1]])
cv_2 = covars %>% filter(stn %in% group_chunks[[2]])
cv_3 = covars %>% filter(stn %in% group_chunks[[3]])
cv_4 = covars %>% filter(stn %in% group_chunks[[4]])

qm_1 = quantile_models %>% filter(stn %in% group_chunks[[1]])
qm_2 = quantile_models %>% filter(stn %in% group_chunks[[2]])
qm_3 = quantile_models %>% filter(stn %in% group_chunks[[3]])
qm_4 = quantile_models %>% filter(stn %in% group_chunks[[4]])




process_station_data <- function(obs_data_fun, quantile_models_fun, covars_fun, i, batch_size) {

bts_rng = seq(i,(i+batch_size-1))

for(b in bts_rng){
  
  print(b)
  bts_data <- rlang::duplicate(obs_data_fun, shallow = FALSE) # make a deep copy
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
  
  expanded_data = expanded_data %>%
    mutate(year = lubridate::year(date)) %>%
    left_join(covars_fun)
  
  print("Starting expanded_data computations")
    
  expanded_data = expanded_data %>%
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
        .x$maxtp_3 = NA
        
        for(i in seq(nrow(.x))){
          
          #if the proba associated with the element is above the exceedance threshold, then compute the tp0 using the fact 
          #we know the data follows a gpd distribution 
          #else use the tau_to_temp function which is defined on the body of the distribution 
          
          # --- model 0
          if(!is.na(my_lambda[i]) && !is.na(.x$unif_0[i]) && .x$unif_0[i] > (1-my_lambda[i])){
            
            proba = (1 + (.x$unif_0[i]-1)/(my_lambda[i]))
            
            if (proba>=1 | proba <=0) {
              .x$maxtp_0[i] = NA
            }
            else{
              .x$maxtp_0[i] =  evd::qgpd(p =(1 + (.x$unif_0[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_0[i],shape=shpe_0[i])
            }
            
          }
          else{
            
            if (!is.na(.x$unif_0[i])) {
              .x$maxtp_0[i] = apply_tau_to_temp(quantile_models_fun, .x$year[i], .x$stn[i], .x$unif_0[i])
            }
            else{
              .x$maxtp_0[i] = NA
            }
            
          }
          
          # --- model 1
          if(!is.na(my_lambda[i]) && !is.na(.x$unif_1[i]) && .x$unif_1[i] > (1-my_lambda[i])){
            
            proba = (1 + (.x$unif_1[i]-1)/(my_lambda[i]))
            
            if (proba>=1 | proba <=0) {
              .x$maxtp_1[i] = NA
            }
            else{
              .x$maxtp_1[i] =  evd::qgpd( p =(1 + (.x$unif_1[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_1[i],shape=shpe_1[i])
            }
            
          }
          else{
            if (!is.na(.x$unif_1[i])) {
              .x$maxtp_1[i] = apply_tau_to_temp(quantile_models_fun, .x$year[i], .x$stn[i], .x$unif_1[i])
            }
            else{
              .x$maxtp_1[i] = NA
            }
          }
          
          # --- model 2
          if(!is.na(my_lambda[i]) && !is.na(.x$unif_2[i]) && .x$unif_2[i] > (1-my_lambda[i])){
            
            proba = (1 + (.x$unif_2[i]-1)/(my_lambda[i]))
            
            if (proba>=1 | proba <=0) {
              .x$maxtp_2[i] = NA
            }
            else{
              .x$maxtp_2[i] =  evd::qgpd( p =(1 + (.x$unif_2[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_2[i],shape=shpe_2[i])
              
            }
            
          }else{
            
            if (!is.na(.x$unif_2[i])) {
              .x$maxtp_2[i] = apply_tau_to_temp(quantile_models_fun, .x$year[i], .x$stn[i], .x$unif_2[i])
            }
            else{
              .x$maxtp_2[i] = NA
            }
            
          }
          
          # --- model 3
          if(!is.na(my_lambda[i]) && !is.na(.x$unif_3[i]) && .x$unif_3[i] > (1-my_lambda[i])){
            
            proba = (1 + (.x$unif_3[i]-1)/(my_lambda[i]))
            
            if (proba>=1 | proba <=0) {
              .x$maxtp_3[i] = NA #temp solution after fixing stuff shouldn't be used anymore 
            }
            else{
              
              .x$maxtp_3[i] =  evd::qgpd( p =(1 + (.x$unif_3[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_3[i],shape=shpe_3[i])
              
            }
            
            
          }else{
            
            if (!is.na(.x$unif_3[i])) {
              .x$maxtp_3[i] = apply_tau_to_temp(quantile_models_fun, .x$year[i], .x$stn[i], .x$unif_3[i])
            }
            else{
              .x$maxtp_3[i] = NA
            }
            
          }
        }
        .x
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()%>%
      dplyr::select(stn, date, scale_9, threshold_9, maxtp_0, maxtp_1, maxtp_2, maxtp_3) %>%
      saveRDS(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_",b))
  }
}

  
i = 1
batch_size = 1


  job::job({
    
    process_station_data(df_1, qm_1, cv_1, i, batch_size)
  
  }, import = c("process_station_data", "apply_tau_to_temp", "df_1", "qm_1", "num_quantiles", "cv_1", "i", "batch_size", "bootstrap_dates", "sites_by_date"))
  
  job::job({
    
    process_station_data(df_2, qm_2, cv_2, i, batch_size)
    
  }, import = c("process_station_data", "apply_tau_to_temp", "df_2", "qm_2", "num_quantiles", "cv_2", "i", "batch_size", "bootstrap_dates", "sites_by_date"))
  
  job::job({
    
    process_station_data(df_3, qm_3, cv_3, i, batch_size)
    
  }, import = c("process_station_data", "apply_tau_to_temp", "df_3", "qm_3", "num_quantiles", "cv_3", "i", "batch_size", "bootstrap_dates", "sites_by_date"))
  
  job::job({
    
   process_station_data(df_4, qm_4, cv_4, i, batch_size)
    
  }, import = c("process_station_data", "apply_tau_to_temp", "df_4", "qm_4", "num_quantiles", "cv_4", "i", "batch_size", "bootstrap_dates", "sites_by_date")) 
