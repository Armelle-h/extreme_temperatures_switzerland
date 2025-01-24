
gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#if the above has already been run

num_quantiles = 30

obs_data = readRDS(paste0("Data/processed/plain_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

legend = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")%>%
  select(stn, altitude)%>%
  unique()

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year") %>% 
  left_join(legend, by = "stn")

# ---- get covariates for prediction
temporal_covariates = obs_data %>%  #issue, be wary of, glob anom is defined for June July and August so each year is associated with 3 different values
  dplyr::select(year, glob_anom) %>% 
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # by default, they should've already been computed

if(fit_clim_quants){
  
  file.remove(paste0("Data/processed/plain_glob_anom_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
  # Loop over each quantile and fit the quantile regression model
  for(q in seq_along(quantiles_to_estimate_bulk)){
    
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    print(paste0("Fitting  tau = ", zeta))
    
    #in prep data bulk model we had associated to each data a quantile climate. 
    #regression only in terms of the climate quantile covariate
    
    # Set the 'value' column for current quantile
    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    
    if(q==1){
      inits_coeff = c(30, 0.2, 4, -4)
    }else {
      #to ensure smoothness in the estimated coefficients
      inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2],quantile_model_fit$location$coefficients[3])
    }
    
    
    # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
    quantile_model_fit <- evgam(maxtp ~ value + glob_anom + log(altitude), obs_data_for_quant_reg,
                                family = "ald", ald.args = list(tau = zeta))
    
    # Save the fitted parameter estimates for each quantile
    tibble(tau = zeta,
           beta_0 = quantile_model_fit$location$coefficients[1],
           beta_1 = quantile_model_fit$location$coefficients[2],
           beta_2 = quantile_model_fit$location$coefficients[3],
           beta_3 = quantile_model_fit$location$coefficients[4]) %>%
      write_csv(paste0("Data/processed/plain_glob_anom_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("Data/processed/plain_glob_anom_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3'))

# # --- creates a tibble with each station and its quantile model
# For each station, interpolate the quantiles and predict quantile values for each year
obs_smoothed_quantiles = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    # --- get the climate quantile estimates closest to current station
    clim_vals <<- obs_data %>%
      filter(stn == .x$stn[1]) %>%
      dplyr::select(quantile, value) %>%
      unique() %>%
      pull(value) %>%
      unlist()
    
    # predict (observed) quantile for each year and site
    quant_reg_pars = quant_reg_model_pars %>%
      arrange(tau)
    
    altitude_df = obs_data %>% filter(stn == .x$stn[1])%>%pull(altitude)
    altitude = altitude_df[[1]]
    
    res = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      # Calculate the predicted quantile value using the regression parameters
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates$year,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2) *(temporal_covariates$glob_anom)+(qpars$beta_3) *altitude))
    }
    
    print(paste0("Interpolating quantile estimates for ", .x$stn[1]))
    
    # interpolate quantiles over tau for each year
    res %>%
      group_by(year) %>%
      group_map(~{
        tibble(year = .x$year[1],
               #creates a spline function to go from tau to temp or temp to tau. 
               #the function depends on the year
               tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
               temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      mutate(stn = .x$stn[1])
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


# save quantile models
obs_smoothed_quantiles %>% saveRDS(paste0("output/plain_glob_anom_alt_quant_models_num_quantiles_",num_quantiles,".csv"))

#------------------------------------------------------------------------------------------------------------------

#obs_smoothed_quantiles = readRDS(paste0("output/plain_glob_anom_alt_quant_models_num_quantiles_",num_quantiles,".csv"))

# Calculate exceedance probability (lambda) for a threshold (threshold_9)
lambda_thresh_ex = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    # Calculate exceedance probability for threshold_9 
    #(threshold_9 using a regression model is an estimate of the 0.9 obs quantile)
    
    #using the temp_to_tau function computed above to find the tau associated with threshold_9
    thresh_exceedance_9 = obs_smoothed_quantiles%>%
      filter(stn == .x$stn[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_9[1], x))
    
    
    tibble(stn = .x$stn[1],
           year = temporal_covariates$year,
           thresh_exceedance_9 = 1-thresh_exceedance_9)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


lambda_thresh_ex %>%
  write_csv(paste0("Data/processed/plain_glob_anom_alt_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))