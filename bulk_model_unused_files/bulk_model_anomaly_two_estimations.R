gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

L = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")
L_alt = L %>% 
  select(c("stn", "Altitude.m."))%>%
  rename(altitude = Altitude.m.)

obs_data = obs_data %>%
  left_join(L_alt, by="stn")

obs_data_1 = obs_data %>%
  filter(altitude<1500)
  
obs_data_2 = obs_data  %>%
  anti_join(obs_data_1)

obs_data = obs_data_1

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

# ---- get covariates for prediction
temporal_covariates = obs_data %>%  #issue, be wary of, glob anom is defined for June July and August so each year is associated with 3 different values
  dplyr::select(year, glob_anom) %>% 
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # by default, they should've already been computed

#inside the if loop depends only on climate data, not the observed data
if(fit_clim_quants){
  
  file.remove(paste0("Data/processed/debug_1_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
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
      inits_coeff = c(-2, 1, 3)
    }
    if (q ==21){
      inits_coeff = c(0, 1, 3)
    }
    if (q %in% c(24, 25, 27, 28) ){
      inits_coeff = c(0.3, 1, 3)
    }
    if (!(q %in% c(1, 21, 24, 25, 27, 28))) {
      #to ensure smoothness in the estimated coefficients
      inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2],quantile_model_fit$location$coefficients[3])
    }
    
    # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
    quantile_model_fit <- evgam(maxtp ~ value + glob_anom, obs_data_for_quant_reg,
                                family = "ald", inits = inits_coeff ,ald.args = list(tau = zeta))
    
    # Save the fitted parameter estimates for each quantile
    tibble(tau = zeta,
           beta_0 = quantile_model_fit$location$coefficients[1],
           beta_1 = quantile_model_fit$location$coefficients[2],
           beta_2 = quantile_model_fit$location$coefficients[3]) %>%
      write_csv(paste0("Data/processed/debug_1_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("Data/processed/debug_1_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))




obs_data = obs_data_2

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

if(fit_clim_quants){
  
  file.remove(paste0("Data/processed/v2_debug_2_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
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
      inits_coeff = c(-2, 1, 3)
    }else {
      #to ensure smoothness in the estimated coefficients
      inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2],quantile_model_fit$location$coefficients[3])
    }
    
    # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
    quantile_model_fit <- evgam(maxtp ~ value + glob_anom, obs_data_for_quant_reg,
                                family = "ald", inits = inits_coeff ,ald.args = list(tau = zeta))
    
    # Save the fitted parameter estimates for each quantile
    tibble(tau = zeta,
           beta_0 = quantile_model_fit$location$coefficients[1],
           beta_1 = quantile_model_fit$location$coefficients[2],
           beta_2 = quantile_model_fit$location$coefficients[3]) %>%
      write_csv(paste0("Data/processed/v2_debug_2_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed




gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

# ---- get covariates for prediction
temporal_covariates = obs_data %>%  #issue, be wary of, glob anom is defined for June July and August so each year is associated with 3 different values
  dplyr::select(year, glob_anom) %>% 
  unique() %>%
  arrange(year)

L =read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>% 
  select(c("stn", "Altitude.m."))%>%
  rename(altitude = Altitude.m.)

stations_below_2000 = L%>% 
  filter(altitude<2000)%>% 
  select(stn)



quant_reg_model_pars_1 = read_csv(paste0("Data/processed/debug_1_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))
# --- read in fitted quantile regression coefficients
quant_reg_model_pars_2 = read_csv(paste0("Data/processed/v2_debug_2_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))















quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

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
    if (.x$stn[1] %in% stations_below_2000$stn){
      quant_reg_pars = quant_reg_model_pars_1 %>%
        arrange(tau)
      
    } else {
      quant_reg_pars = quant_reg_model_pars_2 %>%
        arrange(tau)      
    }
    
    res = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      # Calculate the predicted quantile value using the regression parameters
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates$year,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2) *(temporal_covariates$glob_anom)))
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
obs_smoothed_quantiles %>% saveRDS(paste0("output/debug_glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

obs_smoothed_quantiles = readRDS(paste0("output/debug_glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

#I have negative proba :( !!! NOt a good thing, needs to be fixed

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









quant_reg_model_pars = quant_reg_model_pars[-((nrow(quant_reg_model_pars) - 5):nrow(quant_reg_model_pars)), ]

write.csv(quant_reg_model_pars, paste0("Data/processed/debug_1_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), row.names = FALSE)
