
gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/processed/plain_1971_2022_JJA_obs_data_bulk_model.csv")

num_quantiles = 30

clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

obs_data = obs_data %>%
  left_join(clim_quantiles_subset)

obs_data %>%
  saveRDS(paste0("Data/processed/plain_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

#obs_data = readRDS(paste0("Data/processed/plain_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

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

if(fit_clim_quants){
  
  file.remove(paste0("Data/processed/plain_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
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
      write_csv(paste0("Data/processed/plain_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("Data/processed/plain_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

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
obs_smoothed_quantiles %>% saveRDS(paste0("output/plain_glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

#------------------------------------------------------------------------------------------------------------------

#obs_smoothed_quantiles = readRDS(paste0("output/plain_glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))

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
  write_csv(paste0("Data/processed/plain_glob_anomaly_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))


# ------------ get splines on clim scale

gc()
rm(list = ls())
library(tidyverse)
library(evgam)
library(data.table)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

temporal_covariates = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)%>% 
  filter(year %in% c(1971, 2022))%>% 
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

quant_reg_model_pars = read_csv(paste0("Data/processed/plain_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

#need to create
clim_grid = read_csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")%>%
  filter(id %% 10 == 0)


clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  filter(id %% 10 == 0)

#to have scale_9 and quantile and value in the same dataframe
clim_grid = clim_grid %>%
  left_join(clim_quantiles_subset, by = "id")


#climate data with quantile modified
clim_date_w_quantile_mod = c()

for(i in seq(nrow(clim_grid))){
  
  if (i%%100 == 0){print(i)}
  
  # Create data for prediction by merging with temporal covariates
  this_data_for_pred = tibble(year = c(1971, 2022)) %>%
    #one is before observed range, we have the extremes and one is after the observed range
    left_join(temporal_covariates, by = "year") %>%
    mutate(id = clim_grid[i,]$id,
           quantile = clim_grid[i,]$quantile,
           value = clim_grid[i,]$value) 
  
  # predict quantile for each year and site
  quant_reg_pars = quant_reg_model_pars %>%
    arrange(tau)
  
  clim_vals = clim_grid[i,]$value[[1]]
  res = c()
  
  # Loop through each quantile to estimate predicted values
  for(q in seq_along(quantiles_to_estimate_bulk)){
    qpars = quant_reg_pars[q,]
    res = rbind(res,
                tibble(quantile =  qpars$tau,
                       year = this_data_for_pred$year,
                       quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates$glob_anom)   ))
  }
  
  # Fit splines to interpolate between quantile values
  #in order to find relationship between tau and temp
  res = res %>%
    group_by(year) %>%
    group_map(~{
      tibble(year = .x$year[1],
             tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
             temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble() 
  
  # Merge predictions with the grid data
  #quantile mod as now they are estimated using a regression model (no longer empirical quantile)
  clim_date_w_quantile_mod = rbind(clim_date_w_quantile_mod, 
                                   this_data_for_pred %>% left_join(res, by = "year"))
  
}

clim_date_w_quantile_mod %>%
  saveRDS(paste0("output/plain_glob_anomaly_quant_models_clim_num_quantiles_",num_quantiles,".csv"))

quantile_model_fit_9 = readRDS("output/plain_threshold_model_9")
clim_date_w_quantile_mod = readRDS(paste0("output/plain_quant_models_clim_num_quantiles_",num_quantiles,".csv"))

#loop over my climate data files

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])%>%
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

clim_thresh = do.call(rbind, clim_thresh_values_list)

# Add the threshold value to the model data
clim_date_w_quantile_mod = clim_date_w_quantile_mod %>% left_join(clim_thresh)
clim_date_w_quantile_mod$threshold_9= predict(quantile_model_fit_9, clim_date_w_quantile_mod)$location

# Calculate lambda (exceedance probability) for the threshold
# (difference compared to before is that now we're working with the regression estimated climate quantile)
lambda_thresh_ex = clim_date_w_quantile_mod %>%
  group_by(id) %>%
  group_map(~{
    
    thresh_exceedance_9 = clim_date_w_quantile_mod%>%
      filter(id == .x$id[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_9[1], x))# Apply the threshold function
    
    # Create a tibble with exceedance probabilities
    tibble(id = .x$id[1],
           year = c(1971, 2022),
           thresh_exceedance_9 = 1-thresh_exceedance_9)# Inverse exceedance probability  
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

#saving result
lambda_thresh_ex %>% 
  write_csv(paste0("Data/processed/plain_climate_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))

# Merge lambda results with the main model data and save
clim_date_w_quantile_mod %>% 
  left_join(lambda_thresh_ex) %>%
  saveRDS(paste0("output/plain_quant_models_clim_num_quantiles_",num_quantiles,".csv"))