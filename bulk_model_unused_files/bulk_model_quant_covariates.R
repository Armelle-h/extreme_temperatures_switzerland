# This script fits and saves quantile regression model and lambda estimates

#until spline models for climate data, takes 12 minutes to run

gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

# ---- get covariates for prediction
temporal_covariates = obs_data %>%
  dplyr::select(year) %>% #keeping only year
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # bool, re-estimate clim quantiles?

if(fit_clim_quants){
  
  file.remove(paste0("bulk_model_unused_files/processed/quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
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
    
    if (q==1 | q==7 | q==8){
      init_coeff = c(-1, 1)
    } else {
      init_coeff = c(quantile_model_fit$location$coefficients[1], quantile_model_fit$location$coefficients[2])
    }
    
    # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
    quantile_model_fit <- evgam(maxtp ~ value, obs_data_for_quant_reg,
                                family = "ald", ald.args = list(tau = zeta), init_coeff)
    
    # Save the fitted parameter estimates for each quantile
    tibble(tau = zeta,
           beta_0 = quantile_model_fit$location$coefficients[1],
           beta_1 = quantile_model_fit$location$coefficients[2]) %>%
      write_csv(paste0("bulk_model_unused_files/processed/quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("bulk_model_unused_files/processed/quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1'))

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
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q]))
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
#obs_smoothed_quantiles %>% saveRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
#obs_smoothed_quantiles = readRDS("output/quant_models")


# Calculate exceedance probability (lambda) for a threshold (threshold_9)
lambda_thresh_ex = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    print(.x$stn[1])
    
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
  write_csv(paste0("Data/processed/thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))


# ------------ get splines on clim scale

#need to create
clim_grid = read_csv("Data/id_lon_lat_correspondance.csv")

#depending on what is needed, might need to do a for loop
clim_dat_full = read_csv("Data/Climate_data/to_complete.csv")

clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

#don't see the point, unless it is to have more flexibility about which localisation data we're considering
clim_grid = clim_grid %>%
  left_join(clim_quantiles_subset)


#climate data with quantile modified
clim_date_w_quantile_mod = c()

for(i in seq(nrow(clim_grid))){
  print(clim_grid[i,]) #for debugging
  
  # Create data for prediction by merging with temporal covariates
  this_data_for_pred = tibble(year = c(1960, 1971, 2023, 2024)) %>%
    #one is before observed range, we have the extremes and one is after the observed range
    left_join(temporal_covariates) %>%
    mutate(longitude = clim_grid[i,]$longitude,
           latitude = clim_grid[i,]$latitude,
           id = clim_grid[i,]$id,
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
                       quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] ))
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
                                   this_data_for_pred %>% left_join(res))
  
}

clim_date_w_quantile_mod %>%
  saveRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))

#can stop here as the file is read right below

# threshold at these points 

quantile_model_fit_9 = readRDS("output/threshold_model_9")
clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))

#loop over my climate data files

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

for (i in seq_along(files)){

  clim_data = fread(files[[i]])
  
  # Calculate the climate threshold values (0.9 quantile of 'maxtp')
  sing_clim_thresh = clim_dat_full %>%
    group_by(id) %>%
    summarise(clim_thresh_value_9 = quantile(maxtp, 0.9, na.rm = TRUE))
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
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
           year = c(1960, 1971, 2023, 2024),
           thresh_exceedance_9 = 1-thresh_exceedance_9)# Inverse exceedance probability  
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

#saving result
lambda_thresh_ex %>% 
  write_csv(paste0("Data/processed/climate_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))

# Merge lambda results with the main model data and save
clim_date_w_quantile_mod %>% 
  left_join(lambda_thresh_ex) %>%
  saveRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))