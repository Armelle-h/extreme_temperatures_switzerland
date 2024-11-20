gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


sigmoid= function(altitude, alt_threshold){
  result = 1/(1+exp(-200*(altitude-alt_threshold)))
  
  return(result)
}

num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>%
  rename(altitude = Altitude.m.)

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, altitude), by="stn")

#joining obs data with global anomaly

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

temporal_covariates_year = temporal_covariates$year
temporal_covariates_glob_anom = temporal_covariates$glob_anom

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # bool, re-estimate clim quantiles?

if(fit_clim_quants){
  
  file.remove(paste0("Data/processed/debug_glob_anom_indicator_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  
  # Loop over each quantile and fit the quantile regression model
  for(q in seq_along(quantiles_to_estimate_bulk)){
    
    if(q==1 ){
      inits_coeff = c(-3, 1, 3, -1)
    } else{
      inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2],quantile_model_fit$location$coefficients[3])
    }
    
    if (q == 23){
      inits_coeff = c(1, 2, -6, 1)
    }
    
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    print(paste0("Fitting  tau = ", zeta))
    
    #in prep data bulk model we had associated to each data a quantile climate. 
    #regression only in terms of the climate quantile covariate
    
    # Set the 'value' column for current quantile
    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    
    # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
    quantile_model_fit <- evgam(maxtp ~ (1-sigmoid(altitude, 1500))*value + glob_anom + sigmoid(altitude, 1500)*value, obs_data_for_quant_reg,
                                family = "ald", inits = inits_coeff, ald.args = list(tau = zeta))
    
    # Save the fitted parameter estimates for each quantile
    tibble(tau = zeta,
           beta_0 = quantile_model_fit$location$coefficients[1],
           beta_1 = quantile_model_fit$location$coefficients[2],
           beta_2 = quantile_model_fit$location$coefficients[3],
           beta_3 = quantile_model_fit$location$coefficients[4]
           ) %>%
      write_csv(paste0("Data/processed/debug_glob_anom_indicator_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

#I can run only until here for starters, the line below is especially used if the quantile
#regression coeff were not already computed

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("Data/processed/debug_glob_anom_indicator_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3'))



stn_name = "AAR"

clim_vals <<- obs_data %>%
  filter(stn == stn_name) %>%
  dplyr::select(quantile, value) %>%
  unique() %>%
  pull(value) %>%
  unlist()

# predict (observed) quantile for each year and site
quant_reg_pars = quant_reg_model_pars %>%
  arrange(tau)


alt_df = legend_data%>%filter(stn == stn_name) %>%select(altitude)
alt = alt_df$altitude

res = c()
for(q in seq_along(quantiles_to_estimate_bulk)){
  qpars = quant_reg_pars[q,]
  # Calculate the predicted quantile value using the regression parameters
  res = rbind(res,
              tibble(quantile =  qpars$tau,
                     year = temporal_covariates_year,
                     quant_value = qpars$beta_0 + qpars$beta_1*(1-sigmoid(alt, 1500))*clim_vals[q] + (qpars$beta_2)*(temporal_covariates_glob_anom) 
                     + (qpars$beta_3) *sigmoid(alt, 1500)*clim_vals[q] ) )
}




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
    
    
    alt_df = legend_data%>%filter(stn == .x$stn[1]) %>%select(altitude)
    alt = alt_df$altitude
    
    res = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      # Calculate the predicted quantile value using the regression parameters
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates_year,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates_glob_anom) 
                         + (qpars$beta_3) *sigmoid(alt, 1500)*clim_vals[q] ) )
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


num_quantiles = 100 
quantiles_to_estimate = seq(0.001,0.99,length.out = num_quantiles)

obs_data_quantile <- obs_data %>%
  select(-id) %>%
  mutate(year = year(date)) %>%
  select(-date) %>%
  group_by(stn, year) %>%
  summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))

extract_quantile <- function(df){
  df_pred_quantiles = df %>%
    mutate(quantiles = map(tau_to_temp, ~ sapply(quantiles_to_estimate, .x))) %>%
    select(year, stn, quantiles)
}

pred = extract_quantile(obs_smoothed_quantiles)

calculate_rmse <- function(list1, list2) {
  if (length(list1) == length(list2)) {
    return(sqrt(mean((list1 - list2)^2, na.rm = TRUE)))  # Calculate RMSE
  } else {
    return(NA)  # Return NA if the lengths are different
  }
}

total_rmse <- function(df){
  RMSE_df = df%>%
    full_join(obs_data_quantile, by = c("stn", "year"), suffix = c("_pred", "_empirical"))
  
  # Calculate the RMSE for corresponding quantile lists
  results <- RMSE_df %>%
    mutate(rmse = map2_dbl(quantiles_pred, quantiles_empirical, calculate_rmse))
  
  # Sum the RMSE over all years and stations
  total_rmse_df <- results %>%
    summarise(total_rmse = sqrt(sum(rmse^2, na.rm = TRUE) / sum(!is.na(rmse)))) 
  
  total_rmse_value <- total_rmse_df$total_rmse
}

rmse_val = total_rmse(pred)

print(rmse_val)
























obs_smoothed_quantiles %>% saveRDS(paste0("output/glob_anom_indicator_alt_quant_models_num_quantiles_",num_quantiles,".csv"))


lambda_thresh_ex = obs_data %>%  #five elements have proba less than 0. 
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
           year = temporal_covariates_year,
           thresh_exceedance_9 = 1-thresh_exceedance_9)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()





































num_quantiles = 100 
quantiles_to_estimate = seq(0.001,0.99,length.out = num_quantiles)

obs_data_quantile <- obs_data %>%
  select(-id) %>%
  mutate(year = year(date)) %>%
  select(-date) %>%
  group_by(stn, year) %>%
  summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))

extract_quantile <- function(df){
  df_pred_quantiles = df %>%
    mutate(quantiles = map(tau_to_temp, ~ sapply(quantiles_to_estimate, .x))) %>%
    select(year, stn, quantiles)
}

pred = extract_quantile(obs_smoothed_quantiles)

calculate_rmse <- function(list1, list2) {
  if (length(list1) == length(list2)) {
    return(sqrt(mean((list1 - list2)^2, na.rm = TRUE)))  # Calculate RMSE
  } else {
    return(NA)  # Return NA if the lengths are different
  }
}

total_rmse <- function(df){
  RMSE_df = df%>%
    full_join(obs_data_quantile, by = c("stn", "year"), suffix = c("_pred", "_empirical"))
  
  # Calculate the RMSE for corresponding quantile lists
  results <- RMSE_df %>%
    mutate(rmse = map2_dbl(quantiles_pred, quantiles_empirical, calculate_rmse))
  
  # Sum the RMSE over all years and stations
  total_rmse_df <- results %>%
    summarise(total_rmse = sqrt(sum(rmse^2, na.rm = TRUE) / sum(!is.na(rmse)))) 
  
  total_rmse_value <- total_rmse_df$total_rmse
}

rmse_val = total_rmse(pred)

print(rmse_val)