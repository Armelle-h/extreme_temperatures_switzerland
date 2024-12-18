#the difference between this file for the plain and the whole country is that the model is more stable for the plain so we need less regularising conditions

library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#putting obs_data and obs_data_quantile as arguments
extract_temperature <- function(df, quantiles_to_estimate){
  df_pred_quantiles = df %>%
    mutate(quantiles = map(tau_to_temp, ~ sapply(quantiles_to_estimate, .x))) %>%
    select(year, stn, quantiles)
  
  return(df_pred_quantiles)
}

calculate_rmse <- function(list1, list2) {
  if (length(list1) == length(list2)) {
    return(sqrt(mean((list1 - list2)^2, na.rm = TRUE)))  # Calculate RMSE
  } else {
    return(NA)  # Return NA if the lengths are different
  }
}

total_rmse <- function(df, obs_data_quantile){
  RMSE_df = df%>%
    full_join(obs_data_quantile, by = c("stn", "year"), suffix = c("_pred", "_empirical"))
  
  # Calculate the RMSE for corresponding quantile lists
  results <- RMSE_df %>%
    mutate(rmse = map2_dbl(quantiles_pred, quantiles_empirical, calculate_rmse))
  
  return(mean(results$rmse, na.rm = TRUE))
}

calculate_mae <- function(list1, list2) {
  if (length(list1) == length(list2)) {
    return((mean(abs(list1 - list2), na.rm = TRUE)))  # Calculate RMSE
  } else {
    return(NA)  # Return NA if the lengths are different
  }
}

total_mae <- function(df, obs_data_quantile){
  MAE_df = df%>%
    full_join(obs_data_quantile, by = c("stn", "year"), suffix = c("_pred", "_empirical"))
  
  # Calculate the RMSE for corresponding quantile lists
  results <- MAE_df %>%
    mutate(mae = map2_dbl(quantiles_pred, quantiles_empirical, calculate_mae))
  
  return(mean(results$mae, na.rm = TRUE)) #NOT OPTIONAL AS EMPIRICAL ARE NOT DEFINED FOR ALL STN YEAR PAIR
}

total_mae_quantile <- function(df, quantiles_to_estimate){
  
  results <- df %>%
    rowwise() %>%
    mutate(
      mae = calculate_mae(unlist(quantile_estimate), quantiles_to_estimate)
    ) %>%
    ungroup()
  
  return(mean(results$mae, na.rm = TRUE)) #NOT OPTIONAL AS EMPIRICAL ARE NOT DEFINED FOR ALL STN YEAR PAIR
}

total_rmse_quantile <- function(df, quantiles_to_estimate){
  
  results <- df %>%
    rowwise() %>%
    mutate(
      rmse = calculate_rmse(unlist(quantile_estimate), quantiles_to_estimate)
    ) %>%
    ungroup()
  
  return(mean(results$rmse, na.rm = TRUE)) #NOT OPTIONAL AS EMPIRICAL ARE NOT DEFINED FOR ALL STN YEAR PAIR
}

regression <- function(qpars, q, clim_vals, model_name, temporal_covariates, altitude){
  if (!(model_name %in% c("no_covariate", "quant", "quant_glob_anom", "quant_log_alt"))) {
    print("model name is not accepted")
    break
  }
  
  if (model_name == "no_covariate"){
    result = qpars$beta_0
  }
  
  if (model_name == "quant"){
    result = qpars$beta_0 + qpars$beta_1*clim_vals[q]
  }
  if (model_name == "quant_glob_anom"){
    result = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates$glob_anom)
  }
  if (model_name == "quant_log_alt"){
    
    result = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates$glob_anom) + (qpars$beta_3) *log(altitude)
  }
  
  return(result)
}


estimate_parameters <- function(quantiles_to_estimate_bulk, obs_data, model_name){ #model_name describes which bulk model are we computing results for
  
  if (!(model_name %in% c("no_covariate", "quant", "quant_glob_anom", "quant_log_alt"))) {
    print("model name is not accepted")
    break
  }
  
  results_df = tibble() #initialising empty dataframe
  
  for(q in seq_along(quantiles_to_estimate_bulk)){
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    #print(paste0("Fitting  tau = ", zeta))
    
    #in prep data bulk model we had associated to each data a quantile climate. 
    #regression only in terms of the climate quantile covariate
    
    # Set the 'value' column for current quantile
    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    
    if (model_name == "no_covariate"){
      
      # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
      quantile_model_fit <- evgam(maxtp ~ 1, obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1])
    }
    
    if (model_name == "quant"){
      # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
      if(q==1){
        inits_coeff = c(-1, 1)
      }else {
        #to ensure smoothness in the estimated coefficients
        inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2])
      }
      
      quantile_model_fit <- evgam(maxtp ~ value, obs_data_for_quant_reg,
                                  family = "ald", inits = inits_coeff ,ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1],
                      beta_1 = quantile_model_fit$location$coefficients[2])
    }
    
    if (model_name == "quant_glob_anom"){
      
      # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
      if(q==1){
        inits_coeff = c(-2, 1, 3)
      }else {
        #to ensure smoothness in the estimated coefficients
        inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2],quantile_model_fit$location$coefficients[3])
      }
      
      quantile_model_fit <- fit_with_warning_handling(obs_data_for_quant_reg, zeta, inits_coeff)
      
      #evgam(maxtp ~ value + glob_anom, obs_data_for_quant_reg,
      #                          family = "ald", inits = inits_coeff ,ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1],
                      beta_1 = quantile_model_fit$location$coefficients[2],
                      beta_2 = quantile_model_fit$location$coefficients[3])
    }
    
    if (model_name == "quant_log_alt"){
      # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
      
      if(q==1){
        inits_coeff = c(1,1,1,1)
      }else {
        #to ensure smoothness in the estimated coefficients
        inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2],quantile_model_fit$location$coefficients[3])
      }
      
      
      quantile_model_fit <- evgam(maxtp ~ value + glob_anom + log(altitude), obs_data_for_quant_reg,
                                  family = "ald", inits = inits_coeff ,ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1],
                      beta_1 = quantile_model_fit$location$coefficients[2],
                      beta_2 = quantile_model_fit$location$coefficients[3],
                      beta_3 = quantile_model_fit$location$coefficients[4])
    }
    results_df = bind_rows(results_df, result) 
  }
  
  quant_reg_model_pars = results_df
  
  #defining temporal_covariates
  
  if (model_name %in% c("no_covariate", "quant")){
    temporal_covariates = obs_data %>%
      dplyr::select(year) %>% 
      unique() %>%
      arrange(year)
  }
  
  if (model_name %in% c("quant_glob_anom", "quant_log_alt")){
    temporal_covariates = obs_data %>%
      dplyr::select(year, glob_anom) %>% 
      unique() %>%
      arrange(year)
  }
  
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
      
      if (model_name == "quant_log_alt"){
        altitude_df = obs_data %>% filter(stn == .x$stn[1])%>%pull(altitude)
        altitude = altitude_df[[1]]
      }
      
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
                           quant_value = regression(qpars, q, clim_vals, model_name, temporal_covariates, altitude)
                    ))
      }
      
      
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
  
  return(obs_smoothed_quantiles)
}

estimate_RMSE <- function (fitting_quantiles, obs_data, model_name, quantiles_to_estimate, test_data){
  obs_smoothed_quantiles = estimate_parameters (fitting_quantiles, obs_data, model_name)
  
  #compute the quantiles of test_data here
  test_data_quantile <- test_data %>%
    select(-id) %>%
    mutate(year = year(date)) %>%
    select(-date) %>%
    group_by(stn, year) %>%
    summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))
  
  pred = extract_temperature(obs_smoothed_quantiles, quantiles_to_estimate) #predicting the quantiles
  
  rmse_val = total_rmse(pred, test_data_quantile)
  
  return(rmse_val)
}

extract_quantile = function(df, obs_smoothed_quantiles){
  
  result <- df %>%
    inner_join(obs_smoothed_quantiles, by = c("year", "stn")) %>%
    rowwise() %>%
    mutate(
      quantile_estimate = list(match.fun(temp_to_tau)(unlist(quantiles)))
    ) %>%
    select(stn, year, quantile_estimate)
  
  return(result)
  
}

estimate <- function (fitting_quantiles, obs_data, model_name, quantiles_to_estimate, test_data){
  obs_smoothed_quantiles = estimate_parameters (fitting_quantiles, obs_data, model_name)
  
  #compute the quantiles of test_data here
  test_data_quantile <- test_data %>%
    select(-id) %>%
    mutate(year = year(date)) %>%
    select(-date) %>%
    group_by(stn, year) %>%
    summarise(quantiles = list(as.numeric(quantile(maxtp, probs = quantiles_to_estimate, na.rm = TRUE))))%>%
    select(c(year, stn, quantiles))
  
  pred = extract_temperature(obs_smoothed_quantiles, quantiles_to_estimate) #predicting the quantiles
  
  pred_quantile = extract_quantile(test_data_quantile, obs_smoothed_quantiles)
  
  mae_val_tau_to_temp = total_mae(pred, test_data_quantile)
  
  rmse_val_tau_to_temp = total_rmse(pred, test_data_quantile)
  
  mae_val_temp_to_tau = total_mae_quantile(pred_quantile, quantiles_to_estimate)
  
  rmse_val_temp_to_tau = total_rmse_quantile(pred_quantile, quantiles_to_estimate)
  
  return (c(mae_val_tau_to_temp, rmse_val_tau_to_temp, mae_val_temp_to_tau, rmse_val_temp_to_tau)) #the higher the rmse/the mae, the worst the fit
}

