#in this file, I describe all the functions needed do to a cross validation on bulk models 

#I need functions that return the quantiles estimated by the wanted bulk model 


#no covariate bulk model 

library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#putting obs_data and obs_data_quantile as arguments
extract_quantile <- function(df, quantiles_to_estimate){
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
  
  # Sum the RMSE over all years and stations
  total_rmse_df <- results %>%
    summarise(total_rmse = sqrt(sum(rmse^2, na.rm = TRUE) / sum(!is.na(rmse)))) 
  
  total_rmse_value <- total_rmse_df$total_rmse
  
  return(total_rmse_value)
}

regression <- function(qpars, q, clim_vals, model_name, temporal_covariates){
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
    result = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2) *log(temporal_covariates$altitude)
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
    print(paste0("Fitting  tau = ", zeta))
    
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
      quantile_model_fit <- evgam(maxtp ~ value, obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1],
                      beta_1 = quantile_model_fit$location$coefficients[2])
    }
    
    if (model_name == "quant_glob_anom"){
      # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
      quantile_model_fit <- evgam(maxtp ~ value + glob_anom, obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1],
                      beta_1 = quantile_model_fit$location$coefficients[2],
                      beta_2 = quantile_model_fit$location$coefficients[3])
    }
    
    if (model_name == "quant_log_alt"){
      # Fit the quantile regression model using EVGAM package with asymmetric Laplace distribution
      quantile_model_fit <- evgam(maxtp ~ value + log(altitude), obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(tau = zeta,
                      beta_0 = quantile_model_fit$location$coefficients[1],
                      beta_1 = quantile_model_fit$location$coefficients[2],
                      beta_2 = quantile_model_fit$location$coefficients[3])
    }

    results_df <- bind_rows(results_df, result) 
  }
  return(results_df)
}

estimate_RMSE <- function (fitting_quantiles, obs_data, model_name, quantiles_to_estimate, obs_data_quantile, temporal_covariates){
  
  quant_reg_model_pars = estimate_parameters(fitting_quantiles, obs_data, model_name)
  
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
      for(q in seq_along(fitting_quantiles)){
        qpars = quant_reg_pars[q,]
        # Calculate the predicted quantile value using the regression parameters
        res = rbind(res,
                    tibble(quantile =  qpars$tau,
                           year = temporal_covariates$year,
                           quant_value = regression(qpars, q, clim_vals, model_name, temporal_covariates)
                           ))
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
  
  pred = extract_quantile(obs_smoothed_quantiles, quantiles_to_estimate)
  rmse_val = total_rmse(pred, obs_data_quantile)
  
  return(rmse_val)
}

