library(tidyverse)
library(evd)

gpd_survival <- function(t, sigma, xi) {
  ifelse(
    xi == 0,
    exp(-t / sigma),
    (1 + xi * t / sigma)^(-1 / xi)
  )
}

# - - - - - - - - - - MODEL 1

# Function to compute the negative log-likelihood for Model 1
ngll_1 = function(par){
  
  # Check that the second parameter is non-negative
  #if(par[2] < 0) return(2^30) # Return a large penalty if invalid
  
  # Estimate scale parameter with additional climatic factor
  scale_est_uncensored = exp(par[1] + par[2]*clim_scale_uncensored + par[3]*glob_anom_uncensored)
  scale_est_censored = exp(par[1] + par[2]*clim_scale_censored + par[3]*glob_anom_censored)
  
  # Shape parameter
  shape_est_uncensored = par[4]
  shape_est_censored = par[4]
  
  # Check for valid scale estimates
  if(any(scale_est_uncensored <= 0)) return(2^30)
  if(any(scale_est_censored <= 0)) return(2^30)
  if(any(scale_est_uncensored > -1/shape_est_uncensored)) return(2^30)
  if(any(scale_est_censored > -1/shape_est_censored)) return(2^30)
  if(any((1+shape_est_uncensored*excess_dat_uncensored/scale_est_uncensored)< 0)) return(2^30)
  if(any((1+shape_est_censored*excess_dat_censored/scale_est_censored)< 0)) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat_uncensored, loc=0, scale = scale_est_uncensored, shape=shape_est_uncensored, log=T)) - sum(log(gpd_survival(excess_dat_censored, scale_est_censored, shape_est_censored)))
}

# Function to fit Model 1
fit_mod_1 = function(this_dat_uncensored, this_clim_scale_uncensored, this_glob_anom_uncensored, this_dat_censored, this_clim_scale_censored, this_glob_anom_censored, initial_pars = c(0.6, 0.1, 0.5, -0.1)){ #used to be c(0.3510713,  0.7598344, 0.3735851, -0.1429355)
  # Set the excess data and other variables globally
  excess_dat_uncensored <<- this_dat_uncensored
  clim_scale_uncensored <<- log(this_clim_scale_uncensored)
  glob_anom_uncensored <<- this_glob_anom_uncensored
  excess_dat_censored <<- this_dat_censored
  clim_scale_censored <<- log(this_clim_scale_censored)
  glob_anom_censored <<- this_glob_anom_censored
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_1)$par
}

# Function to predict scale and shape from fitted parameters for Model 1
my_predict_1 = function(estimates_pars, this_clim_scale, this_glob_anom){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*this_glob_anom),
         shape = estimates_pars[4]) # Return predictions as a tibble
}