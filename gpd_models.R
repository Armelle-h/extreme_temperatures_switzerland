library(tidyverse)
library(evd)

# - - - - - - - - - MODEL 0

#in the first part, the shape parameter has to be estimated, in the second, it is assumed to be known

# Function to compute the negative log-likelihood for Model 0
ngll_0 = function(par){
  # Estimate scale parameter using a log-linear model
  scale_est = exp(par[1] + par[2]*clim_scale)
  
  # Shape parameter from the parameter vector
  shape_est = par[3]  
  
  # Check if the scale estimates are valid (positive)
  if(any(scale_est <= 0)) return(2^30) # Return a large penalty if invalid
  
  # Calculate the log-likelihood using the generalized Pareto distribution
  sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 0
fit_mod_0 = function(this_dat, this_clim_scale, initial_pars = c(0.95, -0.3, -0.05)){ #-0.21 is the value of the climate scale, 21 is the value of beta0 in the no covariate quantile regression 
  # Set the excess data globally
  excess_dat <<- this_dat
  # Log-transform the climatic scale
  clim_scale <<- log(this_clim_scale) 
  # Optimize parameters using ngll_0
  optim(par = initial_pars, fn = ngll_0, control = list(fnscale = -1))$par
}

# Function to compute the negative log-likelihood for Model 0 with fixed shape
ngll_0_fix_shape = function(par){
  # Estimate scale parameter using a log-linear model
  scale_est = exp(par[1] + par[2]*clim_scale)
  
  # Check for valid scale estimates
  if(any(scale_est <= 0)) return(-2^30) # Return a large penalty if invalid
  if(any(scale_est > -1/shape_est)) return(-2^30) # Check against the shape
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(-2^30) # Ensure validity
  
  # Calculate the log-likelihood
  sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T)) 
}

# Function to fit Model 0 with fixed shape
fit_mod_0_fix_shape  = function(this_dat, this_clim_scale, this_shape_est, initial_pars = c(0.95, -0.3)){
  print(initial_pars) # Print initial parameters for debugging
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  shape_est <<- this_shape_est
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_0_fix_shape, control = list(fnscale = -1))$par
}

# Function to predict scale and shape from fitted parameters for Model 0
my_predict_0 = function(estimates_pars, this_clim_scale){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale)),
         shape = estimates_pars[3]) # Return predictions as a tibble
}

# Function to calculate return levels for Model 0
rl_mod_0 = function(estimates_pars, rl_quantile, thresh, this_clim_scale){
  # Estimate scale and shape
  estimated_scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale))
  estimated_shape = estimates_pars[3]
  # Calculate return level based on threshold
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
}

# - - - - - - - - - - MODEL 1

# Function to compute the negative log-likelihood for Model 1
ngll_1 = function(par){
  # Check that the second parameter is non-negative
  if(par[2] < 0) return(2^30) # Return a large penalty if invalid
  
  # Estimate scale parameter with additional climatic factor
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*glob_anom)
  
  # Shape parameter
  shape_est = par[4]
  
  # Check for valid scale estimates
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 1
fit_mod_1 = function(this_dat, this_clim_scale, this_glob_anom, initial_pars = c(0.3510713,  0.7598344, 0.3735851, -0.1429355)){
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale)
  glob_anom <<- this_glob_anom #setting it up as a general parameter so that it can be accessed by the inside of the function
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_1)$par
}

# Function to compute the negative log-likelihood for Model 1 with fixed shape
ngll_1_fix_shape = function(par){
  # Check that the second parameter is non-negative
  if(par[2] < 0) return(2^30) 
  
  # Estimate scale parameter with additional climatic factor
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*glob_anom)
  
  # Check for valid scale estimates
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 1 with fixed shape
fit_mod_1_fix_shape  = function(this_dat, this_clim_scale, this_glob_anom, this_shape_est, initial_pars = c(0.3510713,  0.7598344, 0.3735851)){
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  glob_anom <<- this_glob_anom
  shape_est <<- this_shape_est
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_1_fix_shape)$par
}

# Function to predict scale and shape from fitted parameters for Model 1
my_predict_1 = function(estimates_pars, this_clim_scale, this_glob_anom){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*this_glob_anom),
         shape = estimates_pars[4]) # Return predictions as a tibble
}

# Function to calculate return levels for Model 1
rl_mod_1 = function(estimates_pars, rl_quantile, thresh, this_clim_scale, this_glob_anom){
  # Estimate scale and shape
  estimated_scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*this_glob_anom)
  estimated_shape = estimates_pars[4]
  # Calculate return level based on threshold
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
}

# - - - - - - - - - - MODEL 2

# Function to compute the negative log-likelihood for Model 2
ngll_2 = function(par){
  # Estimate scale parameter with additional climatic and geographic factors
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*dist_sea + par[4]*loess_temp_anom + par[5]*dist_sea*loess_temp_anom) 
  shape_est = par[6]
  
  # Check for valid scale estimates
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 2
fit_mod_2 = function(this_dat, this_clim_scale, this_loess_temp_anom, this_dist_sea, initial_pars = c( 0.0713, 0.700, 0.0513, 0.123, 0.00117, -0.15)){
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale)
  loess_temp_anom <<- this_loess_temp_anom
  dist_sea <<- this_dist_sea
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_2)$par
}

# Function to predict scale and shape from fitted parameters for Model 2
my_predict_2 = function(estimates_pars, this_clim_scale, this_loess_temp_anom, this_dist_sea){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*this_dist_sea + estimates_pars[4]*this_loess_temp_anom + estimates_pars[5]*this_dist_sea*this_loess_temp_anom),
         shape = estimates_pars[6]) # Return predictions as a tibble
}

# Function to compute the negative log-likelihood for Model 2 with fixed shape
ngll_2_fix_shape = function(par){
  # Estimate scale parameter with additional climatic and geographic factors
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*dist_sea + par[4]*loess_temp_anom + par[5]*dist_sea*loess_temp_anom) 
  # Check for valid scale estimates
  if(any(scale_est <= 0)) return(2^30)
  if(any(scale_est > -1/shape_est)) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0)) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 2 with fixed shape
fit_mod_2_fix_shape  = function(this_dat, this_clim_scale, this_loess_temp_anom, this_dist_sea, this_shape_est, initial_pars = c( 0.0713, 0.700, 0.0513, 0.123, 0.00117)){
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  loess_temp_anom <<- this_loess_temp_anom
  dist_sea <<- this_dist_sea
  shape_est <<- this_shape_est
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_2_fix_shape)$par
}

# Function to calculate return levels for Model 2
rl_mod_2 = function(estimates_pars, rl_quantile, thresh, this_clim_scale, this_loess_temp_anom, this_dist_sea){
  # Estimate scale and shape
  estimated_scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*this_dist_sea + estimates_pars[4]*this_loess_temp_anom + estimates_pars[5]*this_dist_sea*this_loess_temp_anom)
  estimated_shape = estimates_pars[6]
  # Calculate return level based on threshold
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
}