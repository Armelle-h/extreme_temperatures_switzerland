gc()
rm(list = ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

library(tidyverse)
library(spatialsample)
library(scoringRules)
library(job)
library(sf)
library(rnaturalearth) #library for map of switzerland

ngll_2 = function(par){
  # Estimate scale parameter with additional climatic and geographic factors
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*altitude + par[4]*glob_anom + par[5]*altitude*glob_anom) 
  shape_est = par[6]
  
  # Check for valid scale estimates
  if(any(scale_est <= 0) || any(is.na(scale_est))) return(2^30)
  if(any(scale_est > -1/shape_est) || any(is.na(scale_est))) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0) || any(is.na(scale_est))) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 2
fit_mod_2 = function(this_dat, this_clim_scale, this_glob_anom, this_altitude, initial_pars = c( 0.0713, 0.700, 0.0513, 0.123, 0.00117, -0.15)){
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale)
  glob_anom <<- this_glob_anom
  altitude <<- log(this_altitude)
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_2)$par
}

# Function to predict scale and shape from fitted parameters for Model 2
my_predict_2 = function(estimates_pars, this_clim_scale, this_glob_anom, this_altitude){
  tibble(scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*log(this_altitude) + estimates_pars[4]*this_glob_anom + estimates_pars[5]*log(this_altitude)*this_glob_anom),
         shape = estimates_pars[6]) # Return predictions as a tibble
}

# Function to compute the negative log-likelihood for Model 2 with fixed shape
ngll_2_fix_shape = function(par){ #altitude is set globally to log(this altitude)-> no need to add the log
  # Estimate scale parameter with additional climatic and geographic factors
  scale_est = exp(par[1] + par[2]*clim_scale + par[3]*altitude + par[4]*glob_anom + par[5]*altitude*glob_anom) 
  # Check for valid scale estimates
  if(any(scale_est <= 0) || any(is.na(scale_est))) return(2^30)
  if(any(scale_est > -1/shape_est) || any(is.na(scale_est))) return(2^30)
  if(any((1+shape_est*excess_dat/scale_est)< 0) || any(is.na(scale_est))) return(2^30)
  
  # Calculate the log-likelihood
  -sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}

# Function to fit Model 2 with fixed shape
fit_mod_2_fix_shape  = function(this_dat, this_clim_scale, this_glob_anom, this_altitude, this_shape_est, initial_pars = c( 0.0713, 0.700, 0.0513, 0.123, 0.00117)){
  # Set the excess data and other variables globally
  excess_dat <<- this_dat
  clim_scale <<- log(this_clim_scale) 
  glob_anom <<- this_glob_anom
  altitude <<- log(this_altitude)
  shape_est <<- this_shape_est
  # Optimize parameters
  optim(par = initial_pars, fn = ngll_2_fix_shape)$par
}

# Function to calculate return levels for Model 2
rl_mod_2 = function(estimates_pars, rl_quantile, thresh, this_clim_scale, this_glob_anom, this_altitude){
  # Estimate scale and shape
  estimated_scale = exp(estimates_pars[1] + estimates_pars[2]*log(this_clim_scale) + estimates_pars[3]*log(this_altitude) + estimates_pars[4]*this_glob_anom + estimates_pars[5]*log(this_altitude)*this_glob_anom)
  estimated_shape = estimates_pars[6]
  # Calculate return level based on threshold
  return(thresh + estimated_scale * ((1-rl_quantile)^(-estimated_shape) - 1)/estimated_shape)
}


obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 
obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

threshold_9_df = vroom::vroom("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

obs_data = obs_data %>% 
  mutate(year = year(date), week=week(date)) %>% #week is used for the temporal cross validation
  left_join(glob_anomaly_reshaped, by = "year")%>%
  left_join(threshold_9_df, by="id")

# onlykeep full years for cv
obs_data <- obs_data %>%
  group_by(year, stn) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n == 92) %>% 
  left_join(obs_data) %>%
  ungroup()


# -------- CREATE CV Folds
set.seed(51964)

obs_sites = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv") %>%
  select(stn, longitude, latitude)

obs_sites_sf = st_as_sf(obs_sites, coords = c("longitude", "latitude"), crs = 4326)
#obs_sites_sf <- st_set_geometry(obs_sites_sf, "coord_pts")

#splits data based on coordinates
num_spatial_folds  = 5
clustered = spatial_clustering_cv(obs_sites_sf, v = num_spatial_folds)

spatial_folds = c()
for(i in seq(num_spatial_folds)){
  spatial_folds = rbind(spatial_folds, 
                        assessment(clustered$splits[i][[1]]) %>% 
                          mutate(longitude = st_coordinates(geometry)[, 1],
                                 latitude = st_coordinates(geometry)[, 2]) %>%
                          dplyr::select(stn, longitude, latitude) %>% 
                          st_drop_geometry %>% 
                          unique %>%
                          mutate(spatial_fold = i))
}

week_chunks = list(c(22,25,28,31,34),
                   c(23,26,29,32,35),
                   c(24,27,30,33)) #22 to 35 correspond to summer weeks

obs_data$temporal_fold = 123456789
obs_data[obs_data$week %in% week_chunks[[1]],]$temporal_fold=1
obs_data[obs_data$week %in% week_chunks[[2]],]$temporal_fold=2
obs_data[obs_data$week %in% week_chunks[[3]],]$temporal_fold=3

get_metrics = function(orig_dat, excess, quant, threshold,  pred_scale, pred_shape){
  likelihood_vals = evd::dgpd(x = excess, loc = 0, scale = pred_scale, shape = pred_shape, log = T)
  ll = sum(likelihood_vals[!is.infinite(likelihood_vals)])
  ll_standardised = sum(likelihood_vals[!is.infinite(likelihood_vals)])/length(likelihood_vals[!is.infinite(likelihood_vals)])
  
  # - RMSE
  x = orig_dat
  y = evd::qgpd(p = quant, loc = threshold, scale = pred_scale, shape = pred_shape[1])
  
  my_rmse = sqrt(mean((x-y)^2))
  scoring=0
  scoring = crps(y = excess, family = "gpd", location = 0, scale = pred_scale, shape = pred_shape[1], mass = 0) %>% mean
  return(paste0(ll, ",", ll_standardised, ",",my_rmse, ",", scoring))
}


debug = function (thresh_qnt, obs_data, spatial_folds, get_metrics, num_spatial_folds = "NA", week_chunks="NA"){
  
  extreme_data = obs_data %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0) %>%
    left_join(spatial_folds, by="stn") %>%
    group_by(year, stn) %>%
    group_map(~{
      
      .x %>%
        mutate(quant = rank(excess)/(length(excess)+1)) #defining quant 
      
    },.keep=T) %>%
    plyr::rbind.fill() %>%
    as_tibble()
  
  
  for(s in seq(num_spatial_folds)){
    for(t in week_chunks){
      
      cat(s, "and", t)
      
      test = extreme_data %>% 
        filter(spatial_fold == s) %>%
        filter(week %in% t)
      
      train = extreme_data %>%
        anti_join(test)
    
      this_fit_mod_2 = fit_mod_2(train$excess, train$scale_9, train$glob_anom, train$altitude)
      
      # calculate scale parameter on climate grid
      pred_2 = my_predict_2(this_fit_mod_2, test$scale_9, test$glob_anom, test$altitude)
      
      get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_2$scale, pred_2$shape[1])
    }
  }
}

debug(0.9, obs_data, spatial_folds, get_metrics, num_spatial_folds, week_chunks)
