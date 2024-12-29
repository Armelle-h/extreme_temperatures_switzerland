gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(data.table)
library(sf)
library(rnaturalearth) #library for map of switzerland

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


my_pal = c( 
  '#062c30', # extra dark 
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)

ngll = function(par){
  if(par <= 0) return(2^30)
  if(par > -1/shape_param) return(2^30)
  if(any((1+shape_param*this.dat/par)< 0)) return(2^30)
  
  -sum(evd::dgpd(x = this.dat, loc=0, scale = par, shape=shape_param, log=T))
}

estimate_scale_fixed_shape = function(x,shape_c, initial_pars = c(0)){
  this.dat <<- x
  shape_param <<- shape_c
  optim(par = initial_pars, fn = ngll, method = 'Brent', lower=0, upper = 5)
}


ngll_free_xi = function(par){
  if(par[1] <= 0) return(2^30)
  if(par[1] > -1/par[2]) return(2^30)
  if(any((1+par[2]*this.dat/par[1])< 0)) return(2^30)
  
  -sum(evd::dgpd(x = this.dat, loc=0, scale = par[1], shape=par[2], log=T))
}

estimate_scale_free_shape = function(x, initial_pars = c(1.5,-0.15)){
  this.dat <<- x
  optim(initial_pars, fn = ngll_free_xi)
}



# 3a. ---- Climate model

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]]) %>%
    filter(id %% 5 == 0)
  
  clim_data$maxtp = as.integer(clim_data$maxtp)
  
  sing_clim_thresh_values = clim_data %>%
    group_by(id) %>%
    mutate(threshold = quantile(maxtp, 0.9),
           excess = maxtp - threshold) %>%
    filter(excess > 0) %>%
    ungroup()
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_thresh_values = do.call(rbind, clim_thresh_values_list)

rm(clim_thresh_values_list)
gc()

clim_data_extreme_9 = clim_thresh_values

optimal_shape_9 = -0.2019049 #has been changed


grid = read_csv("Data/Climate_data/clim_scale_grid_gpd_model.csv") %>%
  filter(id %% 5 == 0)

scales_fixed_xi = c()
loglik_fixed_xi = c()

scales_free_xi = c()
loglik_free_xi = c()

index = 1

for(i in grid$id){
  
  if (index %% 500 == 0){
    print(index)
  }
  
  index = index + 1
  
  this_clim_extm_irel = clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)
  
  # --- fixed shape
  model_fit = estimate_scale_fixed_shape(this_clim_extm_irel, optimal_shape_9, c(0))  #could change to c(2.4)
  scales_fixed_xi = c(scales_fixed_xi, model_fit$par)
  loglik_fixed_xi = c(loglik_fixed_xi, model_fit$value)
  
  # --- free shape
  model_fit = estimate_scale_free_shape(this_clim_extm_irel, c(2.5,-0.3)) #c(1.5,-0.15)   could be c(2.8, -0.34)
  scales_free_xi = c(scales_free_xi, model_fit$par)
  loglik_free_xi = c(loglik_free_xi, (evd::fpot(this_clim_extm_irel, threshold = 0, std.err=FALSE) %>% logLik %>% as.numeric()))
}


grid$LL_null = loglik_fixed_xi
grid$LL_alt = -loglik_free_xi
grid$LL_ratio = -2*log(grid$LL_alt/grid$LL_null)

id_lon_lat = read.csv("Data/id_lon_lat_correspondance.csv")

grid_lonlat = grid %>%
  left_join(id_lon_lat, by = "id") 


plt_clim = grid_lonlat %>%  #the weird position of the points comes from the fact that we're selecting only every 10 points to be plot.
  ggplot()+
  geom_point(aes(longitude, latitude, col = LL_ratio))+
  coord_map()+
  theme_minimal()+
  geom_sf(data = switzerland, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks= -c(10, 8, 6))+
  theme_minimal(12)+
  labs(col = '2ln(LR)', x = '', y = '')


obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv") 

threshold_9_df = vroom::vroom("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, stn)%>%
  unique()

obs_dat = obs_data %>% 
  left_join(threshold_9_df, by="stn")%>%
  mutate(excess = maxtp - threshold_9) %>%
  filter(excess>0)

loglik_sum = c()

#investigate which range seems good before running the whole thing !!!!

for(potential_shape in seq(-0.26, -0.21, length.out = 100)){ 
  print(potential_shape)
  
  loglik = c()
  
  for(i in (obs_dat$stn %>% unique())){
    this_extm = obs_dat %>% filter(stn == i) %>% pull(excess)
    
    model_fit = estimate_scale_fixed_shape(this_extm, potential_shape)
    
    loglik = c(loglik, model_fit$value)
    
  }
  loglik_sum = c(loglik_sum, sum(loglik))
}

print(sort(loglik_sum))

optimal_shape = seq(-0.26, -0.21, length.out = 100)[which.min(loglik_sum)]

loglik_fixed_xi = c()
loglik_free_xi = c()
std_errors_scale = c()
std_errors_shape = c()

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

obs_grid = obs_dat %>% left_join(legend_data, by = "stn") %>% dplyr::select(longitude, latitude, stn) %>% unique

for(i in obs_grid$stn){
  print(i)
  this_extm_irel = obs_dat %>% filter(stn == i) %>% pull(excess)
  
  # --- fixed shape
  model_fit = estimate_scale_fixed_shape(this_extm_irel, optimal_shape) 
  loglik_fixed_xi = c(loglik_fixed_xi, model_fit$value)
  
  if (i %in% c("WSLLEB", "NAS", "ROM", "MAE", "INNEBI", "KAGOD", "PERGEM", "GEN")){
    # --- free shape
    loglik_free_xi = c(loglik_free_xi, (evd::fpot(this_extm_irel, threshold = 0, std.err = FALSE) %>% logLik %>% as.numeric()))
    std_errors_scale = c(std_errors_scale, NA)
    std_errors_shape = c(std_errors_shape, NA)
  }
  else{
    # --- free shape
    model = evd::fpot(this_extm_irel, threshold = 0)
    loglik_free_xi = c(loglik_free_xi, model %>% logLik %>% as.numeric())
    std_errors_scale = c(std_errors_scale, model$std.err[[1]])
    std_errors_shape = c(std_errors_shape, model$std.err[[2]])
  }
}

obs_grid$LL_null = loglik_fixed_xi
obs_grid$LL_alt = -loglik_free_xi
obs_grid$LL_ratio = -2*log(obs_grid$LL_alt/obs_grid$LL_null)
obs_grid$LL_alt_std_error_scale = std_errors_scale
obs_grid$LL_alt_std_error_shape = std_errors_shape


obs_grid %>%
  ggplot()+
  geom_point(aes(longitude, latitude, col = LL_ratio), alpha = 0.5, size = 2.5)+ #size = 2.5
  coord_map()+
  theme_minimal()+
  geom_sf(data = switzerland, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks= -c(10, 8, 6))+
  theme_minimal(12)+
  labs(col = '2ln(LR)', x = '', y = '')


plt = gridExtra::grid.arrange(plt_clim, plt_obs, nrow = 1, widths = c(1, 0.95))

