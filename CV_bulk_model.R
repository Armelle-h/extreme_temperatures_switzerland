gc()
rm(list = ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('gpd_models.R')

library(tidyverse)
library(spatialsample)
library(scoringRules)
library(job)
library(sf)
library(rnaturalearth) #library for map of switzerland

# colour palatte 
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


# map of switzerland sf
switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")


# ========  ========  Global parameters ========

# ========  ========  Read in data  ========  ========
# Observational data with covariates

#new code
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
               names_to = "month", 
               values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

threshold_9_df = vroom::vroom("Data/Observed_data/1971_2023_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

obs_data = obs_data %>% 
  mutate(year = year(date), month = month(date), week=week(date)) %>% #week is used for the temporal cross validation
  left_join(glob_anomaly_reshaped, by = c("year", "month"))%>%
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

obs_sites = read.csv("Data/Observed_data/1971_2023_JJA_obs_legend.csv") %>%
  select(stn, longitude, latitude)

obs_sites_sf = st_as_sf(obs_sites, coords = c("longitude", "latitude"), crs = 4326)
#obs_sites_sf <- st_set_geometry(obs_sites_sf, "coord_pts")

#splits data based on coordinates
num_spatial_folds  = 30
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

get_metrics = function(orig_dat, ){ #will need to add other arguments
  
  # - RMSE
  x = orig_dat
  y = evd::qgpd(p = quant, loc = threshold, scale = pred_scale, shape = pred_shape[1]) #will need to be changed
  
  my_rmse = sqrt(mean((x-y)^2))
  
  return(my_rmse) #only computing the rmse
}



run_cv = function(cv_method, thresh_qnt, obs_data, spatial_folds, get_metrics, num_spatial_folds = "NA", week_chunks="NA"){
  
  if(cv_method == "spatial-temporal"){
    
    extreme_data = obs_data %>%
      mutate(excess = maxtp - threshold_9) %>%
      filter(excess > 0) %>%
      left_join(spatial_folds) %>%
      group_by(year, stn) %>%
      group_map(~{
        
        .x %>%
          mutate(quant = rank(excess)/(length(excess)+1)) #defining quant 
        
      },.keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    for(s in seq(num_spatial_folds)){
      for(t in week_chunks){
        #to complete
      }
    }
  }
  
  if(cv_method == "10fold"){
    
    extreme_data = obs_data %>%
      mutate(excess = maxtp - threshold_9) %>%
      filter(excess > 0) %>%
      group_by(year, stn) %>%
      group_map(~{
        
        .x %>%
          mutate(quant = rank(excess)/(length(excess)+1))
        
      },.keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    #cv_method = '10fold'
    # ---------- define random folds
    num_random_folds = 90
    set.seed(1234567)
    extreme_data$random_fold = sample(seq(num_random_folds), size = nrow(extreme_data), replace = T)
    
    
    for(i in seq(num_random_folds)){
      #to complete
    }
  }
}