setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')

library(tidyverse)
library(spatialsample)
library(scoringRules)
library(job)
library(sf)

# ========  ========  Global parameters ========

# ========  ========  Read in data  ========  ========
# Observational data with covariates

#new code
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
num_spatial_folds  = 4
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


list_par = list(c(0.8, -0.05, 0.001, -0.1), c(0.5, 0.01, 0.01, -0.1), c(0.1, 0.01, 0.1, -0.1))

for ( par in list_par){
  par = c(0.5, -0.01, 0.04, -0.2)
  opt_par = fit_mod_3(extreme_data$excess, extreme_data$scale_9, extreme_data$altitude, par)
  print(opt_par)
  print(ngll_3(opt_par))
}