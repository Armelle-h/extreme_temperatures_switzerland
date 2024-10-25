#run spatial temporal !!!! Should take 1 hour 10.


gc()
rm(list = ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('bulk_models_function_cv.R')

library(tidyverse)
library(spatialsample)
library(scoringRules)
library(job)
library(sf)
# ========  ========  Global parameters ========

# ========  ========  Read in data  ========  ========
# Observational data with covariates

#new code
num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

#removing clim_thresh_value_9 and threshold_9 as they are unused in this code 

obs_data <- obs_data %>% 
  select(-c(clim_thresh_value_9, threshold_9))

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
               names_to = "month", 
               values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

obs_data = obs_data %>% 
  mutate(year = year(date), month = month(date), week=week(date)) %>% #week is used for the temporal cross validation
  left_join(glob_anomaly_reshaped, by = c("year", "month"))

# onlykeep full years for cv
obs_data <- obs_data %>%
  group_by(year, stn) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n == 92) %>% 
  left_join(obs_data) %>%
  ungroup()

#remove column n which contains only 92

obs_data <- obs_data %>% 
  select(-n)

num_quantiles_to_estimate = 50
quantiles_to_estimate = seq(0.001,0.99,length.out = num_quantiles_to_estimate)

num_quantiles_to_fit = 30
fitting_quantiles = seq(0.001,0.99,length.out = num_quantiles_to_fit)


# -------- CREATE CV Folds
set.seed(51964)

obs_sites = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv") %>%
  select(stn, longitude, latitude)

obs_sites_sf = st_as_sf(obs_sites, coords = c("longitude", "latitude"), crs = 4326)
#obs_sites_sf <- st_set_geometry(obs_sites_sf, "coord_pts")

#splits data based on coordinates
num_spatial_folds  = 3
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

obs_data = obs_data %>% 
  left_join(spatial_folds, by="stn")

obs_data$week <- as.numeric(obs_data$week)

week_chunks = list(c(22,24,26,28,30,32,34),
                   c(23,25,27,29,31,33,35)) #22 to 35 correspond to summer weeks
                    

obs_data$temporal_fold = 123456789
obs_data[obs_data$week %in% week_chunks[[1]],]$temporal_fold=1
obs_data[obs_data$week %in% week_chunks[[2]],]$temporal_fold=2

run_cv = function(cv_method, obs_data, fitting_quantiles, model_name, quantiles_to_estimate, num_spatial_folds = "NA", week_chunks="NA"){
  
  if(cv_method == "spatial-temporal"){
    
    for(s in seq(num_spatial_folds)){
      for(t in week_chunks){
        
        test_data = obs_data %>% 
          filter(spatial_fold == s) %>%
          filter(week %in% t)
        
        train_data = obs_data %>%
          anti_join(test_data)
        
        rmse_result = estimate_RMSE(fitting_quantiles, train_data, model_name, quantiles_to_estimate, test_data)
        
        write.table(paste0(model_name, ",", rmse_result, ",", mean(t)), 
                    file = paste0("output/cv_bulk_models/spatio_temporal_cv.csv"), 
                    sep=",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE     
                    )
      }
    }
  }
  
  if(cv_method == "10fold"){
    
    #cv_method = '10fold'
    # ---------- define random folds
    num_random_folds = 6
    set.seed(1234567)
    obs_data$random_fold = sample(seq(num_random_folds), size = nrow(obs_data), replace = T)
    
    
    for(i in seq(num_random_folds)){
      test_data = obs_data %>% filter(random_fold == i)
      train_data = obs_data %>% filter(random_fold != i)
      
      rmse_result = estimate_RMSE(fitting_quantiles, train_data, model_name, quantiles_to_estimate, test_data)
      
      write.table(paste0(model_name, ",", rmse_result), 
                  file = paste0("output/cv_bulk_models/six_fold_cv.csv"), 
                  sep=",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE     
      )
      
    }
  }
}

for(cv_method in c("10fold", "spatial-temporal")){ #not beautifully written but the idea is that I'm lauching 4 jobs 
                                                   #and then when they are finished the 4 next
  
  if(cv_method == "10fold" ){
    job::job({run_cv("10fold", obs_data%>% select(-c(altitude, glob_anom, week, temporal_fold)), fitting_quantiles, "no_covariate", quantiles_to_estimate)}) #takes 41 minutes
    job::job({run_cv("10fold", obs_data%>% select(-c(altitude, glob_anom, week, temporal_fold)), fitting_quantiles, "quant", quantiles_to_estimate)}) #takes 51 minutes
    job::job({run_cv("10fold", obs_data%>% select(-c(altitude, week, temporal_fold)), fitting_quantiles, "quant_glob_anom", quantiles_to_estimate)}) #takes 53 minutes
    job::job({run_cv("10fold", obs_data%>% select(-c(glob_anom, week, temporal_fold)), fitting_quantiles, "quant_log_alt", quantiles_to_estimate)}) #takes 1 hour 11 minutes
    #in total, the 4 jobs take 1h10 to run
  }
  else{
    job::job({run_cv("spatial-temporal", obs_data%>% select(-c(altitude, glob_anom)), fitting_quantiles, "no_covariate", quantiles_to_estimate, num_spatial_folds, week_chunks)})
    job::job({run_cv("spatial-temporal", obs_data%>% select(-c(altitude, glob_anom)), fitting_quantiles, "quant", quantiles_to_estimate, num_spatial_folds, week_chunks)})
    job::job({run_cv("spatial-temporal", obs_data%>% select(-c(altitude)), fitting_quantiles, "quant_glob_anom", quantiles_to_estimate, num_spatial_folds, week_chunks)})
    job::job({run_cv("spatial-temporal", obs_data%>% select(-c(glob_anom)), fitting_quantiles, "quant_log_alt", quantiles_to_estimate, num_spatial_folds, week_chunks)})
  }

}