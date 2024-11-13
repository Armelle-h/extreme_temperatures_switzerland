#takes 


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
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date), week=week(date)) %>% #week is used for the temporal cross validation
  left_join(glob_anomaly_reshaped, by = "year")

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

obs_data = obs_data %>% 
  left_join(spatial_folds, by="stn")

obs_data$week <- as.numeric(obs_data$week)

week_chunks = list(c(22,25,28,31,34),
                   c(23,26,29,32,35),
                   c(24,27,30,33)) #22 to 35 correspond to summer weeks
                    

obs_data$temporal_fold = 123456789
obs_data[obs_data$week %in% week_chunks[[1]],]$temporal_fold=1
obs_data[obs_data$week %in% week_chunks[[2]],]$temporal_fold=2
obs_data[obs_data$week %in% week_chunks[[3]],]$temporal_fold=3

run_cv = function(cv_method, obs_data, fitting_quantiles, model_name, quantiles_to_estimate, num_spatial_folds = "NA", week_chunks="NA"){
  
  if(cv_method == "spatial-temporal"){
    
    for(s in seq(num_spatial_folds)){
      for(t in week_chunks){
        
        test_data = obs_data %>% 
          filter(spatial_fold == s) %>%
          filter(week %in% t)
        
        train_data = obs_data %>%
          anti_join(test_data)
        
        result = estimate(fitting_quantiles, train_data, model_name, quantiles_to_estimate, test_data)
        
        write.table(paste0(model_name, ",", result[1], ",", mean(t)), 
                    file = paste0("output/cv_bulk_models/mae_spatio_temporal_cv.csv"), 
                    sep=",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE     
                    )
        
        write.table(paste0(model_name, ",", result[2], ",", mean(t)), 
                    file = paste0("output/cv_bulk_models/rmse_spatio_temporal_cv.csv"), 
                    sep=",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE     
        )
      }
    }
  }
  
  if(cv_method == "10fold"){
    
    #cv_method = '10fold'
    # ---------- define random folds
    num_random_folds = 12 
    set.seed(1234567)
    obs_data$random_fold = sample(seq(num_random_folds), size = nrow(obs_data), replace = T)
    
    
    for(i in seq(num_random_folds)){
      test_data = obs_data %>% filter(random_fold == i)
      train_data = obs_data %>% filter(random_fold != i)
      
      result = estimate(fitting_quantiles, train_data, model_name, quantiles_to_estimate, test_data)
      
      write.table(paste0(model_name, ",", result[1]), 
                  file = paste0("output/cv_bulk_models/mae_twelve_fold_cv.csv"), 
                  sep=",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE     
      )
      
      write.table(paste0(model_name, ",", result[2]), 
                  file = paste0("output/cv_bulk_models/rmse_twelve_fold_cv.csv"), 
                  sep=",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE     
      )
      
    }
  }
}
  

job::job({run_cv("10fold", obs_data%>% select(-c(altitude, glob_anom, week, temporal_fold)), fitting_quantiles, "no_covariate", quantiles_to_estimate)}) # 1hour 12
job::job({run_cv("10fold", obs_data%>% select(-c(altitude, glob_anom, week, temporal_fold)), fitting_quantiles, "quant", quantiles_to_estimate)})  #1h29
job::job({run_cv("10fold", obs_data%>% select(-c(altitude, week, temporal_fold)), fitting_quantiles, "quant_glob_anom", quantiles_to_estimate)}) #1h32
job::job({run_cv("10fold", obs_data%>% select(-c(glob_anom, week, temporal_fold)), fitting_quantiles, "quant_log_alt", quantiles_to_estimate)}) #1h50
    #in total, the 4 jobs take 1h10 to run  --> expected 2h20

job::job({run_cv("spatial-temporal", obs_data%>% select(-c(altitude, glob_anom)), fitting_quantiles, "no_covariate", quantiles_to_estimate, num_spatial_folds, week_chunks)})
job::job({run_cv("spatial-temporal", obs_data%>% select(-c(altitude, glob_anom)), fitting_quantiles, "quant", quantiles_to_estimate, num_spatial_folds, week_chunks)})
job::job({run_cv("spatial-temporal", obs_data%>% select(-c(altitude)), fitting_quantiles, "quant_glob_anom", quantiles_to_estimate, num_spatial_folds, week_chunks)})
job::job({run_cv("spatial-temporal", obs_data%>% select(-c(glob_anom)), fitting_quantiles, "quant_log_alt", quantiles_to_estimate, num_spatial_folds, week_chunks)})


#computing the mean 

file_12fols_mae = read_csv("output/cv_bulk_models/mae_twelve_fold_cv.csv",
                       col_names=c("model_name", "mae"))

file_12fols_rmse = read_csv("output/cv_bulk_models/rmse_twelve_fold_cv.csv",
                       col_names=c("model_name", "rmse"))


mean_12fold_mae = file_12fols_mae %>%
  group_by(model_name) %>%
  summarise(mean_MAE = mean(mae, na.rm = TRUE))

mean_12fold_rmse = file_12fols_rmse %>%
  group_by(model_name) %>%
  summarise(mean_RMSE = mean(rmse, na.rm = TRUE))

mean_12_fold = mean_12fold_mae %>%
  left_join(mean_12fold_rmse, by="model_name") 




