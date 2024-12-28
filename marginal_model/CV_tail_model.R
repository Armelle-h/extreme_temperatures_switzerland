# Description: This script provides comparison between potential marginal GPD
# models using bootstrapped AIC,BIC and log likelihood. Bootstraps are 
# constructed to preserve spatial and temporal dependence.

#takes 3 minutes to run

#used to be done for 90 folds, now done for 15 folds 


gc()
rm(list = ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')

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



run_cv = function(cv_method, thresh_qnt, obs_data, spatial_folds, get_metrics, num_spatial_folds = "NA", week_chunks="NA"){
  
  source('marginal_model/gpd_models.R')
  
  if(cv_method == "spatial-temporal"){
    
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
        
        test = extreme_data %>% 
          filter(spatial_fold == s) %>%
          filter(week %in% t)
        
        train = extreme_data %>%
          anti_join(test)
        
        this_fit_mod_0 = fit_mod_0(train$excess, train$scale_9, c(0.8, -0.05, -0.1))
        this_fit_mod_1 = fit_mod_1(train$excess, train$scale_9, train$glob_anom, c(0.8, -0.05, 0.001, -0.1))
        this_fit_mod_2 = fit_mod_2(train$excess, train$scale_9, train$glob_anom, train$altitude, c(0.03, -0.03,  0.08,  1, -0.08, -0.2))
        this_fit_mod_3 = fit_mod_3(train$excess, train$scale_9, train$altitude, c(0.5, 0.01, 0.01, -0.1))
        
        # calculate scale parameter on climate grid
        pred_0 = my_predict_0(this_fit_mod_0, test$scale_9)
        pred_1 = my_predict_1(this_fit_mod_1, test$scale_9, test$glob_anom)
        pred_2 = my_predict_2(this_fit_mod_2, test$scale_9, test$glob_anom, test$altitude)
        pred_3 = my_predict_3(this_fit_mod_3, test$scale_9, test$altitude)
        
        write.table(paste0("model.0", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_0$scale, pred_0$shape[1]),",",mean(t)),
                    file = paste0('output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
        
        write.table(paste0("model.1", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_1$scale, pred_1$shape[1]),",",mean(t)),
                    file = paste0('output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
        
        write.table(paste0("model.2", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_2$scale, pred_2$shape[1]),",",mean(t)),
                    file = paste0('output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
        
        write.table(paste0("model.3", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_3$scale, pred_3$shape[1]),",",mean(t)),
                    file = paste0('output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
      }
    }
  }
  
  if(cv_method == "12fold"){
    
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
    
    
    #cv_method = '12fold'
    # ---------- define random folds
    num_random_folds = 12
    set.seed(1234567)
    extreme_data$random_fold = sample(seq(num_random_folds), size = nrow(extreme_data), replace = T)
    
    
    for(i in seq(num_random_folds)){
      test = extreme_data %>% filter(random_fold == i)
      train = extreme_data %>% filter(random_fold != i)
      
      this_fit_mod_0 = fit_mod_0(train$excess, train$scale_9)
      this_fit_mod_1 = fit_mod_1(train$excess, train$scale_9, train$glob_anom)
      this_fit_mod_2 = fit_mod_2(train$excess, train$scale_9, train$glob_anom, train$altitude)
      this_fit_mod_3 = fit_mod_3(train$excess, train$scale_9, train$altitude)
      
      # calculate scale parameter on climate grid
      pred_0 = my_predict_0(this_fit_mod_0, test$scale_9)
      pred_1 = my_predict_1(this_fit_mod_1, test$scale_9, test$glob_anom)
      pred_2 = my_predict_2(this_fit_mod_2, test$scale_9, test$glob_anom, test$altitude)
      pred_3 = my_predict_3(this_fit_mod_3, test$scale_9, test$altitude)
      
      write.table(paste0("model.0", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_0$scale, pred_0$shape[1])),
                  file = paste0('output/cv_res/twelve_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    
      write.table(paste0("model.1", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_1$scale, pred_1$shape[1])),
                  file = paste0('output/cv_res/twelve_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)

      write.table(paste0("model.2", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_2$scale, pred_2$shape[1])),
                  file = paste0('output/cv_res/twelve_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
      
      write.table(paste0("model.3", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_3$scale, pred_3$shape[1])),
                  file = paste0('output/cv_res/twelve_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
      
    }
  }
}


job::job({run_cv('12fold', 0.9, obs_data, spatial_folds, get_metrics)}, import = c("obs_data", 'run_cv', "spatial_folds", "get_metrics"))
job::job({run_cv('spatial-temporal', 0.9, obs_data, spatial_folds, get_metrics, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', "spatial_folds",  "get_metrics", "num_spatial_folds", "week_chunks"))


spatio_temporal_metrics = read.csv("output/cv_res/spatio_temporal_cv_test_mod_n9.csv", header=FALSE)
colnames(spatio_temporal_metrics) = c("model", "sum_ll", "stand_sum_ll", "rmse", "crps", "avg_week")

final_metrics = spatio_temporal_metrics %>%
  select(-avg_week) %>%
  group_by(model) %>%
  summarise(
    sum_ll_mean = mean(sum_ll, na.rm = TRUE),
    stand_sum_ll_mean = mean(stand_sum_ll, na.rm = TRUE),
    rmse_mean = mean(rmse, na.rm = TRUE),
    crps_mean = mean(crps, na.rm = TRUE)
  )

twelvefold_metrics = read.csv("output/cv_res/twelve_fold_cv_test_mod_n9.csv", header=FALSE)
colnames(twelvefold_metrics) = c("model", "sum_ll", "stand_sum_ll", "rmse", "crps")

final_metrics_twelvefold = twelvefold_metrics %>%
  group_by(model) %>%
  summarise(
    sum_ll_mean = mean(sum_ll, na.rm = TRUE),
    stand_sum_ll_mean = mean(stand_sum_ll, na.rm = TRUE),
    rmse_mean = mean(rmse, na.rm = TRUE),
    crps_mean = mean(crps, na.rm = TRUE)
  )

