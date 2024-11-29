#in this file, we investigate the choice of the threshold for the climate model and the observed data model 
#by using the property that scale(u) = scale(u_0)+ epsilon(u-u_0) where espilon is the shape parameter, u and u_0 are two threshold values


# I need scale_85, scale_9, scale_95 , threshold_85, threshold_9, threshold_95 and depending on the model chosen, glob anom, altitude


gc()
rm(list = ls())

library(tidyverse)
library(data.table)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#should be app. the same
shape_85 = readRDS("optimal_shape_85.rds")
shape_9 = readRDS("optimal_shape.rds")
shape_95 = readRDS("optimal_shape_95.rds")

print(shape_85)
print(shape_9)
print(shape_95)

shape = (shape_85 + shape_9 + shape_95)/3

#using clim_thresh_ as threshold 

#contains id and scale_9
scale_data_85 = read_csv("Data/Climate_data/clim_scale_grid_85_gpd_model.csv")
scale_data_9 = read_csv("Data/Climate_data/clim_scale_grid_gpd_model.csv")
scale_data_95 = read_csv("Data/Climate_data/clim_scale_grid_95_gpd_model.csv")

#computing the thresholds 

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_thresh_values_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])%>%
    filter(id %%1000 == 0)
  
  clim_data$maxtp = as.integer(clim_data$maxtp)
  
  sing_clim_thresh_values = clim_data %>%
    group_by(id) %>%
    mutate(
    clim_thresh_value_85 = quantile(maxtp, 0.85, na.rm = TRUE),
    clim_thresh_value_9 = quantile(maxtp, 0.9, na.rm = TRUE),
    clim_thresh_value_95 = quantile(maxtp, 0.95, na.rm = TRUE)
    )
  
  clim_thresh_values_list[[i]] = sing_clim_thresh_values
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_data = do.call(rbind, clim_thresh_values_list)%>%
  select(c("id", "clim_thresh_value_85", "clim_thresh_value_9", "clim_thresh_value_95"))%>%
  unique()

rm(clim_thresh_values_list)

clim_data = clim_data %>%
  left_join(scale_data_85, by = "id")%>%
  left_join(scale_data_9, by = "id")%>%
  left_join(scale_data_95, by = "id")

#creating 3 columns, 95-90, 90-85 and 95-85

linearity = function(scale_u, scale_u_0, u, u_0){
  result = abs(scale_u - (scale_u_0+shape*(u-u_0)))
  return(result)
}

clim_data = clim_data %>%
  rowwise() %>%
  mutate(
    res_95_90 = linearity(scale_95, scale_9, clim_thresh_value_95, clim_thresh_value_9),
    res_90_85 = linearity(scale_9, scale_85, clim_thresh_value_9, clim_thresh_value_85),
    res_95_85 = linearity(scale_95, scale_85, clim_thresh_value_95, clim_thresh_value_85)
  ) %>%
  ungroup()
  










#loading the data for the observed model --- for now this technique wouldn't work as I do't have the same shape parameters for the different thresholds

num_quantiles = 30

#for the threshold
obs_data_85 = readRDS(paste0("Data/processed/obs_data_85_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  select(c(stn, year, threshold_85))%>%
  unique()

obs_data_9 = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  select(c(stn, year, threshold_9))%>%
  unique()

obs_data_95 = readRDS(paste0("Data/processed/obs_data_95_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  select(c(stn, year, threshold_95))%>%
  unique()


#for the scale
scale_data_85 = vroom::vroom("Data/Observed_data/obs_data_85_gpd_model.csv")%>%
  mutate(year = lubridate::year(date))%>%
  select(c(stn, year, scale_85))%>%
  unique()

scale_data_9 = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")%>%
  mutate(year = lubridate::year(date))%>%
  select(c(stn, year, scale_9))%>%
  unique()

scale_data_95 = vroom::vroom("Data/Observed_data/obs_data_95_gpd_model.csv")%>%
  mutate(year = lubridate::year(date))%>%
  select(c(stn, year, scale_95))%>%
  unique()

#glob anom
glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

#altitude

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 
obs_data = obs_data_85 %>%
  left_join (obs_data_9, by=c("stn", "year"))%>%
  left_join (obs_data_95, by=c("stn", "year"))%>%
  left_join (scale_data_85, by=c("stn", "year"))%>%
  left_join (scale_data_9, by=c("stn", "year"))%>%
  left_join (scale_data_95, by=c("stn", "year"))%>%
  left_join (glob_anomaly_reshaped, by="year")%>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

#testing for model 0


