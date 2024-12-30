#In this file, we plot the 100 year return level 

gc()
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')
library(tidyverse)
library(sf)
library(rnaturalearth)

# map for plotting
switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)

# colours for plotting
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

marg_mod = "mod_2"
num_quantiles = 30

id_lon_lat = read.csv("Data/id_lon_lat_correspondance.csv")

clim_date_w_quantile_mod = readRDS(paste0("output/glob_anomaly_quant_models_clim_num_quantiles_",num_quantiles,".csv")) %>%
  dplyr::select(-c(quantile, value)) %>% filter(year %in% c(1971, 2022)) #the filter by years not realy needed, already filtered

grid_pred = clim_date_w_quantile_mod %>%
  left_join(id_lon_lat, by = "id") %>%
  dplyr::select(id, longitude, latitude, year, threshold_9, thresh_exceedance_9, glob_anom) %>%
  left_join(read_csv('Data/Climate_data/clim_scale_grid_gpd_model.csv'))%>%
  left_join(read_csv('Data/clim_id_altitude_correspondance.csv'), by = "id")

true_change = grid_pred
true_change$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92)
true_change$rl = rl_mod_2((read_csv("output/gpd_model_fits/model_2_true.csv") %>% unlist %>% as.numeric),
                           true_change$rl_qnt, true_change$threshold_9, true_change$scale_9, true_change$glob_anom, true_change$altitude     )

rl_plot = gridExtra::grid.arrange(true_change %>%
                                    filter(year == 2022) %>%
                                    dplyr::select(longitude, latitude, id, rl) %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = rl), size = 0.7)+
                                    geom_sf(data = switzerland, col = 'black', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal)+
                                    labs(col = expression(paste('°C')))+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                          plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                  true_change %>%
                                    filter(year == 1971) %>%
                                    rename(rl_100_1971 = rl) %>%
                                    dplyr::select(longitude, latitude, id, rl_100_1971) %>%
                                    left_join(true_change %>%
                                                filter(year == 2022) %>%
                                                rename(rl_100_2022 = rl) %>%
                                                dplyr::select(longitude, latitude, id, rl_100_2022)) %>%
                                    mutate(rl_100_diff = rl_100_2022 - rl_100_1971) %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = rl_100_diff), size = 0.7)+
                                    geom_sf(data = switzerland, col = 'black', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal)+
                                    labs(col = expression(paste(nabla, '°C')))+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                          plot.margin = unit(c(0,-0.1,0,-0.1),"cm")))