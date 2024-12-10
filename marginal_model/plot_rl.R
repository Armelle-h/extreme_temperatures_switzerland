#corresponds to Figure 5 in the paper


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


#not doing bootstrap 

model_fits = read_csv("output/gpd_model_fits/model_1_true.csv", 
                      col_names = c('bts', 'b0', 'b1', 'b2', 'xi')) #to change
marg_mod = "mod_1"
num_quantiles = 30

id_lon_lat = read.csv("Data/id_lon_lat_correspondance.csv")

clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_", marg_mod, "_quantiles_",num_quantiles,".csv")) %>%
  dplyr::select(-c(quantile, value)) %>% filter(year %in% c(1971, 2022)) #the filter by years not realy needed, already filtered

grid_pred = clim_date_w_quantile_mod %>%
  left_join(id_lon_lat, by = "id") %>%
  dplyr::select(id, longitude, latitude, year, threshold_9, thresh_exceedance_9, glob_anom) %>%
  left_join(read_csv('Data/Climate_data/clim_scale_grid_gpd_model.csv')) 

#This is "if I have bootstraps "
grid_pred$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92) 

#recheck, I want just the parameters
this_fit_pars = read_csv(paste0("Data/processed/glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                         col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))
this_scale_fit = my_predict_1(this_fit_pars, grid_pred$scale_9, grid_pred$glob_anom)
#recheck def of rl_mod_1 
this_rl_fit = rl_mod_1(this_fit_pars, grid_pred$rl_qnt, grid_pred$threshold_9, grid_pred$scale_9, grid_pred$glob_anom)

#When I don't have bootstraps
grid_pred = clim_date_w_quantile_mod %>%
  left_join(id_lon_lat, by = "id") %>%
  dplyr::select(id, longitude, latitude, year, threshold_9, thresh_exceedance_9, glob_anom) %>%
  left_join(read_csv('Data/Climate_data/clim_scale_grid_gpd_model.csv'))

true_change = grid_pred
true_change$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92)
true_change$rl = rl_mod_1((read_csv("output/gpd_model_fits/model_1_true.csv") %>% unlist %>% as.numeric),
                           true_change$rl_qnt, true_change$threshold_9, true_change$scale_9, true_change$glob_anom)

#corresponds to Figure 5 in the paper
rl_plot = gridExtra::grid.arrange(true_change %>%
                                    filter(year == 2022) %>%
                                    dplyr::select(longitude, latitude, id, rl) %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = rl))+
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
                                    geom_point(aes(longitude, latitude, col = rl_100_diff))+
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