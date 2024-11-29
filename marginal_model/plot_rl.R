
gc()
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('gpd_models.R')
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

# read in corrected fits
bootstrap_model_fits = read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                                col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi'))


# # --- get scale prediction for each fit
bts_res = c()
marg_mod = "mod_2"
num_quantiles = 25

for(i in seq(300)){
  print(i)
  clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",i,".csv")) %>%
    dplyr::select(-c(glob_anom)) %>% filter(year %in% c(1931, 1942, 2020, 2022))
  
  grid_pred = clim_date_w_quantile_mod %>%
    dplyr::select(id, longitude, latitude, year, threshold_9, thresh_exceedance_9, glob_anom) %>%
    left_join(read_csv('data/processed/clim_scale_grid.csv')) %>%
    left_join(read_csv("data/processed/sites_clim_sea_dist.csv"))
  
  grid_pred$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92)
  
  this_fit_pars = bootstrap_model_fits[i,] %>% unlist %>% .[-1] %>% as.numeric()
  this_scale_fit = my_predict_2b(this_fit_pars, grid_pred$scale_9, grid_pred$glob_anom, grid_pred$dist_sea)
  this_rl_fit = rl_mod_2b(this_fit_pars, grid_pred$rl_qnt, grid_pred$threshold_9, grid_pred$scale_9, grid_pred$glob_anom, grid_pred$dist_sea)
  
  rbind(bts_res, grid_pred %>%
          dplyr::select(-c(threshold_9, thresh_exceedance_9, glob_anom, scale_9, rl_qnt)) %>%
          mutate(bts = i, rl = this_rl_fit)) %>%
    write_csv("output/scale_and_rl_for_plotting_mod_2.csv", append = T)
}



num_quantiles = 30
clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv")) %>%
  dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>% filter(year %in% c(1931, 1942, 2020, 2022))

grid_pred = clim_date_w_quantile_mod %>%
  dplyr::select(id, longitude, latitude, year, threshold_9, thresh_exceedance_9, glob_anom) %>%
  left_join(read_csv('Data/processed/clim_scale_grid.csv')) %>%
  left_join(read_csv("Data/processed/sites_clim_sea_dist.csv"))

true_change = grid_pred
true_change$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92)
true_change$rl = rl_mod_2b((read_csv("output/gpd_model_fits/model_2_true.csv") %>% unlist %>% as.numeric),
                           true_change$rl_qnt, true_change$threshold_9, true_change$scale_9, true_change$glob_anom, true_change$dist_sea)

bts_res = read_csv("output/scale_and_rl_for_plotting_mod_2.csv",
                   col_names = c("id", "longitude", "latitude", "year", "bts", "rl"))

rl_plot = gridExtra::grid.arrange(true_change %>%
                                    filter(year == 2020) %>%
                                    dplyr::select(longitude, latitude, id, rl) %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = rl))+
                                    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
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
                                  bts_res %>%
                                    filter(year  == 1942) %>%
                                    group_by(id) %>%
                                    summarise(longitude, latitude, id, lower_1940  = quantile(rl, 0.025)) %>% unique() %>%
                                    left_join(bts_res %>%
                                                filter(year  == 2020) %>%
                                                group_by(id) %>%
                                                summarise(longitude, latitude, id,
                                                          lower_2020  = quantile(rl, 0.025)) %>% unique()) %>%
                                    mutate(lower_change = lower_2020 - lower_1940) %>%
                                    unique() %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = lower_change))+
                                    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
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
                                          plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                  true_change %>%
                                    filter(year == 1942) %>%
                                    rename(rl_100_1942 = rl) %>%
                                    dplyr::select(longitude, latitude, id, rl_100_1942) %>%
                                    left_join(true_change %>%
                                                filter(year == 2020) %>%
                                                rename(rl_100_2020 = rl) %>%
                                                dplyr::select(longitude, latitude, id, rl_100_2020)) %>%
                                    mutate(rl_100_diff = rl_100_2020 - rl_100_1942) %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = rl_100_diff))+
                                    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
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
                                          plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                  bts_res %>%
                                    filter(year %in% c(1942)) %>%
                                    group_by(id) %>%
                                    summarise(longitude, latitude, id,lower_1940  = quantile(rl, 0.975)) %>% unique() %>%
                                    left_join(bts_res %>%
                                                filter(year %in% c(2020)) %>%
                                                group_by(id) %>%
                                                summarise(longitude, latitude, id,
                                                          lower_2020  = quantile(rl, 0.975)) %>% unique()) %>%
                                    mutate(lower_change = lower_2020 - lower_1940) %>%
                                    ggplot()+
                                    geom_point(aes(longitude, latitude, col = lower_change))+
                                    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal, breaks = c(2.5,3.5))+
                                    labs(col = expression(paste(nabla, '°C')))+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(),
                                          axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                          plot.margin = unit(c(0,-0.1,0,-0.1),"cm")), nrow = 1, widths = c(0.01,0.01,0.01016,0.01))

#ggsave(rl_plot, filename ="output/figs/marginal_rl.pdf", height = 2, width = 9)