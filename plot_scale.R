rm(list = ls())
library(tidyverse)
library(rnaturalearth) #library for map of switzerland
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('gpd_models.R')

marg_mod = "mod_0" # --- model to plot  --> will do both model 0 and model 1 to see the difference
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

#adding threshold 9 

threshold_9_df = vroom::vroom("Data/Observed_data/1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom("Data/processed/thresh_exceedance_lambda_num_quantiles_30.csv") 

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(-JJA)%>%
  rename("06" = Jun, "07" = Jul, "08" = Aug)%>%
  pivot_longer(cols = c("06", "07", "08"), 
               names_to = "month", 
               values_to = "glob_anom")%>%
  mutate(month = as.numeric(month))

obs_data = obs_data %>% 
  mutate(year = year(date), month = month(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = c("year", "month")) %>%
  left_join(threshold_9_df, by="id")


# map for plotting
# map of switzerland sf
switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

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

# --- get data for prediction
dat_for_pred = readRDS("output/quant_models_clim_num_quantiles_30.csv") %>% #we're working with model 0 so far
  dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>%
  left_join(read_csv("data/processed/sites_clim_sea_dist.csv"))%>%
  left_join(read_csv("data/processed/clim_scale_grid.csv"))

if(marg_mod == "mod_2"){
  this_fit_mod_2 = read_csv("output/gpd_model_fits/model_2_true.csv") %>% unlist %>% as.numeric
  pred_2 = my_predict_2(this_fit_mod_2, dat_for_pred$scale_9, dat_for_pred$loess_temp_anom, dat_for_pred$dist_sea)
  dat_for_pred$scale = pred_2$scale
  dat_for_pred$shape = pred_2$shape
}else if(marg_mod == 'mod_0'){
  this_fit_mod_0 = read_csv("output/gpd_model_fits/model_0_true.csv") %>% unlist %>% as.numeric
  
  pred_0 = my_predict_0(this_fit_mod_0, dat_for_pred$scale_9)
  dat_for_pred$scale = pred_0$scale
  dat_for_pred$shape = pred_0$shape
}else if(marg_mod == 'mod_1'){
  this_fit_mod_1 = read_csv("output/gpd_model_fits/model_1_true.csv") %>% unlist %>% as.numeric
  
  pred_1 = my_predict_1(this_fit_mod_1, dat_for_pred$scale_9, dat_for_pred$loess_temp_anom)
  dat_for_pred$scale = pred_1$scale
  dat_for_pred$shape = pred_1$shape
}


scale_and_thresh = gridExtra::grid.arrange(dat_for_pred %>%
                                             filter(year == 2020) %>%
                                             ggplot()+
                                             geom_point(aes(Long, Lat, col = threshold_9))+
                                             geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                             scale_x_continuous(breaks = -c(10,8,6))+
                                             
                                             scale_color_gradientn( colours=my_pal)+
                                             labs(col = expression(u[o]),
                                                  x = "", y = "")+
                                             theme_minimal()+
                                             theme(axis.text.x  = element_blank(),
                                                   axis.text.y  = element_blank(),
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                           dat_for_pred %>%
                                             filter(year == 1942) %>%
                                             ggplot()+
                                             geom_point(aes(Long, Lat, col = scale))+
                                             geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                             scale_color_gradientn( colours=my_pal)+
                                             scale_x_continuous(breaks = -c(10,8,6))+
                                             
                                             labs(col = expression(sigma[o]),
                                                  x = "", y = "")+
                                             theme_minimal()+
                                             theme(axis.text.x  = element_blank(),
                                                   axis.text.y  = element_blank(),
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                           dat_for_pred %>%
                                             filter(year == 1942) %>%
                                             mutate(diff = (dat_for_pred %>%
                                                              filter(year == 2020) %>%
                                                              pull(scale))-
                                                      (dat_for_pred %>%
                                                         filter(year == 1942)%>%
                                                         pull(scale))) %>%
                                             ggplot()+
                                             geom_point(aes(Long, Lat, col = diff))+
                                             geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                             scale_x_continuous(breaks = -c(10,8,6))+
                                             scale_color_gradientn( colours=my_pal)+
                                             labs(col = expression(nabla~sigma[o]),
                                                  x = "", y = "")+
                                             theme_minimal()+
                                             theme(axis.text.x  = element_blank(),
                                                   axis.text.y  = element_blank(),
                                                   panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   plot.margin = unit(c(0,-0.1,0,-0.1),"cm")), nrow=1)

ggsave(scale_and_thresh, filename = paste0("output/figs/thresh_scale_plot_",marg_mod,".pdf"), width = 6.75, height = 2)



# ----- figs to compare mod 1 and 2
this_fit_mod_0 = read_csv("output/gpd_model_fits/model_0_true.csv") %>% unlist %>% as.numeric
pred_0 = my_predict_0(this_fit_mod_0, dat_for_pred$scale_9)
dat_for_pred$scale1 = pred_0$scale

this_fit_mod_1 = read_csv("output/gpd_model_fits/model_1_true.csv") %>% unlist %>% as.numeric
pred_1 = my_predict_1(this_fit_mod_1, dat_for_pred$scale_9, dat_for_pred$loess_temp_anom)
dat_for_pred$scale2 = pred_1$scale

compare_sigma_model_0_and_1 = gridExtra::grid.arrange(dat_for_pred %>%
                                                        filter(year == 2020) %>%
                                                        ggplot()+
                                                        geom_point(aes(Long, Lat, col = scale1))+
                                                        geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                                        scale_x_continuous(breaks = -c(10,8,6))+
                                                        
                                                        scale_color_gradientn( colours=my_pal, limits = c(1.6, 2.6))+
                                                        labs(col = expression(sigma[o]),
                                                             x = "", y = "")+
                                                        theme_minimal(12)+
                                                        theme(axis.text.x  = element_blank(),
                                                              axis.text.y  = element_blank(),
                                                              panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                              axis.title.x = element_blank(),
                                                              axis.title.y = element_blank(),
                                                              plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                                      dat_for_pred %>%
                                                        filter(year == 1942) %>%
                                                        ggplot()+
                                                        geom_point(aes(Long, Lat, col = scale2))+
                                                        geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                                        scale_color_gradientn( colours=my_pal, limits = c(1.6, 2.6))+
                                                        scale_x_continuous(breaks = -c(10,8,6))+
                                                        
                                                        labs(col = expression(sigma[o]),
                                                             x = "", y = "")+
                                                        theme_minimal(12)+
                                                        theme(axis.text.x  = element_blank(),
                                                              axis.text.y  = element_blank(),
                                                              panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                              axis.title.x = element_blank(),
                                                              axis.title.y = element_blank(),
                                                              plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                                                      dat_for_pred %>%
                                                        filter(year == 1942) %>%
                                                        mutate(diff = (dat_for_pred %>%
                                                                         filter(year == 2020) %>%
                                                                         pull(scale2))-
                                                                 (dat_for_pred %>%
                                                                    filter(year == 1942)%>%
                                                                    pull(scale2))) %>%
                                                        ggplot()+
                                                        geom_point(aes(Long, Lat, col = diff))+
                                                        geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                                        scale_x_continuous(breaks = -c(10,8,6))+
                                                        scale_color_gradientn( colours=my_pal)+
                                                        labs(col = expression(nabla~sigma[o]),
                                                             x = "", y = "")+
                                                        theme_minimal(12)+
                                                        theme(axis.text.x  = element_blank(),
                                                              axis.text.y  = element_blank(),
                                                              panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                              axis.title.x = element_blank(),
                                                              axis.title.y = element_blank(),
                                                              plot.margin = unit(c(0,-0.1,0,-0.1),"cm")), nrow=1)

ggsave(compare_sigma_model_0_and_1, filename = paste0("output/figs/scale_compare_model_0_and_1.pdf"), width = 6.75, height = 2)