#corresponds to Fig 5 in the paper

gc()
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(doParallel)
library(foreach)
library(sf)
library(rnaturalearth)

source('marginal_model/gpd_models.R')
source('mvPot/simulPareto.R')

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

#general variables 
nu_name = "015"
nu_val = 0.15
robust = TRUE
  
if (robust == TRUE){
  true_folder = "true_robust"
} else{
  true_folder = "true"
}

filter_nb = 10 

# map for plotting
switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)


marg_mod = 'mod_1'
obs_smoothed_quantiles = readRDS(paste0("output/quant_models_clim_num_", marg_mod, "_quantiles_30.csv"))%>% 
  filter(id %%filter_nb == 0)

# --- marginal model names

#clim_grid_simulated_on is the same as id_lon_lat restricted to the columns longitude and latitude for me where longitude and latitude are projected
#clim_grid_simulated_on_not_projected is the same as id_lon_lat restricted to the columns longitude and latitude for me

#filtering mode 10 to start with and then with mod 5.
locs_to_pred = read.csv("Data/plain_id_lon_lat_correspondance.csv")%>% filter(id %%filter_nb == 0) %>% select(longitude_proj, latitude_proj) %>% unique() 

#read_csv("data/processed/clim_data_full.csv") %>% dplyr::select(Long.projected, Lat.projected) %>% unique %>% write_csv("Data/processed/clim_grid_simulated_on.csv")
#locs_to_pred = as.data.frame(read_csv("data/processed/clim_grid_simulated_on.csv"))
#read_csv("data/processed/clim_data_full.csv") %>% dplyr::select(Long, Lat) %>% unique %>% write_csv("data/processed/clim_grid_simulated_on_not_projected.csv")

variogram_model = function(h){
  nu=nu_val 
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}

fit = read_csv(paste0("output/rpareto_model_fits/robust_plain_nu_", nu_name,"_true_rpareto_fits_model_", marg_mod, ".csv"),
               col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()

alpha_var <<- fit[1]
beta_var <<- fit[2]


nlocs = nrow(locs_to_pred) # number of sites
loc.id.pairs = expand.grid(seq(nlocs) ,seq(nlocs))
loc.pairs = cbind(locs_to_pred[loc.id.pairs[,1],], locs_to_pred[loc.id.pairs[,2],]) # all pairs of locations

dstncs = loc.pairs %>%
  apply(MARGIN = 1, FUN = function(x) sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2))

#variogram matrix
vario_matrix = (variogram_model(dstncs)) %>% matrix(nrow = nlocs, byrow = T)

# this returns samples in the unit Pareto margin
simulations <- simulPareto(n = 500,
                          vario_mat = vario_matrix,
                          robust = TRUE) #when the median cost function is used

#most_ex = map(simulations, mean) %>% unlist %>% sort() %>% tail(5)
most_ex = map(simulations, median) %>% unlist %>% sort() %>% tail(5) #should set a general boolean "robust"
#ind_of_ex = which((map(simulations, mean) %>% unlist) %in% most_ex)
ind_of_ex = which((map(simulations, median) %>% unlist) %in% most_ex)

sims = read.csv("Data/plain_id_lon_lat_correspondance.csv")%>% filter(id %%filter_nb ==0) %>% select(id) %>% unique()  %>%
  mutate(sim_1 = simulations[[ind_of_ex[1]]],
         sim_2 = simulations[[ind_of_ex[2]]],
         sim_3 = simulations[[ind_of_ex[3]]],
         sim_2 = simulations[[ind_of_ex[4]]],
         sim_5 = simulations[[ind_of_ex[5]]]) %>%
  pivot_longer(-c(id)) %>%
  rename(pareto = value)

sims = sims %>%
  left_join(read_csv("Data/Climate_data/clim_scale_grid_gpd_model.csv") %>% 
              filter(id %% filter_nb ==0) %>%
              dplyr::select(id, scale_9), by = "id")

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

sims = rbind(sims %>% mutate(year = 2022), sims %>% mutate(year = 1971)) %>%
  left_join(glob_anomaly_reshaped %>%
              filter(year %in% c(2022, 1971)) %>% 
              unique()) %>%
  left_join(obs_smoothed_quantiles %>%
              dplyr::select(id, tau_to_temp, threshold_9, thresh_exceedance_9, year))

this_fit_mod = read_csv("output/gpd_model_fits/model_1_true.csv") %>%
  unlist() %>% as.numeric
sims$shape = this_fit_mod[length(this_fit_mod)]
sims$scale = my_predict_1(this_fit_mod, sims$scale_9, sims$glob_anom)$scale


# min of 5
sims$unif = evd::pgev(sims$pareto, 1,1,1)
extreme_ind = sims$unif > (1-sims$thresh_exceedance_9)
sims$date_scale = NA
sims$date_scale[extreme_ind] = evd::qgpd(p = 1 + (sims$unif[extreme_ind]-1)/sims$thresh_exceedance_9[extreme_ind],
                                         shape = sims$shape[1],
                                         loc = sims$threshold_9[extreme_ind],
                                         scale = sims$scale[extreme_ind])

for(i in seq(nrow(sims))){
  if(!extreme_ind[i]){
    sims$date_scale[i] = sims[i,]$tau_to_temp[[1]](sims[i,]$unif)
  }
}

sims = sims %>%
  rename(lab = name)%>%
  left_join(read.csv("Data/plain_id_lon_lat_correspondance.csv")%>% filter(id %%filter_nb == 0) %>% select(id, longitude, latitude) %>% unique() )


#supposed to have 5 simulations. Something went wrong with one of them

plt = gridExtra::grid.arrange(sims %>%
                                dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                filter(lab %in% c("sim_1", "sim_2"))%>%
                                filter(year == 2022) %>%
                                ggplot()+
                                geom_point(aes(longitude, latitude, col = date_scale), size = 0.5)+
                                facet_grid(year~lab)+
                                geom_sf(data = switzerland, alpha = 0, col = 'black')+
                                ggplot2::scale_color_gradientn(colors = my_pal)+
                                theme_minimal()+
                                labs(col = '°C')+
                                theme(axis.text = element_blank(),
                                      strip.text = element_blank(),
                                      axis.title = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      plot.margin=unit(c(0,0,-0.1,0),"cm")),
                              sims %>%
                                dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                filter(lab %in% c("sim_1", "sim_2"))%>%
                                filter(year == 2022) %>% 
                                rename(temp_new = date_scale) %>%
                                dplyr::select(-year) %>%
                                left_join(sims %>%
                                            dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                            filter(year == 1971) %>% 
                                            rename(temp_old = date_scale)) %>%
                                mutate(diff = temp_new - temp_old) %>% 
                                ggplot()+
                                geom_point(aes(longitude, latitude, col = diff), size = 0.5)+
                                geom_sf(data = switzerland, alpha = 0, col = 'black')+
                                facet_grid(year~lab)+
                                ggplot2::scale_color_gradientn(colors = my_pal)+
                                theme_minimal()+
                                labs(col = expression(paste(nabla, '°C')))+
                                theme(axis.text = element_blank(),
                                      strip.text = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_blank(),
                                      plot.margin=unit(c(-0.1,0,0,0),"cm")), nrow = 2)



plt = gridExtra::grid.arrange(sims %>%
                                dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                filter(lab %in% c("sim_3", "sim_5"))%>%
                                filter(year == 2022) %>%
                                ggplot()+
                                geom_point(aes(longitude, latitude, col = date_scale), size = 0.5)+
                                facet_grid(year~lab)+
                                geom_sf(data = switzerland, alpha = 0, col = 'black')+
                                ggplot2::scale_color_gradientn(colors = my_pal)+
                                theme_minimal()+
                                labs(col = '°C')+
                                theme(axis.text = element_blank(),
                                      strip.text = element_blank(),
                                      axis.title = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      plot.margin=unit(c(0,0,-0.1,0),"cm")),
                              sims %>%
                                dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                filter(lab %in% c("sim_3", "sim_5"))%>%
                                filter(year == 2022) %>% 
                                rename(temp_new = date_scale) %>%
                                dplyr::select(-year) %>%
                                left_join(sims %>%
                                            dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                            filter(year == 1971) %>% 
                                            rename(temp_old = date_scale)) %>%
                                mutate(diff = temp_new - temp_old) %>% 
                                ggplot()+
                                geom_point(aes(longitude, latitude, col = diff), size = 0.5)+
                                geom_sf(data = switzerland, alpha = 0, col = 'black')+
                                facet_grid(year~lab)+
                                ggplot2::scale_color_gradientn(colors = my_pal)+
                                theme_minimal()+
                                labs(col = expression(paste(nabla, '°C')))+
                                theme(axis.text = element_blank(),
                                      strip.text = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_blank(),
                                      plot.margin=unit(c(-0.1,0,0,0),"cm")), nrow = 2)








#8 plots in one plot


plt = gridExtra::grid.arrange(sims %>%
                                dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                filter(year == 2022) %>%
                                ggplot()+
                                geom_point(aes(longitude, latitude, col = date_scale))+
                                facet_grid(year~lab)+
                                geom_sf(data = switzerland, alpha = 0, col = 'black')+
                                ggplot2::scale_color_gradientn(colors = my_pal)+
                                theme_minimal()+
                                labs(col = '°C')+
                                theme(axis.text = element_blank(),
                                      strip.text = element_blank(),
                                      axis.title = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      plot.margin=unit(c(0,0,-0.1,0),"cm")),
                              sims %>%
                                dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                filter(year == 2022) %>% 
                                rename(temp_new = date_scale) %>%
                                dplyr::select(-year) %>%
                                left_join(sims %>%
                                            dplyr::select(longitude, latitude, year, date_scale, lab) %>%
                                            filter(year == 1971) %>% 
                                            rename(temp_old = date_scale)) %>%
                                mutate(diff = temp_new - temp_old) %>% 
                                ggplot()+
                                geom_point(aes(longitude, latitude, col = diff))+
                                geom_sf(data = switzerland, alpha = 0, col = 'black')+
                                facet_grid(year~lab)+
                                ggplot2::scale_color_gradientn(colors = my_pal)+
                                theme_minimal()+
                                labs(col = expression(paste(nabla, '°C')))+
                                theme(axis.text = element_blank(),
                                      strip.text = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_blank(),
                                      plot.margin=unit(c(-0.1,0,0,0),"cm")), nrow = 2)

#Need to plot 4 and 4 (else too small)