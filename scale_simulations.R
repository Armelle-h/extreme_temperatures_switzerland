# scale up simulations

#NEED TO REDO THE SIMULATIONS (even tough I don't think it will change much)


rm(list=ls())
library(tidyverse)
library(evd)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source("marginal_model/gpd_models.R")

#no temporal covariates in this code !!!! Only spatial ones (if needed)

scale_true_sim = function(marg_mod, yr, tmp, nu_name, robust = TRUE){
  
  if (robust == TRUE){
    true_folder = "true_robust"
  } else{
    true_folder = "true"
  }
  
  #extreme_temp_frechet = read_csv(paste0("output/plain_obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv"))
  
  #doing the filtering here to reduce the memory strain
  #defined for all stations, for all observed years
  extreme_temp_pareto = read_csv(paste0("output/plain_obs_sites_extreme_temps_frechet_pareto_scale_",marg_mod,".csv"))%>%
       filter(temp == tmp, year == yr)%>%
       select(-c(frechet_value, temp, -year))
    
  data_sets = list.files(paste0('output/simulations/simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/'))
  data_sets = data_sets[which(grepl(marg_mod, data_sets))]
  
  legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")
  
  threshold_9_df = vroom::vroom("Data/processed/plain_1971_2022_JJA_obs_data_bulk_model.csv")%>%
    dplyr::select(threshold_9, stn)%>%
    unique()
  
  obs_sites = vroom::vroom("Data/Observed_data/plain_obs_data_gpd_model.csv") %>%
    select(stn, scale_9)%>%
    unique()%>%
    left_join(threshold_9_df, by="stn")%>%
    filter(!(stn %in% c("WSLBTB", "WSLHOB")))
  
  grid_simulated = as.data.frame(read_csv("Data/processed/plain_obs_grid_simulated_on.csv")) %>%
    left_join(legend_data %>%select(stn, longitude_proj, latitude_proj), by = c("longitude_proj", "latitude_proj") )%>%
    left_join(obs_sites, by= "stn")
  
  # ---- for each model get scale shape and lambda 
  this_grid = grid_simulated %>%
    left_join(read_csv("Data/processed/plain_glob_anomaly_thresh_exceedance_lambda_num_quantiles_30.csv") %>% filter((year == yr) & !(stn %in% c("WSLBTB", "WSLHOB")))) %>%
    left_join(extreme_temp_pareto )
  
  my_simulations_extremes = c()
  
  #100 simulations of 2500 simulations of Y_r^P where each component is associated with a station
  for(i in seq(1, 100)){ 
    my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/", true_folder,"/nu_",nu_name, "/nu_", nu_name,"_",marg_mod,"_run_",i)))
  }
  
  #y_i^P(s)/r_i
  my_simulations_standardised = list()
  for (i in 1:length(my_simulations_extremes)) {
    if (robust == TRUE){
      this_cost = median(my_simulations_extremes[[i]]) 
    } else{
      this_cost = mean(my_simulations_extremes[[i]])
    }
    my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost 
  }
  
  # --- maximum pareto margin at each of my simulation sites
  
  #computing omega_(m)(s) for each site
  max_at_each_site = c()
  for(s in seq(length(my_simulations_standardised[[1]]))){
    max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
  }
  
  this_grid$scaler = this_grid$pareto_value/max_at_each_site 
  this_grid$scaler[this_grid$scaler == -Inf] = Inf
  br = min(this_grid$scaler) #corresponds to b_T(t)
  
  my_r = evd::rgpd(n=length(my_simulations_extremes), 1,1,1)
  # --- scale simulations back up
  my_simulations_rescaled = my_simulations_standardised
  for (i in 1:length(my_simulations_standardised)) {
    my_simulations_rescaled[[i]] = (my_r[i] * br)*my_simulations_standardised[[i]] #corresponds to what is inside the importance sampling estimator
  }
  
  my_simulations_rescaled %>% saveRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))
}


rm(nu_name)
nu_name="015"
robust=TRUE

job::job({scale_true_sim("mod_1", 1971, 27, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 27, nu_name, robust)})
job::job({scale_true_sim("mod_1", 1971, 28, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 28, nu_name, robust)})
job::job({scale_true_sim("mod_1", 1971, 29, nu_name, robust)})

job::job({scale_true_sim("mod_1", 2022, 29, nu_name, robust)})
job::job({scale_true_sim("mod_1", 1971, 30, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 30, nu_name, robust)})
job::job({scale_true_sim("mod_1", 1971, 31, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 31, nu_name, robust)})

job::job({scale_true_sim("mod_1", 1971, 32, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 32, nu_name, robust)})
job::job({scale_true_sim("mod_1", 1971, 33, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 33, nu_name, robust)})

job::job({scale_true_sim("mod_1", 1971, 34, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 34, nu_name, robust)})
job::job({scale_true_sim("mod_1", 1971, 35, nu_name, robust)})
job::job({scale_true_sim("mod_1", 2022, 35, nu_name, robust)})