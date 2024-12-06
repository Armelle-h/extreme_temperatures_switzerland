rm(list=ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


marg_mod = 'mod_1'
temp_conditioned_on = 28
obs_sites = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
  select("stn", "longitude_proj", "latitude_proj")%>%
  unique()

yr = 2022
nu_name = "015"
robust = TRUE 

if (robust == TRUE){
  true_folder = "true_robust"
} else{
  true_folder = "true"
}

if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",temp_conditioned_on))){
  
  my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",temp_conditioned_on))
  prop_ex_2022 = c()
  for(tmp in c(temp_conditioned_on)){

    pareto_val = read_csv("Data/processed/plain_obs_grid_simulated_on.csv") %>%
      left_join(obs_sites) %>%
      left_join(read_csv(paste0("output/plain_obs_sites_extreme_temps_frechet_pareto_scale_",marg_mod,".csv")) %>% filter(temp == tmp, year == yr)) %>% 
      pull(pareto_value)
    
    pareto_val[pareto_val == -Inf] = Inf
    
    num_exceed_tmp = my_simulations %>%
      map(~{ sum(.x > pareto_val)}) %>%
      unlist()
    
    prop_ex_2022 = c(prop_ex_2022, mean(num_exceed_tmp[num_exceed_tmp>0]/154))
  }
}

yr = 1971
nu_name = "015"
robust = TRUE 

if (robust == TRUE){
  true_folder = "true_robust"
} else{
  true_folder = "true"
}
  
if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",temp_conditioned_on))){  

  my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",temp_conditioned_on))
  prop_ex_1971 = c()
  for(tmp in c(temp_conditioned_on)){
    
    pareto_val = read_csv("Data/processed/plain_obs_grid_simulated_on.csv") %>%
      left_join(obs_sites) %>%
      left_join(read_csv(paste0("output/plain_obs_sites_extreme_temps_frechet_pareto_scale_",marg_mod,".csv")) %>% filter(temp == tmp, year == yr)) %>% 
      pull(pareto_value)
    
    pareto_val[pareto_val == -Inf] = Inf
    
    num_exceed_tmp = my_simulations %>%
      map(~{ sum(.x > pareto_val)}) %>%
      unlist()
    
    prop_ex_1971 = c(prop_ex_1971, mean(num_exceed_tmp[num_exceed_tmp>0]/154))
  }
}
  
tibble(temp = temp_conditioned_on, prop_ex_1971, prop_ex_2022) %>%
  write_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod,"_temp_condtioned_on_",temp_conditioned_on,".csv"))
  
#do a scatter plot, x axis the stations, y axis the exceedance for 1971, 2022 with different colors for 1971 and 2022


