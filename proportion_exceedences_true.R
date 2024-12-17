#corresponds to Figure 8 of the paper

rm(list=ls())
library(tidyverse)
library(scales)
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

proportion = function(temp_conditioned_on){

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
}

for (temp_conditioned_on in seq(27, 35)){
  proportion(temp_conditioned_on)
}


true_dat = rbind(read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_27.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_28.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_29.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_30.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_31.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_32.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_33.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_34.csv")),
                 read_csv(paste0("output/simulations/simulation_summary/nu_", nu_name,"_prop_exceedance_model_", true_folder,"_",marg_mod ,"_temp_condtioned_on_35.csv")))

prob_observing_T = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                            col_names = c('temp', 'p_1971', 'p_2022'))%>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"p_")) 

prob_observing_T$value = prob_observing_T$value

prob_observing_T = rbind(prob_observing_T %>% filter(name == 'p_1971') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name),
                         prob_observing_T %>% filter(name == 'p_2022') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name)) %>%
  mutate(year = as.numeric(year)) %>%
  rename(unconditional_factor = value)


true_dat = true_dat %>%
  pivot_longer(-temp) %>%
  mutate(year = str_remove(name,"prop_ex_")) 
true_dat$year = as.numeric(true_dat$year)
true_dat = true_dat %>%
  left_join(prob_observing_T)

true_dat$unconditioned = true_dat$value*true_dat$unconditional_factor

true_dat = rbind(true_dat %>%
                   dplyr::select(temp, value, year, typ) %>%
                   pivot_wider(names_from = c(year, typ), values_from = value) %>%
                   mutate(lab = 'conditioned'),
                 true_dat %>%
                   dplyr::select(temp, unconditioned, year, typ) %>%
                   pivot_wider(names_from = c(year, typ), values_from = unconditioned) %>%
                   mutate(lab = 'unconditioned'))

library(grid)

plot1 <- ggplot() +
  geom_smooth(data = true_dat %>% filter(lab == 'conditioned'),
              aes(temp, `1971_actual`, col = '1971', linetype = '1971'), se = FALSE) +
  geom_smooth(data = true_dat %>% filter(lab == 'conditioned'),
              aes(temp, `2022_actual`, col = '2022', linetype = '2022'), se = FALSE) +
  scale_x_continuous(limits = c(27, 34.5),
                     breaks = c(28, 30, 32, 34),
                     labels = paste0(c(28, 30, 32, 34), "°C")) +
  labs(x = "Temperature",
       y = expression(E[o]),
       col = "Year",
       linetype = "Year") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.x = element_text(size = 0))

# Plot 2: Unconditioned
plot2 <- ggplot() +
  geom_smooth(data = true_dat %>% filter(lab == 'unconditioned'),
              aes(temp, `1971_actual`, col = '1971', linetype = '1971'), se = FALSE) +
  geom_smooth(data = true_dat %>% filter(lab == 'unconditioned'),
              aes(temp, `2022_actual`, col = '2022', linetype = '2022'), se = FALSE) +
  scale_x_continuous(limits = c(27, 34.5),
                     breaks = c(28, 30, 32, 34),
                     labels = paste0(c(28, 30, 32, 34), "°C")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Temperature",
       y = "",
       col = "Year",
       linetype = "Year") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.x = element_text(size = 0))

# Combine the two plots side by side
plt <- gridExtra::grid.arrange(plot1, plot2, widths = c(1, 1), nrow = 1)


































plt = gridExtra::grid.arrange(geom_smooth(data = true_dat %>% filter(lab == 'conditioned'),
                                            aes(temp, `1971_actual`, col = '1971', linetype = '1971'), se = F)+
                                geom_smooth(data = true_dat %>% filter(lab == 'conditioned'),
                                            aes(temp, `2022_actual`, col = '2022', linetype = '2022'), se = F)+
                                scale_x_continuous(limits = c(27, 34.5),
                                                   breaks = c(28,  30,  32, 34),
                                                   label = paste0(c(28,  30,  32,34),"°C"))+
                                labs(x = "Temperature",
                                     y = expression(E[o]),
                                     col = "Year",
                                     fill = "Year",
                                     linetype = "Year")+
                                theme_minimal(12)+
                                theme(axis.text.x = element_text(angle = 0),
                                      axis.title.y = element_text(angle = 0, vjust = 0.5),
                                      legend.position = 'none',
                                      strip.text.x = element_text(size=0)),
                              geom_smooth(data = true_dat %>% filter(lab == 'unconditioned'),
                                            aes(temp, `1971_actual`, col = '1971', linetype = '1971'), se = F)+
                                geom_smooth(data = true_dat %>% filter(lab == 'unconditioned'),
                                            aes(temp, `2022_actual`, col = '2022', linetype = '2022'), se = F)+
                                scale_x_continuous(limits = c(27, 34.5),
                                                   breaks = c(28,  30,  32, 34),
                                                   label = paste0(c(28,  30,  32,34),"°C"))+
                                labs(x = "Temperature",
                                     y = "",
                                     col = "Year",
                                     fill = "Year",
                                     linetype = "Year")+
                                theme_minimal(12)+
                                theme(axis.text.x = element_text(angle = 0),
                                      axis.title.y = element_text(angle = 0, vjust = 0.5),
                                      legend.position = 'none',
                                      strip.text.x = element_text(size=0))+
                                facet_wrap(~lab, scales = 'free')+
                                scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                              labels = trans_format("log10", math_format(10^.x))),
                              widths = c(1, 1), nrow = 1)

