#corresponds to figure 6 of the paper

gc()
rm(list = ls())
library(tidyverse)

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('marginal_model/gpd_models.R')

num_quantiles = 30

scales = vroom::vroom("Data/Observed_data/plain_obs_data_gpd_model.csv")%>% 
  select(stn, scale_9) %>%
  filter(!(stn %in% c("WSLBTB", "WSLHOB")))%>% 
  unique()

threshold_9_df = readRDS(paste0("Data/processed/plain_obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv")) %>%
  select(stn, threshold_9)  %>%
  filter(!(stn %in% c("WSLBTB", "WSLHOB")))%>%  
  unique()

legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")%>%
  select(stn, longitude_proj, latitude_proj)%>%
  unique()

obs_sites = scales %>%
  left_join(threshold_9_df, by = "stn")%>%
  left_join(legend_data, by = "stn")%>%
  unique()

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

thresh_ex = vroom::vroom(paste0("Data/processed/plain_glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv"))%>%
  filter(!(stn %in% c("WSLBTB", "WSLHOB")))

obs_grid = thresh_ex  %>%
  left_join(obs_sites, by = "stn") %>%
  left_join(glob_anomaly_reshaped, by = "year")%>%
  filter(year %in% c(1971, 2022))

gpd_pars = read_csv("output/gpd_model_fits/model_1_true.csv")
obs_grid$scale = my_predict_1(unlist(gpd_pars), obs_grid$scale_9, obs_grid$glob_anom)$scale
obs_grid$shape = my_predict_1(unlist(gpd_pars), obs_grid$scale_9, obs_grid$glob_anom)$shape

# ---- Read in simulations
my_simulations_extremes = c()
for(i in seq(1, 100)){
  #stuck I don't have that
  my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true_robust/nu_015/nu_015_mod_1_run_",i)))
}

my_simulations_extremes = my_simulations_extremes[seq(25000)]

# ---- standardise simulations to have cost = 1
my_simulations_standardised = list()
for (i in 1:length(my_simulations_extremes)) {
  this_cost = median(my_simulations_extremes[[i]]) #swithced mean to median
  my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost
}

# --- maximum frechet margin at each of my simulation sites
max_at_each_site = c()
for(s in seq(length(my_simulations_standardised[[1]]))){
  max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
}

m = length(my_simulations_standardised)
L = 300

res_1971 = c()
res_2022 = c()

r_thresh = readRDS("Data/spatial_threshold.rds")

for(temp_i_want in seq(26,34)){
  print(temp_i_want)
  
  obs_grid$pareto = 1/(obs_grid$thresh_exceedance_9*(1-evd::pgpd(temp_i_want - obs_grid$threshold_9,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
  T_2022 = obs_grid %>% filter(year == 2022) %>% pull(pareto) 
  T_1971 = obs_grid %>% filter(year == 1971) %>% pull(pareto) 
  
  T_1971[T_1971 == -Inf] = Inf
  T_2022[T_2022 == -Inf] = Inf
  
  b_2022 = min(T_2022/max_at_each_site)
  b_1971 = min(T_1971/max_at_each_site)
  m = length(my_simulations_standardised)
  
  above_1971 = 0
  above_2022 = 0
  for(i in seq(L)){
    print(i)
    this_set_scaled = my_simulations_standardised %>%
      map(~{ 
        .x*evd::rgpd(n=1, loc = 1, scale = 1, shape = 1)
      }) %>%
      map(~{ 
        if(median(.x)>=r_thresh){ #using median for robust simulation
          c(sum(.x*b_2022 > T_2022) > 0,
            sum(.x*b_1971 > T_1971) > 0)
        }
      })
    above_2022 = above_2022 + (lapply(this_set_scaled, "[[", 1) %>% unlist %>% sum)
    above_1971 = above_1971 + (lapply(this_set_scaled, "[[", 2) %>% unlist %>% sum)
  }
  
  tibble(temp = temp_i_want, 
         p_1971 = above_1971 / (length(my_simulations_standardised)*L*b_1971),
         p_2022 = above_2022 / (length(my_simulations_standardised)*L*b_2022)) %>%
    write_csv("output/importance_sampling/prob_extreme_temp_imp_samp_mod_1.csv", append = T)
}

# ----- plot figures

read_csv("output/importance_sampling/prob_extreme_temp_imp_samp_mod_1.csv",
         col_names = c('temp', 'p_1971', 'p_2022')) %>%
  filter(temp<=35, temp>25) %>%
  group_by(temp) %>%
  summarise(p_1971_lower = quantile(p_1971, 0.1, na.rm=T),
            p_1971_upper = quantile(p_1971, 0.9, na.rm=T),
            p_2022_lower = quantile(p_2022, 0.1, na.rm=T),
            p_2022_upper = quantile(p_2022, 0.9, na.rm=T)) %>%
  ggplot()+
  geom_ribbon(aes(temp,  ymin = 1/(92*p_1971_lower), ymax = 1/(92*p_1971_upper), fill = '1971'), alpha = 0.5)+
  geom_ribbon(aes(temp,  ymin = 1/(92*p_2022_lower), ymax = 1/(92*p_2022_upper), fill = '2022'), alpha = 0.5)+
  geom_line(data = read_csv("output/importance_sampling/prob_extreme_temp_imp_samp_mod_1.csv",
                            col_names = c('temp', 'p_1971', 'p_2022')) %>%
              filter(temp<=35, temp>25), aes(temp, 1/(92*p_1971), col = '1971'))+
  geom_line(data = read_csv("output/importance_sampling/prob_extreme_temp_imp_samp_mod_1.csv",
                            col_names = c('temp', 'p_1971', 'p_2022')) %>%
              filter(temp<=35, temp>25), aes(temp, 1/(92*p_2022), col = '2022'), linetype = 'dashed')+
  theme_minimal()+
  labs(x = "Temperature",
       y = "Return period (Years)",
       fill = 'Year', col = "Year")+
  theme_minimal(12)+
  scale_x_continuous(limits = c(26, 34),
                     breaks = c(26, 28,  30,  32, 34),
                     label = paste0(c(26, 28,  30,  32,34),"°C"))+
  scale_y_log10(breaks = c(0.1, 1, 10,  100,  1000),
                labels = c(0.1, 1, 10,  100,  1000))+
  coord_flip(ylim = c(0.075, 1000), xlim = c(26, 34))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none')