
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

#computes the chi dependence static for a given model, year and temperature
calc_chi_true = function(marg_mod, yr, tmp, nu_name, robust=FALSE){
  
  if (robust == TRUE){
    true_folder = "true_robust"
  } else{
    true_folder = "true"
  }
  
  set.seed(123456)
  #number of site pairs to be sampled
  num_samples = 1000 #used to be 8000, but in the report they mentioned they used 500 bootstraps 
  
  sites = read_csv("Data/processed/plain_obs_pairs_with_dist.csv") %>% sample_n(num_samples)
  
  obs_sites = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
    select("stn", "longitude_proj", "latitude_proj")%>%
    unique()
  
  grid_simulated = as.data.frame(read_csv("Data/processed/plain_obs_grid_simulated_on.csv")) %>%
    left_join(obs_sites) %>% as_tibble()
  
  pareto_val = grid_simulated %>%
    left_join(read_csv(paste0("output/plain_obs_sites_extreme_temps_frechet_pareto_scale_",marg_mod,".csv")) %>% filter(temp == tmp, year == yr)) %>% 
    pull(pareto_value)
  
  if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))){
    my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/", true_folder,"/nu_",nu_name,"/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))
    
    for(s in seq(nrow(sites))){ #we have 8000 sites --> takes a really long time to run :(
      
      if((s %% 100) == 0){
        print(s)
      }
      
      id1 = sites[s,]$V1
      id2 = sites[s,]$V2
      
      # ---- to get all sims at loc with id v1 ---> unlist(lapply(my_simulations, "[[", sites[s,]$V1))
      
      #computes chi = number of simulations where id1 and id2 exceed the threshold / nb od simulations where id1 exceeds.
      #conditional proba that site id2 exceeds threshold given that site id1 does.
      #we don't necessarily have values for all the pairs, it depends on the pairs simulated
      
      #we compute P[X^P_o(t,s)>T_P_o(t,s)] where T is tmp and the transformation of T to pareto margin is done in the df pareto_val
      #as we already have the conditioning on t (t=yr), we only need to extract the value associated with the station.
      #Hence the pareto_val[which(grid_simulated$stn==id1)].
      id1_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$stn == id1))) > pareto_val[which(grid_simulated$stn == id1)])
      id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$stn == id2))) > pareto_val[which(grid_simulated$stn == id2)])
      
      tibble(sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
        write_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_",yr, "_conditioned_on_",tmp,".csv"),append = T)
    }
  }
}

#20 minutes for 1000 simulations. A bit more than 2 hours for 8000 !!!

#27, 28, 29 correspond to the 0.8, 0.85, 0.9 quantile of the observed temperatures
job::job({calc_chi_true(marg_mod = "mod_1", yr = 2022, tmp = 27, nu_name = "007", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 2022, tmp = 28, nu_name = "007", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 2022, tmp = 29, nu_name = "007", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 1971, tmp = 27, nu_name = "007", robust = TRUE)})

job::job({calc_chi_true(marg_mod = "mod_1", yr = 1971, tmp = 28, nu_name = "007", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_1", yr = 1971, tmp = 29, nu_name = "007", robust = TRUE)})

#for comparison purpose
job::job({calc_chi_true(marg_mod = "mod_1", yr = 1971, tmp = 36, nu_name = "012", robust = TRUE)}) #pretty high temp to see what happens

job::job({calc_chi_true(marg_mod = "mod_0", yr = 2022, tmp = 27, nu_name = "012", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_0", yr = 2022, tmp = 28, nu_name = "012", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_0", yr = 2022, tmp = 29, nu_name = "012", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_0", yr = 1971, tmp = 27, nu_name = "012", robust = TRUE)})

job::job({calc_chi_true(marg_mod = "mod_0", yr = 1971, tmp = 28, nu_name = "012", robust = TRUE)})
job::job({calc_chi_true(marg_mod = "mod_0", yr = 1971, tmp = 29, nu_name = "012", robust = TRUE)})

# # ------ PLOT MODELS

marg_mod = 'mod_1'

robust = TRUE

if (robust == TRUE){
  true_folder = "true_robust"
} else{
  true_folder = "true"
}

nu_name = "012"

chi_true = rbind(read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_1971_conditioned_on_28.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "27°C", year = '1971'),
                 read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_1971_conditioned_on_29.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "28°C", year = '1971'),
                 read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_1971_conditioned_on_30.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "29°C", year = '1971'),
                 read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_2022_conditioned_on_28.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "27°C", year = '2022'),
                 read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_2022_conditioned_on_29.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "28°C", year = '2022'),
                 read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_2022_conditioned_on_30.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "29°C", year = '2022'),
                 read_csv(paste0("output/simulations/simulation_summary/", true_folder, "_nu_",nu_name,"_chi_data_scale_clim_grid_model_",marg_mod,"_yr_2022_conditioned_on_30.csv"),
                          col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "36°C", year = '2022'))


max_distance <- max(chi_true$distance, na.rm = TRUE)

#summarize is deprecated but this gives me what I want 
chi_true_summary = chi_true %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, max_distance, length.out = 30))) %>% #used to be 20 but in the report they mentioned 30 binned distances
  group_by(dist_bin) %>%
  summarise(year, temp, chi, distance = median(distance))  %>%
  group_by(distance, year, temp) %>%
  drop_na() %>%
  summarise(mean_chi = mean(chi), median_chi = median(chi))


chi_true_violin = chi_true %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, max_distance, length.out = 30))) %>% #used to be 20 but in the report they mentioned 30 binned distances
  group_by(dist_bin) %>%
  summarise(year, temp, chi, distance = median(distance))  %>%
  group_by(distance, year, temp) %>%
  drop_na()


chi_true_summary %>%
  ggplot(aes(x = distance, y = median_chi)) +
  geom_point() +
  facet_wrap(~ temp, scales = "free") +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "Scatter plots of distance vs mean by temperature")

#scatter plot, unbinned distance vs estimate of chi

chi_true %>%
  filter(temp == "27°C")%>%
  ggplot(aes(x = distance, y = chi)) +
  geom_point() +
  ylim(0, 1) +
  theme_minimal()+
  labs(title = "27°C")

chi_true %>%
  filter(temp == "28°C")%>%
  ggplot(aes(x = distance, y = chi)) +
  geom_point() +
  ylim(0, 1) +
  theme_minimal()+
  labs(title = "28°C")

chi_true %>%
  filter(temp == "29°C")%>%
  ggplot(aes(x = distance, y = chi)) +
  geom_point() +
  ylim(0, 1) +
  theme_minimal()+
  labs(title = "29°C")

chi_true %>%
  filter(temp == "36°C")%>%
  ggplot(aes(x = distance, y = chi)) +
  geom_point() +
  ylim(0, 1) +
  theme_minimal()+
  labs(title = "36°C")


#violin plot of chi for each binned temperature


chi_true_violin %>%
  filter(temp == "27°C")%>%
  ggplot(aes(x = factor(distance), y = chi, fill = temp)) +
  geom_violin() +
  facet_wrap(~ temp, scales = "free") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5) # Centers the title
  ) +
  labs(
    title = "27°C",
    x = "Binned Distance",
    y = "Chi"
  )

chi_true_violin %>%
  filter(temp == "28°C")%>%
  ggplot(aes(x = factor(distance), y = chi, fill = temp)) +
  geom_violin() +
  facet_wrap(~ temp, scales = "free") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5) # Centers the title
  ) +
  labs(
    title = "28°C",
    x = "Binned Distance",
    y = "Chi"
  )


chi_true_violin %>%
  filter(temp == "29°C")%>%
  ggplot(aes(x = factor(distance), y = chi, fill = temp)) +
  geom_violin() +
  facet_wrap(~ temp, scales = "free") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5) # Centers the title
  ) +
  labs(
    title = "29°C",
    x = "Binned Distance",
    y = "Chi"
  )

chi_true_violin %>%
  filter(temp == "36°C")%>%
  ggplot(aes(x = factor(distance), y = chi, fill = temp)) +
  geom_violin() +
  facet_wrap(~ temp, scales = "free") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5) # Centers the title
  ) +
  labs(
    title = "36°C",
    x = "Binned Distance",
    y = "Chi"
  )



# # ----- Uncondition
prob_2022_28 = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1971', 'p_2022'))  %>%
  filter(temp == 28) %>% pull(p_2022)

prob_2022_29  = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                         col_names = c('temp', 'p_1971', 'p_2022'))  %>%
  filter(temp == 29) %>% pull(p_2022)

prob_2022_30 = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1971', 'p_2022'))  %>%
  filter(temp == 30) %>% pull(p_2022)

prob_1971_28 = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1971', 'p_2022'))  %>%
  filter(temp == 28) %>% pull(p_1971)

prob_1971_29 = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1971', 'p_2022'))  %>%
  filter(temp == 29) %>% pull(p_1971)

prob_1971_30 = read_csv(paste0("output/importance_sampling/prob_extreme_temp_imp_samp_",marg_mod,".csv"),
                        col_names = c('temp', 'p_1971', 'p_2022'))  %>%
  filter(temp == 30) %>% pull(p_1971)

true_unconditioned = chi_true

true_unconditioned = rbind(true_unconditioned %>% filter(year == 1971, temp == '28°C') %>% mutate(chi = prob_1971_28*chi),
                           true_unconditioned %>% filter(year == 1971, temp == '29°C') %>% mutate(chi = prob_1971_29*chi),
                           true_unconditioned %>% filter(year == 1971, temp == '30°C') %>% mutate(chi = prob_1971_30*chi),
                           true_unconditioned %>% filter(year == 2022, temp == '28°C') %>% mutate(chi = prob_2022_28*chi),
                           true_unconditioned %>% filter(year == 2022, temp == '29°C') %>% mutate(chi = prob_2022_29*chi),
                           true_unconditioned %>% filter(year == 2022, temp == '30°C') %>% mutate(chi = prob_2022_30*chi))

true_data = rbind(chi_true %>%
                    mutate(lab = "conditioned"),
                  true_unconditioned %>%
                    mutate(lab = "unconditioned"))

max_distance <- max(chi_true$distance, na.rm = TRUE)

chi_true_summary = true_data %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, max_distance, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(year, temp, chi,lab, distance = median(distance))  %>%
  group_by(distance, year, temp,lab) %>%
  drop_na() %>%
  summarise(mn = median(chi))

#don't understand what dat_for_rat isused for
dat_for_rat = chi_true_summary %>%
  ungroup() %>%
  filter(distance < 1.16*4/max_distance, distance > 1*4/max_distance) %>%
  filter(lab == 'unconditioned')  %>%
  dplyr::select(year, mn, temp) 

dat_for_rat[dat_for_rat$year == 2022 & dat_for_rat$temp == "28°C",]$mn/dat_for_rat[dat_for_rat$year == 1971 & dat_for_rat$temp == "28°C",]$mn
dat_for_rat[dat_for_rat$year == 2022 & dat_for_rat$temp == "29°C",]$mn/dat_for_rat[dat_for_rat$year == 1971 & dat_for_rat$temp == "29°C",]$mn
dat_for_rat[dat_for_rat$year == 2022 & dat_for_rat$temp == "30°C",]$mn/dat_for_rat[dat_for_rat$year == 1971 & dat_for_rat$temp == "30°C",]$mn

#corresponds to figure 7 in the paper

ggplot(data = chi_true_summary, 
              aes(x = distance, y = mn, col = year, linetype = year)) +
  geom_smooth(se = FALSE) + 
  facet_grid(lab ~ temp, scales = 'free') +
  labs(x = "Distance (km)",
       y = expression(chi[o]),
       col = "Year",
       shape = "Year",
       linetype = "Year") +
  scale_x_continuous(limits = c(0, 375)) +  # Use scale_x_continuous instead of xlim()
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.y = element_text(size = 0))
