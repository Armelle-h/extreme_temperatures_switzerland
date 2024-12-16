rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(data.table)

filter_id = 25 #else too computationally expensive

num_quantiles = 30

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

clim_1 = fread("Data/Climate_data/By_year/1971_1990_JJA_climate_data.csv")%>%
  filter(id %in% I_plain$id)%>%
  filter(id %% 25 == 0)

clim_2 = fread("Data/Climate_data/By_year/1991_2017_JJA_climate_data.csv")%>%
  filter(id %in% I_plain$id)%>%
  filter(id %% 25 == 0)

clim_3 = fread("Data/Climate_data/By_year/2018_2022_JJA_climate_data.csv")%>%
  filter(id %in% I_plain$id)%>%
  filter(id %% 25 == 0)

clim_data = rbind(clim_1, clim_2, clim_3)

rm(clim_1, clim_2, clim_3)
gc()

clim_thresh = clim_data %>%
  group_by(id) %>%
  summarise(clim_thresh_value_9 = quantile(maxtp, 0.9, na.rm = TRUE))
#first converting clim_data to pareto

clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  filter(id %% 25 == 0)%>%
  filter(id %in% I_plain$id)

clim_grid = read_csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")%>%
  filter(id %% 25 == 0)

clim_smooth = clim_quantiles_subset %>%
  group_by(id) %>%
  group_map(~{
    tibble( id = .x$id[[1]],
            tau_to_temp = list(splinefun(unlist(.x$quantile),unlist(.x$value),  method = 'monoH.FC')),
            temp_to_tau = list(splinefun(unlist(.x$value),unlist(.x$quantile),  method = 'monoH.FC')))
  }, .keep = T)%>%
  plyr::rbind.fill() %>%
  as_tibble() 

# Calculate lambda (exceedance probability) for the threshold
# (difference compared to before is that now we're working with the regression estimated climate quantile)
lambda_thresh_ex = clim_thresh %>%
  group_by(id) %>%
  group_map(~{
    
    thresh_exceedance_9 = clim_smooth%>%
      filter(id == .x$id[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$clim_thresh_value_9[1], x))# Apply the threshold function
    
    # Create a tibble with exceedance probabilities
    tibble(id = .x$id[1],
           thresh_exceedance_9 = 1-thresh_exceedance_9)# Inverse exceedance probability  
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

shape = readRDS("Data/clim_data_gpd_model/plain_optimal_shape.rds")

I_date = read.csv("Data/id_date_correspondance.csv")%>%
  mutate(year = lubridate::year(date))

clim_data = clim_data %>%
  left_join(I_date, by = "date_id")%>%
  select(-date_id)%>%
  left_join(lambda_thresh_ex, by = "id")%>%
  left_join(clim_grid, by = "id")%>%
  left_join(clim_thresh, by = "id")%>%
  filter(year %in% c(1971, 2022))

rm(lambda_thresh_ex, clim_thresh, clim_grid)

clim_data_standardised = clim_data %>%  #will have to do a for loop over the clim data files
  group_by(id, year) %>% #the function is just defined with respect to space, not time
  group_map(~{
    data = .x$maxtp
    threshold = .x$clim_thresh_value_9
    res = rep(NA, length(data))
    num_extremes = sum(data>threshold)
    
    if(num_extremes >0){ # if extreme obs in this year at this site
      scle = .x$scale_9
      shpe = shape
      my_lambda = .x$thresh_exceedance_9
      
      res[data > threshold] = 1 - my_lambda[data > threshold]*(1-evd::pgpd((data[data > threshold] - threshold[data > threshold]), loc = 0,scale = scle[data > threshold], shape = shpe[1]))
    }
    
    this_quant_mod = clim_smooth %>%
      filter(id == .x$id[[1]]) %>%
      pull(temp_to_tau) %>%
      .[[1]]
    
    res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
    res[res<0]=0
    res[res>0.999999]=0.999999
    
    .x$unif = res
    #.x$frechet_marg = -1/log(.x$unif)
    .x$pareto_marg = 1/(1-.x$unif)
    .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble() 

#the qq plot ! 
  
clim_data_standardised$ideal = (clim_data_standardised$maxtp %>% rank(ties.method = "random"))/(nrow(clim_data_standardised)+1) 

library(grid)

gridExtra::grid.arrange(clim_data_standardised %>%
                          drop_na() %>%
                          ggplot()+
                          geom_point(aes(sort(ideal), sort(unif)), size = 0.75)+
                          geom_abline(col = 'red',linetype = 'longdash')+
                          geom_vline(xintercept = 0.9)+
                          theme_minimal(12)+
                          theme(panel.grid.minor = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank())
)


rm(clim_smooth)

calc_clim_chi_true = function(yr, tmp, pareto_val){
  
  set.seed(123456)
  #number of site pairs to be sampled
  num_samples = 1000 #turn into 8000 later
  
  sites = read_csv("Data/processed/plain_clim_pairs_with_dist.csv") %>% sample_n(num_samples)
  
  clim_sites = read.csv("Data/plain_id_lon_lat_correspondance.csv")%>%
    select("id", "longitude_proj", "latitude_proj")%>%
    filter(id %% 25 ==0)%>% #else way too heavy
    unique()
  
  #part below will change because we do an empirical estimate and not a simulated one
  
  pareto_val_pulled = pareto_val %>% 
    pull(pareto_marg)
    
  for(s in seq(nrow(sites))){ 
      
      if((s %% 100) == 0){
        print(s)
      }
      
      id1 = sites[s,]$V1
      id2 = sites[s,]$V2

      #used to have my simulation, not sure if it's still valid
      id1_exceeds = (unlist(lapply(pareto_val, "[[", which(pareto_val$id == id1))) > pareto_val_pulled[which(pareto_val$id == id1)])
      id2_exceeds = (unlist(lapply(pareto_val, "[[", which(pareto_val$id == id2))) > pareto_val_pulled[which(pareto_val$id == id2)])
      
      tibble(sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
        write_csv(paste0("output/simulations/clim_simulation_summary/chi_data_scale_clim_grid_model_yr_",yr, "_conditioned_on_",tmp,".csv"),append = T)
  }
}

#20 minutes for 1000 simulations. A bit more than 2 hours for 8000 !!!


pareto_val_2022_28 = clim_data_standardised %>% 
  filter(pareto_marg > 28, year == 2022)%>%
  select(id, pareto_marg)

pareto_val_2022_29 = clim_data_standardised %>% 
  filter(pareto_marg > 29, year == 2022)%>%
  select(id, pareto_marg)

pareto_val_2022_30 = clim_data_standardised %>% 
  filter(pareto_marg > 30, year == 2022)%>%
  select(id, pareto_marg)

pareto_val_1971_28 = clim_data_standardised %>% 
  filter(pareto_marg > 28, year == 1971)%>%
  select(id, pareto_marg)

pareto_val_1971_29 = clim_data_standardised %>% 
  filter(pareto_marg > 29, year == 1971)%>%
  select(id, pareto_marg)

pareto_val_1971_30 = clim_data_standardised %>% 
  filter(pareto_marg > 30, year == 1971)%>%
  select(id, pareto_marg)

#rm(clim_data_standardised)

job::job({
  calc_clim_chi_true(2022, 28, pareto_val_2022_28)
  
}, import = c("pareto_val_2022_28", "calc_clim_chi_true"))

job::job({
  calc_clim_chi_true(2022, 29, pareto_val_2022_29)
  
}, import = c("pareto_val_2022_29", "calc_clim_chi_true"))

job::job({
  calc_clim_chi_true(2022, 30, pareto_val_2022_30)
  
}, import = c("pareto_val_2022_30", "calc_clim_chi_true"))

job::job({
  calc_clim_chi_true(1971, 28, pareto_val_1971_28)
  
}, import = c("pareto_val_1971_28", "calc_clim_chi_true"))

job::job({
  calc_clim_chi_true(1971, 29, pareto_val_1971_29)
  
}, import = c("pareto_val_1971_29", "calc_clim_chi_true"))

job::job({
  calc_clim_chi_true(1971, 30, pareto_val_1971_30)
  
}, import = c("pareto_val_1971_30", "calc_clim_chi_true"))


























