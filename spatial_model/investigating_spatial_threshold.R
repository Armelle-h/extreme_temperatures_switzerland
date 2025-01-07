rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

marg_mod = "mod_1"

robust = TRUE 

if (robust == TRUE){
  true_folder = "true_robust"
} else{
  true_folder = "true"
}

nu_val = 0.15
nu_name = "015"

#getting the observed data on frechet scale

pareto_comparison = function (threshold_quant){

legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")

num_quantiles = 30

source('marginal_model/gpd_models.R')

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

threshold_9_df = vroom::vroom("Data/processed/plain_1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, stn)%>%
  unique()

obs_data = vroom::vroom("Data/Observed_data/plain_obs_data_gpd_model.csv") %>%
  mutate(year = lubridate::year(date))%>%
  left_join(vroom::vroom(paste0("Data/processed/plain_glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv")))  %>%
  left_join(glob_anomaly_reshaped, by = "year")%>%
  left_join(threshold_9_df, by="stn")%>%
  filter(!(stn %in% c("WSLBTB", "WSLHOB")))

obs_smoothed_quantiles = readRDS(paste0("output/plain_glob_anomaly_quant_models_num_quantiles_", num_quantiles, ".csv"))%>%
  filter(!(stn %in% c("WSLBTB", "WSLHOB")))

# pull date which have 'conditioning' site
#doesn't change anything as KOP has observation for the entirety of the temporal range
dates_to_keep = obs_data %>%
  group_by(date) %>%
  summarise(check_obs_for_s = 'KOP' %in% stn) %>%   #defining the temporal range such that we have observed data for 'to complete' to simplify the fitting
  filter(check_obs_for_s) %>%
  pull(date)

obs_data = obs_data %>% filter(date %in% dates_to_keep) 

if(marg_mod == 'mod_0'){
  this_fit_mod = read_csv("output/gpd_model_fits/model_0_true.csv") %>%
    unlist() %>% as.numeric
  pred <<- my_predict_0(this_fit_mod, obs_data$scale_9)
  
}else if(marg_mod == 'mod_1'){
  
  this_fit_mod = read_csv("output/gpd_model_fits/model_1_true.csv") %>%
    unlist() %>% as.numeric
  pred <<- my_predict_1(this_fit_mod, obs_data$scale_9, obs_data$glob_anom)
  
}

obs_data$scale = pred$scale
obs_data$shape = pred$shape

#standardising the observation to a unit Pareto distribution or to a frechet scale
obs_data_standardised = obs_data %>%
  group_by(stn, year) %>%
  group_map(~{
    data = .x$maxtp
    threshold = .x$threshold_9
    res = rep(NA, length(data))
    num_extremes = sum(data>threshold)
    
    if(num_extremes >0){ # if extreme obs in this year at this site
      scle = .x$scale
      shpe = .x$shape
      my_lambda = .x$thresh_exceedance_9
      
      res[data > threshold] = 1 - my_lambda[data > threshold]*(1-evd::pgpd((data[data > threshold] - threshold[data > threshold]), loc = 0,scale = scle[data > threshold], shape = shpe[1]))
    }
    
    this_quant_mod = obs_smoothed_quantiles %>%
      filter(stn == .x$stn[1],
             year == .x$year[1]) %>%
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

extreme_dates = obs_data_standardised %>%
  dplyr::group_by(date) %>%
  #dplyr::summarise(cost = mean(frechet_marg)) %>%
  dplyr::summarise(cost = median(pareto_marg)) %>%  #dplyr::summarise(cost = mean(pareto_marg)) %>%
  ungroup() %>%
  arrange(desc(cost))

spatial_threshold = quantile(extreme_dates$cost, threshold_quant) %>% as.numeric


# get extreme dates
extreme_dates = extreme_dates %>%
  filter(cost > spatial_threshold) %>%
  arrange(desc(cost))

# temporally decluster events
my_data_tmp = extreme_dates %>% arrange(date)

min_dif = my_data_tmp %>% arrange(date) %>% pull(date) %>% diff() %>% min

print(my_data_tmp)
print("declustering events")

#removes events within 7 days of each other and retain only the highest cost event in a cluster
while(min_dif<7){
  i<-1
  while(i < nrow(my_data_tmp)){
    if(diff(c(my_data_tmp[i,]$date, my_data_tmp[i+1,]$date))<7){
      # pick the largest
      if(my_data_tmp[i,]$cost > my_data_tmp[i+1,]$cost){
        my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i+1,]$date,]
      }else{
        my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i,]$date,]
      }
    }
    i<-i+1
  }
  min_dif <- my_data_tmp$date %>% diff() %>% min %>% abs()
}
extreme_dates = my_data_tmp

# tibble of extreme events
exceedances = obs_data_standardised %>%
  filter(date %in% extreme_dates$date) %>%
  arrange(date)

#to remove was empty

obs_pareto = exceedances %>% filter(stn == "KOP") %>% select(pareto_marg)  %>% mutate(pareto_marg = pareto_marg/spatial_threshold) 

read_simulated = readRDS(paste0("output/simulations/simulations_on_obs_grid/", true_folder,"/nu_", nu_name,"/nu_", nu_name,"_",marg_mod,"_run_1"))

simulated_pareto = sapply(read_simulated, function(x) x[[1]])

return (list(obs_pareto$pareto_marg, simulated_pareto))

}

threshold_quant = 0.8

result = pareto_comparison(threshold_quant)
  
obs_pareto = result[[1]]  
simulated_pareto = result[[2]]


qqplot(obs_pareto, simulated_pareto, 
       xlab = "Quantiles of obs pareto", ylab = "Quantiles of simulated pareto",
       xlim = c(0, quantile(obs_pareto, 0.99)),
       ylim = c(0, quantile(simulated_pareto, 0.99)))
abline(0, 1, col = "red", lwd = 2)

ks.test(obs_pareto, simulated_pareto)


length(obs_pareto)
