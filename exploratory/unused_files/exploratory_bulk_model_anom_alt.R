gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

sigmoid= function(altitude, alt_threshold){
  result = 1/(1+exp(-200*(altitude-alt_threshold)))
  
  return(result)
}

num_quantiles = 30

obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
quantiles_to_estimate_bulk = obs_data$quantile[[1]]
index_clim_quantile = seq_along(quantiles_to_estimate_bulk)

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>%
  rename(altitude = Altitude.m.)

#joining obs data with the altitude 

obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, altitude), by="stn")

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

quant_reg_model_pars = read_csv(paste0("Data/processed/glob_anom_indicator_alt_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3'))

temporal_covariates = obs_data %>%
  dplyr::select(year, glob_anom) %>% 
  unique() %>%
  arrange(year)


temporal_covariates_year = temporal_covariates$year
temporal_covariates_glob_anom = temporal_covariates$glob_anom

obs_smoothed_quantiles = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    # --- get the climate quantile estimates closest to the current station
    clim_vals <<- obs_data %>%
      filter(stn == .x$stn[1]) %>%
      dplyr::select(quantile, value) %>%
      unique() %>%
      pull(value) %>%
      unlist()
    
    # predict (observed) quantile for each year and site
    quant_reg_pars = quant_reg_model_pars %>%
      arrange(tau)
    
    alt_df = legend_data%>%filter(stn == .x$stn[1]) %>%select(altitude)
    alt = alt_df$altitude
    
    res = c()
    for(q in index_clim_quantile){
      qpars = quant_reg_pars[q,]
      
      # Calculate the predicted quantile value using the regression parameters
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates_year,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates_glob_anom)
                         + (qpars$beta_3) *sigmoid(alt, 1500),
                         stn = .x$stn[1]))
    }
    # Return the reshaped result with the desired columns
    res
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


empirical_quantile <- obs_data %>%
  select(c("stn", "year", "maxtp")) %>%
  group_by(stn, year) %>%
  reframe(
    quantiles = map_df(
      sort(unique(obs_smoothed_quantiles$quantile)), #need to do this trick otherwise don't recognize the quantiles as the same (although they are)
      ~ tibble(
        quantile = .x,
        quant_value = quantile(maxtp, probs = .x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) %>%
  unnest(cols = quantiles)


compare_quantiles = empirical_quantile %>%
  select(-.groups) %>%
  left_join(obs_smoothed_quantiles, by=c("stn", "year", "quantile"), suffix = c("_emp", "_est")) %>%
  mutate(diff_quant = abs(quant_value_emp-quant_value_est))

print(sum(is.na(compare_quantiles$diff_quant)))

C = compare_quantiles %>%
  filter(!(quantile %in% c(0.001, 0.99)))





C = C %>% left_join(legend_data%>% select(c("stn", "altitude")), by = "stn" )

