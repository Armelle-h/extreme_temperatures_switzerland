gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


num_quantiles = 30

obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
quantiles_to_estimate_bulk = obs_data$quantile[[1]]
index_clim_quantile = seq_along(quantiles_to_estimate_bulk)

L =read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>% 
  select(c("stn", "Altitude.m."))%>%
  rename(altitude = Altitude.m.)
  
stations_below_2000 = L%>% 
  filter(altitude<2000)%>% 
  select(stn)


glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

# ---- get covariates for prediction
temporal_covariates = obs_data %>%  #issue, be wary of, glob anom is defined for June July and August so each year is associated with 3 different values
  dplyr::select(year, glob_anom) %>% 
  unique() %>%
  arrange(year)

quant_reg_model_pars_1 = read_csv(paste0("Data/processed/debug_1_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))
quant_reg_model_pars_2 = read_csv(paste0("Data/processed/debug_2_glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                  col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

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
    
    if (.x$stn[1] %in% stations_below_2000$stn){
      quant_reg_pars = quant_reg_model_pars_1 %>%
        arrange(tau)
      
    } else {
      quant_reg_pars = quant_reg_model_pars_2 %>%
        arrange(tau)      
    }
    
    res = c()
    for(q in index_clim_quantile){
      qpars = quant_reg_pars[q,]
      # Calculate the predicted quantile value using the regression parameters
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates$year,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates$glob_anom),
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
  unnest(cols = quantiles)%>%
  unique()

compare_quantiles = empirical_quantile %>%
  select(-.groups) %>%
  left_join(obs_smoothed_quantiles, by=c("stn", "year", "quantile"), suffix = c("_emp", "_est")) %>%
  mutate(diff_quant = abs(quant_value_emp-quant_value_est))


L = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")
L_alt = L %>% select(c("stn", "Altitude.m."))
C = C %>% left_join(L_alt, by = "stn")

print(sum(is.na(C$diff_quant)))


obs_data_orig = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")%>%
  select(stn)

count_stn = obs_data_orig %>%
  count(stn)

C_fin = C %>%
  left_join(count_stn, by = c("stn"))

plot(C_fin$n, C_fin$diff_quant)