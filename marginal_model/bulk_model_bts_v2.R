gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

source('marginal_model/gpd_models.R') 

fit_quant_reg = T
calc_lam = T

fit_quant_regression = function(bts_range, marg_mod, num_quantiles, zeta_list, obs_data){
  
  quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
  
  # # ------ Bootstrap quant reg models to get CI for each parameter
  glob_anom_df = obs_data %>%
    dplyr::select(year, glob_anom) %>%
    unique()
  
  res = c()
  dat = c()
  for(file_name in bts_range){
    if(marg_mod == "mod_0"){
      dat <- read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_", file_name, ".csv"), header = FALSE) %>%
        setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3")) %>%
        rename(maxtp = maxtp_0) %>%
        select(-c("maxtp_1", "maxtp_2", "maxtp_3"))
    }
    
    if(marg_mod == "mod_1"){
      dat <- read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_", file_name, ".csv"), header = FALSE) %>%
        setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3")) %>%
        rename(maxtp = maxtp_1) %>%
        select(-c("maxtp_0", "maxtp_2", "maxtp_3"))
    }
    
    if(marg_mod == "mod_2"){
      dat <- read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_", file_name, ".csv"), header = FALSE) %>%
        setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3")) %>%
        rename(maxtp = maxtp_2)%>%
        select(-c("maxtp_0", "maxtp_1", "maxtp_3"))
    }
    
    if(marg_mod == "mod_3"){
      dat <- read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_", file_name, ".csv"), header = FALSE) %>%
        setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3")) %>%
        rename(maxtp = maxtp_3) %>%
        select(-c("maxtp_0", "maxtp_1", "maxtp_2"))
    }
    
    dat$date = as.Date(dat$date)
    
    dat = dat %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(glob_anom_df, by ="year")
    
    #print(paste0("Bootstrap number ", file_name))
    for(q in seq_along(quantiles_to_estimate_bulk)){
      
      zeta = zeta_list[q] # quantile to estimate
      
      #print(paste0("quantile ", zeta))
      
      # spatial covariate [q^tau_c(s)]
      obs_data$thisval = obs_data$value %>% lapply(`[[`, q) %>% unlist
      
      dat = dat %>%
        dplyr::select(date,stn,maxtp,year,glob_anom) %>%
        left_join(obs_data %>% dplyr::select(stn, value = thisval) %>% unique, by = "stn") %>% drop_na()
      
      #issue when the hessian is positive negative definite. Enforcing regularisation conditions
      if(q == 1){
        quantile_model_fit <- evgam(maxtp ~ value + glob_anom, dat,
                                    family = "ald", ald.args = list(tau = zeta))
        
        new_coeff = c(quantile_model_fit$location$coefficients[1], quantile_model_fit$location$coefficients[2], quantile_model_fit$location$coefficients[3])
      }
      
      if (q>1){
        quantile_model_fit <- evgam(maxtp ~ value + glob_anom, dat,
                                    family = "ald", inits = prev_coeff ,ald.args = list(tau = zeta))
        
        new_coeff = c(quantile_model_fit$location$coefficients[1], quantile_model_fit$location$coefficients[2], quantile_model_fit$location$coefficients[3])
        
        if (max(abs(prev_coeff-new_coeff))>2){
          new_coeff = prev_coeff
        }
      }
      
      # Save the fitted parameter estimates for each quantile
      result = tibble(bts = file_name,
                      tau = zeta,
                      beta_0 = new_coeff[1],
                      beta_1 = new_coeff[2],
                      beta_2 = new_coeff[3]) %>%
        write_csv(paste0("output/bts_quant_reg/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"), append = T)
      prev_coeff = new_coeff
    }
  }
}

marg_mod = "mod_2"

num_quantiles = 30 

obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

zeta_list = obs_data$quantile[[1]]

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  select(-c("id", "clim_thresh_value_9", "threshold_9", "quantile"))

rm('glob_anomaly', 'glob_anomaly_reshaped')

#then  seq(1, 32, 8) seq(33, 64, 8) seq(72, 103, 8)
#each of them take 50 minutes
for(i in seq(72, 103, 8)){ 
  seq_ = seq(i, i+7)
  job::job({fit_quant_regression(seq_,  'mod_2', 30, zeta_list, obs_data)}, import = c("fit_quant_regression", "obs_data", "zeta_list", "seq_")) 
}

calc_lambda_bts = function(bts_range, marg_mod, num_quantiles, obs_data, zeta_list){
  
  quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
  
  # # # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year, glob_anom) %>%
    unique() %>%
    arrange(year)
  
  #already has colnames
  bootstrapped_models = read.csv(paste0("output/bts_quant_reg/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"), header = FALSE)
  colnames(bootstrapped_models) <- c("bts", "tau", "beta_0", "beta_1", "beta_2")
  bootstrapped_models = unique(bootstrapped_models) #removing dupplicates
  
  for(bts_num in bts_range){
    print(paste0("Calculating lambda for bts num ", bts_num))
    this_est_of_qnt_mods = bootstrapped_models %>%
      filter(bts == bts_num)
    
    # # --- creates a tibble with each station and its quantile model
    obs_smoothed_quantiles = obs_data %>%
      group_by(stn) %>%
      group_map(~{
        
        # --- get the climate quantile estimates closest to current station
        clim_vals = obs_data %>%
          filter(stn == .x$stn[1]) %>%
          dplyr::select(value) %>%
          unique() %>%
          pull(value) %>%
          unlist()
        
        # predict quantile for each year and site
        quant_reg_pars = this_est_of_qnt_mods %>%
          arrange(tau)
        
        res = c()
        for(q in seq_along(quantiles_to_estimate_bulk)){
          
          qpars = quant_reg_pars[q,]
          res = rbind(res,
                      tibble(quantile =  quantiles_to_estimate_bulk[q],
                             year = temporal_covariates$year,
                             glob_anom = temporal_covariates$glob_anom,
                             quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + qpars$beta_2*temporal_covariates$glob_anom))
        }
        
        # interpolate quantiles over tau for each year
        res %>%
          group_by(year) %>%
          group_map(~{
            tibble(year = .x$year[1],
                   quant_spline = list(splinefun(.x$quant_value, .x$quantile,  method = 'monoH.FC')))
            # tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
            # temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
          }, .keep = T) %>%
          plyr::rbind.fill() %>%
          as_tibble() %>%
          mutate(stn = .x$stn[1])
        
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    # save this 
    obs_smoothed_quantiles %>% saveRDS(paste0("output/bulk_model_fits/", marg_mod, "_num_quantiles_",num_quantiles, "_bts_", bts_num))
    
    # look at threshold exceedance for this model
    obs_data %>%
      group_by(stn) %>%
      group_map(~{
        
        thresh_exceedance = obs_smoothed_quantiles %>%
          filter(stn == .x$stn[1]) %>%
          pull(quant_spline) %>%
          sapply(function(x) sapply(.x$threshold_9[1], x))
        
        tibble(bts_num = bts_num,
               stn = .x$stn[1],
               year = temporal_covariates$year,
               thresh_exceedance = 1-thresh_exceedance)  %>%
          write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts_num,".csv"), append = T)
      }, .keep = T)
  }
}


num_quantiles = 30 

obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

zeta_list = obs_data$quantile[[1]]

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  select(-c("id", "clim_thresh_value_9", "quantile")) #keeping threshold_9

rm('glob_anomaly', 'glob_anomaly_reshaped')


#seq(1, 32, 8) seq(33, 64, 8) seq(72, 103, 8)
#each take 15 min
for (i in seq(33, 64, 8)){
  
  seq_ = seq(i,i+7)
  
  job::job({calc_lambda_bts(seq_,  'mod_2', 30, obs_data, zeta_list)} , import = c("calc_lambda_bts", "obs_data", "zeta_list", "seq_"))
  
}



#bootstrap fitting the gpd model 

bts_files = list.files("data/processed/bootstrap_data/bts_under_gpd_models/")

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)
  
  # ----- fit all models
  for(file_name in bts_files){
    
    bts_value <- sub(".*_bts_(\\d+)\\.csv", "\\1", file_name)
    bts_value <- as.numeric(bts_value)
    
    if (bts_value < 65 | bts_value > 71) {next}
    
    print(paste0("fitting to bootstrap ", file_name))
    dat = read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models/", file_name), header = FALSE)%>%
      setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3")) %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
      rename(altitude = Altitude.m.) %>%
      left_join(glob_anomaly_reshaped, by = "year")
    
    # --- remove observations poorly interpolated by quantile model (i.e. <0)
    dat = dat %>% select(-c(maxtp_0, maxtp_1)) %>% filter(maxtp_2 > 0)
    
    obs_data_to_pred = dat %>%
      mutate(excess = maxtp_2 - threshold_9) %>%
      filter(excess > 0)
    this_fit_mod_2 = fit_mod_2(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
                               obs_data_to_pred$glob_anom,
                               obs_data_to_pred$altitude,
                               initial_pars = c(0.03, -0.03,  0.08,  1, -0.08, -0.2))
    
    bts_value <- sub(".*_bts_(\\d+)\\.csv", "\\1", file_name)
    bts_value <- as.numeric(bts_value)
    
    c(bts_value, this_fit_mod_2)  %>% matrix() %>% t() %>% as.data.frame() %>%
      write_csv("output/gpd_model_fits/bts/model_2_bts.csv", append = T)
  }













# ------------------- PLOTS 

marg_mod = 'mod_2'
num_quantiles = 30

true = read_csv(paste0("Data/processed/glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

files = list.files("output/bts_thresh_ex_lambda/")
files = files[grepl(marg_mod, files) & grepl(paste0("num_quantiles_",num_quantiles), files)]

all_bts_dat = c()
for(f in files){
  print(f)
  all_bts_dat = rbind(all_bts_dat,read_csv(paste0("output/bts_thresh_ex_lambda/", f),
                                           col_names = c("bts", "stn", "year", "thresh_exceedance")))
}


all_bts_dat = rbind(all_bts_dat,read_csv(paste0("output/bts_thresh_ex_lambda/", f),
                                         col_names = c("bts", "stn", "year", "thresh_exceedance")))

true_res = read_csv(paste0("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv")) %>%
  group_by(year) %>%
  summarise(thresh_exceedance = mean(thresh_exceedance_9))

res = read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles,".csv"), skip = 1,
               col_names = c('bts', 'tau', 'b0', 'b1', 'b2')) 

# --- read in fitted quantile regression coefficients
true_pars = read_csv(paste0("Data/processed/glob_anomaly_quantile_model_fit_pars_num_quantiles_", num_quantiles,".csv"),
                     col_names = c('tau', 'beta_0', 'beta_1', 'beta_2')) #need to check if need to rename columns 


plts=gridExtra::grid.arrange(res %>%
                               group_by(tau) %>%
                               summarise(upper = quantile(b2, 0.975),
                                         lower = quantile(b2, 0.025)) %>% #used to be 0.030
                               ggplot()+
                               geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
                               geom_line(data = true_pars, aes(tau, beta_2))+
                               xlim(c(0.08, 0.9))+
                               ylim(c(-1.5, 3))+
                               theme_minimal(12)+
                               labs(y = expression(beta[3]),
                                    x = expression(tau))+
                               theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                             #second plot
                             all_bts_dat %>% 
                               group_by(bts, year) %>%
                               summarise(thresh_exceedance = mean(thresh_exceedance))  %>%
                               group_by(year) %>%
                               summarise(upper = quantile(thresh_exceedance, 0.975),
                                         lower = quantile(thresh_exceedance, 0.025)) %>%
                               ungroup() %>%
                               ggplot()+
                               geom_ribbon(aes(year, ymin = lower, ymax = upper), alpha = 0.3) +
                               geom_line(data = true_res, aes(year, thresh_exceedance))+
                               theme_minimal(12)+
                               labs(y = expression(lambda[t]),
                                    x = "Year")+
                               theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
                               xlim(c(1971, 2022))+
                               ylim(c(0.05, 0.19)), nrow = 1)