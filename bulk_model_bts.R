gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

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
    
    
    dat = dat %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(glob_anom_df)
    
    print(paste0("Bootstrap number ", file_name))
    for(q in seq_along(quantiles_to_estimate_bulk)){
      
      zeta = zeta_list[q] # quantile to estimate
      
      print(paste0("quantile ", zeta))
      
      # spatial covariate [q^tau_c(s)]
      obs_data$thisval = obs_data$value %>% lapply(`[[`, q) %>% unlist
      
      dat = dat %>%
        dplyr::select(date,stn,maxtp,year,glob_anom) %>%
        left_join(obs_data %>% dplyr::select(stn, value = thisval) %>% unique) %>% drop_na()
      
      
      quantile_model_fit <- evgam(maxtp ~ value + glob_anom, dat,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # save parameter estimates w CI
      tibble(bts = file_name,
             tau = zeta,
             beta_0 = quantile_model_fit$location$coefficients[1],
             beta_1 = quantile_model_fit$location$coefficients[2],
             beta_2 = quantile_model_fit$location$coefficients[3])%>%
        write_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"), append = T)
    }
  }
}

num_quantiles = 40 

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

#takes   time 

for(i in seq(15, 20, 2)){ #should take ~20 minutes
  seq_ = seq(i, i+1)
  job::job({fit_quant_regression(seq_,  'mod_0', 40, zeta_list, obs_data)}, import = c("fit_quant_regression", "obs_data", "zeta_list", "seq_")) 
}

job::job({fit_quant_regression(seq(1,1),  'mod_0', 40)}, import = c("fit_quant_regression", "obs_data")) 
job::job({fit_quant_regression(seq(1,20),  'mod_1', 40)}) 
job::job({fit_quant_regression(seq(1,20),  'mod_2', 40)}) 
job::job({fit_quant_regression(seq(1,20),  'mod_3', 40)}) 



# job::job({fit_quant_regression(seq(1,20),  'mod_0', 30)}) 
# job::job({fit_quant_regression(seq(21,40), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(41,60), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(61,80), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(81,100),'mod_0', 30)})
# job::job({fit_quant_regression(seq(101,120), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(121,140), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(141,160), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(161,180), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(181,200), 'mod_0', 30)})

# job::job({fit_quant_regression(seq(1,20),  'mod_1', 30)}) 
# job::job({fit_quant_regression(seq(21,40), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(41,60), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(61,80), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(81,100),'mod_1', 30)})
# job::job({fit_quant_regression(seq(101,120), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(121,140), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(141,160), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(161,180), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(181,200), 'mod_1', 30)})

# job::job({fit_quant_regression(seq(1,20),    'mod_2', 30)}) 
# job::job({fit_quant_regression(seq(21,40),   'mod_2', 30)})
# job::job({fit_quant_regression(seq(41,60),   'mod_2', 30)})
# job::job({fit_quant_regression(seq(61,80),   'mod_2', 30)})
# job::job({fit_quant_regression(seq(81,100),  'mod_2', 30)})
# job::job({fit_quant_regression(seq(101,120), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(121,140), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(141,160), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(161,180), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(181,200), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(201,220), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(221,240), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(241,260), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(261,280), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(281,300), 'mod_2', 30)})



calc_lambda_bts = function(bts_range, marg_mod, num_quantiles, obs_data, zeta_list){
  
  quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
  
  # # # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year, glob_anom) %>%
    unique() %>%
    arrange(year)
  
  #already has colnames
  bootstrapped_models = read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"))
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


num_quantiles = 40 

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


#Takes <20 minutes
for (i in seq(1,20, 5)){ #20 dividable by 5 so no need for extra precautions
  
  seq_ = seq(i,i+4)
  
  job::job({calc_lambda_bts(seq_,  'mod_0', 40, obs_data, zeta_list)} , import = c("calc_lambda_bts", "obs_data", "zeta_list", "seq_"))
  
}


job::job({calc_lambda_bts(seq(1,20),  'mod_0', 40)})
job::job({calc_lambda_bts(seq(1,20),  'mod_1', 40)})
job::job({calc_lambda_bts(seq(1,20),  'mod_2', 40)})
job::job({calc_lambda_bts(seq(1,20),  'mod_3', 40)})


# job::job({calc_lambda_bts(seq(1,100),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(101,200),  'mod_0', 30)})
# job::job({calc_lambda_bts(seq(201,300),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(301,400),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(401,500),  'mod_0', 30)})

# job::job({calc_lambda_bts(seq(1,100),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(21,40), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(41,60), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(61,80), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(81,100),'mod_0', 30)})
# job::job({calc_lambda_bts(seq(101,120), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(121,140), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(141,160), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(161,180), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(181,200), 'mod_0', 30)})

# job::job({calc_lambda_bts(seq(1,100),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(101,200),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(201,300),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(301,400),  'mod_1', 30)})
# job::job({calc_lambda_bts(seq(401,500),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(1,20),  'mod_1', 30)})
# job::job({calc_lambda_bts(seq(21,40), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(41,60), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(61,80), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(81,100),'mod_1', 30)})
# job::job({calc_lambda_bts(seq(101,120), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(121,140), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(141,160), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(161,180), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(181,200), 'mod_1', 30)})

# job::job({calc_lambda_bts(seq(1,60),    'mod_2', 30)})
# job::job({calc_lambda_bts(seq(61, 120),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(121,180),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(181,240),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(241,300),  'mod_2', 30)})
# job::job({calc_lambda_bts(seq(301,360),    'mod_2', 30)})
# job::job({calc_lambda_bts(seq(361, 420),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(421,480),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(481,500),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(201,220), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(221,240), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(241,260), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(261,280), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(281,300), 'mod_2', 30)})


# --- on clim grid
get_bulk_bts_on_clim_grid = function(bts_range, marg_mod, num_quantiles, temporal_covariates){
  
  clim_grid = read_csv("Data/Climate_data/clim_scale_grid_gpd_model.csv")%>%
    filter(id %% 10 == 0)
    
  # estimate empiricle quantiles for climate data, was already done and saved in a previous file
  clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))  %>%
    filter(id %% 10 == 0)
    
  #done on one index out of 10 because otherwise too memory expensive (!) and too slow(?)
    
  for(this_bts in bts_range){
    
    quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
    
    #columns are already named
    quant_reg_model_pars =  read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv")) %>%
      filter(bts == this_bts)
    
    clim_grid = clim_grid %>%
      left_join(clim_quantiles_subset)
    
    clim_date_w_quantile_mod = c()
    
    
    # here
    for(i in seq(nrow(clim_grid))){
      cat(i, "over", nrow(clim_grid))
      
      this_data_for_pred = tibble(year = c(1971, 2022)) %>% #was (1960, 1971, 2022, 2024) but was unsured how to handle 1960 and 2024 so changed to data for which I do have results
        left_join(temporal_covariates) %>%
        mutate(id = clim_grid[i,]$id,
               # threshold = clim_grid[i,]$threshold,
               # clim_scale = clim_grid[i,]$clim_scale,
               quantile = clim_grid[i,]$quantile,
               value = clim_grid[i,]$value) 
      
      # predict quantile for each year and site
      quant_reg_pars = quant_reg_model_pars %>%
        arrange(tau)
      
      clim_vals = clim_grid[i,]$value[[1]]
      res = c()
      for(q in seq_along(quantiles_to_estimate_bulk)){
        qpars = quant_reg_pars[q,]
        res = rbind(res,
                    tibble(quantile =  qpars$tau,
                           year = this_data_for_pred$year,
                           quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(this_data_for_pred$glob_anom)))
      }
      
      res = res %>%
        group_by(year) %>%
        group_map(~{
          tibble(year = .x$year[1],
                 tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
                 temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
        }, .keep = T) %>%
        plyr::rbind.fill() %>%
        as_tibble() 
      
      clim_date_w_quantile_mod = rbind(clim_date_w_quantile_mod, 
                                       this_data_for_pred %>% left_join(res))
      
    }
    
    clim_date_w_quantile_mod %>%
      saveRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
    
    # thresold at these points 
    
    quantile_model_fit_9 = readRDS("output/threshold_model_9")
    clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
    
    # get clim threshold values
    clim_thresh = clim_dat_full %>%
      group_by(id) %>%
      summarise(clim_thresh_value_9 = quantile(maxtp, 0.9))
    
    clim_date_w_quantile_mod = clim_date_w_quantile_mod %>% left_join(clim_thresh)
    clim_date_w_quantile_mod$threshold_9= predict(quantile_model_fit_9, clim_date_w_quantile_mod)$location
    
    #---- get lambda for climate model
    # Calculate lambda
    lambda_thresh_ex = clim_date_w_quantile_mod %>%
      group_by(id) %>%
      group_map(~{
        
        thresh_exceedance_9 = clim_date_w_quantile_mod%>%
          filter(id == .x$id[1]) %>%
          pull(temp_to_tau) %>%
          sapply(function(x) sapply(.x$threshold_9[1], x))
        
        
        tibble(id = .x$id[1],
               year = c(1971, 2022), #in the paper they never talk about the results from the untrained years so not doing it for untrained years 
               thresh_exceedance_9 = 1-thresh_exceedance_9)  
        
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    lambda_thresh_ex %>% 
      write_csv(paste0("Data/processed/climate_thresh_exceedance_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
    
    clim_date_w_quantile_mod %>% 
      left_join(lambda_thresh_ex) %>%
      saveRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
  }
}



glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

temporal_covariates = glob_anomaly_reshaped %>%
  dplyr::select(year, glob_anom) %>%
  unique() %>%
  arrange(year)

rm('glob_anomaly', 'glob_anomaly_reshaped')

#one iteration takes >20 minutes, maybe 30 minutes. Check if doing paralleization would be possible.

job::job({get_bulk_bts_on_clim_grid(seq(1,1),  'mod_0', 40, temporal_covariates)}, import=c("get_bulk_bts_on_clim_grid", "temporal_covariates"))


# job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_2', 40)}) #doing it on 40 quantiles
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,220), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(221,240), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(241,260), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(261,280), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(281,300), 'mod_2', 30)})


#weird ...........
# job::job({get_bulk_bts_on_clim_grid(seq(1,100),  'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,200), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,300), 'mod_0', 30)}
 
 
 job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_0', 40)}) #doing it on 40 quantiles 
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_0', 30)})

#weird ............
# job::job({get_bulk_bts_on_clim_grid(seq(1,100),  'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,200), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,300), 'mod_1', 30)}
 
 job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_1', 40)}) #doing it on 40 quantiles
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_1', 30)})



# ------------------- PLOTS 
marg_mod = 'mod_0'
num_quantiles = 40

true = read_csv(paste0("Data/processed/quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))

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

true_res = read_csv("Data/processed/thresh_exceedance_lambda_num_quantiles_05.csv") %>%
  group_by(year) %>%
  summarise(thresh_exceedance = mean(thresh_exceedance_9))

res = read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles,".csv"),
               col_names = c('bts', 'tau', 'b0', 'b1', 'b2')) #need to check if need to name the columns

# --- read in fitted quantile regression coefficients
true_pars = read_csv(paste0("data/processed/quantile_model_fit_pars_num_quantiles_", num_quantiles,".csv"),
                     col_names = c('tau', 'beta_0', 'beta_1', 'beta_2')) #need to check if need to rename columns 


plts=gridExtra::grid.arrange(res %>%
                               group_by(tau) %>%
                               summarise(upper = quantile(b2, 0.975),
                                         lower = quantile(b2, 0.030)) %>%
                               ggplot()+
                               geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
                               geom_line(data = true_pars, aes(tau, beta_1))+
                               xlim(c(0.08, 0.9))+
                               ylim(c(-1.5, 3))+
                               theme_minimal(12)+
                               labs(y = expression(beta[2]),
                                    x = expression(tau))+
                               theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                             all_bts_dat %>% 
                               group_by(bts, year) %>%
                               summarise(thresh_exceedance = mean(thresh_exceedance))  %>%
                               group_by(year) %>%
                               summarise(upper = quantile(thresh_exceedance, 0.975),
                                         lower = quantile(thresh_exceedance, 0.030)) %>%
                               ungroup() %>%
                               ggplot()+
                               geom_ribbon(aes(year, ymin = lower, ymax = upper), alpha = 0.3) +
                               geom_line(data = true_res, aes(year, thresh_exceedance))+
                               theme_minimal(12)+
                               labs(y = expression(lambda[t]),
                                    x = "Year")+
                               theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
                               xlim(c(1931, 2020))+
                               ylim(c(0.05, 0.19)), nrow = 1)

ggsave(plts,filename = "output/figs/bulk_quantile_regression.pdf", height = 3, width = 7.5)