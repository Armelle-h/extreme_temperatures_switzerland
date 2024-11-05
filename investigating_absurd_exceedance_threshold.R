library(tidyverse)
library(evgam)
library(Monospline)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


L_absurd = lambda_thresh_ex %>% filter(thresh_exceedance_9>1 | thresh_exceedance_9<0)

#clim_thresh_value_9 and threshold_9 are defined spatially and not temporally

L_absurd = L_absurd %>% left_join (obs_data %>% select("stn", "clim_thresh_value_9", "threshold_9"), by= "stn" )

L_absurd = unique(L_absurd)

#investigating if any threshold exceedance is associated to observed values 


for (i in seq(1, nrow(L_absurd))){
  stn_val = L_absurd$stn[[i]] #prints with repetition, all good !
  year_val = L_absurd$year[[i]]
  
  nb_row = nrow(obs_data %>% filter(stn==stn_val & year==year_val))
  
  if (nb_row != 0){
    cat(i, ", ", stn_val, ", ", year_val, ", ", nb_row, "\n")
  }
  
}


temp_to_tau_function = (obs_smoothed_quantiles %>% filter(stn== "AIR" & year== 1971) )$temp_to_tau[[1]]

list_temp = c(
2.191140, 7.440990, 8.851104, 9.816488, 10.722202, 16.672270, 12.173902, 18.400830, 13.268828, 13.827848, 14.263268, 
14.647108, 14.991075, 15.296208, 15.618966, 15.943421, 16.253215, 16.615179, 16.927505, 17.284826, 17.662910, 18.070279,
18.504375, 18.910435, 19.365001, 19.853072, 20.514873, 24.530417, 29.906171, 24.495306)

for (tp in list_temp){
  cat(tp, ":", temp_to_tau_function(tp), "\n")
}

list_temp = seq(1,30)

for (tp in list_temp){
  cat(tp, ":", temp_to_tau_function(tp), "\n")
}



obs_smoothed_quantiles_test = obs_data %>%
  group_by(stn) %>%
  group_map(~{
    
    if(.x$stn[1] == "AIR"){
      
      # --- get the climate quantile estimates closest to current station
      clim_vals <<- obs_data %>%
        filter(stn == .x$stn[1]) %>%
        dplyr::select(quantile, value) %>%
        unique() %>%
        pull(value) %>%
        unlist()
      
      # predict (observed) quantile for each year and site
      quant_reg_pars = quant_reg_model_pars %>%
        arrange(tau)
      
      res = c()
      for(q in seq_along(quantiles_to_estimate_bulk)){
        qpars = quant_reg_pars[q,]
        # Calculate the predicted quantile value using the regression parameters
        res = rbind(res,
                    tibble(quantile =  qpars$tau,
                           year = temporal_covariates$year,
                           quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2) *(temporal_covariates$glob_anom)))
        
        if(.x$stn[1]=="AIR"){cat("interp_quantile", qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2) *(temporal_covariates$glob_anom), "\n", "beta_1", qpars$beta_1, "\n", "clim_quant", clim_vals[q], "\n")}
      }
      
      print(paste0("Interpolating quantile estimates for ", .x$stn[1]))
      
      # interpolate quantiles over tau for each year
      res %>%
        group_by(year) %>%
        group_map(~{
          tau_val = .x$quantile
          temp_val = .x$quant_value
          tibble(year = .x$year[1],
                 #creates a spline function to go from tau to temp or temp to tau. 
                 #the function depends on the year
                 
                 tau_to_temp = list(scam(temp_val ~ s(tau_val, bs = "mpi"))),
                 
                 temp_to_tau = list(scam(tau_val ~ s(temp_val, bs = "mpi"))))
        }, .keep = T) %>%
        plyr::rbind.fill() %>%
        as_tibble() %>%
        mutate(stn = .x$stn[1])
    }  
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


#one option but I don't like it so much


temp_to_tau_function_test = (obs_smoothed_quantiles_test %>% filter(stn== "AIR" & year== 1971) )$temp_to_tau[[1]]



list_temp = seq(1,30)

for (tp in list_temp){
  cat(tp, ":", predict(temp_to_tau_function_test, tp ), "\n")
}

list_temp = c(2.5, 10.3, 12.7, 28.9)

print(predict(temp_to_tau_function_test, new_data=data.frame(temp_val = list_temp)))

#neither works, nor is increasing :( 

