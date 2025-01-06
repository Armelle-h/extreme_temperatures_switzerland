gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

quant_reg_bts = read.csv(paste0("output/bts_quant_reg/bts_quant_reg_mod_2_num_quantiles_30.csv"), header = FALSE)
colnames(quant_reg_bts) = c("bts", "tau", "beta_0", "beta_1", "beta_2")

tau_list = sort(unique(quant_reg_bts$tau))

#histogram for the quantile regression 

for (tau_val in tau_list){
  file_path = paste0("pictures/beta_2, tau =",tau_val)
  
  png(filename = file_path, width = 443, height = 290)
  
  df = quant_reg_bts %>%
    filter(tau == tau_val)
  
  hist( df$beta_2, breaks = 100,
        main = paste0("beta_2, tau =", tau_val), 
        xlab = "beta_2", 
        ylab = "Frequency")
  
  dev.off()
}


exceedance_proba_bts = data.frame(bts = numeric(), stn = character(), year = numeric(), ex_proba = numeric())

for (bts_num in seq(1, 100)){
  exceedance_proba_bts_num = read.csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_mod_2_num_quantiles_30_bts_",bts_num,".csv"), header=FALSE)
  
  exceedance_proba_bts = rbind(exceedance_proba_bts_num, exceedance_proba_bts)
}

#histogram for threshold exceedance (one histogram per station)


fit_gpd_bts = read.csv("output/gpd_model_fits/bts/model_2_bts.csv", header = FALSE) %>%
  unique()
colnames(fit_gpd_bts) = c("bts", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4", "xi" )

#histogram for the gpd parameters 
hist( fit_gpd_bts$xi, breaks = 100,
     main = "xi", 
     xlab = "xi", 
     ylab = "Frequency")