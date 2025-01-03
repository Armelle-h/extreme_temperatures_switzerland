gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30

for(file_name in seq(1, 100)) {
  
  dat <- read.csv(paste0("Data/processed/bootstrap_data/bts_under_gpd_models_debug/num_quantiles_",num_quantiles,"_bts_", file_name, ".csv"), header = FALSE)%>%
    setNames(c("stn", "date", "scale_9", "threshold_9", "maxtp_0", "maxtp_1", "maxtp_2", "maxtp_3"))
  
  if (dat$stn[[1]] == "stn"){
    print(file_name)
    dat = dat[-1, ]
    colnames(dat) <- dat[1, ]
    dat = dat[-1, ]
    
    write.csv(dat, (paste0("Data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_", file_name, ".csv")), row.names = FALSE)
  }
}