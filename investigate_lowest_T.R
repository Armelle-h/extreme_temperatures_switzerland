rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

marg_mod = "mod_1"

standardised_data = read.csv(paste0("Data/processed/plain_obs_data_pareto_frechet_scale_", marg_mod,".csv"))%>%
  select(stn, date, pareto_marg, maxtp)

spatial_threshold = readRDS(paste0("Data/spatial_threshold_", marg_mod,".rds"))

standardised_data_filtered = standardised_data%>%
  filter(pareto_marg>spatial_threshold)