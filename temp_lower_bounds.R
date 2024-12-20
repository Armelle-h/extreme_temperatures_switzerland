
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

marg_mod = "mod_1"

data_for_rpareto = readRDS(paste0("Data/processed/data_for_rpareto/true/robust_plain_data_for_rpareto_", marg_mod))

extreme_dates = data_for_rpareto$extreme_dates

obs_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_data_loc_id.csv")%>%
  group_by(date) %>%
  summarise(maxtp_date = max(maxtp, na.rm = TRUE))


obs_data_not_extreme_dates = obs_data %>%
  filter(!(date %in% extreme_dates$date))

obs_data_extreme_dates = obs_data %>%
  filter(date %in% extreme_dates$date)
