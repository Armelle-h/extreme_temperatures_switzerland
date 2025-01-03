gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")%>%
  mutate(year = lubridate::year(date))

obs_data_extreme = obs_data %>%
  group_by(year, stn)%>%
  summarise(max_tp_year = quantile(maxtp, 0.9, na.rm = TRUE), .groups = "drop")


final_df = obs_data_extreme %>%
  group_by(stn)%>%
  summarize(diff = max(max_tp_year)-min(max_tp_year), .groups = "drop")