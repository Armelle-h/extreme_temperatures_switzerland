#computing the percentage of stations with more than 30 years of data associated 

gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

obs_data = obs_data %>%
  mutate(year = lubridate::year(date))

obs_data = obs_data %>%
  select(c(stn, year)) %>%
  unique()

count = obs_data %>%
  count(stn, name = "nb_years")

count_filter = count %>%
  filter(nb_years >=30)

print(nrow(count_filter)/nrow(count))
