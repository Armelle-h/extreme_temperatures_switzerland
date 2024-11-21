gc()
rm(list = ls())
library(tidyverse)
library(evgam)
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


# investigating the evolution in the yearly max tp for each stations 

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv") 

obs_data_max_year = obs_data %>%
  mutate(year = lubridate::year(date))%>%
  group_by(stn, year)%>%
  summarise(max_y_tp = max(maxtp, na.rm = TRUE), .groups = "drop")

extract_coeff = function(col1, col2){
  model = lm(col1 ~ col2)
  
  return (coef(model)[2])
}

count = obs_data_max_year %>%
  count(stn, name = "nb_years")

obs_data_coeff = obs_data_max_year %>%
  group_by(stn) %>%
  summarise(coeff = extract_coeff(max_y_tp, year), .groups = "drop")%>%
  left_join(count, by = "stn")



