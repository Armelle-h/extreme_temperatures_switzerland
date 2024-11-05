gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

n_stn = length(unique(obs_data$stn))

s_stn = unique(obs_data$stn)[seq(1, n_stn, 12)]

obs_data_plot = obs_data %>%
  filter(stn %in% s_stn)
  

ggplot(obs_data_plot, aes(x = maxtp)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_wrap(~ stn, scales = "free_y") +
  labs(
    title = "Histograms of Distribution for selected stations",
    x = "Temperature",
    y = "Frequency"
  ) +
  theme_minimal()

s_stn = unique(obs_data$stn)[seq(1, n_stn, 48)]

obs_data_plot = obs_data %>%
  filter(stn %in% s_stn)


ggplot(obs_data_plot, aes(x = maxtp)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_wrap(~ stn, scales = "free_y") +
  labs(
    title = "Histograms of Distribution for selected stations",
    x = "Temperature",
    y = "Frequency"
  ) +
  theme_minimal()