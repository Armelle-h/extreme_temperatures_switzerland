gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")%>%
  mutate(year = lubridate::year(date))

num_quantiles = 30

threshold_9_df = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv")) %>%
  select(c(stn, threshold_9))%>% #enough as threshold_9 varies only spatially 
  unique()

obs_data = obs_data %>%
  left_join(threshold_9_df, by = "stn")

nb_extreme = c()

seq_years = seq(1971, 2022)

for (year_val in seq_years){
  obs_data_year = obs_data %>%
    filter(year==year_val)
  
  obs_data_year_extreme = obs_data_year %>%
    group_by(stn)%>%
    filter(maxtp>threshold_9)%>%
    ungroup()
  
  nb_extreme = c(nrow(obs_data_year_extreme), nb_extreme)
}

plot(seq_years, nb_extreme)
plot(seq_years, coef(model)[[1]] + seq_years*coef(model)[[2]])

model = lm(nb_extreme ~ seq_years)

plot(seq_years, nb_extreme, 
     main = "Scatter Plot with Linear Fit",
     xlab = "Years", 
     ylab = "Number of Extremes", 
     pch = 16, col = "blue")

# Add the fitted line from the linear model
abline(model, col = "red", lwd = 2)

#doing it for each station

gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")%>%
  mutate(year = lubridate::year(date))

num_quantiles = 30

threshold_9_df = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv")) %>%
  select(c(stn, threshold_9))%>% #enough as threshold_9 varies only spatially 
  unique()

obs_data = obs_data %>%
  left_join(threshold_9_df, by = "stn")

seq_years = seq(1971, 2022)

nb_extreme_df = data.frame(stn = string, nb_extreme = numeric())

seq_years = seq(1971, 2022)

for (year_val in seq_years){
  obs_data_year = obs_data %>%
    filter(year==year_val)
  
  obs_data_year = obs_data_year %>%
    summarize(nb_extreme = sum(maxtp>threshold_9))
  
  nb_extreme_df = rbind(nb_extreme_df, obs_data_year)
}

plot(seq_years, nb_extreme)
plot(seq_years, coef(model)[[1]] + seq_years*coef(model)[[2]])

model = lm(nb_extreme ~ seq_years)

plot(seq_years, nb_extreme, 
     main = "Scatter Plot with Linear Fit",
     xlab = "Years", 
     ylab = "Number of Extremes", 
     pch = 16, col = "blue")

# Add the fitted line from the linear model
abline(model, col = "red", lwd = 2)