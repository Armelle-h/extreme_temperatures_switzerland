#In this file, we visualise the temperature associated with each station on the plain on a Pareto margin 


gc()
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(sf)
library(rnaturalearth)

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

legend = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")%>%
  select(stn, longitude, latitude)

obs_data_standardised = read.csv("Data/processed/plain_obs_data_pareto_frechet_scale.csv") %>%
  filter(date %in% c("1972-06-01", "1982-06-01", "1992-06-01", "2002-06-01", "2012-06-01", "2022-06-01") )%>%
  select(stn, pareto_marg, date, year)%>%
  left_join(legend, by = "stn")


for (date_ in c("1972-06-01", "1982-06-01", "1992-06-01", "2002-06-01", "2012-06-01", "2022-06-01")){
  
  points_sf <- st_as_sf(obs_data_standardised %>% filter(date == date_) , coords = c("longitude", "latitude"), crs = 4326)

  plot = ggplot(data = switzerland) +
    geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
    geom_sf(data = points_sf, aes(color = pareto_marg), size = 2, alpha=0.7)+ 
    ggtitle(date_) + 
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude")+
    scale_color_gradient(low = "blue", high = "red")
  
  print(plot)
  
}

