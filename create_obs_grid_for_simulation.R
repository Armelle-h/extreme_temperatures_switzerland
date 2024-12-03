
gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

# remove sites that are too close to simualte
obs_sites_1 = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
  select("stn", "longitude", "latitude")%>%
  unique()

obs_sites_2 = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
  select("stn", "longitude", "latitude")%>%
  unique()

dist_h <- function(long1, lat1, long2, lat2) {
  
  point1 <- st_point(c(long1, lat1)) # Replace with actual values
  point2 <- st_point(c(long2, lat2)) # Replace with actual values
  
  sf_points <- st_sfc(point1, point2, crs = 4326)
  
  # Calculate the distance
  distance <- st_distance(sf_points[1], sf_points[2])
  
  return(as.numeric(distance))
}

# itterate through all stations and find climate grid point
closest_irel = c()
dist_smallest = c()
for(i in seq(nrow(obs_sites_1))){
  smallest_dist = 9999999999
  id_of_smallest_dist_irel = 9999999999
  
  this_dist = 10000
  for(j in seq(nrow(obs_sites_2))){
    this_dist = dist_h(obs_sites_1[i,]$longitude,
                       obs_sites_1[i,]$latitude,
                       obs_sites_2[j,]$longitude,
                       obs_sites_2[j,]$latitude)
    
    if((this_dist < smallest_dist) & this_dist > 0){
      smallest_dist = this_dist
      id_of_smallest_dist_irel = obs_sites_2$stn[j]
    }
  }
  dist_smallest = c(dist_smallest,smallest_dist)
  
  closest_irel = c(closest_irel, id_of_smallest_dist_irel)
}


sites_to_simulate = tibble(s1 = obs_sites_1$stn, s2 = closest_irel, dist_smallest) %>%
  arrange(desc(dist_smallest)) %>%
  head(86) 

sites_to_keep = c(sites_to_simulate$s1, sites_to_simulate$s2) %>% unique


sites_to_sim = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
  select("stn", "longitude", "latitude")%>%
  unique() %>%
  filter(stn %in% sites_to_keep)

locs_to_pred = sites_to_sim %>%
  dplyr::select(longitude, latitude)%>% 
  as.data.frame()


list(sites_to_sim = sites_to_sim, locs_to_pred = locs_to_pred) %>% 
  saveRDS("Data/processed/plain_obs_grid_to_simulate")