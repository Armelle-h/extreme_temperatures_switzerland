# --- data used for validating r-paretp process
gc()
rm(list = ls())
library(tidyverse)
library(sf)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
  select("stn", "longitude", "latitude")%>%
  unique()

site_pairs = as.data.frame(t(combn(legend_data$stn , 2 )))

dist_h <- function(long1, lat1, long2, lat2) {
  
  point1 <- st_point(c(long1, lat1)) # Replace with actual values
  point2 <- st_point(c(long2, lat2)) # Replace with actual values
  
  sf_points <- st_sfc(point1, point2, crs = 4326)
  
  # Calculate the distance
  distance <- st_distance(sf_points[1], sf_points[2])
  
  return(distance)
}

dists = c()
for(s in seq(nrow(site_pairs))){
  
  site_1 = legend_data[legend_data$stn == site_pairs[s,]$V1,]
  site_2 = legend_data[legend_data$stn == site_pairs[s,]$V2,]
  
  dists = c(dists, dist_h(site_1$longitude[1],
                          site_1$latitude[1],
                          site_2$longitude[1],
                          site_2$latitude[1]))
}

site_pairs$dist = dists

site_pairs %>%
  write_csv("Data/processed/plain_obs_pairs_with_dist.csv") #the distance is in meters