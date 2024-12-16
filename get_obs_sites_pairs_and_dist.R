# --- data used for validating r-paretp process
gc()
rm(list = ls())
library(tidyverse)
library(sf)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#for consistency, computing the distance using the projected coordinates
#projected coordinates are in km so distance will be in km

legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>%
  select("stn", "longitude_proj", "latitude_proj")%>%
  unique()

site_pairs = as.data.frame(t(combn(legend_data$stn , 2 )))

dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

dists = c()
for(s in seq(nrow(site_pairs))){
  
  site_1 = legend_data[legend_data$stn == site_pairs[s,]$V1,]
  site_2 = legend_data[legend_data$stn == site_pairs[s,]$V2,]
  
  dists = c(dists, dist_h(site_1$longitude_proj[1],
                          site_1$latitude_proj[1],
                          site_2$longitude_proj[1],
                          site_2$latitude_proj[1]))
}

site_pairs$dist = dists

site_pairs %>%
  write_csv("Data/processed/plain_obs_pairs_with_dist.csv") #the distance is in kilometers

#same thing for the climate data 

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")%>%
  select("id", "longitude_proj", "latitude_proj")%>%
  filter(id %% 25 ==0)%>% #else way too heavy
  unique()

site_pairs = as.data.frame(t(combn(I_plain$id , 2 )))

dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

dists = c()
for(s in seq(nrow(site_pairs))){
  
  if (s%%100000 == 0){print(s)}
  
  site_1 = I_plain[I_plain$id == site_pairs[s,]$V1,]
  site_2 = I_plain[I_plain$id == site_pairs[s,]$V2,]
  
  dists = c(dists, dist_h(site_1$longitude_proj[1],
                          site_1$latitude_proj[1],
                          site_2$longitude_proj[1],
                          site_2$latitude_proj[1]))
}

site_pairs$dist = dists

site_pairs %>%
  write_csv("Data/processed/plain_clim_pairs_with_dist.csv") #the distance is in kilometers