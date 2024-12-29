setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(sf)


I = read.csv("Data/id_lon_lat_correspondance.csv")

my_coords <- I %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)  # WGS84 (EPSG:4326)

# Transform the coordinates to UTM Zone 29N (EPSG:32629)
proj_cords <- st_transform(my_coords, crs = 32629)

# Extract the transformed coordinates, expressed in meters
proj_coords <- st_coordinates(proj_cords) 

# Add projected coordinates and an ID column to clim_data, expressed in kilometers
I$longitude_proj <- proj_coords[, 1] / 1000
I$latitude_proj <- proj_coords[, 2] / 1000


legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv") %>%
  select(c("stn", "longitude", "latitude"))

my_coords <- legend %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)  # WGS84 (EPSG:4326)

# Transform the coordinates to UTM Zone 29N (EPSG:32629)
proj_cords <- st_transform(my_coords, crs = 32629)

# Extract the transformed coordinates, expressed in meters
proj_coords <- st_coordinates(proj_cords) 

# Add projected coordinates and an ID column to clim_data, expressed in kilometers
legend$longitude_proj <- proj_coords[, 1] / 1000
legend$latitude_proj <- proj_coords[, 2] / 1000

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv") %>%
  left_join(legend %>% select(stn, longitude_proj, latitude_proj), by = "stn") %>%
  select(stn, id, longitude_proj, latitude_proj) %>%
  unique()

#creating column where the location id will be saved
obs_data$id_proj <- NA

#finding the correct id by associating the id of the pair latitude-longitude that is closest
#according to the euclidean distance

# Function to compute the Euclidean distance between two points
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Loop over each row in obs_loc
for (i in 1:nrow(obs_data)) {
  # Extract the current pair from obs_loc
  x1 <- obs_data$longitude_proj[i]
  y1 <- obs_data$latitude_proj[i]
  
  # Compute the distances between this pair and all pairs in id_loc
  distances <- mapply(euclidean_distance, x1, y1, I$longitude_proj, I$latitude_proj)
  
  # Find the index of the minimum distance
  min_index <- which.min(distances)
  
  # Assign the corresponding value from id_loc$id to obs_loc
  obs_data$id_proj[i] <- I$id[min_index]
}