gc()
rm(list = ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(sf)
library(raster)

I = read.csv("Data/id_lon_lat_correspondance.csv")%>%
  mutate(long_trunc = trunc(longitude * 1000) / 1000)%>%
  mutate(lat_trunc = trunc(latitude * 1000) / 1000)

# Specify the path to the .nc file
file_path <- "Data/topography.nc"

# Load the .nc file as a raster
altitude_raster <- raster(file_path)

# Convert raster to a dataframe with coordinates
altitude_df <- as.data.frame(altitude_raster, xy = TRUE, na.rm = FALSE)

# Rename the columns
colnames(altitude_df) <- c("longitude", "latitude", "altitude")

altitude_df = altitude_df %>%
  filter(longitude > min(I$longitude) & longitude <max(I$longitude) & latitude > min(I$latitude) & latitude <max(I$latitude))

altitude_df = altitude_df %>%
  mutate(long_trunc = trunc(longitude * 1000) / 1000)%>%
  mutate(lat_trunc = trunc(latitude * 1000) / 1000)

I_alt = I %>%
  left_join(altitude_df, by =c("long_trunc", "lat_trunc"))%>%
  select(id, altitude)

I_alt_final = I_alt %>%
  group_by(id)%>%
  summarise(altitude = mean(altitude, na.rm = TRUE), .groups = "drop")

write.csv(I_alt_final, "Data/clim_id_altitude_corres.csv")
