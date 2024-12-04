
library(tidyverse)
library(sf)
library(rnaturalearth)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>%
  rename(altitude = Altitude.m.)%>%
  select(stn, Nom, altitude, longitude, latitude)

obs_data = obs_data %>%
  left_join(legend, by="stn")

obs_data_plain = obs_data %>%
  filter(altitude<1200)%>%
  filter(!(longitude>=6.8 & latitude<=46.4))%>%
  filter(!(longitude>=9.5 & latitude<=47.1))%>%
  filter(!(longitude>=8.5 & latitude<=46.7))%>%
  filter(!(stn %in% c("GST", "KALTB", "GWA", "ENG", "GTT", "INT", "BEA", "MMBIS", "ELM", "VAE", "ENNESF", "ILZ", "MMLEN", "FRU", "MMGTT", "MER", "MMLIN", "TFD", "ALT", "INNESF")))

#combining some stations that don't have the same name but are at exactly the same place, and record temperatures at non intersecting times

obs_data_plain <- obs_data_plain %>%
  mutate(stn = if_else(stn == "AGAAR", "AAR", stn))%>%
  mutate(stn = if_else(stn == "BZN", "BEZ", stn))%>%
  mutate(stn = if_else(stn == "WSLLAE", "NABLAE", stn))%>%
  mutate(stn = if_else(stn == "TGFRA", "FRF", stn))


#The stations WSLBTF and WSLBTB are overlapping, WSLBTF has more measurements, keeping WSLBTF, deleting WSLBTB
#The stations WSLHOC and WSLHOB are overlapping, WSLHOC has more measurements, keeping WSLHOC, deleting WSLHOB.

obs_data_plain = obs_data_plain %>%
  filter( !(stn %in% c("WSLBTB", "WSLHOB")) )

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

# Step 3: Convert the dataframe to an sf object, assuming coordinates in WGS84 (EPSG:4326)
pts_sf <- st_as_sf(obs_data_plain, coords = c("longitude", "latitude"), crs = 4326)

# Step 4: Ensure the CRS for both datasets match (WGS84)
switzerland <- st_transform(switzerland, crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = pts_sf) + 
  ggtitle("plain vs mountain") +
  theme_minimal()

legend_filtered <- legend %>%
  filter(stn %in% obs_data_plain$stn)

write.csv(legend_filtered, "Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv", row.names = FALSE)

write.csv(obs_data_plain%>%select(stn, date, maxtp, id), "Data/Observed_data/plain_1971_2022_JJA_obs_data_loc_id.csv", row.names=FALSE)


library(data.table)
library(rnaturalearth)
library(sf)

I = read.csv("Data/id_lon_lat_correspondance.csv")
L = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")

sf1 <- st_as_sf(L, coords = c("longitude", "latitude"), crs = 4326)
sf2 <- st_as_sf(I, coords = c("longitude", "latitude"), crs = 4326)

# Transform to a projected CRS for distance calculation
sf1_proj <- st_transform(sf1, 3857)
sf2_proj <- st_transform(sf2, 3857)

# Find points within 5 km
within_5km <- st_is_within_distance(sf2_proj, sf1_proj, dist = 5000)

# Filter df2 points within 5 km of df1 points
I_within_5_km <- I[apply(within_5km, 1, any), ]

I_filtered = I %>%
  rowwise() %>%
  filter(any(I_within_5_km$longitude > longitude & I_within_5_km$latitude < latitude) | longitude<6.5  ) %>%
  ungroup()



files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_data_extreme_9_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])%>%
    filter(date_id == 130 & id %in% I_filtered$id)
  
  
  clim_data_extreme_9_list[[i]] = clim_data
  
  #to free memory
  rm(clim_data)
  gc()
}

clim_data_extreme_9 =  do.call(rbind, clim_data_extreme_9_list)

clim_data_extreme_9 = clim_data_extreme_9%>%
  left_join(I_second_try, by="id")

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)


clim_data_extreme_9 %>%  #the weird position of the points comes from the fact that we're selecting only every 10 points to be plot.
  ggplot()+
  geom_point(aes(longitude, latitude))+
  coord_map()+
  theme_minimal()+
  geom_sf(data = switzerland, alpha = 0, col = 'black')+
  theme_minimal()

write.csv(I_second_try, "Data/plain_id_lon_lat_correspondance.csv", row.names = FALSE)


#adding a column to legend_data with latitude and longitude projected, is useful for the rPareto process

library(sf)
library(tidyverse)

legend_data = read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv")


my_coords <- legend_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)  # WGS84 (EPSG:4326)

# Transform the coordinates to UTM Zone 29N (EPSG:32629)
proj_cords <- st_transform(my_coords, crs = 32629)

# Extract the transformed coordinates, expressed in meters
proj_coords <- st_coordinates(proj_cords) 

# Add projected coordinates and an ID column to clim_data, expressed in kilometers
legend_data$longitude_proj <- proj_coords[, 1] / 1000
legend_data$latitude_proj <- proj_coords[, 2] / 1000

write.csv(legend_data, "Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv", row.names = FALSE) #distance will be in kilometers
