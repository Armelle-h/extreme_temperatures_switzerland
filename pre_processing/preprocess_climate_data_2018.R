#a fortiori, is successful :) 
setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

rm(list =ls())
library(raster) # package for netcdf manipulation
library(tidyverse)
library(rnaturalearth)
library(magrittr)
library(raster)

file.loc = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/2018_9km/"
file.save = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/2018_JJA_climate_data.csv"

my_files = list.files(file.loc)

lonlat_raster_file = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/TmaxD_ch01r.swisscors_201901010000_201912310000.nc"

lonlat_raster = raster::brick(lonlat_raster_file)
crs(lonlat_raster) = 21781

index=0

last_date_id = 4324 #insider knowledge

date_mapping <- data.frame(date = character(), date_id = integer(), stringsAsFactors = FALSE)

for(f in paste0(file.loc, my_files)){
  
  index <- index+1
  
  my_dat = raster::brick(f) 
  
  #preprocessing the data so that we have the max per day
  total_layers <- nlayers(my_dat)
  
  result_layers <- list()
  
  for (i in seq(1, total_layers, by = 13)) {
    
    # Subset the current group of 13 layers
    layer_group <- my_dat[[i:(i + 12)]]
    
    # Calculate the maximum value for the group of 13 layers
    max_layer <- calc(layer_group, fun = max)
    
    # Store the result in the list
    result_layers[[length(result_layers) + 1]] <- max_layer
  }
  
  result_brick <- stack(result_layers)
  
  crs(result_brick) = 4326
  
  result_brick_proj = projectRaster(result_brick, crs = 21781)
  
  #converting from 9km grid to 1 km grid
  
  raster_1km_interpolated <- resample(result_brick_proj, lonlat_raster, method="bilinear")
  
  my_data_swiss_c = raster::coordinates(lonlat_raster) %>% as_tibble() #swiss coordinates
  names(my_data_swiss_c) = c('longitude', 'latitude')
  
  my_sf_data_swiss_c <- st_as_sf(my_data_swiss_c, coords = c("longitude", "latitude"), crs = 21781)  # CH1903 / LV95
  
  # Transform the coordinates to WGS84
  my_sf_data <- st_transform(my_sf_data_swiss_c, crs = 4326)  # EPSG 4326 is WGS84
  
  # Convert back to a dataframe if needed
  my_data <- st_coordinates(my_sf_data) %>% as_tibble()
  names(my_data) <- c('longitude', 'latitude')
  
  my_data = my_data %>%
    group_by(longitude, latitude) %>%
    dplyr::mutate(id = cur_group_id())
  
  #combine the data with their associated location and unique id
  #id will be used to filter out data not in Switzerland
  all_data = raster_1km_interpolated %>%
    values() %>%
    as_tibble()
  all_data$id = my_data$id
  names(all_data) = names(all_data) %>% str_remove_all("X")
  
  #each column corresponds to a date
  
  #not generalizable, quick fix
  if (index==1){
    time_date="2018-06-%02d"
    new_names <- sprintf(time_date, 1:30)
    names(all_data)[1:30] <- new_names
  }
  if (index==2){
    time_date="2018-07-%02d"
    new_names <- sprintf(time_date, 1:31) #to change
    names(all_data)[1:31] <- new_names
  }
  if (index==3){
    time_date="2018-08-%02d"
    new_names <- sprintf(time_date, 1:31) #to change
    names(all_data)[1:31] <- new_names
  }
  
  switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")
  
  points_sf <- sf::st_as_sf(my_data, coords = c("longitude", "latitude"), crs = 4326)
  points_within_switzerland <- points_sf[sf::st_intersects(points_sf, switzerland, sparse = FALSE), ]
  
  locs = tibble(id = points_within_switzerland$id,
                longitude = sf::st_coordinates(points_within_switzerland)[,1],
                latitude = sf::st_coordinates(points_within_switzerland)[,2])
  
  locs = locs %>% distinct()
  
  #restricts to points inside Switzerland  
  all_data <- all_data %>%
    right_join(locs, by = 'id') %>%
    dplyr::select(-id) %>%
    tidyr::pivot_longer(c(-longitude, -latitude), names_to = "date", values_to = "maxtp") %>%
    mutate(date = lubridate::ymd(date))%>%
    mutate(maxtp = maxtp - 273.15) #converting from Kalvin to Celsius degrees
  
  #doing a date mapping 
  
  num_dates <- length(unique(all_data$date))
  current_dates <- unique(all_data$date)
  new_ids <- seq(last_date_id +1, last_date_id + num_dates)  
  
  new_mapping <- data.frame(date = current_dates, date_id = new_ids, stringsAsFactors = FALSE)
  date_mapping <- rbind(date_mapping, new_mapping)
  
  last_date_id <- max(new_ids)
  
  all_data_with_id <- all_data %>%
    left_join(new_mapping, by = "date") %>%  # Add the date ID
    select(-date) %>%  # Remove the original date column
    group_by(longitude, latitude) %>%
    mutate(id = cur_group_id()) %>% #guarantees correpondance id - lon lat is unique across files
    ungroup() %>%
    mutate(maxtp = round(maxtp, 2))
  
  if (!file.exists("Data/id_lon_lat_correspondance_2018.csv")){
    #to ensure each row is unique
    id_loc <- all_data_with_id %>%
      dplyr::select(longitude, latitude, id) %>%
      distinct()
    write.csv(id_loc, "Data/id_lon_lat_correspondance_2018.csv", row.names = FALSE)
  }
  
  #Write only the date ID column in a CSV
  all_data_with_id %>%
    select(date_id, maxtp, id) %>%  # Choose the columns to save
    write.table(file.save,
                sep = ",",
                col.names = !file.exists(file.save),
                append = TRUE,
                row.names = FALSE)
  
  print(nrow(all_data_with_id)) #sanity check
}

write.table(date_mapping, "Data/Climate_data/id_date_correspondance.csv", sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)