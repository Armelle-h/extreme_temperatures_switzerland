#a fortiori, is successful :) 
setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

rm(list =ls())
library(raster) # package for netcdf manipulation
library(tidyverse)
library(rnaturalearth)
library(magrittr)
library(raster)

file.loc = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/2018_9km/"
file.save = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/temp_2018_JJA_climate_data.csv"

my_files = list.files(file.loc)

#for loop takes a while, will have to parallelize

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
  
  #extracting the latitude-longitude from result_brick (not explicitely given)
  my_data = raster::coordinates(result_brick) %>% as_tibble()
  names(my_data) = c('longitude', 'latitude')
  
  my_data = my_data %>% #creates an id column such that each pair "latitude-longitude"
    mutate(longitude = signif(longitude, 4), #is uniquely associated to an id (from 1 to 1127)
           latitude = signif(latitude, 4)) %>%
    group_by(longitude, latitude) %>%
    dplyr::mutate(id = cur_group_id())
  
  #combine the data with their associated location and unique id
  #id will be used to filter out data not in Switzerland
  all_data = result_brick %>%
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
    mutate(id_2018 = cur_group_id()) %>% #guarantees correpondance id - lon lat is unique across files
    ungroup() %>%
    mutate(maxtp = round(maxtp, 2)) 
  
  if (!file.exists("Data/Climate_data/id_lon_lat_correspondance_2018.csv")){
    #to ensure each row is unique
    id_loc <- all_data_with_id %>%
      select(longitude, latitude, id_2018) %>%
      distinct()
    write.csv(id_loc, "Data/Climate_data/id_lon_lat_correspondance_2018.csv", row.names = FALSE)
  }
  
  #Write only the date ID column in a CSV
  all_data_with_id %>%
    select(date_id, maxtp, id_2018) %>%  # Choose the columns to save
    write.table(file.save,
                sep = ",",
                col.names = !file.exists(file.save),
                append = TRUE,
                row.names = FALSE)
}

#updating the id column such that the id is the same as for the data on 1km grid 

id_loc_2018 = read.csv("Data/Climate_data/id_lon_lat_correspondance_2018.csv")

id_loc = read.csv("Data/id_lon_lat_correspondance.csv")

euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Create an empty vector to store the closest index for each row in Dataframe2
closest_indices <- integer(nrow(id_loc))

# Loop over each row in Dataframe2
for (i in 1:nrow(id_loc)) {
  # Calculate the distance between the current row in Dataframe2 and all rows in Dataframe1
  distances <- euclidean_distance(id_loc_2018$longitude, id_loc_2018$latitude,
                                  id_loc$longitude[i], id_loc$latitude[i])
  
  # Find the index of the minimum distance
  closest_indices[i] <- which.min(distances)
}

# Now join column5 from Dataframe2 to Dataframe1 using the closest indices
id_loc_2018$id <- id_loc$id[closest_indices]

id_loc_2018 = id_loc_2018 %>%
  select(-latitude, -longitude)

climate_data_2018 = read.csv("Data/Climate_data/2018_JJA_climate_data.csv")

final_data_2018 = climate_data_2018 %>%
  left_join(id_loc_2018, by = "id_2018") %>%
  select(-id_2018)  # Remove id_2018

write.csv(final_data_2018, "Data/Climate_data/2018_JJA_climate_data.csv", row.names = FALSE)



#---------------------------------------------------

id_loc_2018 = read.csv("Data/Climate_data/id_lon_lat_correspondance_2018.csv")

id_loc = read.csv("Data/id_lon_lat_correspondance.csv")

id_loc_2018$id <- NA

euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Loop over each row in id_loc_2018
for (i in 1:nrow(id_loc_2018)) {
  # Extract the current pair from id_loc_2018
  x1 <- id_loc_2018$longitude[i]
  y1 <- id_loc_2018$latitude[i]
  
  # Compute the distances between this pair and all pairs in id_loc
  distances <- mapply(euclidean_distance, x1, y1, id_loc$longitude, id_loc$latitude)
  
  # Find the index of the minimum distance
  min_index <- which.min(distances)
  
  # Assign the corresponding value from id_loc$id to obs_loc
  id_loc_2018$id[i] <- id_loc$id[min_index]
}

id_loc_2018 = id_loc_2018 %>%
  select(-latitude, -longitude)

climate_data_2018 = read.csv("Data/Climate_data/temp_2018_JJA_climate_data.csv")

final_data_2018 = climate_data_2018 %>%
  left_join(id_loc_2018, by = "id_2018") %>%
  select(-id_2018)  # Remove id_2018

write.csv(final_data_2018, "Data/Climate_data/2018_JJA_climate_data.csv", row.names = FALSE)

