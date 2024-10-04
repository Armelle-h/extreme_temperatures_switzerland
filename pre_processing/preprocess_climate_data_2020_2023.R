setwd("C:/Users/HOURS/Desktop/PDM/Code_R")


#so far, saves only id and maxtp (not even date). I want to do an id column
#for date

#testing preprocessing and file saving of a netcdf file 
#recall: they were doing things in parallel, iterating through the file, 
#writing the data in a csv file 

#for single month, exactly the same except I don't have the month extraction

rm(list =ls())
library(raster) # package for netcdf manipulation
library(tidyverse)
library(rnaturalearth)
library(magrittr)
library(ncdf4) #needed for loading netcdf data 
library(raster)
library(sf)

file.loc = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/TmaxD_2020_2023_ch01r.swisscors.zip"
file.save = "C:/Users/HOURS/Desktop/PDM/Code_R/Data/Climate_data/2020_2023_JJA_climate_data.csv"

file_list = unzip(file.loc, list = TRUE)

my_files = file_list$Name

#bisextile years
bisextile <-seq(1972, 2020, 4)

#to store date and date id
date_mapping <- data.frame(date = character(), date_id = integer(), stringsAsFactors = FALSE)

last_date_id = 4508 #insider knowledge

nb_file = 0

#function to rename columns under the appropriate date format if needed
rename_columns <- function(df, year, month) {
  
  if (month == 06){
    tot_length = 30}
  else {
    tot_length = 31}
  
  # Create a sequence of dates for the specified month and year
  dates <- seq(as.Date(paste(year, month, "01", sep = "-")),
               by = "day",
               length.out = tot_length)
  
  # Format dates to "yyyy-mm-dd"
  new_colnames <- format(dates, "%Y-%m-%d")
  
  # Keep the last column name unchanged
  new_colnames <- c(new_colnames, colnames(df)[ncol(df)])
  
  # Rename columns
  colnames(df) <- new_colnames
  
  return(df)
}

for(f in my_files){
  
  nb_file = nb_file +1
  
  f_unzip = unzip(file.loc, files = f, exdir=tempdir())
  
  my_dat = raster::brick(f_unzip)
  
  #no need to extract months, already extracted
  
  #extracting geographical coordinates 
  my_data_swiss_c = raster::coordinates(my_dat) %>% as_tibble() #swiss coordinates
  names(my_data_swiss_c) = c('longitude', 'latitude')
  
  #starting 2022, going from old to new coordinate system
  
  year <- as.numeric(substr(f, 23, 26))
  
  if (year %in% c(2020, 2021)){
    
    crs_ = 21781
    
  } else {
    
    crs_ = 2056
  }
  
  my_sf_data_swiss_c <- st_as_sf(my_data_swiss_c, coords = c("longitude", "latitude"), crs = crs_)  # CH1903 / LV95
  
  # Transform the coordinates to WGS84
  my_sf_data <- st_transform(my_sf_data_swiss_c, crs = 4326)  # EPSG 4326 is WGS84
  
  # Convert back to a dataframe if needed
  my_data <- st_coordinates(my_sf_data) %>% as_tibble()
  names(my_data) <- c('longitude', 'latitude')
  
  my_data = my_data %>% #creates an id column such that each pair "latitude-longitude"
    mutate(longitude = signif(longitude, 4), #is uniquely associated to an id (from 1 to 1127)
           latitude = signif(latitude, 4)) %>%
    group_by(longitude, latitude) %>%
    dplyr::mutate(id = cur_group_id())
  
  #combine the data with their associated location and unique id
  #id will be used to filter out data not in Switzerland
  all_data = my_dat %>%
    values() %>%
    as_tibble()
  all_data$id = my_data$id
  names(all_data) = names(all_data) %>% str_remove_all("X")
  
  #in 2020 and 2021, the columns do not have a date structure name
  
  year = as.numeric(substr(f, 23, 26))
  month = as.numeric(substr(f, 27, 28))
  
  if (year %in% c("2020", "2021")){
    
    all_data  = rename_columns(all_data, year, month)
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
    mutate(date = lubridate::ymd(date))
  
  #doing a date mapping 
  
  num_dates <- length(unique(all_data$date))
  current_dates <- unique(all_data$date)
  new_ids <- seq(last_date_id +1, last_date_id + num_dates)
  
  new_mapping <- data.frame(date = current_dates, date_id = new_ids, stringsAsFactors = FALSE)
  date_mapping <- rbind(date_mapping, new_mapping)
  
  last_date_id <- max(new_ids)
  
  #adding the date index to the dataframe
  
  all_data_with_id <- all_data %>%
    left_join(new_mapping, by = "date") %>%  # Add the date ID
    select(-date) %>%  # Remove the original date column
    group_by(longitude, latitude) %>%
    mutate(id = cur_group_id()) %>% #guarantees correpondance id - lon lat is unique across files
    ungroup() %>%
    mutate(maxtp = round(maxtp, 2))
  
  if (!file.exists("Data/id_lon_lat_correspondance.csv")){
    #to ensure each row is unique
    id_loc <- all_data_with_id %>%
      select(longitude, latitude, id) %>%
      distinct()
    write.csv(id_loc, "Data/id_lon_lat_correspondance.csv", row.names = FALSE)
  }
  
  # Write only the date ID column in a CSV
  all_data_with_id %>%
    select(date_id, maxtp, id) %>%  # Choose the columns to save
    write.table(file.save,
                sep = ",",
                col.names = !file.exists(file.save),
                append = TRUE,
                row.names = FALSE)
}

#write.table(date_mapping, "Data/Climate_data/id_date_correspondance.csv", sep = ",", col.names = !file.exists(file.save), col.names = FALSE, row.names = FALSE, append = TRUE)