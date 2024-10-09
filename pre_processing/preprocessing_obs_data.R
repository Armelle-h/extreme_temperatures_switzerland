setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

data <- read_csv("Data/Observed_data/1940_2023_data.csv")
data$tp <- as.integer(data$tp)

#renaming columns
colnames(data)[2:3] <- c("date", "maxtp")

# turning date into date format
data$date <- as.Date(as.character(data$date), format = "%Y%m%d")

#removing months not in June July August and years below 1971
data <- data %>%
  filter(
    (format(date, "%m") %in% c("06", "07", "08")) & (as.numeric(format(date, "%Y")) >= 1971)
  )

#adding a latitude and longitude column

legend <- read_csv("Data/Observed_data/1940_2023_legend.csv") 

data <- data %>%
  left_join(legend %>% select(stn, longitude, latitude), by = "stn")

# Removing stations in legend not appearing in data

legend <- legend %>%
  filter(stn %in% data$stn)

write.csv(data, "1971_2023_JJA_data.csv", row.names = FALSE)
write.csv(legend, "1971_2023_JJA_legend.csv", row.names = FALSE)

#from here ------------------------------------------------------------------

setwd("C:/Users/HOURS/Desktop/PDM/Code_R")
library(tidyverse)

#joining the location id as defined in climate data

obs_file <- read.csv("Data/Observed_data/1971_2023_JJA_obs_data.csv")

#loading the location id correspondance file
id_loc = read.csv("Data/id_lon_lat_correspondance.csv")

obs_loc = unique(obs_file[, c("longitude", "latitude")])

#creating column where the location id will be saved
obs_loc$id <- NA

#finding the correct id by associating the id of the pair latitude-longitude that is closest
#according to the euclidean distance

# Function to compute the Euclidean distance between two points
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Loop over each row in obs_loc
for (i in 1:nrow(obs_loc)) {
  # Extract the current pair from obs_loc
  x1 <- obs_loc$longitude[i]
  y1 <- obs_loc$latitude[i]
  
  # Compute the distances between this pair and all pairs in id_loc
  distances <- mapply(euclidean_distance, x1, y1, id_loc$longitude, id_loc$latitude)
  
  # Find the index of the minimum distance
  min_index <- which.min(distances)
  
  # Assign the corresponding value from id_loc$id to obs_loc
  obs_loc$id[i] <- id_loc$id[min_index]
}

obs_file_with_id <- merge(obs_file, obs_loc[, c("longitude", "latitude", "id")], 
                          by = c("longitude", "latitude"), 
                          all.x = TRUE)

#removing the columns associated with the longitude and latitude
obs_file_with_id <- obs_file_with_id %>%
  select(-longitude, -latitude)

#saving in a csv file
write.csv(obs_file_with_id, "Data/Observed_data/1971_2023_JJA_obs_data_loc_id.csv", row.names=FALSE)
