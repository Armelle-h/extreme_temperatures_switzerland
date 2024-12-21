setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(rnaturalearth) #library for map of switzerland
library(sf)
library(tidyverse)

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

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

#joining the location id as defined in climate data

obs_file <- read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")#read.csv("Data/Observed_data/1971_2023_JJA_obs_data.csv")
legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv") %>%
  select(c("stn", "longitude", "latitude"))

obs_file = obs_file %>%
  select(-id)%>%
  left_join(legend, by = "stn")


#loading the location id correspondance file
id_loc = read.csv("Data/id_lon_lat_correspondance.csv")

obs_loc = unique(obs_file[, c("longitude", "latitude")])

#creating column where the location id will be saved
obs_loc$id <- NA
obs_loc$dist <- NA

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
  obs_loc$dist[i] <- min(distances)
}

obs_file_with_id <- merge(obs_file, obs_loc[, c("longitude", "latitude", "id", "dist")], 
                          by = c("longitude", "latitude"), 
                          all.x = TRUE)

obs_file_with_id = obs_file_with_id %>%
  filter(dist < 0.0089) #keeping only stations for which we have climate data less than 1 km close.

#removing the columns associated with the longitude and latitude
obs_file_with_id <- obs_file_with_id %>%
  select(-longitude, -latitude)

#saving in a csv file
write.csv(obs_file_with_id, "Data/Observed_data/1971_2023_JJA_obs_data_loc_id.csv", row.names=FALSE)

#year 2023 is removed in the file "pre_processing/removing_2023.R"

#detecting outliers 

obs_data = read.csv("Data/Observed_data/Unfiltered/1971_2022_JJA_obs_data_loc_id.csv")

obs_data = unique(obs_data)

#print(nrow(obs_data))

#need to inspect by hand the below to investigate if the value is absurd or not 
obs_data_absurd <- obs_data %>%
  group_by(stn) %>%
  filter(maxtp > quantile(maxtp, 0.75) + 3 * IQR(maxtp) |
           maxtp < quantile(maxtp, 0.25) - 3 * IQR(maxtp))

#remove absurd values 

absurd_stations = c("ATT", "CMA", "DIA", "DUB", "TICOM")

obs_data_absurd = obs_data_absurd %>%
  filter(stn %in% absurd_stations)

obs_data_filtered = anti_join(obs_data, obs_data_absurd)

#removing year 2016 of TICAB 

obs_data_filtered = obs_data_filtered %>%
  mutate(year = lubridate::year(date))%>%
  filter(!(stn == "TICAB" & year==2016))

#removing the station TIT 

obs_data_filtered = obs_data_filtered %>%
  filter(stn != "TIT")

#removing June 2015 for TICOM, remove TICAB on the 2015-08-13, remove BEKAP in 2011-07-04, after observation from the qqplot seemed unrealistic

obs_data_filtered = obs_data_filtered %>%
  mutate(month = lubridate::month(date))%>%
  filter(!(stn == "TICOM" & month == 06))%>%
  filter(!(stn == "TICAB" & date == "2015-08-13"))%>%
  filter(!(stn == "BEKAP" & date == "2011-07-04"))

cat("Number of data removed:", nrow(obs_data_filtered), " over ", nrow(obs_data), "equivalently ", (1-nrow(obs_data_filtered)/nrow(obs_data))*100, "%")

obs_data_filtered = obs_data_filtered %>%
  select(-c(year, month))

write.csv(obs_data_filtered, "Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv", row.names=FALSE)

#updating the legend file 

legend = read.csv("Data/Observed_data/Unfiltered/1971_2022_JJA_obs_legend.csv")

legend_filtered <- legend %>%
  filter(stn %in% obs_data_filtered$stn)

write.csv(legend_filtered, "Data/Observed_data/1971_2022_JJA_obs_legend.csv", row.names = FALSE)


#removing stations outside of Switzerland 

legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

# Step 3: Convert the dataframe to an sf object, assuming coordinates in WGS84 (EPSG:4326)
legend_sf <- st_as_sf(legend, coords = c("longitude", "latitude"), crs = 4326)

# Step 4: Ensure the CRS for both datasets match (WGS84)
switzerland <- st_transform(switzerland, crs = 4326)

# Step 5: Check if points are within Switzerland
points_in_switzerland <- st_intersects(legend_sf, switzerland, sparse = FALSE)

# Step 6: Filter the dataframe to keep only the rows with points inside Switzerland
legend_within_switzerland <- legend[points_in_switzerland, ]

filtered = legend %>%
  anti_join(legend_within_switzerland)

write.csv(legend_within_switzerland, "Data/Observed_data/1971_2022_JJA_obs_legend.csv", row.names = FALSE)

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

obs_data_filtered = obs_data %>%
  filter(stn %in% legend_within_switzerland$stn)

write.csv(obs_data_filtered, "Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv", row.names=FALSE)

#combining similar stations

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

obs_data_filtered <- obs_data %>%
  mutate(stn = if_else(stn == "AGAAR", "AAR", stn))%>%
  mutate(stn = if_else(stn == "BZN", "BEZ", stn))%>%
  mutate(stn = if_else(stn == "WSLLAE", "NABLAE", stn))%>%
  mutate(stn = if_else(stn == "TGFRA", "FRF", stn))%>%
  mutate(stn = if_else(stn == "MMHIR", "HIR", stn))

#deleting stations
obs_data_filtered = obs_data_filtered %>%
  filter( !(stn %in% c("WSLBTB", "WSLHOB", "WSLBAB", "WSLCLB", "NABDAV", "WSLISB", "WSLISF")) )


legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

legend_filtered <- legend %>%
  filter(stn %in% obs_data_filtered$stn)

write.csv(legend_filtered, "Data/Observed_data/1971_2022_JJA_obs_legend.csv", row.names = FALSE)

write.csv(obs_data_filtered, "Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv", row.names = FALSE)



