setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

library(ggplot2) #library for visualisation
library(rnaturalearth) #library for map of switzerland
library(tidyverse)
library(sf)

raw_data<- read.csv("Data_climate/1940_2023_data.csv", header=TRUE)
colnames(raw_data)[3] <- "tp" #renaming last column
data <- raw_data %>%
  filter(tp != "-")#filtering out rows with missing values
data$tp <- as.integer(data$tp) #converting last column elements as integers

legend<- read.csv("Data_climate/1940_2023_legend.csv", header=TRUE)

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

#count nb of measurement associated with each station
nb_meas_df <- data %>%
  group_by(stn) %>%
  summarise(count = n())

merged_data <- legend %>%
  left_join(nb_meas_df, by = "stn")

points_sf <- st_as_sf(merged_data, coords = c("longitude", "latitude"), crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = points_sf, aes(size = count, color = count), alpha=0.7) +  # Plot points sized by measurements
  scale_size_continuous(range = c(1, 4)) +  # Adjust size range
  scale_color_gradientn(name = "Number of Measurements", 
                        colors = c("blue", "purple", "red", "orange")) +
  ggtitle("Weather Stations in Switzerland") +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")