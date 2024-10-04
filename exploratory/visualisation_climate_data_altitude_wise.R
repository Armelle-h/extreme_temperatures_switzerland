setwd("C:/Users/HOURS/Desktop/PDM/Code_R")

library(ggplot2) #library for visualisation
library(rnaturalearth) #library for map of switzerland
library(tidyverse)
library(sf)

data<- read.csv("Data_climate/1940_2023_data.csv", header=TRUE)

legend<- read.csv("Data_climate/1940_2023_legend.csv", header=TRUE)

switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

points_sf <- st_as_sf(legend, coords = c("longitude", "latitude"), crs = 4326)

ggplot(data = switzerland) +
  geom_sf(fill = "lightblue", color = "black") +  # Plot Switzerland
  geom_sf(data = points_sf, aes(size = Altitude.m. , color = Altitude.m.), alpha=0.7) +  # Plot points sized by measurements
  scale_size_continuous(range = c(1, 4)) +  # Adjust size range
  scale_color_gradientn(name = "Altitude", 
                        colors = c("blue", "purple", "red", "orange")) +
  ggtitle("Weather Stations in Switzerland") +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")