#file to delete

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)

obs = read.csv("pre_processing/Data/1940_2022_data.csv")

# Step 2: Convert 'yyyymmdd' to 'yyyy-mm-dd' format

obs$date <- as.Date(obs$date)

# Step 3: Remove rows where the year is 2023
obs <- obs %>%
  filter(format(date, "%Y") > "1970")

write.csv(obs, "pre_processing/Data/1971_2022_data.csv", row.names=FALSE)