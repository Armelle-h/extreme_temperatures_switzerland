#in this file, I am counting the number of stations from which we have temperature measurements for each year and represent the result in a histogram.
#I remove rows associated with years <=1970

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)


data<- read.csv("pre_processing/Data/1940_2022_data.csv", header=TRUE)

data$date <- as.Date(data$date)

#count nb of stations
unique_count <- length(unique(data[[1]]))

#count nb of measurement associated with each station
count_df <- data %>%
  group_by(stn) %>%
  summarise(count = n()) #no need to incorporate in the code "oh I want to see in increasing order"

nb_stat_per_year <- data %>%
  mutate(year = format(date, "%Y")) %>%          #extract year
  group_by(year) %>%                                # Group by year
  summarise(unique_count = n_distinct(stn)) %>%  # Count unique occurrences in column2
  ungroup()  

ggplot(nb_stat_per_year, aes(x = as.numeric(year), y = unique_count)) +
  geom_bar(stat = "identity", fill = "blue") +      # Use bars for the histogram
  labs(title = "Nb of stations per year", 
       x = "year", 
       y = "nb stations") +
  scale_x_continuous(breaks = seq(1950, 2010, by = 10)) +
  theme_minimal() 

obs <- obs %>%
  filter(format(date, "%Y") > "1970")
write.csv(obs, "pre_processing/Data/1971_2022_data.csv", row.names=FALSE)