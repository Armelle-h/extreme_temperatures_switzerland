gc() 
rm(list = ls())

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(lubridate)  # For date manipulation

data<- read.csv("pre_processing/Data/1971_2022_data.csv", header=TRUE)

data$date <- as.Date(data$date)

# Calculate 99th quantile and filter the dataframe -- should be correct
data_quantiles <- data %>%
  group_by(stn) %>%
  summarize(q99 = quantile(tp, 0.99, na.rm = TRUE), .groups = 'drop') %>%
  inner_join(data, by = "stn") %>%
  filter(tp > q99)

#Extract month and count occurrences
station_monthly_count <- data_quantiles %>%
  mutate(month = month(date)) %>%
  group_by(stn, month) %>%
  summarise(count = n()) %>%
  arrange(stn, month)

monthly_count <- data_quantiles %>%
  mutate(month = month(date)) %>%
  group_by(month) %>%
  summarise(quant_count = n()) %>%
  arrange(month)

total_count <- nrow(data_quantiles)

monthly_perc <- monthly_count %>%
  mutate(percentage = (quant_count / total_count)) %>%
  dplyr::select(month, quant_count, percentage)

# Plot the histogram
ggplot(monthly_perc, aes(x = month, y = percentage)) +
  geom_col(fill = "skyblue") +
  labs(title = "99th quantile", 
       x = "Month", 
       y = "Proportion of extreme events") +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11),       
                     labels = c("January", "March", "May", "July", "September", "November")) +  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#percentage of exceedances occuring in the summer months
#percentage of exceedances occuring in the summer months
print( sum(monthly_perc[6:8, "quant_count", na.rm=TRUE])*100/sum(monthly_perc$quant_count)  ) #prints 97.79076


