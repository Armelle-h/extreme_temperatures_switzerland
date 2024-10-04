gc() 
rm(list = ls())


setwd("C:/Users/HOURS/Desktop/PDM/Code_R")
library(dplyr)
library(lubridate)  # For date manipulation
library(ggplot2)

raw_data<- read.csv("Data_climate/1940_2023_data.csv", header=TRUE)
colnames(raw_data)[3] <- "tp" #renaming last column
data <- raw_data %>%
  filter(tp != "-")#filtering out rows with missing values
data$tp <- as.integer(data$tp) #converting last column elements as integers

# Calculate 99th quantile and filter the dataframe -- should be correct
data_quantiles <- data %>%
  group_by(stn) %>%
  summarize(q99 = quantile(tp, 0.99, na.rm = TRUE), .groups = 'drop') %>%
  inner_join(data, by = "stn") %>%
  filter(tp > q99)

#Extract month and count occurrences
station_monthly_count <- data_quantiles %>%
  mutate(month = month(as.Date(as.character(time), format = "%Y%m%d"))) %>%
  group_by(stn, month) %>%
  summarise(count = n()) %>%
  arrange(stn, month)

monthly_count <- data_quantiles %>%
  mutate(month = month(as.Date(as.character(time), format = "%Y%m%d"))) %>%
  group_by(month) %>%
  summarise(quant_count = n()) %>%
  arrange(month)

# Merge monthly sums with count data to get total number of observations per month
#-- now is correct
#total_monthly_count <- data_quantiles %>%
#  mutate(month = month(as.Date(as.character(time), format = "%Y%m%d"))) %>%
#  group_by(month) %>%
#  summarise(total_count=n()) %>%
#  arrange(month)

total_count <- nrow(data_quantiles)

monthly_perc <- monthly_count %>%
  #left_join(total_monthly_count, by = "month") %>%
  mutate(percentage = (quant_count / total_count)) %>%
  #select(month, quant_count, total_count, percentage)
  dplyr::select(month, quant_count, percentage)

# Plot the histogram
ggplot(monthly_perc, aes(x = month, y = percentage)) +
  geom_col(fill = "skyblue") +
  labs(title = "Proportion extreme event, 99 quantile", 
       x = "Month", 
       y = "Proportions") +
  theme_minimal()

#need to fix the labels of the x axis-- to do later