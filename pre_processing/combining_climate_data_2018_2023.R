setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

library(tidyverse)

df1 = read.csv("Data/Climate_data/2018_JJA_climate_data.csv")

df2 = read.csv("Data/Climate_data/2019_JJA_climate_data.csv")

df3 = read.csv("Data/Climate_data/2020_2023_JJA_climate_data.csv")

combined_df = bind_rows(df1, df2, df3)

write.csv(combined_df, "Data/Climate_data/2018_2023_JJA_climate_data.csv", row.names = FALSE)