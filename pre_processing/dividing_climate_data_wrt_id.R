setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

rm(list = ls())
gc()
library(tidyverse)
library(data.table)

files = list.files(path = "Data/Climate_data", full.names = TRUE)
files = files[-length(files)]

col_names = TRUE

for (f in files) {

  df = fread(f)
  
  threshold_1 = 13773
  threshold_2 = 27547

  df1 <- df %>%
    filter(id <= threshold_1)

  df2 <- df %>%
    filter(id > threshold_1 & id<= threshold_2)

  df3 <- df %>%
    filter(id > threshold_2)

  write.table(df1, "Data/Climate_data/1_13773_JJA_climate_data.csv", sep = ",", col.names = col_names, row.names = FALSE, append = TRUE)

  write.table(df2, "Data/Climate_data/13774_27547_JJA_climate_data.csv", sep = ",", col.names = col_names, row.names = FALSE, append = TRUE)

  write.table(df3, "Data/Climate_data/27548_41319_JJA_climate_data.csv", sep = ",", col.names = col_names, row.names = FALSE, append = TRUE)

  col_names = FALSE
  
  rm(df, df1, df2, df3) #removing from memory to lessen the strain in R
  gc()

}