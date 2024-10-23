gc()
rm(list = ls())
library(data.table)
library(tidyverse)

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

id_date = read.csv("Data_2023/id_date_correspondance.csv")
#after visual observation, realize we're keeping only rows where date_id <=4784

id_date_2022 = id_date %>%
  filter(format(as.Date(date), "%Y") != "2023")

write.csv(id_date_2022, "Data/id_date_correspondance.csv", row.names = FALSE)

obs_data = read.csv("Data_2023/Observed_data/1971_2023_JJA_obs_data_loc_id.csv")

obs_data_2022 = obs_data %>%
  filter(format(as.Date(date), "%Y") != "2023")

write.csv(obs_data_2022, "Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv", row.names = FALSE)

obs_data_2022 = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

legend_2023 = read.csv("2023/Data_2023/Observed_data/1971_2023_JJA_obs_legend.csv")

legend_2022 = legend_2023 %>%
  semi_join(obs_data_2022, by = "stn")

write.csv(legend_2022, "Data/Observed_data/1971_2023_JJA_obs_legend.csv", row.names = FALSE)

#moving on to the climate data 

climate_data_2023 = read.csv("Data_2023/Climate_data/2018_2023_JJA_climate_data.csv")

#computationaly faster than doing a join with elements in id_date_2022

climate_data_2022 = climate_data_2023 %>%
  filter(date_id <= 4784)

write.csv(climate_data_2022, "Data/Climate_data/2018_2022_JJA_climate_data.csv", row.names = FALSE)

#climate data by id 

climate_data_2023_part_1 = fread("Data_2023/Climate_data/By_id/1_13773_JJA_climate_data.csv")

climate_data_2022_part_1 = climate_data_2023_part_1%>%
  filter(date_id <= 4784)

fwrite(climate_data_2022_part_1, "Data/Climate_data/By_id/1_13773_JJA_climate_data.csv", row.names = FALSE)

climate_data_2023_part_2 = fread("Data_2023/Climate_data/By_id/13774_27547_JJA_climate_data.csv")

climate_data_2022_part_2 = climate_data_2023_part_2%>%
  filter(date_id <= 4784)

fwrite(climate_data_2022_part_2, "Data/Climate_data/By_id/13774_27547_JJA_climate_data.csv", row.names = FALSE)


climate_data_2023_part_3 = fread("Data_2023/Climate_data/By_id/27548_41319_JJA_climate_data.csv")

climate_data_2022_part_3 = climate_data_2023_part_3%>%
  filter(date_id <= 4784)

fwrite(climate_data_2022_part_3, "Data/Climate_data/By_id/27548_41319_JJA_climate_data.csv", row.names = FALSE)


glob_anom_2023 = read.csv("2023/Data_2023/global_tp_anomaly_JJA.csv")

glob_anom_2022 = glob_anom_2023%>%
  filter(year != "2023")

write.csv(glob_anom_2022, "Data/global_tp_anomaly_JJA.csv", row.names=FALSE)

















