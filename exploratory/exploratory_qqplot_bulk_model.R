gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")%>%
  mutate(year=lubridate::year(date))

obs_smoothed_quantiles = readRDS(paste0("output/glob_anom_indicator_alt_quant_models_num_quantiles_",num_quantiles,".csv"))%>%
  select(-tau_to_temp)

obs_smoothed_quantiles_2 = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))%>%
  select(-tau_to_temp)

#doing qq plot on both body and tail we ignore the threshold thing and the fact they don't follow the same distribution for now.
obs_data$ideal = (obs_data$maxtp %>% rank(ties.method = "random"))/(nrow(obs_data)+1)

unif_body = obs_data %>%
  group_by(stn, year)%>%
  group_map(~{
  this_quant_mod = obs_smoothed_quantiles %>%
    filter(stn == .x$stn[1],
           year == .x$year[1]) %>%
    pull(temp_to_tau) %>%
    .[[1]]
  .x$unif = this_quant_mod(.x$maxtp) 
  .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble() 

unif_body_2 = obs_data %>%
  group_by(stn, year)%>%
  group_map(~{
    this_quant_mod = obs_smoothed_quantiles_2 %>%
      filter(stn == .x$stn[1],
             year == .x$year[1]) %>%
      pull(temp_to_tau) %>%
      .[[1]]
    .x$unif = this_quant_mod(.x$maxtp) 
    .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble()

unif_body_2 = unif_body_2 %>%
  rename(ideal_2 = ideal, unif_2 = unif)

combined_body = bind_cols(unif_body%>%select(c("stn", "year", "ideal", "unif")), unif_body_2%>%select(c("stn", "year", "ideal_2", "unif_2")))

combined_body_reduced = combined_body[seq(1, nrow(unif_body), by = 200), ]



ggplot(combined_body_reduced %>% drop_na()) +
  geom_point(
    aes(sort(ideal), sort(unif)),
    size = 0.75, col="blue", alpha = 0.75
  ) +
  geom_point(
    aes(sort(ideal_2), sort(unif_2)),
    size = 0.75, col = "green", alpha = 0.75
  ) +
  geom_abline(col = 'red', linetype = 'longdash') +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )




ggplot(unif_body_reduced %>% drop_na()) +
  geom_point(
    aes(sort(ideal), sort(unif)),
    size = 0.75
  ) +
  geom_abline(col = 'red', linetype = 'longdash') +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )






ggplot(combined_data_reduced, aes(x=ideal, y=unif, color = dataset))+
                                geom_point(size=0.8, alpha= 0.6)+
                                geom_abline(col = 'red',linetype = 'longdash')+
                                theme_minimal(12)+
                                theme(panel.grid.minor = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.title.x = element_blank())






plt = gridExtra::grid.arrange(combined_data_reduced %>%
                                drop_na() %>%
                                ggplot()+
                                geom_point(aes(sort(ideal), sort(unif), color = dataset), size = 0.9, alpha = 0.3)+
                                geom_abline(col = 'red',linetype = 'longdash')+
                                theme_minimal(12)+
                                theme(panel.grid.minor = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.title.x = element_blank())
)












plt = gridExtra::grid.arrange(unif_body %>%
                                drop_na() %>%
                                ggplot()+
                                geom_point(aes(sort(ideal), sort(unif)), size = 0.75)+
                                geom_abline(col = 'red',linetype = 'longdash')+
                                theme_minimal(12)+
                                theme(panel.grid.minor = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.title.x = element_blank())
)