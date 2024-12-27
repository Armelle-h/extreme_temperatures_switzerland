gc()
rm(list=ls()) 
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(gridExtra)

#loading the custom functions
source('gpd_models_censoring.R') #for now, might need to do a bit of relocating

num_quantiles = 30

fit_uncorrected_models = T
fit_true_models = T

#observed data that has exceeded the threshold
obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

obs_data_censored = obs_data %>%
  select(stn, date, maxtp) %>%
  group_by(stn)%>%
  mutate(q_0.999 = quantile(maxtp, 0.999))%>%
  filter(maxtp>q_0.999)%>%
  ungroup

legend_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

#joining obs data with the altitude 
obs_data = obs_data %>%
  left_join(legend_data %>% select(stn, Altitude.m.), by="stn")%>%
  rename(altitude = Altitude.m.)

#adding threshold 9 

threshold_9_df = vroom::vroom("Data/processed/1971_2022_JJA_obs_data_bulk_model.csv")%>%
  dplyr::select(threshold_9, id)%>%
  unique()

#lambda associated with observed data, vary the quantile model
lambda_df = vroom::vroom(paste0("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv") )

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(lambda_df, by=c("stn", "year")) %>%
  left_join(glob_anomaly_reshaped, by = "year") %>%
  left_join(threshold_9_df, by="id")

#loading the covariates
covars = obs_data %>%
  dplyr::select(stn, year, scale_9, glob_anom, thresh_exceedance_9, threshold_9, altitude) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  
  #extracting the excess data
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)
  
  extreme_dat_true_uncensor = extreme_dat_true %>%
    anti_join(obs_data_censored%>%select(stn, date, maxtp), by = c("stn", "date"))
  
  extreme_dat_true_censor = extreme_dat_true %>%
    semi_join(obs_data_censored, by = c("stn", "date"))%>% #to only keep rows that will be censored
    left_join(obs_data_censored, by= c("stn", "date", "maxtp"))%>% #to add the quant_0.999 column
    mutate(excess_q_0.999 = q_0.999-threshold_9)
  
  #only doing for model 1
  
  fit_mod_1(extreme_dat_true_uncensor$excess, extreme_dat_true_uncensor$scale_9,
            extreme_dat_true_uncensor$glob_anom, extreme_dat_true_censor$excess_q_0.999, extreme_dat_true_censor$scale_9,
            extreme_dat_true_censor$glob_anom, c(0.6, -0.1, 0.5, -0.2))  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/censor_model_1_true.csv")
}

#DOING THE QQPLOT 

obs_data = obs_data %>%
  group_by(stn) %>%
  mutate(
    quantile_999 = quantile(maxtp, 0.999, na.rm = TRUE),  # Calculate the 0.999 quantile
    maxtp = pmin(maxtp, quantile_999)  # Replace values above the quantile
  ) %>%
  ungroup %>%
  select(-quantile_999)  # Remove the temporary column
  

this_fit_mod = read_csv("output/gpd_model_fits/censor_model_1_true.csv") %>%
  unlist() %>% as.numeric

pred <- my_predict_1(this_fit_mod, obs_data$scale_9, obs_data$glob_anom)

obs_data$scale = pred$scale
obs_data$shape = pred$shape

standardised_qq = obs_data %>%
  filter(maxtp>=threshold_9) %>%
  mutate(unif = evd::pgpd(q = (maxtp - threshold_9), loc = 0, scale = scale, shape = shape[1])) %>%
  mutate(exp = -log(1 - unif)) %>%
  mutate(rank = seq(nrow(.))/(nrow(.)+1)) %>%
  mutate(rank = (-log(1-rank)))

num_reps = nrow(standardised_qq)
exp_ci = c()
for(i in seq(250)){
  exp_ci = rbind(exp_ci, sort(rexp(n = num_reps)))
}

lower_ci = exp_ci %>% apply(MARGIN = 2, FUN = quantile, 0.02)
upper_ci = exp_ci %>% apply(MARGIN = 2, FUN = quantile, 0.98)

library(scales)

#!! if there are dupplicate values in maxtp (which there are) they will have the same rank and ideal won't be unique. To fix that, 
#the argument "random" randomly assigns an order to similar maxtp values
obs_data$ideal = (obs_data$maxtp %>% rank(ties.method = "random"))/(nrow(obs_data)+1) 
obs_data$exp_ideal = -log(1-obs_data$ideal)

library(grid)

plt = gridExtra::grid.arrange(standardised_qq %>%
                                ggplot()+
                                geom_ribbon(data = tibble(qnt = standardised_qq$rank, (lower_ci), (upper_ci)), aes(x = qnt, ymin = sort(lower_ci), ymax = sort(upper_ci)), alpha = 0.3)+
                                geom_point(aes(rank, sort(exp)), size = 0.75)+
                                geom_abline(col = 'red',linetype = 'longdash')+
                                theme_minimal(12)+
                                theme(axis.title.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      panel.grid.minor = element_blank()) + 
                                scale_x_continuous(breaks = c(0,2,4,6,8,10))+ 
                                scale_y_continuous(breaks = c(0,2,4,6,8,10,12)), nrow = 1)
