gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


num_quantiles = 30
obs_data = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>% 
  filter(stn == "AAR", year==1971)

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)%>%
  filter(year == 2000)

obs_data = obs_data %>% 
  left_join(glob_anomaly_reshaped, by = "year")

# ---- get covariates for prediction
temporal_covariates = glob_anomaly_reshaped %>% 
  dplyr::select(year, glob_anom) %>% 
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

quant_reg_model_pars = read_csv(paste0("Data/processed/glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

    
# --- get the climate quantile estimates closest to current station
clim_vals <<- obs_data %>%
  dplyr::select(quantile, value) %>%
  unique() %>%
  pull(value) %>%
  unlist()
    
# predict (observed) quantile for each year and site
quant_reg_pars = quant_reg_model_pars %>%
  arrange(tau)
    
res = c()
for(q in seq_along(quantiles_to_estimate_bulk)){
  qpars = quant_reg_pars[q,]
  # Calculate the predicted quantile value using the regression parameters
  res = rbind(res,
              tibble(quantile =  qpars$tau,
                      year = temporal_covariates$year,
                      quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2) *(temporal_covariates$glob_anom)))
}

obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_",num_quantiles,".csv"))%>% 
  filter(stn == "AAR", year==2000)

fun = obs_smoothed_quantiles%>%
  pull(tau_to_temp)

fun = fun[[1]]

quantile_seq =  seq(0.001,0.99,length.out = 150)


par(mar = c(5, 6, 4, 2))

plot(res$quantile, res$quant_value, 
     xlab = expression(tau),
     ylab = expression(q[o]^{tau}),
     cex=1.5, 
     col="red",
     cex.lab = 1.25,
     pch=19
)

points(quantile_seq, fun(quantile_seq), col = "blue", pch = 19, cex=0.75)

legend("topleft",                    # Position of the legend
       legend = c("interpolating quantiles", "spline function"),  # Labels for the data
       col = c("red", "blue"),        # Colors of the points/lines
       pch = c(19, 19),                  # Shape of the points (1 = circle, 19 = solid circle)
       bty = "n",                       # Remove border around the legend
       cex = 1.2)                       # Size of the legend text


gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

L = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")

num_quantiles = 30

obs_data = read.csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv") %>%
  mutate(year=lubridate::year(date))

threshold_9_df = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv")) %>%
  select(stn, threshold_9) %>%
  filter(stn %in% c("MAG", "JUN", "RIE", "LSN"))%>%
  unique()

threshold_MAG = threshold_9_df[threshold_9_df$stn == "MAG", ]$threshold_9[[1]]
threshold_JUN = threshold_9_df[threshold_9_df$stn == "JUN", ]$threshold_9[[1]]
threshold_RIE = threshold_9_df[threshold_9_df$stn == "RIE", ]$threshold_9[[1]]
threshold_LSN = threshold_9_df[threshold_9_df$stn == "LSN", ]$threshold_9[[1]]

MAG_df <- obs_data[obs_data$stn == "MAG", ]
JUN_df <- obs_data[obs_data$stn == "JUN", ]
RIE_df <- obs_data[obs_data$stn == "RIE", ]
LSN_df <- obs_data[obs_data$stn == "LSN", ]


# Plotting histogram of col3 in the filtered data
hist(MAG_df$maxtp,
     breaks = 30,
     main = "Magadino/Cadenazzo, 203 m",
     xlab = "max tp (°C)",
     col = "skyblue",
     border = "black")



