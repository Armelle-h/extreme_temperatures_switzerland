gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

quant_reg_bts = read.csv(paste0("output/bts_quant_reg/bts_quant_reg_mod_2_num_quantiles_30.csv"), header = FALSE)
colnames(quant_reg_bts) = c("bts", "tau", "beta_0", "beta_1", "beta_2")

tau_list = sort(unique(quant_reg_bts$tau))

#histogram for the quantile regression 

for (tau_val in tau_list){
  file_path = paste0("uncertainty_pictures/beta_2_tau_", tau_val, ".png")
  
  png(filename = file_path, width = 443, height = 290)
  
  df = quant_reg_bts %>%
    filter(tau == tau_val)
  
  hist( df$beta_2, breaks = 25,
        main = paste0("beta_2, tau =", tau_val), 
        xlab = "beta_2", 
        ylab = "Frequency")
  
  dev.off()
}

#uncertainty

quant_reg_model_pars = read_csv(paste0("Data/processed/glob_anomaly_quantile_model_fit_pars_num_quantiles_30.csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

columns_to_compute <- c('beta_0', 'beta_1', 'beta_2')

quantile_0025 <- quant_reg_bts %>%
  group_by(tau) %>%
  summarise(across(all_of(columns_to_compute), 
                   ~ mean(., na.rm = TRUE) + sd(., na.rm = TRUE) * qnorm(0.025), 
                   .names = "quantile_0025_{col}"))

quantile_0975 <- quant_reg_bts %>%
  group_by(tau) %>%
  summarise(across(all_of(columns_to_compute), 
                   ~ mean(., na.rm = TRUE) + sd(., na.rm = TRUE) * qnorm(0.975), 
                   .names = "quantile_0975_{col}"))

#debugging, else issue with the left join
quantile_0025$tau = quant_reg_model_pars$tau
quantile_0975$tau = quant_reg_model_pars$tau

final_df = quant_reg_model_pars %>%
  left_join(quantile_0025, by = "tau")%>%
  left_join(quantile_0975, by = "tau")

ggplot(final_df, aes(x = tau)) +
  geom_line(aes(y = beta_0, color = "Beta_0")) +
  geom_ribbon(aes(ymin = quantile_0025_beta_0, ymax = quantile_0975_beta_0), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = beta_1, color = "Beta_1")) +
  geom_ribbon(aes(ymin = quantile_0025_beta_1, ymax = quantile_0975_beta_1), fill = "red", alpha = 0.2) +
  geom_line(aes(y = beta_2, color = "Beta_2")) +
  geom_ribbon(aes(ymin = quantile_0025_beta_2, ymax = quantile_0975_beta_2), fill = "green", alpha = 0.2) +
  scale_color_manual(values = c("Beta_0" = "blue", "Beta_1" = "red", "Beta_2" = "green")) +
  labs(color = "Beta Curves") +
  xlab("tau") +  
  ylab("values") +
  theme_minimal()

#exceedance proba

exceedance_proba_bts = data.frame(bts = numeric(), stn = character(), year = numeric(), ex_proba = numeric())

for (bts_num in seq(1, 100)){
  exceedance_proba_bts_num = read.csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_mod_2_num_quantiles_30_bts_",bts_num,".csv"), header=FALSE)
  
  exceedance_proba_bts = rbind(exceedance_proba_bts_num, exceedance_proba_bts)
}

colnames(exceedance_proba_bts) = c("bts", "stn", "year", "proba")

#threshold exceedance varies with year and stn 
stn_list = sort(unique(exceedance_proba_bts$stn))
year_list = sort(unique(exceedance_proba_bts$year))

set.seed(123)
selected_stn_list <- sample(stn_list, size = 10, replace = FALSE)
selected_year_list <- sample(year_list, size = 10, replace = FALSE)


for (stn_val in selected_stn_list) {
  for (year_val in selected_year_list) {

  file_path = paste0("uncertainty_pictures/exceedance_proba_", stn_val, "_",year_val,".png")
  
  png(filename = file_path, width = 443, height = 290)
  
  df = exceedance_proba_bts %>%
    filter(stn == stn_val, year == year_val)
  
  hist( df$proba, breaks = 25,
        main = paste0(stn_val,"_", year_val), 
        xlab = "exceedance_proba", 
        ylab = "Frequency")
  
  dev.off()
  }
}

#plotting the threshold uncertainty for a random sample of 30 years and 10 stations 

set.seed(456)
stn_list = sort(unique(exceedance_proba_bts$stn))
selected_stn_list <- sample(stn_list, size = 10, replace = FALSE)

legend = read.csv("Data/Observed_data/1971_2022_JJA_obs_legend.csv")%>%
  select(stn, Nom)%>%
  unique()

for (stn_val in selected_stn_list){

exceedance_proba_bts_filter = exceedance_proba_bts %>%
  filter(stn == stn_val)

thresh_exceedance = read_csv("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_30.csv")%>%
  filter(stn == stn_val & year %in% exceedance_proba_bts_filter$year )

columns_to_compute = c("proba")

quantile_0025 <- exceedance_proba_bts_filter %>%
  group_by(year, stn) %>%
  summarise(across(all_of(columns_to_compute), 
                   ~ mean(., na.rm = TRUE) + sd(., na.rm = TRUE) * qnorm(0.025), 
                   .names = "quantile_0025"))

quantile_0975 <- exceedance_proba_bts_filter %>%
  group_by(year, stn) %>%
  summarise(across(all_of(columns_to_compute), 
                   ~ mean(., na.rm = TRUE) + sd(., na.rm = TRUE) * qnorm(0.975), 
                   .names = "quantile_0975"))

final_df = thresh_exceedance %>%
  left_join(quantile_0025, by = c("stn", "year"))%>%
  left_join(quantile_0975, by = c("stn", "year"))

name_stn = legend %>%
  filter(stn == stn_val)%>%
  pull(Nom)

p = ggplot(final_df, aes(x = year)) +
  geom_line(aes(y = thresh_exceedance_9)) +  # Line plot for thresh_exceedance_9
  geom_ribbon(aes(ymin = `quantile_0025`, ymax = `quantile_0975`), alpha = 0.2) +  # Ribbon for quantiles
  labs(x = "Year", y = "Exceedance Probability", title = paste0(name_stn)) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(paste0("plot/exceedance_plot_", stn_val, ".png"), plot = p, width = 4, height = 3.5, dpi = 300, bg = "white")
}

#gpd parameters

compute_quantiles <- function(column) {
  mu <- mean(column)   
  sigma <- sd(column) 
  p <- c(0.025, 0.975)   
  quantiles <- mu + sigma * qnorm(p)
  return(quantiles)
}

fit_gpd_bts = read.csv("output/gpd_model_fits/bts/model_2_bts.csv", header = FALSE) %>%
  unique()
colnames(fit_gpd_bts) = c("bts", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4", "xi" )

file_path = paste0("uncertainty_pictures/gpd_xi.png")

png(filename = file_path, width = 443, height = 290)

#histogram for the gpd parameters 
hist( fit_gpd_bts$xi, breaks = 25,
     main = "xi", 
     xlab = "xi", 
     ylab = "Frequency")

dev.off()

#uncertainty gpd parameters beta_0, beta_1, beta_2, beta_3, beta_4, xi

print(compute_quantiles(fit_gpd_bts$beta_0))
print(compute_quantiles(fit_gpd_bts$beta_1))
print(compute_quantiles(fit_gpd_bts$beta_2))
print(compute_quantiles(fit_gpd_bts$beta_3))
print(compute_quantiles(fit_gpd_bts$beta_4))
print(compute_quantiles(fit_gpd_bts$xi))


