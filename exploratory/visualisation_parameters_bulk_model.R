gc()
rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

num_quantiles = 30

quant_reg_model_pars = read_csv(paste0("Data/processed/glob_anomaly_quantile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

long_data <- quant_reg_model_pars %>%
  pivot_longer(cols = starts_with("beta_0"):starts_with("beta_2"), 
               names_to = "variable", values_to = "value")

# Create the plot
ggplot(long_data, aes(x = tau, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(x = "tau", y = "Values") +
  theme_minimal()

















plot(quant_reg_model_pars$tau, quant_reg_model_pars$beta_0, type = "l", col = "red", lwd = 2, 
     xlab = "tau", ylab = "Values")

lines(quant_reg_model_pars$tau, quant_reg_model_pars$beta_1, col = "blue", lwd = 2)

lines(quant_reg_model_pars$tau, quant_reg_model_pars$beta_2, col = "green", lwd = 2)

# Add a legend
legend("topright", legend = c("beta_0", "beta_1", "beta_2"), 
       col = c("red", "blue", "green"), lwd = 2)