
gc()
rm(list = ls())

library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

source('gpd_models.R')

num_quantiles = 30

obs_data = vroom::vroom("Data/Observed_data/obs_data_gpd_model.csv")

glob_anomaly = read.csv("Data/global_tp_anomaly_JJA.csv")

glob_anomaly_reshaped = glob_anomaly %>%
  select(c("year", "JJA"))%>%
  rename(glob_anom = JJA)

threshold_9_df = readRDS(paste0("Data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv")) %>%
  select(c(id, threshold_9))%>% #enough as threshold_9 varies only spatially 
  unique()
  
lambda_df = vroom::vroom(paste0("Data/processed/glob_anomaly_thresh_exceedance_lambda_num_quantiles_", num_quantiles, ".csv"))

#"stn"   "date"    "maxtp"  "id"  "scale_9"   "year" "thresh_exceedance_9" "threshold_9" "glob_anom"
obs_data = obs_data %>% 
  mutate(year = year(date)) %>%
  left_join(lambda_df, by = c("stn", "year")) %>%
  left_join(threshold_9_df, by = c("id"))%>%
  left_join(glob_anomaly_reshaped, by = "year")

# "year"        "tau_to_temp" "temp_to_tau" "stn"
obs_smoothed_quantiles = readRDS(paste0("output/glob_anomaly_quant_models_num_quantiles_", num_quantiles, ".csv"))%>%
  select(- tau_to_temp)

marg_mod = 'mod_1'

this_fit_mod = read_csv("output/gpd_model_fits/model_1_true.csv") %>%
  unlist() %>% as.numeric

pred <- my_predict_1(this_fit_mod, obs_data$scale_9, obs_data$glob_anom)


obs_data$scale = pred$scale
obs_data$shape = pred$shape

obs_data %>% group_by(stn) %>% summarise(count = n()) %>% arrange(count)
standardised_qq = obs_data %>% 
  dplyr::select(stn, year, maxtp, scale, shape, threshold_9) %>%
  filter(maxtp>=threshold_9) %>%
  mutate(unif = evd::pgpd(q = (maxtp - threshold_9), loc = 0, scale = scale, shape = shape[1])) %>%
  group_by(stn) %>%
  group_map(~{
    
    .x$exp = sort(-log(1 - .x$unif))
    .x$rank = seq(nrow(.x))/(nrow(.x)+1)
    .x$rank = sort(-log(1-.x$rank))
    
    # ---- tolerence band
    num_reps = nrow(.x)
    exp_ci = matrix(nrow = 1000, ncol = num_reps)
    for(i in seq(1000)){
      exp_ci[i,] = sort(-log(1-runif(num_reps)))
    }
    
    .x$lower_ci = exp_ci %>% apply(MARGIN = 2, FUN = quantile, 0.02)
    .x$upper_ci = exp_ci %>% apply(MARGIN = 2, FUN = quantile, 0.98)
    .x
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble() 



standardised_qq = obs_data %>%
  filter(maxtp>=threshold_9) %>%
  mutate(unif = evd::pgpd(q = (maxtp - threshold_9), loc = 0, scale = scale, shape = shape[1])) %>%
  mutate(exp = -log(1 - unif)) %>%
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


# ------ Standardised QQ plots body
standardised_body = rbind(obs_data %>%
                            filter(exp_ideal>4.5),
                          obs_data %>%
                            filter(exp_ideal<=4.5) %>%
                            sample_n(5000)) %>%
  filter(maxtp<threshold_9) %>%
  group_by(stn, year) %>%
  group_map(~{
    data = .x$maxtp
    threshold = .x$threshold_9
    res = rep(NA, length(data))
    
    this_quant_mod = obs_smoothed_quantiles %>%
      filter(stn == .x$stn[1],
             year == .x$year[1]) %>%
      pull(temp_to_tau) %>%
      .[[1]]
    
    res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
    res[data <= 0] = 0
    .x$unif = res
    .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble()

#### ---- bulk + tail
bulk_and_tail = rbind(obs_data %>%
                        filter(exp_ideal>4.5),
                      obs_data %>%
                        filter(exp_ideal<=4.5) %>%
                        sample_n(5000)) %>%
  group_by(stn, year) %>%
  group_map(~{
    
    data = .x$maxtp
    threshold = .x$threshold_9
    res = rep(NA, length(data))
    num_extremes = sum(data>threshold)
    
    if(num_extremes >0){ # if extreme obs in this year at this site
      scle = .x$scale
      shpe = .x$shape
      my_lambda = .x$thresh_exceedance_9
      
      res[data > threshold] = 1 - my_lambda[data > threshold]*(1-evd::pgpd((data[data > threshold] - threshold[data > threshold]), loc = 0,scale = scle[data > threshold], shape = shpe[1]))
    }
    
    this_quant_mod = obs_smoothed_quantiles %>%
      filter(stn == .x$stn[1],
             year == .x$year[1]) %>%
      pull(temp_to_tau) %>%
      .[[1]]
    
    res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
    res[res<0]=0
    
    .x$unif = res
    .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble() 


library(grid)

plt = gridExtra::grid.arrange(bulk_and_tail %>%
                                drop_na() %>%
                                ggplot()+
                                geom_point(aes(sort(ideal), sort(unif)), size = 0.75)+
                                geom_abline(col = 'red',linetype = 'longdash')+
                                geom_vline(xintercept = 0.9)+
                                theme_minimal(12)+
                                theme(panel.grid.minor = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.title.x = element_blank()),
                              standardised_qq %>%
                                ggplot()+
                                geom_ribbon(data = tibble(qnt = standardised_qq$rank, (lower_ci), (upper_ci)), aes(x = qnt, ymin = sort(lower_ci), ymax = sort(upper_ci)), alpha = 0.3)+
                                geom_point(aes(rank, sort(exp)), size = 0.75)+
                                geom_abline(col = 'red',linetype = 'longdash')+
                                theme_minimal(12)+
                                theme(axis.title.y = element_blank(),
                                      axis.title.x = element_blank(),
                                      panel.grid.minor = element_blank()) + 
                                scale_x_continuous(breaks = c(0,2,4,6,8,10))+ 
                                scale_y_continuous(breaks = c(0,2,4,6,8,10,12)), nrow = 1,
                              left = textGrob("Estimated\nquantiles", rot = 0, vjust = 0.5, 
                                              gp = gpar(col = "black", fontsize = 12)),
                              bottom = textGrob("Ideal quantiles", rot = 0, vjust = 0, 
                                                gp = gpar(col = "black", fontsize = 12)))
#ggsave(plt, filename = 'output/figs/qqplot_pooledsites.png', height = 3, width = 8)