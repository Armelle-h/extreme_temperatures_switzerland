rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

site_pairs = read_csv("Data/processed/plain_obs_pairs_with_dist.csv")
site_pairs$V2 %>% unique %>% length

num_quantiles = 30

my_qnt=0.8
  
  site_pairs = read_csv("Data/processed/plain_obs_pairs_with_dist.csv")
    
  #might not need all the above
  standardised_data = readRDS(paste0("Data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles))%>%
    select(stn, year, date, unif_1)%>%
    rename(unif = unif_1)
    
  dates_to_keep = standardised_data %>%
    group_by(date) %>%
    summarise(check_obs_for_s = 'KOP' %in% stn) %>%
    filter(check_obs_for_s) %>%
    pull(date)
    
  this_data <<- standardised_data %>% filter(date %in% dates_to_keep)
  this_data$order = seq(nrow(this_data)) # remember the original order of the data
    
  chi.emp=function(U,V,q ){
    sum((U>=q)&(V>=q))/sum((U>=q))
  }
  
    empirical_chi = site_pairs %>%
      mutate(bin = ntile(dist, 25)) %>% # used to be 25
      group_by(bin) %>%
      group_map(~{
        # Print current bin
        print(.x$bin[1])
        
        # Initialize vectors for weighted values
        x1_weighted = c()
        x2_weighted = c()
        
        # Loop over each row in the current bin
        for (s in seq(nrow(.x))) {
          s1 = standardised_data %>% filter(stn == .x[s,]$V1)
          s2 = standardised_data %>% filter(stn == .x[s,]$V2)
          
          # Extract overlapping dates and pull uniform values
          x1_weighted = c(x1_weighted, s1 %>% filter(date %in% s2$date) %>% arrange(date) %>% pull(unif))
          x2_weighted = c(x2_weighted, s2 %>% filter(date %in% s1$date) %>% arrange(date) %>% pull(unif))
        }
        
        # Compute the output tibble for the current group
        tibble(
          bin = .x$bin[1],
          dist = median(.x$dist),  # Median distance for the bin
          chi = chi.emp(x1_weighted, x2_weighted, my_qnt),  # Chi value
          my_qnt = my_qnt  # Add my_qnt as a column
        )
      }, .keep = TRUE) %>%
      bind_rows()  # Combine all grouped tibbles into one

#chekc how to plot the empirical chi from what I just did (could check if it's worth it to save the result)
    
empirical_chi %>%
  ggplot(aes(x = dist, y = chi)) +
  geom_point() +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "p = 0.8")