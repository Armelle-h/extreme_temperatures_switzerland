rm(list = ls())
library(tidyverse)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

site_pairs = read_csv("Data/processed/plain_obs_pairs_with_dist.csv")
site_pairs$V2 %>% unique %>% length

num_quantiles = 30

my_qnt=0.8

marg_mod = "mod_1"

empirical_chi_scale = function(my_qnt, marg_mod){
  
  site_pairs = read_csv("Data/processed/plain_obs_pairs_with_dist.csv")

  standardised_data = read.csv(paste0("Data/processed/plain_obs_data_pareto_frechet_scale_", marg_mod,".csv"))%>%
    select(stn, year, date, unif)
    
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
}

#the block takes 7 min 
job::job({
  emp_chi_mod_0_08 = empirical_chi_scale(0.8, "mod_0")
  })
job::job({
  emp_chi_mod_0_085 = empirical_chi_scale(0.85, "mod_0")
})
job::job({
  emp_chi_mod_0_09 = empirical_chi_scale(0.9, "mod_0")
})


job::job({
  emp_chi_mod_1_08 = empirical_chi_scale(0.8, "mod_1")
})
job::job({
  emp_chi_mod_1_085 = empirical_chi_scale(0.85, "mod_1")
})
job::job({
  emp_chi_mod_1_09 = empirical_chi_scale(0.9, "mod_1")
})


emp_chi_mod_0_08 %>%
  ggplot(aes(x = dist, y = chi)) +
  geom_point(size = 3) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "p = 0.8") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),  # Increase title size
    axis.title = element_text(size = 18),# Increase axis label size
    axis.text = element_text(size = 16)# Increase axis number size
  )