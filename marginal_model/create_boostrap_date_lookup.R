# Description: This script samples blocks of time and stores them in a list

rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

#not for the plain, for the whole of switzerland
obs_data = read_csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")

is.subset = function(A, B) all(A %in% B)

sites_by_date = obs_data %>% 
  dplyr::select(stn, date) %>% 
  group_by(date)%>%
  group_map(~{
    tibble(date = .x$date[1], stns = list(.x$stn))
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

sites_by_date %>% 
  saveRDS(paste0("Data/processed/bootstrap_data/sites_by_date"))


run_bts_lookup = function(bts_id){
  # i. ========  ========  Global parameters ========
  obs_data = read_csv("Data/Observed_data/1971_2022_JJA_obs_data_loc_id.csv")
  
  is.subset = function(A, B) all(A %in% B)
  
  sites_by_date = readRDS("Data/processed/bootstrap_data/sites_by_date")
  
  # groups of observed sites and dates they were observed
  site_groups_and_dates = sites_by_date %>% 
    group_by(stns) %>%
    group_map(~{
      tibble(stns = .x$stns[1], dates = list(.x$date))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble()
  
  
  site_groups_and_dates_subsets = c()
  for(i in seq(nrow(site_groups_and_dates))){
    
    group_has_subset = map(.x = site_groups_and_dates$stns, 
                           .f = ~is.subset(A = site_groups_and_dates[i,]$stns[[1]], 
                                           .x)) %>% unlist()
    
    all_dates_which_sites_observed = site_groups_and_dates[group_has_subset,]$dates %>% unlist %>% unique() %>% lubridate::as_date()
    site_groups_and_dates_subsets = rbind(site_groups_and_dates_subsets,
                                          tibble(stns = site_groups_and_dates[i,]$stns,
                                                 dates = list(all_dates_which_sites_observed)))
  }
  
  
  dates_i_can_sample = function(sites_i_need){
    dts = site_groups_and_dates_subsets[(map(.x = site_groups_and_dates_subsets$stns, 
                                             .f = ~is.subset(A = sites_i_need, 
                                                             .x)) %>% unlist()),] %>%
      pull(dates) %>% unlist() %>% unique()
    
    
    if(is.null(dts)) NULL
    else dts %>% lubridate::as_date()
  }
  
  
  
  all_bootstraps = c()
  for(nnn in bts_id){
    print(paste0("on bootstrap number ", nnn))
    i=1
    
    dates_samples = c()
    while(length(dates_samples) <= nrow(sites_by_date)){
      
      # sample event size
      block_size = rgeom(n=1, prob = 0.2)
      sites_i_need = sites_by_date[c(i:(i+block_size)),]$stns %>% unlist() %>% unique()
      dts_to_samp = dates_i_can_sample(sites_i_need) 
      
      if(!(is.null(dts_to_samp) | (length(dts_to_samp) < block_size))){
        
        if(length(dts_to_samp) == block_size){
          dates_samples = c(dates_samples,dts_to_samp)
          
        }else{
          position_of_start_samp = sample(seq(length(dts_to_samp) - block_size), size = 1)
          dates_samples = c(dates_samples,dts_to_samp[position_of_start_samp:(position_of_start_samp+block_size)])
        }
        
        i=length(dates_samples)+1
      }
    }
    all_bootstraps = c(all_bootstraps, list(lubridate::as_date(dates_samples)[1:nrow(sites_by_date)]))
  }
  
  saveRDS(all_bootstraps, paste0("Data/processed/bootstrap_data/temp/bootstrapped_dates_",floor(runif(1)*100000)))
}

# --- takes 15  min
job::job({run_bts_lookup(seq(1,10))})
job::job({run_bts_lookup(seq(11,20))})
job::job({run_bts_lookup(seq(21,30))})
job::job({run_bts_lookup(seq(31,40))})
job::job({run_bts_lookup(seq(41,50))})
job::job({run_bts_lookup(seq(51,60))})
job::job({run_bts_lookup(seq(61,70))})
job::job({run_bts_lookup(seq(71,80))})
job::job({run_bts_lookup(seq(81,90))})
job::job({run_bts_lookup(seq(91,100))})

list.files("data/processed/bootstrap_data/temp/") %>%
  map(~{
    paste0('data/processed/bootstrap_data/temp/',.x) %>% print()
  })

all_bts = c()
files_to_read = list.files("data/processed/bootstrap_data/temp/") 

for(fff in files_to_read){
  all_bts <- c(all_bts, readRDS(paste0("data/processed/bootstrap_data/temp/", fff)))
}
saveRDS(all_bts, "data/processed/bootstrap_data/bootstrapped_dates")