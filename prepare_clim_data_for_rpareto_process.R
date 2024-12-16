gc()
rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(data.table)

filter_id = 25 #else too computationally expensive

num_quantiles = 30

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

clim_1 = fread("Data/Climate_data/By_year/1971_1990_JJA_climate_data.csv")%>%
  filter(id %in% I_plain$id)%>%
  filter(id %% 25 == 0)

clim_2 = fread("Data/Climate_data/By_year/1991_2017_JJA_climate_data.csv")%>%
  filter(id %in% I_plain$id)%>%
  filter(id %% 25 == 0)

clim_3 = fread("Data/Climate_data/By_year/2018_2022_JJA_climate_data.csv")%>%
  filter(id %in% I_plain$id)%>%
  filter(id %% 25 == 0)

clim_data = rbind(clim_1, clim_2, clim_3)

rm(clim_1, clim_2, clim_3)
gc()

clim_data$ideal = (clim_data$maxtp %>% rank(ties.method = "random"))/(nrow(clim_data)+1) 

clim_thresh = clim_data %>%
  group_by(id) %>%
  summarise(clim_thresh_value_9 = quantile(maxtp, 0.9, na.rm = TRUE))
#first converting clim_data to pareto

clim_quantiles_subset = readRDS(paste0("Data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))%>%
  filter(id %% 25 == 0)%>%
  filter(id %in% I_plain$id)

clim_grid = read_csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")%>%
  filter(id %% 25 == 0)

clim_smooth = clim_quantiles_subset %>%
  group_by(id) %>%
  group_map(~{
    tibble( id = .x$id[[1]],
            tau_to_temp = list(splinefun(unlist(.x$quantile),unlist(.x$value),  method = 'monoH.FC')),
            temp_to_tau = list(splinefun(unlist(.x$value),unlist(.x$quantile),  method = 'monoH.FC')))
  }, .keep = T)%>%
  plyr::rbind.fill() %>%
  as_tibble() 

# Calculate lambda (exceedance probability) for the threshold
# (difference compared to before is that now we're working with the regression estimated climate quantile)
lambda_thresh_ex = clim_thresh %>%
  group_by(id) %>%
  group_map(~{
    
    thresh_exceedance_9 = clim_smooth%>%
      filter(id == .x$id[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$clim_thresh_value_9[1], x))# Apply the threshold function
    
    # Create a tibble with exceedance probabilities
    tibble(id = .x$id[1],
           thresh_exceedance_9 = 1-thresh_exceedance_9)# Inverse exceedance probability  
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

shape = readRDS("Data/clim_data_gpd_model/plain_optimal_shape.rds")

I_date = read.csv("Data/id_date_correspondance.csv")%>%
  mutate(year = lubridate::year(date))

clim_data = clim_data %>%
  left_join(I_date, by = "date_id")%>%
  select(-date_id)%>%
  left_join(lambda_thresh_ex, by = "id")%>%
  left_join(clim_grid, by = "id")%>%
  left_join(clim_thresh, by = "id")

rm(lambda_thresh_ex, clim_thresh, clim_grid)

clim_data_standardised = clim_data %>%  #will have to do a for loop over the clim data files
  group_by(id, year) %>% #the function is just defined with respect to space, not time
  group_map(~{
    data = .x$maxtp
    threshold = .x$clim_thresh_value_9
    res = rep(NA, length(data))
    num_extremes = sum(data>threshold)
    
    if(num_extremes >0){ # if extreme temperat in this year at this site
      scle = .x$scale_9
      shpe = shape
      my_lambda = .x$thresh_exceedance_9
      
      res[data > threshold] = 1 - my_lambda[data > threshold]*(1-evd::pgpd((data[data > threshold] - threshold[data > threshold]), loc = 0,scale = scle[data > threshold], shape = shpe[1]))
    }
    
    this_quant_mod = clim_smooth %>%
      filter(id == .x$id[[1]]) %>%
      pull(temp_to_tau) %>%
      .[[1]]
    
    res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
    res[res<0]=0
    res[res>0.999999]=0.999999
    
    .x$unif = res
    #.x$frechet_marg = -1/log(.x$unif)
    .x$pareto_marg = 1/(1-.x$unif)
    .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble() 


# -------- MULTIVARIATE EVENTS

# get cost of each event
extreme_dates = clim_data_standardised %>%
  dplyr::group_by(date) %>%
  #dplyr::summarise(cost = mean(frechet_marg)) %>%
  dplyr::summarise(cost = median(pareto_marg)) %>%  #dplyr::summarise(cost = mean(pareto_marg)) %>%
  ungroup() %>%
  arrange(desc(cost))

tot_dates =nrow(extreme_dates)

# get cost threshold
threshold = quantile(extreme_dates$cost, 0.8) %>% as.numeric   #robust to outliers

# get extreme dates
extreme_dates = extreme_dates %>%
  filter(cost > threshold) %>%
  arrange(desc(cost))

cat("cr =", nrow(extreme_dates)*threshold/tot_dates)

# temporally decluster events
my_data_tmp = extreme_dates %>% arrange(date)

my_data_tmp$date <- as.Date(my_data_tmp$date)
min_dif = my_data_tmp %>% arrange(date) %>% pull(date) %>% diff() %>% min

print(my_data_tmp)
print("declustering events")

#removes events within 7 days of each other and retain only the highest cost event in a cluster
while(min_dif<7){
  i<-1
  while(i < nrow(my_data_tmp)){
    if(diff(c(my_data_tmp[i,]$date, my_data_tmp[i+1,]$date))<7){
      # pick the largest
      if(my_data_tmp[i,]$cost > my_data_tmp[i+1,]$cost){
        my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i+1,]$date,]
      }else{
        my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i,]$date,]
      }
    }
    i<-i+1
  }
  min_dif <- my_data_tmp$date %>% diff() %>% min %>% abs()
}
extreme_dates = my_data_tmp

# tibble of extreme events
clim_data_standardised = clim_data_standardised %>%
  filter(date %in% extreme_dates$date) %>%
  arrange(date)

print("getting data in correct format")

# get observations in list
exceedances = clim_data_standardised %>%
  group_by(date) %>%
  group_map(~{
    
    # our "conditional" site at the top of the list
    if(13125 %in% .x$id){ #station closest to KOP in clim data
      c(.x %>% filter(id == 13125) %>% pull(pareto_marg),
        .x %>% filter(id != 13125) %>% pull(pareto_marg))
      #c(.x %>% filter(id == 13117) %>% pull(frechet_marg),
      #  .x %>% filter(id != 13117) %>% pull(frechet_marg))
    }
  })

# remove observatiosn that dont have valentia_observatory
to_remove = exceedances %>% map(is.null) %>% unlist %>% which()

exceedances_locs = clim_data_standardised %>%
  left_join(I_plain %>% select(id, longitude_proj, latitude_proj), by = "id" )%>%
  group_by(date) %>%
  group_map(~{
    if(13125 %in% .x$id){
      rbind(.x %>% filter(id == 13125) %>% dplyr::select(longitude_proj, latitude_proj) %>% as.matrix(),  #used to be select(id) and not select(longitude, latitude)
            .x %>% filter(id != 13125) %>% dplyr::select(longitude_proj, latitude_proj ) %>% as.matrix())
    }
  })

if(!is.na(to_remove[1])){
  exceedances = exceedances[-to_remove]
  exceedances_locs = exceedances_locs[-to_remove]
}

dates_to_rem = extreme_dates %>% arrange(date) %>% .[to_remove,] %>% pull(date)

extreme_dates = extreme_dates %>% arrange(date)

list(exceedances = exceedances,
     exceedances_locs = exceedances_locs,
     thresh = threshold,
     extreme_dates = extreme_dates) %>%
  saveRDS("Data/processed/clim_data_for_rpareto/true/robust_plain_data_for_rpareto")

#doing the qq plot 

library(grid)

gridExtra::grid.arrange(clim_data_standardised %>%
                          drop_na() %>%
                          ggplot()+
                          geom_point(aes(sort(ideal), sort(unif)), size = 0.75)+
                          geom_abline(col = 'red',linetype = 'longdash')+
                          geom_vline(xintercept = 0.9)+
                          theme_minimal(12)+
                          theme(panel.grid.minor = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank())
)