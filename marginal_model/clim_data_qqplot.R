gc()
rm(list = ls())

library(tidyverse)
library(data.table)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")


files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_data_extreme_9_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])
  
  #keeping only temperatures above the associated threshold 0.9 quantile
  sing_clim_data_extreme_9 = clim_data %>%
    group_by(id) %>%
    mutate(threshold = quantile(maxtp, 0.9),
           excess = maxtp - threshold) %>%
    filter(excess > 0) %>%
    ungroup()
  
  clim_data_extreme_9_list[[i]] = sing_clim_data_extreme_9
  
  #to free memory
  rm(clim_data)
  gc()
}

#date_id, maxtp, id, threshold, excess
clim_data_extreme_9 =  do.call(rbind, clim_data_extreme_9_list)

rm(clim_data_extreme_9_list)
rm(sing_clim_data_extreme_9)

#to have the scale parameter associated
scales_9_df = read.csv("Data/Climate_data/clim_scale_grid_gpd_model.csv")

shape_9 = readRDS("Data/clim_data_gpd_model/optimal_shape.rds")

clim_data_extreme_9 = clim_data_extreme_9 %>%
  left_join(scales_9_df, by="id")

clim_data_extreme_9$shape = shape_9

#reducing the number of data else computation time would explose 

clim_data_extreme_9_reduced = clim_data_extreme_9 %>%
  filter(id %% 100 == 0)

standardised_qq = clim_data_extreme_9_reduced %>%
  mutate(unif = evd::pgpd(q = (excess), loc = 0, scale = scale_9, shape = shape_9)) %>%
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
clim_data_extreme_9_reduced$ideal = (clim_data_extreme_9_reduced$maxtp %>% rank(ties.method = "random"))/(nrow(clim_data_extreme_9_reduced)+1) 
clim_data_extreme_9_reduced$exp_ideal = -log(1-clim_data_extreme_9_reduced$ideal)


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


library(sf)
library(rnaturalearth)

I = read.csv("Data/id_lon_lat_correspondance.csv")

my_pal = c( 
  '#062c30', # extra dark 
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')


switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)

T = standardised_qq %>%select(id, excess, exp)%>%group_by(id)%>% mutate(max_exp= max(exp) )%>% ungroup %>% inner_join(I, by="id")

#T = standardised_qq %>%select(id, excess, exp) %>% inner_join(I, by="id")

T %>%  #the weird position of the points comes from the fact that we're selecting only every 10 points to be plot.
  ggplot()+
  geom_point(aes(longitude, latitude, col = max_exp))+
  coord_map()+
  theme_minimal()+
  geom_sf(data = switzerland, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks= -c(10, 8, 6))+
  theme_minimal(12)+
  labs(col = 'exp', x = '', y = '')



#SAME THING BUT FOR THE PLAIN 
gc()
rm(list = ls())

library(tidyverse)
library(data.table)
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")

I_plain = read.csv("Data/plain_id_lon_lat_correspondance.csv")

files = list.files(path = "Data/Climate_data/By_id", full.names = TRUE)
clim_data_extreme_9_list = list()

for (i in seq_along(files)){
  
  clim_data = fread(files[[i]])%>%
    filter(id %in% I_plain$id)
  
  #keeping only temperatures above the associated threshold 0.9 quantile
  sing_clim_data_extreme_9 = clim_data %>%
    group_by(id) %>%
    mutate(threshold = quantile(maxtp, 0.9),
           excess = maxtp - threshold) %>%
    filter(excess > 0) %>%
    ungroup()
  
  clim_data_extreme_9_list[[i]] = sing_clim_data_extreme_9
  
  #to free memory
  rm(clim_data)
  gc()
}

#date_id, maxtp, id, threshold, excess
clim_data_extreme_9 =  do.call(rbind, clim_data_extreme_9_list)

rm(clim_data_extreme_9_list)
rm(sing_clim_data_extreme_9)

#to have the scale parameter associated
scales_9_df = read.csv("Data/Climate_data/plain_clim_scale_grid_gpd_model.csv")

shape_9 = readRDS("Data/clim_data_gpd_model/plain_optimal_shape.rds")

clim_data_extreme_9 = clim_data_extreme_9 %>%
  left_join(scales_9_df, by="id")

clim_data_extreme_9$shape = shape_9

#reducing the number of data else computation time would explose 

clim_data_extreme_9_reduced = clim_data_extreme_9 %>%
  filter(id %% 50 == 0)

standardised_qq = clim_data_extreme_9_reduced %>%
  mutate(unif = evd::pgpd(q = (excess), loc = 0, scale = scale_9, shape = shape_9)) %>%
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
clim_data_extreme_9_reduced$ideal = (clim_data_extreme_9_reduced$maxtp %>% rank(ties.method = "random"))/(nrow(clim_data_extreme_9_reduced)+1) 
clim_data_extreme_9_reduced$exp_ideal = -log(1-clim_data_extreme_9_reduced$ideal)


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


#plotting Switzerland with the color depending on the value of exp 

library(sf)
library(rnaturalearth)

my_pal = c( 
  '#062c30', # extra dark 
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')


switzerland <- ne_countries(country = "Switzerland", scale = "medium", returnclass = "sf")

switzerland <- st_transform(switzerland, crs = 4326)

#T = standardised_qq %>%select(id, excess, exp)%>% inner_join(I_plain, by="id")

T = standardised_qq %>%select(id, excess, exp)%>%group_by(id)%>% mutate(max_exp= max(exp) )%>% ungroup %>% inner_join(I_plain, by="id")

T %>%  #the weird position of the points comes from the fact that we're selecting only every 10 points to be plot.
  ggplot()+
  geom_point(aes(longitude, latitude, col = max_exp))+
  coord_map()+
  theme_minimal()+
  geom_sf(data = switzerland, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks= -c(10, 8, 6))+
  theme_minimal(12)+
  labs(col = 'exp', x = '', y = '')

