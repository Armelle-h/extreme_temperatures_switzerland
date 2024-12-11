# simulate under true model
# ---- similate on clim scale grid

rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(parallel)
source('mvPot/simulPareto.R') #mvPot is no compatible with my current version of R 

# --- marginal model names
marg_mod = 'mod_1'

#nu

nu_val = 0.15
nu_name = "015"

read.csv("Data/Observed_data/plain_1971_2022_JJA_obs_legend.csv") %>% 
  dplyr::select(longitude_proj, latitude_proj) %>% 
  unique %>% 
  write_csv("Data/processed/plain_obs_grid_simulated_on.csv")

locs_to_pred = as.data.frame(read_csv("Data/processed/plain_obs_grid_simulated_on.csv"))

fit = read_csv(paste0("output/rpareto_model_fits/robust_plain_nu_", nu_name,"_true_rpareto_fits_model_",marg_mod, ".csv"),
               col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()


alpha_var <<- fit[1]
beta_var <<- fit[2] 

variogram_model <- function(h){
  nu = nu_val 
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}

nlocs = nrow(locs_to_pred) # number of sites
loc.id.pairs = expand.grid(seq(nlocs) ,seq(nlocs))
loc.pairs = cbind(locs_to_pred[loc.id.pairs[,1],], locs_to_pred[loc.id.pairs[,2],]) # all pairs of locations

dstncs = loc.pairs %>%
  apply(MARGIN = 1, FUN = function(x) sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2))

#variogram matrix
vario_matrix = (variogram_model(dstncs)) %>% matrix(nrow = nlocs, byrow = T)

for (index in c(1, 26, 51, 76)) {
  job::job (
    {
      file_name = paste0("output/simulations/simulations_on_obs_grid/true/nu_", nu_name,"_",marg_mod,"_run_")
      
      for(i in seq(index, index+24)){
        print(i)
      
        # this returns samples in the unit frechet margin
        simulation_score <- simulPareto(n = 2500,
                                      vario_mat = vario_matrix,
                                      nCores = nCores,
                                      cl = cl,
                                      robust = TRUE)
      
        simulation_score %>%
          saveRDS(paste0(file_name,i))
      }
    
    }
  )
}


#old version


#again it doesn't like the cluster, better to do the parallelization based on jobs.

simulate = function(file_name){
  for(i in seq(1, 100)){
    print(i)
    
    # this returns samples in the unit frechet margin
    simulation_score <- simulPareto(n = 2500,
                                    vario_mat = vario_matrix,
                                    nCores = nCores,
                                    cl = cl,
                                    robust = TRUE)
    
    simulation_score %>%
      saveRDS(paste0(file_name,i))
  }
}


#simulate(paste0("output/simulations/simulations_on_obs_grid/true/",marg_mod,"_run_"))