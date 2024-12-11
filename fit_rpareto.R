#Model fitting using a likelihood-based approach


#Carefull !! Sometimes the covariance matrix is singular because I was considering stations too close to each other. 
#Should I do some additional preprocessing?

#Should also differenciate in the file name if we fit with robust or not !!!!

rm(list=ls())
setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
library(tidyverse)
library(doParallel)
library(optimx)
library(foreach)
library(sf)

marg_mod = "mod_1"

#defining the likelihood function based on pairwise relationships between locations 
#computes the loglikelihood for a single dataset/realization of the data
#dt is the dataset, lcs the location coordinates, vr the function of rmodeling pairwise relationships (variogram)
#conditioned.site: the index of the site used for conditioning

HRD_ll = function(dt, lcs, vr, conditioned.site = 1, zeta, gamma, eta){
  
  #pairwise location calculations
  nlocs = nrow(lcs) # number of sites
  loc.id.pairs = expand.grid(seq(nlocs) ,seq(nlocs))
  loc.pairs = cbind(lcs[loc.id.pairs[,1],], lcs[loc.id.pairs[,2],]) # all pairs of locations
  
  dstncs = loc.pairs %>%
    apply(MARGIN = 1, FUN = function(x) sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2))
  
  #variogram matrix
  Lambda = (vr(dstncs)) %>% matrix(nrow = nlocs, byrow = T)
  conditioned.site.lambda = Lambda[conditioned.site,]
  Lambda = Lambda[-conditioned.site, -conditioned.site]
  
  psi = outer(conditioned.site.lambda[-1], conditioned.site.lambda[-1], "+") - Lambda # "covar" matrix
  detPsi = det(psi)
  
  if (detPsi == 0) {
    warning("Covariance matrix is singular, adding regularization.")
    epsilon <- 1e-6  # Small regularization parameter
    psi <- psi + epsilon * diag(nrow(psi))  # Regularize matrix
  }
  
  inversePsi = psi %>% solve()  # inverse of "covariance" matrix
  
  detPsi = determinant(psi)$modulus[1] # calculates log det by default
  
  #loglikelihood contribution
  omegas = log(t(sapply(dt, function(x) x[-conditioned.site]))/sapply(dt, "[[", conditioned.site)) %>% sweep(2, conditioned.site.lambda[-1], "+")
  
  if(nlocs == 2){
    omegas = t(omegas)
  }
  
  summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
  length(dt)*detPsi+ summ # penalty term + sum of the LL
}

# negative likelihood function for model fitting for the entire dataset
ngll = function(pars){
  
  print(pars)
  if( pars[2] < 0) return (2^100)
  if( pars[1] < 0) return (2^100)
  
  LL <- foreach(i= 1:length(exceedances), .combine = 'c') %dopar% {   #for i == 223, issue with a singular matrix
    my.vario <- function(h){
      if(stat){
        alpha = pars[1]
        beta = pars[2] 
      }else{
        alpha = pars[1] + pars[3]*temporal_covar[[1]][i]
        beta = pars[2] 
      }
      
      if(alpha <= 0){
        return (2^100)
      } 
      
      if( beta <= 0){
        return (2^100)
      } 
      
      nu = 0.01 #used to be 0.2
      res = rep(NA, length(h))
      res[h == 0] = 0
      res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
      res
    }
    HRD_ll(lcs = exceedances_locs[i][[1]], dt = exceedances[i], vr = my.vario, zeta = 1, gamma = 1, eta = 1)
  }
  LL = sum(LL)
  print(LL)
  LL
}


data_for_rpareto = readRDS(paste0("Data/processed/data_for_rpareto/true/robust_plain_data_for_rpareto_", marg_mod))
.GlobalEnv$exceedances <- data_for_rpareto$exceedances
.GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs #this is a list of location id (just the id, not the longitude-latitude pair)
.GlobalEnv$HRD_ll <- HRD_ll

#defining the temporal covariates, not great to use the year as covariate, for now dropping, might come back to it later

#.GlobalEnv$temporal_covar = data_for_rpareto$extreme_dates%>%#exceedances and exceedances_locs are also aranged by date so all good :)
#  mutate(year = lubridate::year(date))%>%
#  mutate(year_index = year-1970)%>%
#  select(year)

.GlobalEnv$stat <- TRUE #do we want alpha and beta to vary with time or not ? 

fit <<- optimx::optimx(par = c(0.7, 600), fn = ngll ,hessian = F, method = 'Nelder-Mead') #used to have par = c(0.5,30)

# method = 'L-BFGS-B'
tibble(fit$p1, fit$p2) %>%
  write_csv(paste0("output/rpareto_model_fits/robust_plain_true_rpareto_fits_model_", marg_mod, ".csv"), append = T)


tibble(0.4644155, 45.313428) %>%
  write_csv(paste0("output/rpareto_model_fits/robust_plain_nu_015_true_rpareto_fits_model_", marg_mod, ".csv"), append = T)














num_cores <- detectCores() - 1
cl <- parallel::makeCluster(num_cores, type = "PSOCK")
clusterExport(cl, list("HRD_ll", "exceedances_locs"))
clusterEvalQ(cl, {library(tidyverse); library(sf)})
doParallel::registerDoParallel(cl)
fit <<- optimx::optimx(par = c(1,1), fn = ngll ,hessian = F, method = 'Nelder-Mead') # Variance + range constant
parallel::stopCluster(cl)

#debug 
num_cores <- detectCores() - 1
cl <- parallel::makeCluster(num_cores, type = "PSOCK")
clusterExport(cl, list("HRD_ll", "exceedances_locs"))
clusterEvalQ(cl, {library(tidyverse); library(sf)})
doParallel::registerDoParallel(cl)
ngll(c(1,1))
parallel::stopCluster(cl)