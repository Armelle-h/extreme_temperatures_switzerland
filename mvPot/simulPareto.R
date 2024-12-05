
library(evd)


simulPareto <- function(n, vario_mat, nCores = 1, cl = NULL, robust = FALSE){
  
  gamma = vario_mat
  
  dim <- nrow(gamma) #number of locations considered
  
  k <- sample(1:dim, n, replace = TRUE) #selects randomly with replacements n integer elements in the interval [1, dim]. 
  d <- nrow(gamma)
  
  #it uses location 1 as reference
  cov.mat <- (outer(gamma[-1,1],gamma[-1,1], "+") - (gamma[-1,-1])) #defining the covariance matrix
  chol.mat <- chol(cov.mat) #cholesky decomposition of the covariance matrix (which we will use to generate the log gaussian stochastic process)
  
  #defines a matrix where each column is a realization of the simulation of Wr(s), the first line is 0 as ...
  proc <- matrix(0,d,n) 
  proc[-1,] <- t(chol.mat)%*%matrix(rnorm((d- 1)*n), ncol=n) #matrix(rnorm((d-1)*n), ncol=n) is a standard normal random matrix.
  #it defines a gaussian process with coavriance matrix cov
  
  sims <- lapply(1:n, function(i){
    buffer <- exp(proc[,i] - proc[k[i],i] - gamma[,k[i]]) #log-gaussian process 
    if (robust == TRUE){
      proc[,i] <-evd::rgpd(1, loc=1, scale=1, shape=1) * buffer / median(buffer)
    } else{
      proc[,i] <-evd::rgpd(1, loc=1, scale=1, shape=1) * buffer / mean(buffer) #from the decomposition Y_r = R W_r ,
      #R unit pareto, W_r is buffer normalized such that it has risk 1, the risk function is defined as the mean.
    }
    
    proc[proc == 0] = .Machine$double.xmin
    proc[,i]})
  return(sims)
}