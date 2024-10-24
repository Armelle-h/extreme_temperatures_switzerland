#testing if bulk_models_function_cv works well

setwd("C:/Users/HOURS/Desktop/PDM/extreme_temperatures_switzerland")
source('bulk_models_function_cv.R')

library(tidyverse)



process_id = function(i, optimal_shape){
  this_clim_extm_irel_9 <- id_clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)
  # Estimate scale parameter with the optimal shape parameter
  model_fit_9 <- estimate_scale_fixed_shape(this_clim_extm_irel_9, optimal_shape)
  return(list(par = model_fit_9$par, loglik= model_fit_9$value))
}

job_process_id = function (indices, optimal_shape){
  scales = c()
  loglik_sum = 0
  for (i in indices){
    fun_output = process_id(i, optimal_shape)
    scales = c(scales, fun_output$par)
    loglik_sum = loglik_sum + fun_output$value
  }

  return (list(scales_=scales, loglik_sum_=loglik_sum))
}


result_1 = job_process_id(c(10, 11), optimal_shape)

result_2 = job_process_id(c(12, 13), optimal_shape)

result_total = c(result_1$scales_, result_2$scales_)