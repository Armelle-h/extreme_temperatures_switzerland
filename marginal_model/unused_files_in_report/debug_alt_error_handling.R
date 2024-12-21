library(evgam)

fit_with_warning_handling <- function(obs_data_for_quant_reg, zeta, inits_coeff) {
  # Define a function to fit the model and capture warnings
  fit_model <- function(start_params) {
    # Capture warnings
    fit <- tryCatch({
      withCallingHandlers(
        evgam(maxtp ~ value + glob_anom, obs_data_for_quant_reg,
              family = "ald", inits = inits_coeff, ald.args = list(tau = zeta)),
        warning = function(w) {
          assign("last_warning", conditionMessage(w), envir = .GlobalEnv)
          invokeRestart("muffleWarning")
        }
      )
    }, error = function(e) e)
    
    return(fit)
  }
  
  # Fit the model initially
  assign("last_warning", NULL, envir = .GlobalEnv)
  model <- fit_model(initial_start)
  
  # Check if the warning matches the specific warning
  if (!is.null(last_warning) && grepl("Final Hessian of negative penalized log-likelihood not numerically positive definite", last_warning)) {
    message("Warning encountered: ", last_warning)
    message("Retrying with modified initial parameters...")
    
    # Modify initial parameters (example: random initialization)
    new_start <- runif(3, min = -1, max = 1)
    
    model <- fit_model(new_start)
  }
  
  return(model)
}

q = 1
  obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
  zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
  obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist

  if(q==1){
    inits_coeff = c(-1, 1, 0.5)
  }else {
    #to ensure smoothness in the estimated coefficients
    inits_coeff = c(quantile_model_fit$location$coefficients[1],quantile_model_fit$location$coefficients[2])
  }
  assign("last_warning", NULL, envir = .GlobalEnv)
  quantile_model_fit <- fit_with_warning_handling(obs_data_for_quant_reg, zeta, c(-1,1,0.5))
  
  
  
  
  