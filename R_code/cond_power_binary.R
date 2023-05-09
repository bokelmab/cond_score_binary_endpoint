# conditional power calculation

cond_power_binary <- function(t1, # Z_1 value first stage
                       lambda, # true or estimated parameter
                       n_1, # interim sample size per group
                       n_2, # initial incremental sample size for stage 2 per group
                       n_2new, # recalculated sample size for stage 2 per group
                       alpha_loc # local significance level
){
  
  w_1 <- sqrt(n_1)
  w_2 <- sqrt(n_2)
  cp <- 1 - pnorm((sqrt(w_1^2+w_2^2)/w_2*qnorm(1-alpha_loc) - w_1/w_2*t1 - lambda*sqrt(n_2new/2))/sqrt(1-1/4*lambda^2))
  
  
  ## special case: observed lambda=0
  cp[lambda==0] <- 0
  
  return(cp)
}