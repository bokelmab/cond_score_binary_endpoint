## function to calculate fixed design sample size
calculate_Nfix <- function(p_lambda, p_pow_targ, p_alpha_glob){
  return((sqrt(2)*qnorm(1-p_alpha_glob)/p_lambda + qnorm(p_pow_targ)*sqrt(2/p_lambda^2 - 1/2))^2)
}

## function for fixed design power
calculate_fix_power <- function(p_lambda, p_Nfix, p_alpha_glob){
  return(1- pnorm((qnorm(1-p_alpha_glob)-p_lambda*sqrt(p_Nfix/2)))/sqrt(1-1/4*p_lambda^2))
}