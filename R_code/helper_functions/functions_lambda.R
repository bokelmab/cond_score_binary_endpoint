## helper function to calculate lambda
calculate_lambda <- function(p_pC, p_pI){
  return((p_pI - p_pC)/ sqrt((p_pI+p_pC)/2 * (1-(p_pI+p_pC)/2)))
}

## helper function to calculate treatment effect from pC and lambda
calculate_treatment <- function(p_lambda, p_pC){
  
  sum1 <- -p_pC*p_lambda^2 + p_lambda^2 + 4 * p_pC
  sum2 <- p_lambda^4 - 16 * (p_lambda^2) * (p_pC^2) + 16 * (p_lambda^2) * p_pC
  sum2 <- sqrt(sum2)
  denom <- p_lambda^2 + 4
  return((sum1 + sum2)/denom - p_pC)
}
