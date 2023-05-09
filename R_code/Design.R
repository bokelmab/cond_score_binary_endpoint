# Storing design parameters and caluclation of local signficance level according
# to Pocock

library(mvtnorm)

design <- function(
  n_1, # interim sample size per group
  alpha_glob, # Global significance level (one-sided)
  n_2, # initial incremental sample size for stage 2 per group
  alpha_0, # futility stop bound
  f #factor for maximal sample size
){
  
  #----------------------------------------------------------------------------#
  # specify correlation between interim and final test statistic 
  #----------------------------------------------------------------------------#
  w_1 <- sqrt(n_1)
  w_2 <- sqrt(n_2)
  r <- matrix(c(1, w_1 / sqrt(w_1^2 + w_2^2), w_1 / sqrt(w_1^2 + w_2^2), 1), ncol = 2)
  
  #--------------------------------------------------------------#
  # Calculate local significance levels according to Pocock      #
  #--------------------------------------------------------------#
  glob <- 1
  alpha_loc <- alpha_glob / 2
  while(glob > (1 - alpha_glob)){
    alpha_loc <- alpha_loc + 0.00001
    glob <- pmvnorm(lower = c(-Inf, -Inf), 
                  upper = c(qnorm(1 - alpha_loc), qnorm(1 - alpha_loc)), 
                  mean = c(0, 0), sigma = r)
    glob <- glob[1]
  }
  alpha_loc <- alpha_loc - 0.00001  
  
  res <- list(
    n_1 = n_1,
    alpha_1 = alpha_loc,
    n_2 = n_2,
    alpha_0 = alpha_0,
    f = f,
    alpha_glob = alpha_glob
  )
  class(res) <- c("TwoStageDesign", class(res))
  return(res) 
}