## function to simulate n_max x n_sim observations of a binary variable with specified proportion
simulate_bin_table <- function(p_pi1, p_nmax, p_nsim){
  
  bin_table <- NULL
  for(i_n in 1:p_nmax){
    bin_table %<>% cbind(rbinom(n = p_nsim,size = 1, prob = p_pi1))
  }
  return(as.data.table(bin_table))
}

## calculate normal approximation test statistic from two empirical proportions
norm_app <- function(p_meanI, p_meanC, p_n){
  
  return(sqrt(p_n/2)*(p_meanI-p_meanC)/sqrt((p_meanI+p_meanC)/2*(1-(p_meanI+p_meanC)/2)))
}

## apply normal approximation test to tables of binary variables
normapp_from_bin_table <- function(p_table_I, p_table_C, p_n){
  
  ## reduce to chosen sample size per group
  p_table_I <- p_table_I[,1:p_n]
  p_table_C <- p_table_C[,1:p_n]
  
  ## calculate group averages
  av_I <- apply(p_table_I, 1, mean)
  av_C <- apply(p_table_C, 1, mean)
  
  ## apply normal approximation test
  combined_av <- cbind(av_I, av_C)
  res <- mapply(norm_app, combined_av[,1], combined_av[,2], MoreArgs = list(p_n = p_n))
  
  ## return test statistic results
  return(res)
}

