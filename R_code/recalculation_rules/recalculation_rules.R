## this function calculates values and functions which are relevant for all or most recalculation rules
## each recalculation rule applies this function at the beginning
rec_base <- function(design, t1){
  
  ## local variables
  base_env <- environment()
  lambda_obs <- t1 * sqrt(2 / design$n_1)   
  w_1 <- sqrt(design$n_1)
  w_2 <- sqrt(design$n_2)
  n_ocp <- design$n_1*((sqrt(w_1^2+w_2^2)/w_2*qnorm(1 - design$alpha_1)-qnorm(0.2)*sqrt(1-1/2/design$n_1*t1^2))/t1-w_1/w_2)^2
  n_ocp <- ceiling(n_ocp)
  effic <- (t1 >= qnorm(1 - design$alpha_1))
  fut <- (t1 < qnorm(1 - design$alpha_0))
  Nrec_max <- design$n_1 * design$f - design$n_1
  
  ## observed conditional power function for a chosen second stage sample size
  obs_cp_func <- function(p_Nrec){
    
    return(cond_power_binary(t1 = t1, lambda = lambda_obs, n_1 = design$n_1, n_2 = design$n_2, p_Nrec, design$alpha_1))
  }
  
  return(list(lambda_obs = lambda_obs, n_ocp = n_ocp, effic = effic, fut = fut, Nrec_max = Nrec_max, obs_cp_func = obs_cp_func))
  
}

rec_OCP <-function(design, t1){
  
  base_res <- rec_base(design, t1)
  
  ## limit to n_max
  n_recalc <- base_res %$% ifelse(n_ocp > Nrec_max, Nrec_max, n_ocp)
  
  ## if no stop at interim, return Nrec
  n_recalc[base_res$effic | base_res$fut] <- 0
  return(n_recalc)
  
}

rec_restrOCP <- function(design, t1){
  
  base_res <- rec_base(design, t1)
  
  ## observed conditional power for maximum recalculated sample size
  obs_cp <- base_res$obs_cp_func(base_res$Nrec_max)
  
  ## restriction
  n_recalc <- base_res %$% ifelse(obs_cp <= 0.6, 0, ifelse(n_ocp > Nrec_max, Nrec_max, n_ocp))
  
  ## if no stop at interim, return Nrec
  n_recalc[base_res$effic | base_res$fut] <- 0
  return(n_recalc)
}

rec_Promising <- function(design, t1){
  
  base_res <- rec_base(design, t1)
  
  ## cp with originally planned second stage
  obs_cp_ini <- base_res$obs_cp_func(design$n_2)
  
  ## recalculation
  n_recalc <- base_res %$% ifelse(obs_cp_ini < 0.36, design$n_2,
                                  ifelse(obs_cp_ini >= 0.8, design$n_2,
                                         ifelse(n_ocp > Nrec_max, Nrec_max, n_ocp)))
  
  ## if no stop at interim, return Nrec
  n_recalc[base_res$effic | base_res$fut] <- 0
  return(n_recalc)
}

rec_OptFunc <- function(design, t1){

  base_res <- rec_base(design, t1)
  gamma <- 0.005/4 ## for comparison binary/normaly distr

  table_cp <- NULL
  for(i_n in design$n_2:base_res$Nrec_max){
    table_cp %<>% cbind(base_res$obs_cp_func(i_n))
  }

  table_n <- NULL
  for(i_n in design$n_2:base_res$Nrec_max){
    table_n %<>% cbind(rep(i_n - design$n_2, nrow(table_cp)))
  }

  table_of <- table_cp - gamma * table_n

  ## the index i of the optimum of the function corresponds to Nrec = i + design$n_2 - 1
  n_recalc <- apply(table_of, 1, which.max) + design$n_2 - 1

  ## if no stop at interim, return Nrec
  n_recalc[base_res$effic | base_res$fut] <- 0
  return(n_recalc)

}

rec_classicGS <- function(design, t1){
  
  base_res <- rec_base(design, t1)
  
  n_recalc <- rep(design$n_1, length(t1))
  ## if no stop at interim, return Nrec
  n_recalc[base_res$effic | base_res$fut] <- 0
  return(n_recalc)
}








