## required libraries
source('R_Code/helper_functions/functions_power_sample_size.R')

calculate_conditional_performance_score <- function(p_meanN, p_varN, p_meanCP, p_varCP, p_Nfix, p_n1, p_Nmax, p_alpha, p_beta, p_lambda){
  
  ## calculate target values
  N_target <- ifelse(p_Nfix > p_Nmax, p_n1, p_Nfix)
  CP_target <- ifelse(p_Nfix > p_Nmax | p_lambda == 0, p_alpha, 1-p_beta)
  
  # location component n:
  e_n <- 1 - abs(p_meanN - N_target) / (p_Nmax - p_n1)
  
  # variation component n:
  var_n_max <- ((p_Nmax - p_n1)^2)/4
  v_n <- 1 - sqrt(p_varN / var_n_max)
  # total n subscore:
  score_n <- 0.5*e_n + 0.5*v_n
  # location component cp:
  e_cp <- 1 - abs(p_meanCP - CP_target) / (1 - p_alpha)
  # variation component cp:
  var_cp_max <- 0.25
  v_cp <- 1 - sqrt(p_varCP / var_cp_max)
  # total cp subscore:
  score_cp <- 0.5*e_cp + 0.5*v_cp
  
  return(list(
    e_cp = e_cp, # location component conditional pwer
    v_cp = v_cp, # variation component conditional power
    score_cp = score_cp, # sub-score conditional power
    e_n = e_n, # location component conditional sample size
    v_n = v_n, # variation component conditional sample size
    score_n = score_n, # sub-score conditional sample size
    score_cond = 0.5*score_n + 0.5*score_cp, # conditional score
    cp_target = CP_target
  ))
}

