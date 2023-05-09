## required libraries
library(data.table)
library(plyr)
library(dplyr)
library(magrittr)

## required helper functions
source('R_Code/recalculation_rules/recalculation_rules.R')
source("R_Code/Design.r")
source("R_Code/cond_power_binary.R")
source('R_code/helper_functions/functions_lambda.R')

## read lambda table and specify n1 values
lambda_table <- readRDS('bin_tables/lambda_table.RDS')
n1_values <- c(10, 20, 50)

## table to save results for the chosen recalcucation rule
results_rec_rule <- NULL

## simulation fixed settings
n_sim <- 10000
beta <- 0.2
alpha_glob <- 0.025
design <- design(n_1 = 50, alpha_glob = alpha_glob, n_2 = 50, alpha_0 = 0.5, f = 4)
w_1 <- sqrt(design$n_1)
w_2 <- sqrt(design$n_2)
list_rec_rules <- c('rec_OCP', 'rec_restrOCP', 'rec_Promising', 'rec_OptFunc', 'rec_classicGS')

## for loop over recalculation rules
for(i_rule in 1:length(list_rec_rules)){
  
  rec_rule <- list_rec_rules[i_rule]
  
  ## for loop over n1
  for(i_n1 in 1:length(n1_values)){
    n_1 <- n1_values[i_n1]
    
    ## for loop over combinations of pI and pC
    for(i_lambda in 1:nrow(lambda_table)){
      
      for(i_pC in 1:ncol(lambda_table)){
        
        pC <- as.numeric(names(lambda_table)[i_pC])
        pI <- lambda_table[i_lambda,i_pC] + pC
        lambda <- as.numeric(row.names(lambda_table))[i_lambda]
        
        ## helper function to decide about rejection
        reject <- function(p_Z, p_crit){
          return(p_Z >= p_crit)
        }
        
        ## helper function to impute NA's of test statistic by zero
        impute_normapprox <- function(p_Zvalues){
          p_Zvalues[is.na(p_Zvalues)] <- 0
          return(p_Zvalues)
        }
        
        ## distribution of Z_1
        distr_Z1 <- fread(paste0('approx_test_tables/int', pI, '_con', pC, '.csv'), select = paste0('sample_size_', n_1))[[1]]
        distr_Z1 <- impute_normapprox(distr_Z1)
        rejection_Z1 <- sapply(distr_Z1, reject, p_crit = qnorm(1-design$alpha_1))
        
        ## distribution of N_rec
        distr_Nrec <- do.call(what = rec_rule, args = list(design = design, t1 = distr_Z1))
        sim_table <- data.table(Z1 = distr_Z1, lambda_obs = distr_Z1*sqrt(2/n_1), lambda_real = lambda, n_1 = n_1, Nrec = distr_Nrec)
        
        ## rejection information
        sim_table[, q_alpha0 := qnorm(1-design$alpha_0)]
        sim_table[, q_alpha1 := qnorm(1-design$alpha_1)]
        sim_table[, q_alpha12 := qnorm(1-design$alpha_1)]
        sim_table[, crit2 := (qnorm(1-design$alpha_1) * sqrt(w_1^2+w_2^2) - w_1 * Z1)/w_2]
        sim_table[, effic := Z1 >= qnorm(1-design$alpha_1)]
        sim_table[, fut := Z1 < q_alpha0]
        sim_table[, continue := !(effic | fut)]
        sim_table[, asymp_cp := cond_power_binary(Z1, lambda, n_1, n_1, n_2new = Nrec, alpha_loc = 0.01469)]
        sim_table[, obs_cp := cond_power_binary(Z1, lambda_obs, n_1, n_1, n_2new = Nrec, alpha_loc = 0.01469)]
        sim_table[, true_cp := NA]
        
        ## second stage test statistic distribution
        distr_Z2 <- fread(paste0('approx_test_tables/int', pI, '_con', pC, '.csv'))
        
        ## calculate conditional power for each simulation
        for(i_sim in 1:nrow(sim_table)){
          
          ## special case: futility for restrOCP
          if(sim_table$Nrec[i_sim]==0){
            distr_Z2_Nrec <- rep(0, n_sim)
          }else{
            distr_Z2_Nrec <- distr_Z2[, paste0('sample_size_', sim_table$Nrec[i_sim]), with = F][[1]]
            distr_Z2_Nrec <- impute_normapprox(distr_Z2_Nrec)
          }
          true_cp <- mean(distr_Z2_Nrec >= sim_table$crit2[i_sim])
          sim_table$true_cp[i_sim] <- true_cp
        }
        
        ## obtain test statistics
        options(scipen=999)
        calculate_summary_statistics <- function(p_sim_table){
          
          global_stats <- p_sim_table %$% c(efficacy = mean(effic), futility = mean(fut), power = mean(effic) + mean(continue*true_cp))
          cond_stats <- p_sim_table[continue ==T,] %$% c(mean_true_cp = mean(true_cp), var_true_cp = var(true_cp), mean_obs_cp = mean(obs_cp), var_obs_cp = var(obs_cp), mean_n = mean(n_1 + Nrec), var_n = var(Nrec))
          target_values <- c(lambda = lambda, n_1 = n_1, pI = pI, pC = pC, pdiff = pI-pC, n_fix = (sqrt(2)*qnorm(0.975)/lambda + qnorm(0.8)*sqrt(2/lambda^2 - 1/2))^2, beta = beta, alpha_glob = alpha_glob, Nmax = design$f*design$n_1)
          return(c(target_values, global_stats, cond_stats))
        }
        
        results_rec_rule %<>% rbind(round(calculate_summary_statistics(sim_table),4))
      }
      print(paste0('lambda ', lambda, ' simulated'))
    }
    
    ## save results
    fwrite(as.data.table(results_rec_rule), paste0('results/score_results/', rec_rule, '_n1_', n_1, '.csv'))
    results_rec_rule <- NULL
  }
  
}







