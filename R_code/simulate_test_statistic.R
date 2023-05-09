## required libraries
library(data.table)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)

## requirer helper functions
source('R_Code/helper_functions/functions_test_statistic.R')
source('R_code/helper_functions/functions_lambda.R')

## simulation settings
n_sim <- 10000
n_max <- 200

## set seed
set.seed(11052022)

## calculate lambda table
lambda_values <- (0:12)/20
pC_values <- c(0, 0.02, (1:10)/20)
lambda_table <- NULL
for(i_pC in 1:length(pC_values)){
  lambda_table %<>% cbind(calculate_treatment(p_lambda = lambda_values, p_pC = pC_values[i_pC]))
}
lambda_table %<>% as.data.frame()
names(lambda_table) <- pC_values
row.names(lambda_table) <- lambda_values
lambda_table %<>% round(4)

## simulation of endpoint in control group
for(i_prop in 1:length(pC_values)){
  
  ## simulate binary tables for intervention and control group
  bin_table_C <- simulate_bin_table(pC_values[i_prop], n_max, n_sim)
  
  ## save binary tables
  fwrite(bin_table_C, paste0('bin_tables/con', pC_values[i_prop], '.csv'))
  
  ## print progress
  print(paste0('Control group ', i_prop, ' simuluted.'))
}

## simulation of endpoint in intervention group
pI_table <- lambda_table
for(i_col in 1:ncol(pI_table)){
  pI_table[,i_col] <- pI_table[,i_col] + pC_values[i_col]
}
pI_values <- as.vector(as.matrix(pI_table)) %>% unique
for(i_prop in 1:length(pI_values)){
  
  ## simulate binary tables for intervention and control group
  bin_table_I <- simulate_bin_table(pI_values[i_prop], n_max, n_sim)
  
  ## save binary tables
  fwrite(bin_table_I, paste0('bin_tables/int', pI_values[i_prop], '.csv'))
  
  ## print progress
  print(paste0('intervention group ', i_prop, ' simuluted.'))
}

for(i_C in 1:ncol(pI_table)){
  
  for(i_I in 1:nrow(pI_table)){
    
    prop_I <- pI_table[i_I,i_C]
    prop_C <- pC_values[i_C]
    table_I <- fread(paste0('bin_tables/int', prop_I, '.csv'))
    table_C <- fread(paste0('bin_tables/con', prop_C, '.csv'))
    Z_values <- NULL
    for(i_samp in 1:n_max){
      Z_values %<>% cbind(normapp_from_bin_table(table_I , table_C, i_samp))
    }
    Z_values %<>% as.data.table
    names(Z_values) <- paste0('sample_size_', 1:n_max)
    fwrite(Z_values, paste0('approx_test_tables/int', prop_I, '_con', prop_C, '.csv'))
    
    ## print progress
    print(paste0('Approximation test table ', (i_C-1)*nrow(pI_table) + i_I, ' simuluted.'))
  }
}

## save lambda table
saveRDS(lambda_table, 'bin_tables/lambda_table.RDS')



