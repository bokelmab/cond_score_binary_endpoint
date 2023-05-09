## required libraries
library(data.table)
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)

## required helper functions
source('R_code/helper_functions/functions_score.R')

list_rec_rules <- c('rec_OCP', 'rec_restrOCP', 'rec_Promising', 'rec_OptFunc', 'rec_classicGS')
results <- fread(paste0('results/score_results/rec_OCP_n1_50.csv'))
pC_values <- results$pC %>% unique 
n1_values <- c(10, 20, 50)
n_1 <- 50

## Calculation of performance measures for n1=50, pC=0.3
pC_sel <- 0.3
power_results <- NULL
sample_size_results <- NULL
score_results <- NULL
NPow_tradeoff_results <- NULL
lambda_as <- 0.3 ## assumed treatment effect for fixed design
for(i_rule in 1:length(list_rec_rules)){
  
  rec_rule <- list_rec_rules[i_rule]
  results <- fread(paste0('results/score_results/', rec_rule, '_n1_', n_1, '.csv'))
  
  ## power 
  power_results %<>% cbind(results[pC == pC_sel,]$power)
  
  ## sample size plot
  unc_sample_size <- results[pC == pC_sel,]$mean_n*(1-results[pC == pC_sel,]$futility-results[pC == pC_sel,]$efficacy) + 
    n_1*(results[pC == pC_sel,]$futility+results[pC == pC_sel,]$efficacy)
  sample_size_results %<>% cbind(unc_sample_size)
  
  ## conditional score results
  cond_results <- results[pC == pC_sel,] %$% calculate_conditional_performance_score(p_meanN = mean_n, p_varN = var_n,
                                                          p_meanCP = mean_obs_cp, p_varCP = var_obs_cp,
                                                          p_Nfix = n_fix, p_n1 = n_1, p_Nmax = Nmax,
                                                          p_alpha = alpha_glob, p_beta = beta, p_lambda = lambda)
  score_results %<>% cbind(cond_results$score_cond)
  
}

## save power results
power_results %<>% as.data.frame
names(power_results) <- str_remove(list_rec_rules, 'rec_')
row.names(power_results) <- results[pC == pC_sel,]$lambda
saveRDS(power_results, paste0('results/for_paper/power_results_n1_', n_1, '_pC_', pC_sel,'.RDS'))

## save sample size results
sample_size_results %<>% as.data.frame
names(sample_size_results) <- str_remove(list_rec_rules, 'rec_')
row.names(sample_size_results) <- results[pC == pC_sel,]$lambda
saveRDS(sample_size_results, paste0('results/for_paper/sample_size_results_n1_', n_1, '_pC_', pC_sel,'.RDS'))

## save score results
score_results %<>% as.data.frame
names(score_results) <- str_remove(list_rec_rules, 'rec_')
row.names(score_results) <- results[pC == pC_sel,]$lambda
saveRDS(score_results, paste0('results/for_paper/score_results_n1_', n_1, '_pC_', pC_sel,'.RDS'))

## statistics fixed design
fix_power <- c()
for(i_lambda in 1:nrow(results[pC == pC_sel,])){
  Nfix_new <- calculate_Nfix(lambda_as, 0.8, 0.025)
  fix_pow_new <- calculate_fix_power(results[pC == pC_sel,]$lambda[i_lambda], Nfix_new, 0.025)
  fix_power <- c(fix_power, fix_pow_new)
}
stats_fixed <- data.frame(lambda = results[pC == pC_sel,]$lambda, fix_power = fix_power, Nfix = calculate_Nfix(lambda_as, 0.8, 0.025))

## Global Power
data_plot_power <- cbind(lambda = as.numeric(row.names(power_results)), power_results)
data_plot_power <-  pivot_longer(data = data_plot_power, cols = !lambda, names_to = 'rec_rule', values_to = 'power')
pl_power <- ggplot(data=data_plot_power, aes(x=lambda, y=power, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  geom_line(data=stats_fixed,aes(x=lambda,y=fix_power, color = "fixed"), linetype="dashed", size = 1.1)+
  ylim(0, 1)+
  labs(x = expression(lambda))+
  labs(y = 'power')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle("(Global) Power")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = "rec rules", values = c("fixed" = "red", "classicGS"="#F8766D", "OCP" = "#CD9600", 
                                                    "OptFunc" = "#00BE67", "Promising"= "#00A9FF", "restrOCP"="#C77CFF"))

## Global mean sample size
data_plot_sample_size <- cbind(lambda = as.numeric(row.names(sample_size_results)), sample_size_results)
data_plot_sample_size <-  pivot_longer(data = data_plot_sample_size, cols = !lambda, names_to = 'rec_rule', values_to = 'sample_size')
pl_sample <- ggplot(data=data_plot_sample_size, aes(x=lambda, y=sample_size, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  geom_line(data=stats_fixed,aes(x=lambda,y=Nfix, color = "fixed"), linetype="dashed", size = 1.1)+
  ylim(0, 200)+
  labs(x = expression(lambda))+
  labs(y = 'mean sample size')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle("")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = "rec rules", values = c("fixed" = "red", "classicGS"="#F8766D", "OCP" = "#CD9600", 
                                                    "OptFunc" = "#00BE67", "Promising"= "#00A9FF", "restrOCP"="#C77CFF"))
ggarrange(pl_power, pl_sample, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

## Conditional performance score for binary endpoints
data_plot_score <- cbind(lambda = as.numeric(row.names(score_results)), score_results)
data_plot_score <- data_plot_score[!data_plot_score$lambda %in% c(0.05, 0.15, 0.25, 0.45, 0.55),] ## plot at the same points as for normally distr. endpoints
data_plot_score <-  pivot_longer(data = data_plot_score, cols = !lambda, names_to = 'rec_rule', values_to = 'score')
pl_score_bin <- ggplot(data=data_plot_score, aes(x=lambda, y=score, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  ylim(0, 1)+
  labs(x = expression(lambda))+
  labs(y = expression("S("*lambda*")"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle("")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))

## all three performance measures for different values of pc
three_perf_meas <- ggarrange(pl_sample, pl_power, pl_score_bin, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
annotate_figure(three_perf_meas, top = text_grob(bquote(paste(p[c], " = ", .(pC_sel))), 
                                      color = "black", face = "bold", size = 16))

## Conditional Performance score for normally distributed endpoints (results from Hermann et al. (2020))
norm_res_0 <- data.frame(delta = 0, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.477, 0.615, 0.652, 0.463, 0.778))
norm_res_01 <- data.frame(delta = 0.1, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.431, 0.541, 0.593, 0.399, 0.743))
norm_res_02 <- data.frame(delta = 0.2, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.398, 0.480, 0.547, 0.351, 0.711))
norm_res_03 <- data.frame(delta = 0.3, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.624, 0.390, 0.527, 0.599, 0.611))
norm_res_035 <- data.frame(delta = 0.35, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.584, 0.493, 0.620, 0.553, 0.698))
norm_res_04 <- data.frame(delta = 0.4, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.555, 0.544, 0.619, 0.525, 0.758))
norm_res_05 <- data.frame(delta = 0.5, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.543, 0.527, 0.597, 0.518, 0.722))
norm_res_06 <- data.frame(delta = 0.6, rec_rule = c('OCP', "restrOCP", "Promising", "OptFunc", "classicGS"),
                         score = c(0.557, 0.534, 0.595, 0.533, 0.714))
data_plot_norm_score <- rbind(norm_res_0, norm_res_01) 
data_plot_norm_score %<>% rbind(norm_res_02)
data_plot_norm_score %<>% rbind(norm_res_03)
data_plot_norm_score %<>% rbind(norm_res_035)
data_plot_norm_score %<>% rbind(norm_res_04)
data_plot_norm_score %<>% rbind(norm_res_05)
data_plot_norm_score %<>% rbind(norm_res_06)
pl_score_norm <- ggplot(data=data_plot_norm_score, aes(x=delta, y=score, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  ylim(0, 1)+
  labs(x = expression(Delta))+
  labs(y = expression("S("*Delta*")"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  ggtitle("CPS for norm. distr. EP")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))
ggarrange(pl_score_norm, pl_score_bin, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

## calculate type 1 error rate
for(i_n1 in 1:length(n1_values)){
  
  n_1 <- n1_values[i_n1]
  type1_error <- NULL
  for(i_pC in 1:length(pC_values)){
    
    pC_sel <- pC_values[i_pC]
    type1_error_row <- c()
    for(i_rule in 1:length(list_rec_rules)){
      
      rec_rule <- list_rec_rules[i_rule]
      results <- fread(paste0('results/score_results/', rec_rule, '_n1_', n_1, '.csv'))
      
      ## type 1 error rate
      type1_error_row <- c(type1_error_row, results[pC == pC_sel & lambda == 0,]$power)
    }
    type1_error %<>% rbind(type1_error_row)
  }
  type1_error %<>% as.data.frame
  names(type1_error) <- str_remove(list_rec_rules, 'rec_')
  row.names(type1_error) <- pC_values
  saveRDS(type1_error, paste0('results/for_paper/type1_error_n1_', n_1,'.RDS'))
  
}

## generate graphics for type one error rate
for(i_n1 in 1:length(n1_values)){
  n_1 <- n1_values[i_n1]
  type1_error <- readRDS(paste0('results/for_paper/type1_error_n1_', n_1,'.RDS'))
  data_plot_type1 <- cbind(pC = as.numeric(row.names(type1_error)), type1_error)
  data_plot_type1 <-  pivot_longer(data = data_plot_type1, cols = !pC, names_to = 'rec_rule', values_to = 'type1_error')
  new_plot <- (ggplot(data=data_plot_type1, aes(x=pC, y=type1_error, group=rec_rule)) +
    geom_line(aes(color=rec_rule), size = 1.1)+
    geom_hline(yintercept = 0.025, linetype="dashed", color = "red", size = 1.1)+
    ylim(0, 0.05)+
    labs(x = expression(p[c]))+
    labs(y = ifelse(i_n1 == 1,'type 1 error rate', ''))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16,face="bold"))+
    ggtitle(paste0("n1 = ", n_1))+
    theme(plot.title = element_text(size=22))+
    theme(plot.title = element_text(hjust = 0.5))) 
  assign(paste0('t1e_', i_n1), new_plot)
}

ggarrange(t1e_1, t1e_2, t1e_3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")


######################### APSAC trial example ###################################################################
n_1 <- 90

## Calculation of performance measures for n1=50, pC=0.3
pC_sel <- 0.05
power_results_apsac <- NULL
sample_size_results_apsac <- NULL
score_results_apsac <- NULL
lambda_as <- 0.295 ## assumed treatment effect for fixed design
for(i_rule in 1:length(list_rec_rules)){
  
  rec_rule <- list_rec_rules[i_rule]
  results <- fread(paste0('results/apsac_results/', rec_rule, '_n1_', n_1, '.csv'))
  
  ## power 
  power_results_apsac %<>% cbind(results[pC == pC_sel,]$power)
  
  ## sample size plot
  unc_sample_size <- results[pC == pC_sel,]$mean_n*(1-results[pC == pC_sel,]$futility-results[pC == pC_sel,]$efficacy) + 
    n_1*(results[pC == pC_sel,]$futility+results[pC == pC_sel,]$efficacy)
  sample_size_results_apsac %<>% cbind(unc_sample_size)
  
  ## conditional score results
  cond_results <- results[pC == pC_sel,] %$% calculate_conditional_performance_score(p_meanN = mean_n, p_varN = var_n,
                                                                                     p_meanCP = mean_obs_cp, p_varCP = var_obs_cp,
                                                                                     p_Nfix = n_fix, p_n1 = n_1, p_Nmax = Nmax,
                                                                                     p_alpha = alpha_glob, p_beta = beta, p_lambda = lambda)
  score_results_apsac %<>% cbind(cond_results$score_cond)
  
}

## save power results
power_results_apsac %<>% as.data.frame
names(power_results_apsac) <- str_remove(list_rec_rules, 'rec_')
row.names(power_results_apsac) <- results[pC == pC_sel,]$lambda
saveRDS(power_results_apsac, paste0('results/for_paper/power_results_apsac_n1_', n_1, '_pC_', pC_sel,'.RDS'))

## save sample size results
sample_size_results_apsac %<>% as.data.frame
names(sample_size_results_apsac) <- str_remove(list_rec_rules, 'rec_')
row.names(sample_size_results_apsac) <- results[pC == pC_sel,]$lambda
saveRDS(sample_size_results_apsac, paste0('results/for_paper/sample_size_results_apsac_n1_', n_1, '_pC_', pC_sel,'.RDS'))

## save score results
score_results_apsac %<>% as.data.frame
names(score_results_apsac) <- str_remove(list_rec_rules, 'rec_')
row.names(score_results_apsac) <- results[pC == pC_sel,]$lambda
saveRDS(score_results_apsac, paste0('results/for_paper/score_results_apsac_n1_', n_1, '_pC_', pC_sel,'.RDS'))

## statistics fixed design
fix_power <- c()
for(i_lambda in 1:nrow(results[pC == pC_sel,])){
  Nfix_new <- calculate_Nfix(lambda_as, 0.8, 0.025)
  fix_pow_new <- calculate_fix_power(results[pC == pC_sel,]$lambda[i_lambda], Nfix_new, 0.025)
  fix_power <- c(fix_power, fix_pow_new)
}
stats_fixed <- data.frame(lambda = results[pC == pC_sel,]$lambda, fix_power = fix_power, Nfix = calculate_Nfix(lambda_as, 0.8, 0.025))

## Global Power
data_plot_power <- cbind(lambda = as.numeric(row.names(power_results_apsac)), power_results_apsac)
data_plot_power <-  pivot_longer(data = data_plot_power, cols = !lambda, names_to = 'rec_rule', values_to = 'power')
pl_power <- ggplot(data=data_plot_power, aes(x=lambda, y=power, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  geom_line(data=stats_fixed,aes(x=lambda,y=fix_power, color = "fixed"), linetype="dashed", size = 1.1)+
  ylim(0, 1)+
  labs(x = expression(lambda))+
  labs(y = 'power')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle("(Global) Power")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = "rec rules", values = c("fixed" = "red", "classicGS"="#F8766D", "OCP" = "#CD9600", 
                                                    "OptFunc" = "#00BE67", "Promising"= "#00A9FF", "restrOCP"="#C77CFF"))

## Global mean sample size
data_plot_sample_size <- cbind(lambda = as.numeric(row.names(sample_size_results_apsac)), sample_size_results_apsac)
data_plot_sample_size <-  pivot_longer(data = data_plot_sample_size, cols = !lambda, names_to = 'rec_rule', values_to = 'sample_size')
pl_sample <- ggplot(data=data_plot_sample_size, aes(x=lambda, y=sample_size, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  geom_line(data=stats_fixed,aes(x=lambda,y=Nfix, color = "fixed"), linetype="dashed", size = 1.1)+
  ylim(0, 270)+
  labs(x = expression(lambda))+
  labs(y = 'mean sample size')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle("")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = "rec rules", values = c("fixed" = "red", "classicGS"="#F8766D", "OCP" = "#CD9600", 
                                                    "OptFunc" = "#00BE67", "Promising"= "#00A9FF", "restrOCP"="#C77CFF"))
ggarrange(pl_power, pl_sample, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")


## Conditional performance score for binary endpoints
data_plot_score <- cbind(lambda = as.numeric(row.names(score_results_apsac)), score_results_apsac)
data_plot_score <-  pivot_longer(data = data_plot_score, cols = !lambda, names_to = 'rec_rule', values_to = 'score')
pl_score_bin <- ggplot(data=data_plot_score, aes(x=lambda, y=score, group=rec_rule)) +
  geom_line(aes(color=rec_rule), size = 1.1)+
  ylim(0, 1)+
  labs(x = expression(lambda))+
  labs(y = expression("S("*lambda*")"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle("")+
  theme(plot.title = element_text(size=22))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(name = "rec rules", values = c("fixed" = "red", "classicGS"="#F8766D", "OCP" = "#CD9600", 
                                                    "OptFunc" = "#00BE67", "Promising"= "#00A9FF", "restrOCP"="#C77CFF"))

## all three performance measures for different values of pc
three_perf_meas <- ggarrange(pl_sample, pl_power, pl_score_bin, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
annotate_figure(three_perf_meas)
