########################################################################################################
# Goal: Compare the Bayesian Gaussian copula factor (BGCF) approach, with the diagonally
#       weighted least squares (DWLS), and the robust maximum likelihood (MLR) on complete 
#       ordinal data (4 categories) for learning parameters of CFA models. 
#
# Metrics: 1) ARB; 2) ARMSE.
#
# Settings: 1) a correlated 4-factor model with specific parameters where each factor has 4 indicators
#           2) mixed Gaussian and ordinal (4 levels) indicators
#           3) n = {100, 200, 500, 1000}
#
########################################################################################################

#### 0. Load Packages ####

##
library(psych)
library(dplyr)
library(ggplot2)
library(gridExtra)

## help function
data2stat <- function(data){
  
  describe(data) %>%
    mutate(ci = qt(0.975,n-1) * se) %>%
    mutate(Methods = c(rep('MLR',4), rep('DWLS',4), rep('BGCF',4))) %>%
    mutate(Sample.size = as.factor(rep(c(100, 200, 500, 1000), 3))) %>%
    dplyr::select(Methods, Sample.size, mean, ci)
}


#### 1. Read Data ####

## raw data
results = read.table(file = "results/factor_analysis/k4_ordinal_complete.txt")

#### 2.  Plot Figure in PDF ####
pd <- position_dodge(0.9) # move them to the left and right
par(mfrow = c(2,2))
par(mar=c(3.5,4,2,1.5))

#### Figure 1: for interfactor correlations

## plot 1: ARB
stat = data2stat(results[,13:24])

ARB.corr = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) + coord_cartesian(ylim = c(-0.15, 0.15)) +
  geom_line() + 
  geom_point(aes(shape=Methods)) + 
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  geom_hline(yintercept = 0.1, linetype = "dotted") +
  geom_hline(yintercept = -0.1, linetype = "dotted") +
  xlab('sample size') + ylab('ARB') +
  #geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + 
  theme_bw() +
  theme(legend.justification=c(0.99,0.99),legend.position=c(0.99,0.99),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 2: RMSE
stat = data2stat(results[,25:36])

RMSE.corr = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, fill = Methods)) + 
  geom_bar(position=pd, stat="identity") + 
  geom_point(aes(shape=Methods), position=pd) + 
  xlab('sample size') + ylab('RMSE') +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + 
  theme_bw() +
  theme(legend.justification=c(0.01,0.01),legend.position=c(0.01,0.01),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## show Figure
grid.arrange(ARB.corr, RMSE.corr, nrow = 1, ncol = 2)
# ## save Figure
# pdf(file = 'results/factor_analysis/k4_complete_corr.pdf', width = 4.5, height = 1.9)
# grid.arrange(ARB.corr, RMSE.corr, nrow = 1, ncol = 2)
# dev.off()


#### Figure 2: for factor loadings

## plot 1: ARB
stat = data2stat(results[,37:48])

ARB.corr = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) + coord_cartesian(ylim = c(-0.15, 0.15)) +
  geom_line() + 
  geom_point(aes(shape=Methods)) + 
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(yintercept = -0.05, linetype = "dashed") +
  geom_hline(yintercept = 0.1, linetype = "dotted") +
  geom_hline(yintercept = -0.1, linetype = "dotted") +
  xlab('sample size') + ylab('ARB') +
  #geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + 
  theme_bw() +
  theme(legend.justification=c(0.99,0.99),legend.position=c(0.99,0.99),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 2: RMSE
stat = data2stat(results[,49:60])

RMSE.corr = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, fill = Methods)) + 
  geom_bar(position=pd, stat="identity") + 
  geom_point(aes(shape=Methods), position=pd) + 
  xlab('sample size') + ylab('RMSE') +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + 
  theme_bw() +
  theme(legend.justification=c(0.01,0.01),legend.position=c(0.01,0.01),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## show Figure
grid.arrange(ARB.corr, RMSE.corr, nrow = 1, ncol = 2)
# ## save Figure
# pdf(file = 'results/factor_analysis/k4_complete_loadings.pdf', width = 4.5, height = 1.9)
# grid.arrange(ARB.corr, RMSE.corr, nrow = 1, ncol = 2)
# dev.off()
