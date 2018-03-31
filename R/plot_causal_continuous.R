################################################################################
# Goal: Plot the performance of the Copula Factor PC algorithm, the greedy 
#       step-wise method and fully-integrated method respecively.
#
# Metrics: TPR, FPR, and SHD
#
# Settings: fully Gaussian; 
#           all factors have multiple indicators;
#           k = {4, 10};
#           n = {500, 1000, 2000}.
################################################################################

#### 0. Load Packages ####

##
#install.packages("psych")
#install.packages("tikzDevice")
library(psych)
library(tikzDevice)
library(dplyr)
library(psych)
library(ggplot2)
library(gridExtra)

## help function
data2stat <- function(data){
  
  describe(data) %>%
    mutate(ci = qt(0.975,n-1) * se) %>%
    mutate(Methods = c(rep('MBPC',3), rep('GSPC',3), rep('CFPC',3))) %>%
    mutate(Sample.size = as.factor(rep(c(500,1000,2000),3))) %>%
    select(Methods, Sample.size, mean, ci)
}


#### 1. Read Data ####

p4 = read.table("results/causal_discovery/TPR_FPR_SHD_p4_all_conti.txt", header = F)
p10 = read.table("results/causal_discovery/TPR_FPR_SHD_p10_all_conti.txt", header = F)

# each data frame above is 100 * 27 (dimension)
# 100 means repeated experiments times


#### 2.  Plot ####

pd <- position_dodge(0.3) # move them to the left and right
par(mfrow = c(2,2))
par(mar=c(3.5,4,2,1.5))

#### Figure

## plot 1: k = 4, TPR
# statistical values
stat = data2stat(p4[,1:9])
# plot
TPR.k4 = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) + coord_cartesian(ylim = c(0, 1)) +
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('sample size') + ylab('TPR') + ggtitle('TPR (k=4)') + 
  theme_bw() +
  theme(legend.justification=c(0.99,0.01),legend.position=c(0.99,0.01),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 2: k = 4, FPR
# statistical values
stat = data2stat(p4[,10:18])
# plot
FPR.k4 = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) +
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('sample size') + ylab('FPR') + ggtitle('FPR (k=4)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 3: k = 4, SHD
# statistical values
stat = data2stat(p4[,19:27])
# plot
SHD.k4 = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) +
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('sample size') + ylab('SHD') + ggtitle('SHD (k=4)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 4: k = 10, TPR
# statistical values
stat = data2stat(p10[,1:9])
# plot
TPR.k10 = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) + coord_cartesian(ylim = c(0.1, 1)) +
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('sample size') + ylab('TPR') + ggtitle('TPR (k=10)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 5: k = 10, FPR
# statistical values
stat = data2stat(p10[,10:18])
# plot
FPR.k10 = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) +
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('sample size') + ylab('FPR') + ggtitle('FPR (k=10)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 6: k = 10, SHD
# statistical values
stat = data2stat(p10[,19:27])
# plot
SHD.k10 = ggplot(stat, aes(x = Sample.size, y = mean, group = Methods, color = Methods)) +
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('sample size') + ylab('SHD') + ggtitle('SHD (k=10)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))


#### 3. Show and Save the Figure ####

##
grid.arrange(TPR.k4, FPR.k4, SHD.k4, TPR.k10, FPR.k10, SHD.k10, nrow = 2, ncol = 3)

# ##
# pdf(file = 'results/causal_discovery/TFPR_SHD_conti.pdf', width = 7.2, height = 4.05)
# ## arrange
# grid.arrange(TPR.k4, FPR.k4, SHD.k4, TPR.k10, FPR.k10, SHD.k10, nrow = 2, ncol = 3)
# 
# dev.off()

