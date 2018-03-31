###############################################################################
# Goal: Plot the performance of the Copula Factor PC algorithm, the greedy 
#       step-wise method and fully-integrated method respecively.
#
# Metrics: TPR, FPR, and SHD (we plot three figures for the three metrics)
#
# Settings: Case 1, (4 factors where all have more measurements)
#           Case 2, (4 factors where 2 have more measurements while 2 have one)
#           Case 3, (4 factors where all have more measurements)
#           Case 4, (4 factors where 2 have more measurements while 2 have one)
###############################################################################

#### 0. Load Packages ####

##
library(psych)
library(dplyr)
library(ggplot2)
library(gridExtra)

## help function
data2stat <- function(data){
  stat_data = describe(data) %>%
    mutate(ci = qt(0.975,n-1) * se) %>%
    mutate(Methods = c(rep('MBPC',3),rep('GSPC',3),rep('CFPC',3))) %>%
    mutate(sample_size = as.factor(rep(c(500,1000,2000),3))) %>%
    select(Methods, sample_size, mean, ci)
}

#### 1. Read Data ####

p4_all = read.table("results/causal_discovery/TPR_FPR_SHD_p4_all.txt", header = F)
p4_half = read.table("results/causal_discovery/TPR_FPR_SHD_p4_half.txt", header = F)
p10_all = read.table("results/causal_discovery/TPR_FPR_SHD_p10_all.txt", header = F)
p10_half = read.table("results/causal_discovery/TPR_FPR_SHD_p10_half.txt", header = F)


#### 2. Plot ####

#
pd <- position_dodge(0.3) # move them to the left and right

#### Figure 1: TPR

## plot 1: for case 1
# statistical values
stat_p4_all = data2stat(p4_all[,1:9])
# ggplot
TPR_p4_all <- ggplot(stat_p4_all, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) + 
  geom_line(position=pd) + 
  geom_point(aes(shape=Methods), position=pd) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('TPR') +ggtitle('Only multi-indicators (k=4)') + 
  theme_bw() +
  theme(legend.justification=c(0.99,0.01),legend.position=c(0.99,0.01), legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 2: for case 2
# statistical values
stat_p4_half = data2stat(p4_half[,1:9])
# ggplot
TPR_p4_half <- ggplot(stat_p4_half, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) +  
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('TPR') +ggtitle('Also single indicator (k=4)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 3: for case 3
# statistical values
stat_p10_all = data2stat(p10_all[,1:9])
# ggplot
TPR_p10_all <- ggplot(stat_p10_all, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('TPR') +ggtitle('Only multi-indicators (k=10)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 4: for case 4
# statistical values
stat_p10_half = data2stat(p10_half[,1:9])
# ggplot
TPR_p10_half <- ggplot(stat_p10_half, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) +  
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('TPR') +ggtitle('Also single indicator (k=10)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## show Figure
grid.arrange(TPR_p4_all, TPR_p4_half, TPR_p10_all, TPR_p10_half, nrow = 2, ncol = 2)
## save Figure
pdf(file = 'results/causal_discovery/TPR.pdf', width = 5, height = 4.8)
grid.arrange(TPR_p4_all, TPR_p4_half, TPR_p10_all, TPR_p10_half, nrow = 2, ncol = 2)
dev.off()

#### Figure 2: FPR

## plot 1: for case 1
# statistical values
stat_p4_all = data2stat(p4_all[,10:18])
# ggplot
TPR_p4_all <- ggplot(stat_p4_all, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('FPR') +ggtitle('Only multi-indicators (k=4)') + 
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 2: for case 2
# statistical values
stat_p4_half = data2stat(p4_half[,10:18])
# ggplot
TPR_p4_half <- ggplot(stat_p4_half, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('FPR') +ggtitle('Also single indicator (k=4)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 3: for case 3
# statistical values
stat_p10_all = data2stat(p10_all[,10:18])
# ggplot
TPR_p10_all <- ggplot(stat_p10_all, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('FPR') +ggtitle('Only multi-indicators (k=10)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 4: for case 4
# statistical values
stat_p10_half = data2stat(p10_half[,10:18])
# ggplot
TPR_p10_half <- ggplot(stat_p10_half, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('FPR') +ggtitle('Also single indicator (k=10)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## show Figure
grid.arrange(TPR_p4_all, TPR_p4_half, TPR_p10_all, TPR_p10_half, nrow = 2, ncol = 2)
# ## save Figure
# pdf(file = 'results/causal_discovery/FPR.pdf', width = 5, height = 4.8)
# grid.arrange(TPR_p4_all, TPR_p4_half, TPR_p10_all, TPR_p10_half, nrow = 2, ncol = 2)
# dev.off()

#### Figure 3: SHD

## plot 1: for case 1
# statistical values
stat_p4_all = data2stat(p4_all[,19:27])
# ggplot
TPR_p4_all <- ggplot(stat_p4_all, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 3.5)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('SHD') +ggtitle('Only multi-indicators (k=4)') + 
  theme_bw() +
  theme(legend.justification=c(0.01,0.01),legend.position=c(0.01,0.01), legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 2: for case 2
# statistical values
stat_p4_half = data2stat(p4_half[,19:27])
# ggplot
TPR_p4_half <- ggplot(stat_p4_half, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('SHD') +ggtitle('Also single indicator (k=4)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 3: for case 3
# statistical values
stat_p10_all = data2stat(p10_all[,19:27])
# ggplot
TPR_p10_all <- ggplot(stat_p10_all, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('SHD') +ggtitle('Only multi-indicators (k=10)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## plot 4: for case 4
# statistical values
stat_p10_half = data2stat(p10_half[,19:27])
# ggplot
TPR_p10_half <- ggplot(stat_p10_half, aes(x=sample_size, y=mean, group=Methods, color=Methods)) + 
  geom_line(position=pd) +
  geom_point(aes(shape=Methods), position=pd) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) +
  xlab('sample size (n)') + ylab('SHD') +ggtitle('Also single indicator (k=10)') + 
  theme_bw() + theme(legend.position = 'none') +
  theme(plot.title = element_text(hjust = 0.5,  size = 10))

## show Figure
grid.arrange(TPR_p4_all, TPR_p4_half, TPR_p10_all, TPR_p10_half, nrow = 2, ncol = 2)
# ## save Figure
# pdf(file = 'results/causal_discovery/SHD.pdf', width = 5, height = 4.6)
# grid.arrange(TPR_p4_all, TPR_p4_half, TPR_p10_all, TPR_p10_half, nrow = 2, ncol = 2)
# dev.off()
