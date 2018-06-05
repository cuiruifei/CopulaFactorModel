##################################################################################################
## Goal: 1. Learn the parameters of the hypothesised factor structure from the 
#           'HolzingerSwineford1939' data provided in Package 'lavaan';
#        2. Compare BGCF (Bayesian Gaussian copula factor) with MLR (robust maximum likelihood)
#           quantitatively with this dataset.
###################################################################################################

#### 0. Dependent Packages and Parameter Settings ####

## packages
library(caret)
library(lavaan)
library(moments)
source('R/inferCopulaFactorModel.R')
source('R/copulaPredict.R')

## parameters
# No. of latent factors
pl <- 3


#### 1. Load Data ####

## raw data set
data(HolzingerSwineford1939)

## select the 9 variables involved
HS9 <- dplyr::select(HolzingerSwineford1939, x1:x9)

## summary of the data set
# VU: No. of unique values of a variable
# p.value: the p-value of normality test
summary.HS9 <- data.frame(VU = apply(HS9, 2, function(x) length(unique(x))), 
                             p.value = apply(HS9, 2, function(x) shapiro.test(x)$p.value),
                             skewness = skewness(HS9),
                             kurtosis = kurtosis(HS9))
# show the summary statistics
summary.HS9


## model description to be used in Method 1 (MLR)
HS.model <- ' visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9 '

## sparsity pattern representation of factor loadinds: lambda to be used in Method 2 (BGCF)
lambda <- matrix(0, 9, 3)
lambda[1:3, 1] = lambda[4:6, 2] = lambda[7:9, 3] = 1


#### 2. Learn Model Parameters with BGCF ####

## inference
HS9.obj <- inferCopulaFactorModel(Y = HS9, Lambda = lambda, nsamp = 1000)
# get Sigma (the correlation matrix over integrated vector)
Sigma <- apply(HS9.obj$Sigma.psamp[,, 501:1000], c(1,2), mean)
# C: correlation matrix over latent factors
C <- Sigma[1:pl, 1:pl]
# Lambda: matrix of factor loadings
Lambda <- round(t(solve(C) %*% Sigma[1:pl, -(1:pl)]), 3)
# D: residual variance
D <- round(Sigma[-(1:pl), -(1:pl)] - Lambda %*% C %*% t(Lambda), 3)


#### 3. Perform Prediction (as a regression problem) ####


#### Method 1: MLR

## help function
lavaanPredict <- function(lav.obj, ind = 9, newData){
  #
  # args:
  #   lav.obj: fitted model
  #   ind: the index of outcome variable 
  #   newData: new input data
  #
  # output: 
  #   predicted values
  #
  S = fitted(lav.obj)$cov
  mu = fitted(lav.obj)$mean
  p = ncol(S)
  b = solve(S[(1:p)[-ind], (1:p)[-ind]], S[(1:p)[-ind],ind])
  b0 = mu[ind] - sum(b * mu[(1:p)[-ind]])
  as.matrix(newData) %*% matrix(b) + b0
}

## help function
test.ml <- function(a, HS9){
  #
  # args: 
  #   a: the set of index of the train samples
  #   HS9: the full data set
  # 
  # output: 
  #   MSE and R-squared on the test set
  
  ## fit the model to training set
  model.fit = cfa(HS.model, data = HS9[a,], auto.var = T, auto.cov.lv.x = T, std.lv = T, 
                  meanstructure = TRUE)
  ##
  p = ncol(HS9)
  ##
  MSE = numeric(p)
  R2 = numeric(p)
  
  ## take i-th variable as outcome alternately
  for(i in 1:p){
    # predicted values
    values.pre = lavaanPredict(model.fit, ind = i, HS9[-a,-i])
    # true values
    values.true = HS9[-a,i]
    # MSE
    MSE[i] = mean((values.pre - values.true)^2)
    # R squared
    R2[i] = 1 - MSE[i]/var(values.true)
  }
  # return
  c(MSE, R2)
}

## run MLR 
# take the first 200 samples as the training set while the rest as the testing set
result <- test.ml(1:200, HS9 = HS9)
# show result
matrix(result, 9, dimnames = list(paste0('Task', 1:9), c('MSE', 'R-squared')))


#### Method 2

## help function
test.cf <- function(a, HS9){
  #
  # args: 
  #   a: the set of index of the train samples
  #   HS9: the full data set
  # 
  # output: 
  #   MSE and R-squared on the test set
  
  ## train
  cop.fac.obj = inferCopulaFactorModel(HS9[a,], odens = 1, Lambda = lambda, nsamp = 500, verb = F)
  Sigma = apply(cop.fac.obj$Sigma.psamp[,,101:500], c(1,2), mean)
  S = Sigma[4:12, 4:12]
  p = ncol(S)
  MSE = numeric(p)
  R2 = numeric(p)
  ## test
  for(i in 1:p){
    Y = HS9
    Y[-a,i] = NA
    pre.obj = copulaPredict(S, Y, nsamp = 500, verb = F, odens = 1)
    values.pre = apply(pre.obj$Y.impute[-a, i, 101:500], 1, mean)#
    # true values
    values.true = HS9[-a,i]
    # MSE
    MSE[i] = mean((values.pre - values.true)^2)
    # R squared
    R2[i] = 1 - MSE[i]/var(values.true)
  }
  # return
  c(MSE, R2)
}

## run BGCF
# take the first 200 samples as the training set while the rest as the testing set
result <- test.cf(1:200, HS9 = HS9)
# show result
matrix(result, 9, dimnames = list(paste0('Task', 1:9), c('MSE', 'R-squared')))


#### 4. Perform Prediction with Cross-validation ####


seed <- 12345
set.seed(seed)

## No. of times for k-fold cross-validation
iters = 10

## matrix to store results
# dimension: 100 * 18
# columns 1:9, MSE estimates for the 9 tasks
# columns 10:18, R-squared estimates for the 9 tasks
# rows 100 = 10*10, 10 times 10-fold cross-validation
results.ml <- NULL
results.cf <- NULL

#### perform multiple (iters) k-fold cross-validation
for(i in 1:iters){
  
  ## create folds
  cv_splits <- createFolds(1:nrow(HS9), k = 10, returnTrain = TRUE)
  
  ## perform cross validation
  # for method 1
  cv.ml <- lapply(cv_splits, test.ml, HS9 = HS9)
  # for method 2
  cv.cf <- lapply(cv_splits, test.cf, HS9 = HS9)
  
  ## collect estimates
  for(j in 1:length(cv.ml)) results.ml=rbind(results.ml,cv.ml[[j]])
  for(j in 1:length(cv.cf)) results.cf=rbind(results.cf,cv.cf[[j]])
  
  #
  print(c('iters:', i))
}

#### Show the results
## the mean over 100 experiments for MLR
matrix(apply(results.ml, 2, mean), 9, dimnames = list(paste0('Task', 1:9), c('MSE', 'R-squared')))
## the mean over 100 experiments for BGCF
matrix(apply(results.cf, 2, mean), 9, dimnames = list(paste0('Task', 1:9), c('MSE', 'R-squared')))


#### 5. Visualize the Results (Cross-validation) ####

#### read saved results
results = read.table(file = 'results/factor_analysis/MSE_R2_HS.txt')
# 1:9, MSE.ml over the 9 variables
# 10:18, R2.ml
# 19:27, MSE.cf
# 28:36, R2.cf

#### show results in number
## the mean over 100 experiments (10 * 10) for each variables
mean.each = rbind(MSE.ml = apply(results[, 1:9], 2, mean), MSE.cf = apply(results[, 19:27], 2, mean),
                  R2.ml = apply(results[, 10:18], 2, mean), R2.cf = apply(results[, 28:36], 2, mean))
mean.each
## the standard error of each mean estimate above
mean.se = rbind(MSE.ml = apply(results[, 1:9], 2, sd), MSE.cf = apply(results[, 19:27], 2, sd),
                R2.ml = apply(results[, 10:18], 2, sd), R2.cf = apply(results[, 28:36], 2, sd))/sqrt(nrow(results))
mean.se
## the mean over all experiments over all variables
apply(mean.each, 1, mean)

#### plot the results
## statistical results
stat.results = data.frame(outcome = rep(paste0('Y', 1:9), 2), methods = c(rep('MLR', 9), rep('BGCF', 9)), MSE = c(mean.each[1,], mean.each[2,]), 
                          MSE.se = c(mean.se[1,], mean.se[2,]), R2 = c(mean.each[3,], mean.each[4,]), R2.se = c(mean.se[3,], mean.se[4,]))
## for plot
pd <- position_dodge(0.9) 
## MSE
MSE.fig = ggplot(stat.results, aes(x = outcome, y = MSE, group = methods, fill = methods)) + coord_cartesian(ylim = c(0.4, 1.3)) +
  geom_bar(stat = 'identity', position = pd) +
  geom_point(aes(shape=methods), position=pd) + 
  xlab('outcome variable') +
  geom_errorbar(aes(ymin=MSE-MSE.se, ymax=MSE+MSE.se), width=.2, position=pd) + 
  theme_bw() +
  theme(legend.justification=c(0.99,0.99),legend.position=c(0.99,0.99),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=10))
## MSE
R2.fig = ggplot(stat.results, aes(x = outcome, y = R2, group = methods, fill = methods)) + #coord_cartesian(ylim = c(0.4, 1.3)) +
  geom_bar(stat = 'identity', position = pd) +
  geom_point(aes(shape=methods), position=pd) + 
  ylab('R-squared') + xlab('outcome variable') +
  geom_errorbar(aes(ymin=R2-R2.se, ymax=R2+R2.se), width=.2, position=pd) + 
  theme_bw() +
  theme(legend.justification=c(0.99,0.99),legend.position=c(0.99,0.99),legend.title = element_blank(), legend.key.size = unit(0.35, "cm")) +
  theme(text = element_text(size=10))

#### show figures
## MSE
MSE.fig
## R2
R2.fig

# #### save plots into files
# ## MSE
# pdf(file = 'result/factor_analysis/HS_MSE.pdf', width = 5, height = 2)
# MSE.fig
# dev.off()
# ## both MSE and R2
# pdf(file = 'result/factor_analysis/HS_MSE_R2.pdf', width = 5, height = 4)
# gridExtra::grid.arrange(MSE.fig, R2.fig,  nrow = 2, ncol = 1)
# dev.off()
