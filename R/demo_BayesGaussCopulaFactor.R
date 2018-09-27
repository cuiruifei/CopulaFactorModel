###########################################################################################################################
# Goal: This is a demo to show how our BGCF approach works.
# Author: Ruifei Cui
###########################################################################################################################


#### 0. Dependencies and Basic Setting ####

## dependent functions and packages
source('R/rmvCopulaFactorModel.R')
source('R/inferCopulaFactorModel.R')
library(infotheo)
library(mvtnorm)
library(BDgraph)
library(sbgcop)

## parameters
# sample size
n <- 1000
# No. of latent factors
pl <- 4


#### 1. Generate Data ####

## interfactor correlations (correlation matrix over factors)
set.seed(12345)
R <- matrix(runif(pl^2, 0.2, 0.4), ncol=pl)
R <- (R * lower.tri(R)) + t(R * lower.tri(R))
diag(R) <- 1
C0 <- R
## simulate data
data.obj <- rmvCopulaFactorModel(n, sigma = C0, range.indicators = 4, lambda.min = 0.7, lambda.max = 0.7, sd.residual = sqrt(0.51))
## factorloading matrix
Lambda0 <- data.obj$Lambda
## data of response random vector
Z <- data.obj$data
## data of observed random vector (here all are binary; one could do whatever he or she wants.)
Y <- discretize(Z, nbins = 2)


### 2. Run BGCF Approach ####

# burn-in
b.in <- 100
# No. of samples
n.samp <- 500
## inference
cop.fac.obj <- inferCopulaFactorModel(Y, Lambda0, nsamp = n.samp, rand.start = F, odens = 1, verb = T)
# 
Sigma.samp <- cop.fac.obj$Sigma.psamp[,,-(1:b.in)]
# get Sigma (the correlation matrix over integrated vector)
Sigma <- apply(Sigma.samp, c(1,2), mean)
# C: correlation matrix over latent factors
C <- Sigma[1:pl, 1:pl]
# Lambda: matrix of factor loadings
Lambda <- round(t(solve(C) %*% Sigma[1:pl, -(1:pl)]), 3)
# D: residual variance
D <- round(Sigma[-(1:pl), -(1:pl)] - Lambda %*% C %*% t(Lambda), 3)


#### 3. Convergence Diagnostic ####

## Autocorrelation
# for 
acf(Sigma.samp[5,3,], plot = T)

## Potential Scale Reduction Factor


