######################################################################################
# This is a demo to show how the Copula Factor PC algorithm works.
# Author: Ruifei Cui
######################################################################################


#### 0. Dependencies and Parameters ####

## dependent packages
library(pcalg)
library(infotheo)
library(BDgraph)
library(sbgcop)
source('R/rmvCopulaFactorModel.R')
source('R/inferCopulaFactorModel.R')

## parameters
# No. of factors
pl <- 6
# sample size
n <- 1000


#### 1. Simulate Data ####

## simulate a DAG 
# a random DAG, serving as the true DAG in latent space
g <- randomDAG(pl,2/(pl-1),lB = 0.1,uB=1)
# true CPDAG
g.cpdag <- dag2cpdag(g)
## generate data
data.obs <- rmvCopulaFactorModel(n, g = g, range.indicators = 4)
# mapping from latents to response variables
Lambda <- (data.obs$Lambda != 0) *1
# data in response space (fully Gaussian)
Z <- data.obs$data
# data in observed space (mixed continuous and ordinal)
p <- ncol(Z)
Y <- Z
for (i in sample(1:p, round(p/2))) {
  Y[,i] = matrix(unlist(discretize(Z[,i], nbins = sample(2:5,1))), byrow=FALSE, nrow=n )
}


#### 2. Inference ####

## inference
cop.fac.obj <- inferCopulaFactorModel(Y, Lambda = Lambda, nsamp = 1000)
# extract samples of the correlation matrix over latent variables, ignoring the first 500 samples (burn-in)
C.samples <- cop.fac.obj$Sigma.psamp[1:pl, 1:pl, 501:1000]
# the posterior mean
C <- apply(C.samples, c(1,2), mean)
# standard deviations
C.sd <- apply(C.samples, c(1,2), sd)
# effective sample size
C.ess <- ((1-C^2)/C.sd)^2

#### 3. Causal discovery ####

## call the order independent version of the standard PC algorithm
graph.cfpc <- pc(suffStat = list(C = C, n = mean(C.ess[upper.tri(C.ess)])), indepTest = gaussCItest, 
                 alpha = 0.05, p = pl, skel.method = "stable", maj.rule = T, solve.confl = T)

## show results
par(mfrow = c(1,2))
plot(g.cpdag, main = 'True Graph')
plot(graph.cfpc, main = 'Copula Factor PC')
