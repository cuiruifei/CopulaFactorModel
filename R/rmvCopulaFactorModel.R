################################################################################
# Goal: This function aims to generate data from a Gaussian copula factor model
#      given a graph or a covariance matrix in the latent space.
# Author: ...
################################################################################
rmvCopulaFactorModel <-function(n = 1000, g = NULL, sigma = NULL, impurity = 0, range.indicators = 3:10, lambda.min = 0.1, lambda.max = 1, sd.residual = 1){
  #library(infotheo)
  #library(pcalg)
  library(mvtnorm)
  
  #### arguments
  # range.indicators: range of the number of indicators per factor
  # lambda.min, lambda.max: the minimum and maximum of factor loadings
  
  #### 0. Check inputs and Initialization
  ## check input g and sigma
  if(!is.null(g)){
    ## pl: the number of latent factors (or variables)
    pl = length(g@nodes)
    ## sigma: the population covariance matrix
    sigma = trueCov(g)
  }else if(!is.null(sigma))
    pl = ncol(sigma)
  else
    stop('Input g and sigma connot be empty simutanously.')
  
  
  ## choose randomly some factors from pl as factors with more (>1) response variables
  # the number of factors with more (>1) response variables
  pm = round(pl)
  index_more = 1:pl#sample(1:pl, pm)
  ## choose the number of response variables for each element of index_more from unif(3,10)
  #num_res_var = sample(3:10, pm, replace = T)#rep(4, pm)#
  ifelse(length(range.indicators) == 1, num_res_var <-  rep(range.indicators, pm), num_res_var <-  sample(range.indicators, pm, replace = T))
  
  ## pr: the number of response variables.
  pr = sum(num_res_var)
  
  
  #### 1. Generate normal data in latent sapce
  ## d_gauss stands for \eta
  d_gauss <- rmvnorm(n, sigma = sigma)#rmvDAG(n,g)
  
  
  #### 2. Generate response data
  ## initializing
  Y=mat.or.vec(n,pl+pr)
  Y[,1:pl] = d_gauss
  
  ## randomly generate factor loadings (Lamda) from latent factors to response variables
  Lamda = matrix(0,pr,pm)
  index_res_var = mat.or.vec(pm,1)
  for (i in 1:pm) {
    index_res_var[i] = sum(num_res_var[1:i])
  }
  index_res_var = c(0,index_res_var)
  for (i in 1:pm) {
    Lamda[(index_res_var[i]+1):index_res_var[i+1],i] = runif(num_res_var[i], min = lambda.min, max = lambda.max)  
  }
  Lamda.pure = Lamda
  ## noise
  error = matrix(rnorm(pr*n, 0, sd.residual),n,pr)
  ## add some impurities
  if (impurity > 0){
    # the number of type 1 and 2 impurities
    ni.1 = round(impurity)
    ni.2 = impurity - ni.1
    # type 1: to Lamda
    if (ni.1 > 0){
      Lamda[sample(which(Lamda == 0), ni.1)] = runif(ni.1, 0.1, 1)
    }
    # type 2: to error
    if (ni.2 > 0){
      vec = rep(0, pr*(pr-1)/2)
      vec[sample(pr*(pr-1)/2, ni.2)] = runif(ni.2)
      sigma.error = diag(pr)
      sigma.error[upper.tri(sigma.error)] = vec
      sigma.error[lower.tri(sigma.error)] = t(sigma.error)[lower.tri(sigma.error)]
      error = rmvnorm(n, sigma = sigma.error)
    }
  }
  
  
  ## generate responsing data
  Y[,(pl+1):(pl+pr)] = d_gauss[,index_more] %*% t(Lamda) + error
  
  
  #### 3. Return
  list(data = Y[,-(1:pl)], Sigma = sigma, index_more = index_more, Lambda = Lamda.pure)
}