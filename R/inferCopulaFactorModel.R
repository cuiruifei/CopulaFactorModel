inferCopulaFactorModel <- function (Y, Lambda = diag(ncol(Y)), trueSigma = NULL, nsamp = 100, 
                                    odens = max(1, round(nsamp/1000)), impute = any(is.na(Y)), 
                                    plugin.threshold = 20, plugin.marginal = (apply(Y, 2, function(x) {
                                      length(unique(x))
                                    }) > plugin.threshold), verb = TRUE) 
{
  #
  # This is the main function to perform inference for Gaussian copula factor models.
  # 
  # Args:
  #   Y, a n by p data matrix
  #   Lambda, a p by k matrix representing the mapping from factors to observed variables
  #   nsamp, No. of samples
  #   For details about the arguements, refer to function 'sbgcop.mcmc' in R package 'sbgcop'
  #
  # Return: a list that contains samples of the correlation matrix over all variables (factors + observed variables)
  #   
  
  require(BDgraph)
  library(mvtnorm)
  library(sbgcop)
  Y <- as.matrix(Y)
  vnames <- colnames(Y)
  colnames(Y) <- vnames
  # sample size
  n <- dim(Y)[1]
  # No. of observed variables
  p <- dim(Y)[2]
  
  #### handle Lambda and get prior graph G
  Lambda = 1 * (Lambda!=0)
  # No. of factors
  k = ncol(Lambda)
  # index of factors with a single indicator
  index.1 = which(colSums(Lambda) == 1)
  # No. of factors with a single indicator 
  k1 = length(index.1)
  # index of factors with multiple indicators
  index.2 = which(colSums(Lambda) > 1)
  # No. of factors with multiple indicators
  k2 = length(index.2)
  ## get the pior graph G
  G1 = matrix(1,k,k)-diag(k)
  if (k1 == 0) G2 = t(Lambda) else G2 = t(Lambda[-index.1,]) 
  G3 = matrix(0,p-k1,p-k1)
  G = rbind(cbind(G1,G2),cbind(t(G2),G3))
  G[lower.tri(G)] = 0
  ## prior parameters for the G-Wishart distribution
  # degrees of freedom
  n0 = p+k2+1
  # scale matrix
  S0 = diag(p+k2)/n0
  
  #### initialize Z, eta, and S
  R <- NULL
  for (j in 1:p) {
    R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
  }
  Rlevels <- apply(R, 2, max, na.rm = TRUE)
  Ranks <- apply(Y, 2, rank, ties.method = "average", na.last = "keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t(t(Ranks)/(N + 1))
  Z <- qnorm(U)
  # # handle categorical variable
  # Z[, ind.cat] = Y[, ind.cat]
  #
  Zfill <- matrix(rnorm(n * p), n, p)
  Z[is.na(Y)] <- Zfill[is.na(Y)]
  # psuedo data of factors with a single indicator
  Z1 = eta1 = Z[, index.1]
  # psuedo data of response variables
  if (k1 == 0) Z2 = Z else Z2 = Z[, -index.1]
  # psuedo data of factors with multiple indicators
  eta2 = matrix(rnorm(n*k2), n)
  eta = cbind(eta1, eta2)
  X = cbind(eta, Z2)
  S <- cov(X)
  
  ####
  Y.pmean <- Y
  Z.pmean <- Z2
  if (impute) {
    Y.pmean <- matrix(0, nrow = n, ncol = p)
  }
  LPC <- NULL
  C.psamp <- array(dim = c(p+k2, p+k2, floor(nsamp/odens)))
  Y.imp <- NULL
  Z.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  if (impute) {
    Y.imp <- array(dim = c(n, p, floor(nsamp/odens)))
  }
  #dimnames(C.psamp) <- list(colnames(Y), colnames(Y), 1:floor(nsamp/odens))
  
  #### start of Gibbs sampling scheme
  for (ns in 1:nsamp) {
    
    ## sample Z1 (=eta1)
    for (j in index.1) {
      # if (!(j %in% ind.cat)){
        ind.tmp = (1:k)[-j]
        Sjc <- S[j, ind.tmp] %*% solve(S[ind.tmp, ind.tmp])
        sdj <- sqrt(S[j, j] - Sjc %*% S[ind.tmp, j])
        muj <- X[, ind.tmp] %*% t(Sjc)
        if (!plugin.marginal[j]) {
          for (r in 1:Rlevels[j]) {
            ir <- (1:n)[R[, j] == r & !is.na(R[, j])]
            lb <- suppressWarnings(max(X[R[, j] == r - 1, j], na.rm = TRUE))
            ub <- suppressWarnings(min(X[R[, j] == r + 1, j], na.rm = TRUE))
            X[ir, j] <- qnorm(runif(length(ir), pnorm(lb, muj[ir], sdj), pnorm(ub, muj[ir], sdj)), muj[ir], sdj)
          }
        }
        ir <- (1:n)[is.na(R[, j])]
        X[ir, j] <- rnorm(length(ir), muj[ir], sdj)
      # }else{
      #   ir <- (1:n)[is.na(R[, j])]
      #   X[ir, j] = rnorm(length(ir))
      # }
    }
    Z1 = eta1 = X[,index.1]
    
    ## sample Z2
    if (k1 == 0) index.tmp = sample(1:p) else index.tmp = sample((1:p)[-index.1])
    if (k2 > 0 ){
      for(j in index.tmp) {
        q <- which(Lambda[j,]!=0)
        a = S[k2+j, q] / S[q, q]
        sdj <- sqrt( S[k2+j,k2+j] - a * S[q,k2+j] )
        muj <- X[,q] * a
        if (!plugin.marginal[j]) {
          for(r in sort(unique(R[,j]))){
            ir<- (1:n)[R[,j]==r & !is.na(R[, j])]
            lb<-suppressWarnings(max( X[ R[,j]<r,j+k2],na.rm=T))
            ub<-suppressWarnings(min( X[ R[,j]>r,j+k2],na.rm=T))
            X[ir,j+k2]<-qnorm(runif(length(ir),
                                    pnorm(lb,muj[ir],sdj),pnorm(ub,muj[ir],sdj)),muj[ir],sdj)
          }
        }
        ir <- (1:n)[is.na(R[, j])]
        X[ir, j+k2] <- rnorm(length(ir), muj[ir], sdj)
      }
      Z2 = X[, (k+1):(k2+p)]
      
      #
      Z = cbind(Z1, Z2)
      
      ## sample eta2
      ind.tmp = (1:(p+k2))[-index.2]
      A <- S[index.2,ind.tmp] %*% solve(S[ind.tmp,ind.tmp])
      sdj <-  S[index.2,index.2] - A %*% S[ind.tmp,index.2]
      muj <- X[,ind.tmp] %*% t(A)
      X[,index.2] = muj + rmvnorm(n, sigma =sdj)
      eta2 = X[,index.2]
      
      eta = cbind(eta1, eta2)
      
      ## identification condition
      for (j in index.2) {
        X[,j] = X[,j] * sign(cov(X[,j], X[,k2+which(Lambda[,j] != 0)[1]]))
      }
    }
    
    ## relocate the mean to zero
    X = t( (t(X)-apply(X,2,mean)))
    
    ## sample S
    P <- rgwish(n = 1, adj.g = G, b = n+n0, D = S0*n0+crossprod(X))[,,1]
    S <- solve(P)
    S = cov2cor(S)
    
    if (ns%%odens == 0) {
      C <- S#/(sqrt(diag(S)) %*% t(sqrt(diag(S))))
      if (is.null(trueSigma))
        lpc <- ldmvnorm(scale(X), C)
      else
        lpc <- ldmvnorm(scale(X), trueSigma)
      LPC <- c(LPC, lpc)
      C.psamp[, , ns/odens] <- C #cor(X)#
      if (impute) {
        Y.imp.s <- Y
        for (j in 1:p) {
          Y.imp.s[is.na(Y[, j]), j] <- quantile(Y[,j], 
                                                pnorm(Z[is.na(Y[, j]), j], 0, sd(Z[,j])), 
                                                na.rm = TRUE, type = 1)
        }
        Y.imp[, , ns/odens] <- Y.imp.s
        Y.pmean <- ((ns/odens - 1)/(ns/odens)) * Y.pmean + (1/(ns/odens)) * Y.imp.s
      }
      Z.imp[,,ns/odens] = cbind(Z1, Z2)
    }
    if (verb == TRUE & (ns%%(odens * 10)) == 0) {
      cat(round(100 * ns/nsamp), "percent done ", 
          date(), "\n")
    }
  }
  #
  Z.pmean = apply(Z.imp, c(1,2), mean)
  
  # # 
  #G.ps <- list(Sigma.psamp = C.psamp, Y.pmean = Y.pmean, Y.impute = Y.imp, Z.pmean = Z.pmean, Z.impute = Z.imp, LPC = LPC)
  G.ps <- list(Sigma.psamp = C.psamp)
  class(G.ps) <- "psgc"
  return(G.ps)
}
