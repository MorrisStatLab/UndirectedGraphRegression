##load Pathway data
## immune, growth hormone, and angiogenesis pathways plus alpha fetaprotein (AFP)
load('Pathway.rda')
load('Hcc.profile.rda')
load('HepatoScore.rda')
yy = t(hcc.pro[gene.id, ])
xx = cbind(1 - score, score)
x = as.numeric(xx[,2])
# lambdav; bandwidthv
#######################################################################################################################################
# Nonparametric Finite Mixture of Gaussian Graphical Models #
#######################################################################################################################################
lambdav = seq(0.01, 2.01, length.out = 101)
bandwidthv = seq(0.01, 2.01, length.out =101)
bicl = c()
lam = c()
nG = dim(yy)[2]
auc1l = c(); auc2l = c();

# Required packages
library(mclust)
library(mvtnorm)
library(glasso)

log.like <- function(X, MU, COVARIANCE, PI){
  return(dmvnorm(X,MU,COVARIANCE, log=TRUE))
}

sum.log.like <- function(X, z, u, est.mu, est.covariance, est.pi){
  sum.log.like.value <- 0
  for(i in 1:dim(X)[1]){
    topz = ceiling(10*z[i])/10
    botz = floor(10*z[i])/10
    if(topz == botz){
      g=which(u==topz)
      sum.log.like.value = sum.log.like.value + log.like(X[i,], est.mu[,,g], est.covariance[,,g,], est.pi[,g])
    }else{
      ratz = (z[i] - botz)/(topz - botz) # (x - x0)/(x1 - x0)
      topg = which(u == topz)
      botg = which(u == botz)
      top.mu = est.mu[,,topg]
      bot.mu = est.mu[,,botg]
      inpl.mu = bot.mu + (top.mu - bot.mu)*ratz
      top.cov = est.covariance[,,topg,]
      bot.cov = est.covariance[,,botg,]
      inpl.cov = bot.cov + (top.cov - bot.cov)*ratz
      top.pi = est.pi[,topg]
      bot.pi = est.pi[,botg]
      inpl.pi = bot.pi + (top.pi - bot.pi)*ratz
      sum.log.like.value = sum.log.like.value + log.like(X[i, ], inpl.mu, inpl.cov, inpl.pi)
    }
  }
    return(sum.log.like.value)
}

kepan <- function(y){
  0.75*(1-y^2)*(abs(y)<=1) # Epanechinikov kernel
}

df <- function(t, K, est.precision){
  df.comps <- numeric(K)
  for(k in 1:K){
    df.comps[k] <- sum(est.precision[,,t,k][upper.tri(est.precision[,,t,k])] != 0)
  }
  return(sum(df.comps))
}

# Function for BIC score

# Inputs:
# Data: X
# The set of grid points: u (use the set of grid points such that the number of grid points is a divisor of a number of observations)
# Number of mixtures: K 
# Bandwidth: h
# Estimated mixing proportions (K*nu): est.pi
# Estimated covariance matrices (p*p*nu*K): est.covariance
# Estimated precision matrices (p*p*nu*K): est.precision

# Output:
# BIC score

BIC_score <- function(X, z, u, K, h, est.pi, est.covariance, est.precision){
  p <- dim(X)[2] 
  N <- dim(X)[1] 
  nu <- length(u) # number of grid points
  est.mu <- array(0, c(K,p,nu))
  
  grid.df <- numeric(nu)
  for(t in 1:nu){
    grid.df[t] <- df(t, K, est.precision) + K*p
  }
  sum.df <- sum(grid.df)
  avg.df <- (1/nu)*sum.df
  
  sum.log.like.value <- sum.log.like(X, z, u, est.mu, est.covariance, est.pi)
  BICscore <- -2*sum.log.like.value + log(N)*((K - 1) + avg.df)*(2.1153*(1/h)*0.45)
  return(BICscore)
}

# Algorithm 1: Proposed generalized effective EM algorithm
# Inputs:
# Data: X
# Covariate: z
# The set of grid points: u (use the set of grid points such that the number of grid points is a divisor of a number of observations)
# Number of mixtures: K 
# Bandwidth: h
# Penalization parameter: rho
# Number of maximum interations: niter
# Convergence criterion: toll

# Outputs:
# Estimated mixing proportions (K*nu)
# Estimated covariance matrices (p*p*nu*K)
# Estimated precision matrices (p*p*nu*K)

NFMGGM <- function(X, z, u, K, h, rho, niter, toll){
  ## Xn: normal samples
  p <- dim(X)[2] 
  N <- dim(X)[1]
  nu <- length(u) # number of grid points
  gn <- N%/%nu 

  pi.f <- array(0, c(K,nu))
  mu.f <- array(0, c(K,p,nu))
  covariance.f <- array(0, c(p,p,nu,K))
  precision.f <- array(0, c(p,p,nu,K))
  
  iz_sort <- order(z, decreasing=FALSE)
  X_sort <- X[iz_sort,]
  # Set initial values for the EM algorithm using mclust package
  for(i in 1:nu){
    if(i == nu){
      pi.f[,i] <- summary(Mclust(X_sort[((i-1)*gn+1):N,], G=K))$pro
      covariance.f[,,i,] <- summary(Mclust(X_sort[((i-1)*gn+1):N,], G=K))$variance
    }else{
      pi.f[,i] <- summary(Mclust(X_sort[((i-1)*gn+1):(i*gn),], G=K))$pro
      covariance.f[,,i,] <- summary(Mclust(X_sort[((i-1)*gn+1):(i*gn),], G=K))$variance
    }

  }
  
  pi.f.new <- array(0, c(K,nu))
  covariance.f.new <- array(0, c(p,p,nu,K))
  precision.f.new <- array(0, c(p,p,nu,K))
  gamma <- matrix(1,K,N)

  Num.pi.f <- matrix(0,N,nu)
  Den.pi.f <- matrix(0,N,nu)
  Num.A.f <- array(0, c(p,p,N,nu))
  Den.A.f <- matrix(0,N,nu)

  for (m in 1:niter){
    for (t in 1:nu){
      for (k in 1:K){
        for (i in 1:N){
          Num.pi.f[i,t] <- gamma[k,i]*(1/h)*kepan((z[i]-u[t])/h)
          Den.pi.f[i,t] <- (1/h)*kepan((z[i]-u[t])/h)
          Num.A.f[,,i,t] <- gamma[k,i]*(1/h)*kepan((z[i]-u[t])/h)*((X[i,]-mu.f[k,,t])%*%t(X[i,]-mu.f[k,,t]))
          Den.A.f[i,t] <- gamma[k,i]*(1/h)*kepan((z[i]-u[t])/h)
        }
        pi.f.new[k,t] <- sum(Num.pi.f[,t])/sum(Den.pi.f[,t])
        covariance.f.new[,,t,k] <- glasso(apply(Num.A.f[,,,t],c(1,2),sum)/sum(Den.A.f[,t]),rho=rho)$w
        precision.f.new[,,t,k] <- glasso(apply(Num.A.f[,,,t],c(1,2),sum)/sum(Den.A.f[,t]),rho=rho)$wi
      }
    }
    
    conv <- matrix(0, K, nu)
    for (k in 1:K){
      for (t in 1:nu){
        conv[k,t] <- norm((precision.f[,,t,k] - precision.f.new[,,t,k]), "F")
      } 
    }
    conv.rule <- (1/nu)*sum(conv)
  
    if (conv.rule < toll){
      outputs <- list("est.pi" = pi.f.new, "est.covariance" = covariance.f.new, "est.precision" = precision.f.new)
      return(outputs)
    }
    else {
      pi.f <- pi.f.new
      covariance.f <- covariance.f.new
      precision.f <- precision.f.new
    }
  }
}
##################Prepare Input#######################
z = as.numeric(xx[,2])
for(m in 1:length(lambdav)){
  for(n in 1:length(bandwidthv)){
    lambda = lambdav[m] 
    bandwidth = bandwidthv[n] 
    start_time <- Sys.time()
    result = try(NFMGGM(X = yy, z = z, u = (0:10)/10, K = 1, h = bandwidth, rho = lambda, niter = 20, toll = 1e-04))
    if(inherits(result, "try-error")) next
    bicvl = BIC_score(X = yy, z=z, u = (0:10)/10, K = 1, h = bandwidth, est.pi = result$est.pi, est.covariance = result$est.covariance, est.precision = result$est.precision)
    bicl = c(bicl, bicvl)
    lam = rbind(lam, c(lambda, bandwidth))
    ##############evaluate auc by fixing similarity parameter##############
    lc1 = result$est.precision[,,1,1]; lc2 = result$est.precision[,,11,1]; lc3 = result$est.precision[,,6,1]
    l.heatmapR = matrix(0, ncol = 3, nrow = nG*(nG - 1)/2)
    h = 0 
    for(i in 1:(nG-1)){
      for(j in (i+1):nG){
        h = h + 1
        l.heatmapR[h, 1] = sum(lc1[i,j] != 0)
        l.heatmapR[h, 2] = sum(lc2[i,j] != 0)
        l.heatmapR[h, 3] = sum(lc3[i,j] != 0)
      }
    }
  end_time <- Sys.time()
  }
}

ibl = which.min(bicl)
best.lambda = lam[ibl, 1]
#
best.log.like = -1*.Machine$double.xmax
best.bandwidth = bandwidthv[1]
log.like.vc=c()
for(n in 1:length(bandwidthv)){
    nu = length((0:10)/10)
    bandwidth = bandwidthv[n] 
    log.like.v = 0
    est.mu <- array(0, c(1,nG,nu))
    print('start cross validation')
    indexs = sample(rep(1:5, length.out = dim(yy)[1]))

    for(kf in 1:5){
      print(kf)
      yy_train = yy[indexs!=kf,]
      yy_test = yy[indexs==kf,]
      z_train = z[indexs!=kf]
      z_test = z[indexs==kf]
      result = try(NFMGGM(X = yy_train, z = z_train, u = (0:10)/10, K = 1, h = bandwidth, rho = best.lambda, niter = 20, toll = 1e-04))
      if(inherits(result, "try-error")){
        log.like.v=-1*.Machine$double.xmax
        break
      }
      sum.log.like.value <- sum.log.like(yy_test, z_test, (0:10)/10, est.mu, result$est.covariance, result$est.pi)
      log.like.v = log.like.v + sum.log.like.value
    }
    log.like.vc = c(log.like.vc, log.like.v)
    #print(log.like.v)
    if(log.like.v > best.log.like){
        best.log.like = log.like.v
        best.bandwidth = bandwidth
    }
}
print(log.like.vc)
#print(best.bandwidth)
#print(best.log.like)
print(best.lambda)
print(best.bandwidth)
result = NFMGGM(X = yy, z = z, u = (0:10)/10, K = 1, h = best.bandwidth, rho = best.lambda, niter = 20, toll = 1e-04)
lc1 = result$est.precision[,,1,1]; lc2 = result$est.precision[,,11,1]; lc3 = result$est.precision[,,6,1]
l.heatmapR = matrix(0, ncol = 3, nrow = nG*(nG - 1)/2)
h = 0 
for(i in 1:(nG-1)){
  for(j in (i+1):nG){
    h = h + 1
    l.heatmapR[h, 1] = sum(lc1[i,j] != 0)
    l.heatmapR[h, 2] = sum(lc2[i,j] != 0)
    l.heatmapR[h, 3] = sum(lc3[i,j] != 0)
  }
}

save.image(file='nfm_case.rda')


