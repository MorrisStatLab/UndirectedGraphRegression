id = commandArgs(trailingOnly = T)
id = as.character(id[1])
load(paste0('all',id,'.rda'))
rm(list=setdiff(ls(), c("id", "omega.mat", "yn", "ye")))
library(LASICH)
library(msos)
source('roc.R')
########################################################################
## Input:
##  D: distance matrix
##  alpha: positive number
## Output:
##  L: graph Lapacian
########################################################################
D2L <- function(D,alpha=1){
  if(alpha>0){
    W <- 1/D^alpha
    W <- ifelse(abs(W)==Inf,0,W)
    diag(W) <- 0
  }else{
    W <- D
  }
  ds <- colSums(W)
  Ds <- (1/sqrt(ds))%*%t(1/sqrt(ds))
  Ds <- ifelse(abs(Ds)==Inf,0,Ds)
  L <- -W*Ds
  diag(L) <- 1
  return(L)
}

D1L <- function(D,alpha=1){
  if(alpha>0){
    W <- 1/D^alpha
    W <- ifelse(abs(W)==Inf,0,W)
    diag(W) <- 0
  }else{
    W <- D
  }
  ds <- diag((colSums(W))^(-0.5))
  return(diag(2) - ds%*%W%*%ds)
}

clust <- function(xs,mthd=1){
  x.mat <- NULL
  for(i in 1:length(xs)){
    x.mat <- rbind(x.mat,xs[[i]])
  }
  lab <- NULL
  for(i in 1:length(xs)){
    lab <- c(lab,rep(i,times=nrow(xs[[i]])))
  }
  ## distance
  d <- dist(x.mat,p=2)
  ## clustering
  if(mthd==1){
    a <- hclust(d,"complete")
  }else if (mthd==2){
    a <- hclust(d,"average")
  }else if (mthd==3){
    a <- hclust(d,"single")
  }
  estclass <- cutree(a,k=length(xs))
  print(table(lab,estclass))
  ## distance
  d.mat <- dist.clust(d=d,cls=lab,mthd=mthd)
  
}

########################################################################
## Input:
##  d: distance matrix
##  cls: vector of estimated cluster labels
##  mthd: clustering method 1=complete 2=average 3=single
## Output:
##  dmat: distance matrix for estimated clusters
########################################################################
dist.clust <- function(d,cls,mthd=1){
  num.cl <- length(unique(cls))
  d <- as.matrix(d)
  dmat <- diag(num.cl)*0
  for(i in 1:(num.cl-1)){
    for(j in (i+1):num.cl){
      elt1 <- which(cls==i)
      elt2 <- which(cls==j)
      mat <- d[elt1,elt2]
      if(is.null(dim(mat))){
        break
      }
      if(mthd==1){
        dmat[i,j] <- max(mat)
      }else if(mthd==2){
        dmat[i,j] <- sum(mat)/2/dim(mat)[1]/dim(mat)[2]
      }else if(mthd==3){
        dmat[i,j] <- min(mat)
      }
    }
  }
  dmat <- dmat+t(dmat)
  return(dmat)
}
#######a function to calculate the approximate AIC of fgl objects
jbic <- function(ff, yy){
  K = length(yy)
  aicv = 0
  for(k in 1:K){
    htheta <- ff$est[,,k]
    xx <- yy[[k]]
    n <- nrow(xx); p <- ncol(xx)
    S <- cov(xx)*(n-1)/n
    Ek <- (sum(htheta != 0) - p)/2
    tr <- sum(diag(S%*%htheta))
    tmp = n*tr - n*logdet(htheta) + Ek*log(n)
    aicv = aicv + tmp
  }
  return(aicv)
}
#evaluation function
TPR <- function(fit, truth){ #given truth = 1, fit = 1
  tp <- sum((fit == 1)&(truth == 1))
  fn <- sum((fit == 0)&(truth == 1))
  return(tp/(tp+fn))
}
FPR <- function(fit, truth){ #given truth = 0, fit = 1
  fp <- sum((fit == 1)&(truth == 0))
  tn <- sum((fit == 0)&(truth == 0))
  return(fp/(fp+tn))
}

##########################################################################
load(paste0('sim',id,'.rda'))
load(paste0('icl',id,'.rda'))
#y1 = yy[xx[,2] == 0,]; y2 = yy[xx[,2] > 0,]
y1 = scale(yn, center=T, scale=T)
y2 = scale(ye, center=T, scale=T)
# lambda1v; lambda2v
lambda1v = seq(0.01, 50, length.out = 101)
lambda2v = seq(0.01, 50, length.out = 101)
bicl = c()
lam = c()
nG = 20
auc1l = c(); auc2l = c();
###################define lSensitivity and lSpecificity matrix######################
#fused
#1
lSens1 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v)) #lSensitivity matrix 
rownames(lSens1) <- paste("l1",lambda1v,sep='')
colnames(lSens1) <- paste("l2",lambda2v,sep='')
lSpec1 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))
rownames(lSpec1) <- paste("l1",lambda1v,sep='')
colnames(lSpec1) <- paste("l2",lambda2v,sep='')
#2
lSens2 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v)) #lSensitivity matrix 
rownames(lSens2) <- paste("l1",lambda1v,sep='')
colnames(lSens2) <- paste("l2",lambda2v,sep='')
lSpec2 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))
rownames(lSpec2) <- paste("l1",lambda1v,sep='')
colnames(lSpec2) <- paste("l2",lambda2v,sep='')
##################Prepare Input#######################
yy = list(y1,y2)
s1 = cov(y1); s2 = cov(y2)
Sigmas=array(c(s1,s2), dim=c(nG,nG,2))
Sigma0 = array(c(omega.mat[[1]], omega.mat[[2]]), dim=c(nG,nG,2))
Omega0 = array(c(solve(omega.mat[[1]]),solve(omega.mat[[2]])), dim=c(nG,nG,2))
ns = as.vector(c(nrow(y1), nrow(y2)))
D = clust(yy, mthd=1)
CL = D2L(D)
CL = CL+diag(2)*0.001
tol=1e-4; rho=1
###################################################################################
for(m in 1:length(lambda1v)){
  for(n in 1:length(lambda2v)){
    lambda1 = lambda1v[m] 
    lambda2 = lambda2v[n] 
    start_time <- Sys.time()
    ff = try(lasich(Sigmas, Sigma0, Omega0, ns, CL, lambda1, lambda2, tol, rho))
    if(inherits(ff, "try-error")) next
    bicvl = jbic(ff, yy)
    bicl = c(bicl, bicvl)
    lam = rbind(lam, c(lambda1, lambda2))
    ##############evaluate auc by fixing similarity parameter##############
    lc1 = ff$est[,,1]; lc2 = ff$est[,,2]
    l.heatmapR = matrix(0, ncol = 2, nrow = nG*(nG - 1)/2)
    h = 0 
    for(i in 1:(nG-1)){
      for(j in (i+1):nG){
        h = h + 1
        l.heatmapR[h, 1] = sum(lc1[i,j] != 0)
        l.heatmapR[h, 2] = sum(lc2[i,j] != 0)
      }
    }
  #######################
  l1 = l.heatmapR[,1]; l2 = l.heatmapR[,2] 
  lSens1[m, n] = TPR(l1, icl[[1]]); lSpec1[m, n] = 1 - FPR(l1, icl[[1]])
  lSens2[m, n] = TPR(l2, icl[[2]]); lSpec2[m, n] = 1 - FPR(l2, icl[[2]])
  end_time <- Sys.time()
  }
}
  
leauc1 <- EAUC(lSens1, lSpec1)
leauc2 <- EAUC(lSens2, lSpec2)


ibl = which.min(bicl)
laml = lam[ibl, ]
print(bicl)
print(laml)
#
ff = lasich(Sigmas, Sigma0, Omega0, ns, CL, laml[1], laml[2], tol, rho)
lc1 = ff$est[,,1]; lc2 = ff$est[,,2]
l.heatmapR = matrix(0, ncol = 2, nrow = nG*(nG - 1)/2)
h = 0 
for(i in 1:(nG-1)){
  for(j in (i+1):nG){
    h = h + 1
    l.heatmapR[h, 1] = sum(lc1[i,j] != 0)
    l.heatmapR[h, 2] = sum(lc2[i,j] != 0)
  }
}

######################
l1 = l.heatmapR[,1]; l2 = l.heatmapR[,2]
ll <-rbind(
  c(TPR(l1, icl[[1]]), FPR(l1, icl[[1]]), leauc1),
  c(TPR(l2, icl[[2]]), FPR(l2, icl[[2]]), leauc2)
  )
colnames(ll) <- c('TPR', 'FPR','AUC')
rownames(ll) <- c('Group 1', 'Group 2')
print(ll)
save.image(file = paste0('lsi', id, '.rda'))


