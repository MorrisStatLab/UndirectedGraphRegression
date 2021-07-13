args = commandArgs(trailingOnly = T)
aaa = as.character(args[1])
library(MASS)
library(igraph)
options(scipen=999)
#
library(igraph)
library(rags2ridges)
library(matrixcalc)
set.seed(as.numeric(aaa))
dan <- function(x){
  dele <- diag(x)
  for (i in 1:nrow(x)) x[i,i] = 0
  row.sum = 1.5*apply(abs(x), 1, sum)
  x = x/row.sum
  for (i in 1:nrow(x)) x[i,i] = dele[i]
  return((x+t(x))/2)
}


gnn <- function(p, k){
  # generate the nearest network with p and k
  # p is the number of genes and k is for the k-nearest
  pts = c()
  for(ii in 1:p){
    pts = rbind(pts, c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 1)))
  }
  ptd = as.matrix(dist(pts, diag=TRUE, upper=TRUE))
  res = diag(1, nrow=p, ncol=p)
  for(ii in 1:p){
    pd = ptd[ii, ]
    ipd = setdiff(order(pd), ii)[1:k]
    for (jj in ipd){
      res[ii, jj] = res[jj, ii] = runif(1, min = 0.3, max = 0.5)*((-1)^(rbinom(1,2,prob=0.5)))
    }
  }
  return(res)
}

gsc <- function(p){
  # generate the scale-free network 
  g = barabasi.game(p)
  eg = as_edgelist(g, names = TRUE)
  res = diag(1, nrow=p, ncol=p)
  for(eid in 1:dim(eg)[1]){
      ii = eg[eid,1]; jj = eg[eid,2]
      res[ii, jj] = res[jj, ii] = runif(1, min = 0.3, max = 0.5)*((-1)^(rbinom(1,2,prob=0.5)))
  }
  return(res)
}




nG = 20;
N = 200;
Nn = 100;
#find the id of upper triangle
h = 0; hh = c()
for(l in 1:nG){
  for(m in 1:nG){
    h = h + 1
    if(l < m) hh = c(hh, h) 
  }
}
omega.mat = list()
#K1
pcov = gnn(nG, 3)
pcov = dan(pcov)
while(!(is.positive.definite(pcov)==TRUE && 
  max(abs(pcov)) <= 1 && min(abs(pcov[pcov!=0])) > 0.1)){
  pcov = gnn(nG, 3)
  pcov = dan(pcov)
}
ccov = solve(pcov)
yy1 = mvrnorm(N, mu = rep(0, nG), Sigma = ccov)
yn = mvrnorm(Nn, mu = rep(0, nG), Sigma = ccov)
omega.mat[[1]] = pcov
print((sum(pcov != 0) - nG)/2/(nG*(nG-1)/2)) # sparsity
#K2
pcov = gnn(nG, 3)
pcov = dan(pcov)
while(!(is.positive.definite(pcov)==TRUE && 
  max(abs(pcov)) <= 1 && min(abs(pcov[pcov!=0])) > 0.1)){
  pcov = gnn(nG, 3)
  pcov = dan(pcov)
}
ccov = solve(pcov)
yy2 = mvrnorm(N, mu = rep(0, nG), Sigma = ccov)
omega.mat[[2]] = pcov
print((sum(pcov != 0) - nG)/2/(nG*(nG-1)/2)) # sparsity
####################################################################################
nPi <- seq(0.01, 0.99, length.out = N)
ye = log2((2^yy1)*(1-nPi) + (2^yy2)*nPi)
yy = rbind(yn, ye)
x = c(rep(0, Nn), nPi)
xx = cbind(1-x,x)
colnames(xx)=c('x1','x2')
yy = scale(yy, center = T, scale = T)
save(yy,xx,yn,ye,omega.mat, file = 'sim.rda')
## give the edge inclusion for graphs
icl <- list()
for(i in 1:2) {
	tmp = c()
	mat = omega.mat[[i]]
	for(m in 1:(nG-1)){
		for(n in (m+1):nG){
			tmp = c(tmp, sum(mat[m, n] != 0))
		}
	}
	icl[[i]] = tmp
}
save(icl, file = 'icl.rda')

