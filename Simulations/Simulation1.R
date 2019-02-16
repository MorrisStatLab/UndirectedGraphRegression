#########simulate from our model to test MCMC
args = commandArgs(trailingOnly = T)
aaa = as.character(args[1])
library(MASS)
options(scipen=999)
#
set.seed(as.numeric(aaa))
nG = 20;
N = 150;
Nn = 50;
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
pcov = diag(rep(1, nG))
for(i in 1:(nG - 1)){
	for(j in (i+1):nG){
	tmp = 0
	if(j == i + 1) tmp = runif(1, min = 0.3, max = 0.5)*((-1)^(rbinom(1,2,prob=0.5)))
	pcov[i,j] = tmp
	pcov[j,i] = tmp
	}
}
ccov = solve(pcov)
yy1 = mvrnorm(N, mu = rep(0, nG), Sigma = ccov)
yn = mvrnorm(Nn, mu = rep(0, nG), Sigma = ccov)
omega.mat[[1]] = pcov
print((sum(pcov != 0) - nG)/2/(nG*(nG-1)/2)) # sparsity
#K2
pcov = diag(rep(1, nG))
for(i in 1:(nG - 1)){
	for(j in (i+1):nG){
		tmp = 0
		if(j == i + 2) tmp = runif(1, min = 0.3, max = 0.5)*((-1)^(rbinom(1,2,prob=0.5)))
		pcov[i,j] = tmp
		pcov[j,i] = tmp
	}
}
ccov = solve(pcov)
yy2 = mvrnorm(N, mu = rep(0, nG), Sigma = ccov)
omega.mat[[2]] = pcov
print((sum(pcov != 0) - nG)/2/(nG*(nG-1)/2)) # sparsity
###print test result
intersect(which(omega.mat[[1]] != 0), hh)
intersect(which(omega.mat[[2]] != 0), hh)
####################################################################################
nPi <- seq(0.01, 0.99, length.out = N)
ye = log2((2^yy1)*(1-nPi) + (2^yy2)*nPi)
yy = rbind(yn, ye)
x = c(rep(0, Nn), nPi)
xx = cbind(1-x,x)
colnames(xx)=c('x1','x2')
########summarize the final data
yy = scale(yy, center = T, scale = T)
save(yy,xx,file = paste0('sim',aaa,'.rda'))
save.image(paste0('all',aaa,'.rda'))
##first give the edge inclusion for four graphs
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
save(icl, file = paste0('icl',aaa,'.rda'))

