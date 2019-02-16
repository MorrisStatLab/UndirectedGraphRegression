#########simulate from our model to test MCMC
args = commandArgs(trailingOnly = T)
aaa = as.character(args[1])
library(MASS)
options(scipen=999)
#write a function to do a positive definite transformation
dan <- function(x){
  dele <- diag(x)
  for (i in 1:nrow(x)) x[i,i] = 0
  row.sum = 1.5*apply(abs(x), 1, sum)
  x = x/row.sum
  for (i in 1:nrow(x)) x[i,i] = dele[i]
  return((x+t(x))/2)
}
#
set.seed(as.numeric(aaa))
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
pcov = diag(rep(1, nG))
for(i in 1:(nG - 1)){
	for(j in (i+1):nG){
	tmp = 0
	if(j == i + 1) tmp = 0.5
	if(j == i + 2) tmp = 0.4
	pcov[i,j] = tmp
	pcov[j,i] = tmp
	}
}
ccov = solve(pcov)
yy2 = mvrnorm(N, mu = rep(0, nG), Sigma = ccov)
omega.mat[[2]] = pcov
print((sum(pcov != 0) - nG)/2/(nG*(nG-1)/2)) # sparsity
#K2
nzero.id = intersect(which(pcov != 0), hh)
zero.id = intersect(which(pcov == 0), hh)
rnz.id = sample(nzero.id, size = 30, replace = F)
rz.id = sample(zero.id, size = 30, replace = F)
for(i in 1:30){
  iid = (rnz.id[i]-1)%/%nG+1 #find row id 
  jid = (rnz.id[i]-1)%%nG+1 # find column id
  pcov[iid, jid] = 0; pcov[jid, iid] = 0
  iid = (rz.id[i]-1)%/%nG+1 #find row id 
  jid = (rz.id[i]-1)%%nG+1 # find column id
  tmp = runif(1, min = 0.4, max = 0.6)*((-1)^(rbinom(1,2,prob=0.5)))
  pcov[iid, jid] = tmp; pcov[jid, iid] = tmp
}
pcov = dan(pcov)
Mm = (sum(sum(pcov^2)) - sum(diag(pcov^2)))/2
Mm = Mm/(nG*(nG-1)/2); print(Mm)
ccov = solve(pcov)
yy1 = mvrnorm(N, mu = rep(0, nG), Sigma = ccov)
yn = mvrnorm(Nn, mu = rep(0, nG), Sigma = ccov)
omega.mat[[1]] = pcov
print((sum(pcov != 0) - nG)/2/(nG*(nG-1)/2)) # sparsity
###print test result
edge.id1 = intersect(which(omega.mat[[1]] != 0), hh)
edge.id2 = intersect(which(omega.mat[[2]] != 0), hh)
print(length(intersect(edge.id1, edge.id2))/37)
####################################################################################
nPi <- seq(0.01, 0.99, length.out = N)
ye = log2((2^yy1)*(1-nPi) + (2^yy2)*nPi)
yy = rbind(yn, ye)
x = c(rep(0, Nn), nPi)
xx = cbind(1-x,x)
colnames(xx)=c('x1','x2')
########summarize the final data
yy = scale(yy, center = T, scale = T)
save(yy,xx, file = paste0('sim',aaa,'.rda'))
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

