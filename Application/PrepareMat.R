library('MASS')
library('GIGrvg')
library('glasso')
ff <- function(m) class(try(solve(m),silent=T))=="matrix" #judge if it is solvable
Mout <- function(yy, x){
	p = ncol(yy)
	yy.tumor <- yy[which(x >= 0.5), ]
	yy.normal <- yy[which(x < 0.5), ]
	n.t = nrow(yy.tumor)
	n.n = nrow(yy.normal)
	if(!ff(cov(yy.tumor))){
	  prec.tumor <- glasso(cov(yy.tumor), rho = 0.01)$wi
	  prec.tumor <- prec.tumor^2
	}else{
	  prec.tumor <- (solve(cov(yy.tumor)*(n.t - 1)/n.t))^2
	}
	if(!ff(cov(yy.normal))){
	  prec.normal <- glasso(cov(yy.normal), rho = 0.01)$wi
	  prec.normal <- prec.normal^2
	}else{
	  prec.normal <- (solve(cov(yy.normal)*(n.n - 1)/n.n))^2
	}
	M.t <- (sum(sum(prec.tumor)) - sum(diag(prec.tumor)))/2
	M.n <- (sum(sum(prec.normal)) - sum(diag(prec.normal)))/2
	M <- c(M.n, M.t)/(p*(p-1)/2)
	return(M)
}
##load Pathway data
##an union of immune, growth hormone, and angiogenesis pathways plus alpha fetaprotein (AFP)
load('Pathway.rda')
load('Hcc.profile.rda')
##load estimated HepatoScore for all the samples
load('HepatoScore.rda')
yy = t(hcc.pro[gene.id, ])
xx = cbind(1 - score, score)
x = as.numeric(xx[,2])
##calculate the hyperparameter Ms 
Mm <- Mout(yy, x)
library(R.matlab)
filename = 'Profile.mat'
writeMat(filename, xx = xx, yy = yy, Mm = Mm)


