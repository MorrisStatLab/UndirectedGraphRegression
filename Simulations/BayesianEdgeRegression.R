library('MASS')
library('GIGrvg')
##load data from sim_mgg rda
run.id = commandArgs(trailingOnly = T)
run.id = as.numeric(run.id[1])
load(paste0('sim', run.id,'.rda'))
x = as.numeric(xx[,2])
##############################################
# a function to generate Mm
source('M.R')
Mm <- Mout(yy, x)
print(Mm)
############prepare and intialize##############
nG = ncol(yy); N = nrow(yy); nX = ncol(xx)
accept_iter <- rep(0, nX)
idsigmaIJ = matrix(0, nrow = nG, ncol = nG)
h = 0;
for(i in 1:(nG-1)){
  for(j in (i+1):nG){
    h = h + 1; idsigmaIJ[i,j] = h; idsigmaIJ[j,i] = h
  }
}
##initialize
sigmaIJ = matrix(0, nrow = nX, ncol = nG*(nG - 1)/2)
sigmaII = rep(1, nG)
psiIJ = matrix(0.25, nrow = nX, ncol = nG*(nG - 1)/2)
lambdam = rep(1, nX)
#tune parameter
sigma.lambda = c(0.2, 0.25)  
gammaSqrm = 0.5
niteration = 20000 
##finish initialize


#input
##normal gamma to update sigma_ij, sigma_ji as a pair
##sigmaIJ: m* (p*(p-1))/2 matrix; precision off-diagonal element to estimate; initialization begins with 0
##sigmaII: p vector; precision on-diagonal element to estimat; initialization with 1
##yy: n*p matrix of expression; n: sample size; p: node size
##xx: n*m matrix of exogenous variable; m: number of exogenous variable
##psiIJ: m* (p*(p-1))/2 matrix; row size is for m variables; column size count IJ pair: 12, 13, ..., (p-1)p 
##lambdam: m-vector: each global hyper-parameter lambda for psiIJ
##gammaSqrm: m-vector: each global hyper-parameter gamma square for psiIJ
##Mm: m-vector; sum squared effect size of sigamIJ from calculation of least squares
##idsigmaIJ: p*p matrix; an indicator matrix to indicate index of IJ corresponding to column id of sigmaIJ
##sigma.lambda: m-vector; set value for random walk of lambda
### N: sample size; nG: node size; nX: exo-covariate size

#output
##update sigmaIJ; psiIJ; lambdam; gammaSqrm
##rhoIJ: updated partial correlation after calculation from sigmaIJ

#a list to store result
res <- list()

#begin iteration
for(t in 1:niteration){
print(t)
##begin with gibbs sampling
#############################update sigmaIJ(i<j) for each pair (i,j) and (j,i)
for(i in 1:(nG-1)){
	for(j in (i+1):nG){
#gibbs sampler from normal
		h = idsigmaIJ[i, j]
	  psiIJm = solve(diag(psiIJ[,h])) ## m*m matrix
#calculate S1: which is a n*n matrix
		S1 = diag(((yy[,j]^2)/sigmaII[i] + (yy[,i]^2)/sigmaII[j]))
#calculate S2: which is a n*n matrix
		idi = idsigmaIJ[i, -c(i,j)]
		idj = idsigmaIJ[j, -c(i,j)]
		S2 = 2*yy[,i]*yy[,j]
		tmp = diag((xx%*%sigmaIJ[,idi])%*%t(yy[,-c(i,j)]))
	  S2 = S2 + yy[,j]*tmp/sigmaII[i]
	  tmp = diag((xx%*%sigmaIJ[,idj])%*%t(yy[,-c(i,j)]))
		S2 = S2 + yy[,i]*tmp/sigmaII[j]
#sampling parameter
	  S2 = as.matrix(S2)
		Normal.var <- solve(psiIJm + t(xx)%*%S1%*%(xx))
		Normal.mu <- -Normal.var%*%t(xx)%*%S2
#sample
		sigmaIJ[,h] <- mvrnorm(n = 1, mu = Normal.mu, Sigma = Normal.var) #
	}
}

#############################update sigmaII for each partial correlation
print('GIG sigmaII')
gig1.lambda = N/2 + 1
for(i in 1:nG){
	gig1.psi = sum(yy[,i]^2)
	idg = idsigmaIJ[i,-i]
    gig1.chi = sum(diag(xx%*%sigmaIJ[,idg]%*%t(yy[,-i]))^2)
	sigmaII[i] <- rgig(1, chi = gig1.chi, lambda = gig1.lambda, psi = gig1.psi)
	print(c(gig1.lambda, gig1.psi, gig1.chi))
}
			
	
	
#check precision
lambdam.l <- which(lambdam <= 0.5); lln = length(lambdam.l)
if(lln > 0){
	for(l in 1:lln){
		ll = lambdam.l[l]
		tmp = sigmaIJ[ll, ]
		tmp[abs(tmp) <= sqrt(.Machine$double.eps ^ 0.5)] = sign(tmp[abs(tmp) <= sqrt(.Machine$double.eps ^ 0.5)]) * sqrt((.Machine$double.eps ^ 0.5)*2)
		sigmaIJ[ll,] = tmp		  
		}
}	
#############################update psiIJ for all sigmaIJ local shrinkage parameter
print('GIG psiIJ')
for(m in 1:nX){
	gig2.lambda = lambdam[m] - 1/2
    gig2.psi = 1/gammaSqrm
   for(i in 1:(nG-1)){
	  for(j in (i+1):nG){
	      h = idsigmaIJ[i, j]
          gig2.chi = sigmaIJ[m,h]^2
#print(c(gig2.chi, gig2.lambda));
          temp = rgig(1, chi = gig2.chi, lambda = gig2.lambda, psi = gig2.psi)
          ### check precision for psi generated from gig
          if(temp < .Machine$double.eps ^ 0.5) temp = (.Machine$double.eps ^ 0.5)*2
          psiIJ[m, h] = temp
	}
 }
}
  
#############################update lambda and gamma for hyper-parameter above the psiIJ
p = nG*(nG-1)/2
#update lambdam
z.lambda = rnorm(nX, mean = 0, sd = 1)
lambda.pv = exp(sigma.lambda*z.lambda)*lambdam #proposed value
accept = rep(0, nX)
for(l in 1:nX){
  logR.lambda = log(dexp(lambda.pv[l], rate = 1)) - log(dexp(lambdam[l], rate = 1)) + 
    p*(log(gamma(lambdam[l])) - log(gamma(lambda.pv[l]))) + 
    (lambda.pv[l] - lambdam[l])*(-p*log(2*gammaSqrm) + sum(log(psiIJ[l,])))
  logR.lambda = logR.lambda + log(lambda.pv[l]) - log(lambdam[l])
  R.lambda = exp(logR.lambda)
  u = runif(1)
  if(u < R.lambda){
    lambdam[l] = lambda.pv[l]; accept[l] =  1
  }else{
    accept[l] =  0 
  }
}
accept_iter <- accept_iter + accept
#update gammaSqrm
e.star.gamma = 2 + p*sum(lambdam)
f.star.gamma = sum(Mm)/(2*sum(lambdam)) + sum(psiIJ)/2
gammaSqrm = 1/rgamma(1, shape = e.star.gamma, rate = f.star.gamma)

res[[t]] = list(sigmaIJ = sigmaIJ, sigmaII = sigmaII, accept = accept_iter, 
                psiIJ = psiIJ, lambdam = lambdam, gammaSqrm = gammaSqrm)
}
				  
save(res, file = paste0('res',run.id,'.rda') )
