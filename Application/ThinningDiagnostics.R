library(coda)
pv <- function(x) min(pnorm(x), 1- pnorm(x))
# this is used for thinning the chain 
# and check for convergence
## Convert Matlab to R
library(R.matlab)
filename = 'reswn.mat'
input = readMat(filename)
sigmaIJN = input$sigmaIJN
sigmaIJT = input$sigmaIJT
psiIJN = input$psiIJN
psiIJT = input$psiIJT
sigmaII = input$sigmaII
lam = input$lambdam
gamma = input$gammaSqrm
nG = ncol(sigmaII)
########################Check#################
theta = cbind(sigmaIJN, sigmaIJT, sigmaII, psiIJN, psiIJT, lam, gamma)
theta = mcmc(theta)
mh.list = mcmc.list(theta)
testr = geweke.diag(mh.list)
testp = sapply(testr[[1]][[1]], pv)
png('hist.png')
hist(testp, main = 'Histogram of p-value for Geweke\' convergence diagnostic', xlab = 'p-value')
dev.off()
###################plot for autocorrelation#########################
acf(sigmaIJN[,1])
acf(sigmaIJT[,1])
acf(sigmaII[,1])
pacf(sigmaIJN[,1])
pacf(sigmaIJT[,1])
pacf(sigmaII[,1])
i=1;j=1;k=1
plot(sigmaIJT[,i], type = 'l'); i=i+1
plot(sigmaIJN[,j], type = 'l'); j=j+1
plot(sigmaII[,k], type = 'l'); k=k+1
plot(lam[,1], type = 'l')
plot(lam[,2], type = 'l')
plot(gamma[,1], type = 'l')
##########################thinning###################################
use.id <- seq(10, 10000, length.out = 1000)
sigmaIJN = sigmaIJN[use.id, ]
sigmaIJT = sigmaIJT[use.id, ]
sigmaII = sigmaII[use.id, ]
save(sigmaIJN, sigmaIJT, sigmaII, file = 'postm.rda')


