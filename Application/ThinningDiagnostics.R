library(coda)
pv <- function(x) pnorm(abs(x), lower.tail=FALSE)*2
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
theta = cbind(sigmaIJN, sigmaIJT, sigmaII, psiIJN, psiIJT, lam, gamma)
theta = mcmc(theta)
mh.list = mcmc.list(theta)
testr = geweke.diag(mh.list)
testp = sapply(testr[[1]][[1]], pv)
hist(testp, main = 'Histogram of p-value for Geweke\' convergence diagnostic', xlab = 'p-value')
use.id <- seq(10, 10000, length.out = 1000)
sigmaIJN = sigmaIJN[use.id, ]
sigmaIJT = sigmaIJT[use.id, ]
sigmaII = sigmaII[use.id, ]
save(sigmaIJN, sigmaIJT, sigmaII, file = 'postm.rda')
