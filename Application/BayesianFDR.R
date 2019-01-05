## normal and tumor and 0.5 tumor graph
args = commandArgs(trailingOnly = T)
# 1: type 1 error; 2: threshold for edge; 
library('coda')
load('postm.rda')
nG = ncol(sigmaII) ## protein number
##load the data from result
lc <<- 1000
ncount = 0
#type I error = 0.9
sl1 = c(0.05, 0.10, 0.20) # three type I errors
typeI = as.numeric(as.character(args[1]))
typeIerr = sl1[typeI]
#threshold
sl2 = c(0.05, 0.1, 0.15, 0.2) #this order is switched just from this version
thresh <<- sl2[as.numeric(as.character(args[2]))]
#function for local fdr
lfdrf <- function(x) return(sum(abs(x) <= thresh)/lc) #just calculate the probability that below thresh
##set up a matrix to store everything
rhomn <- matrix(0, ncol = (nG-1)*nG/2, nrow = 5) #m; n; graph1, graph2, graph3, graph4 plut two columns for index of edge; store probability of rhomn
pr.heatmapR <- matrix(0, nrow = (nG-1)*nG/2, ncol = 3) #posterior mean of three graphs: normal tumor half
#####first calculate the rho of edge from sigmaij
sigmaij = list()
sigmaij[[1]] = sigmaIJN; sigmaij[[2]] = sigmaIJT

##############################################
for(m in 1:(nG-1)){
  for(n in (m+1):nG){
	 rho.edge = list()
    ncount = ncount + 1
    ###calculate the rho
	  rhomn[1, ncount] = m; rhomn[2, ncount] = n 
    ###calculate the rho.edge
    for(i in 1:2) {
	  tmp.mat = sigmaij[[i]]
	  tmp = as.numeric(tmp.mat[, ncount])
	  tmp = -tmp/sqrt(as.numeric(sigmaII[,m])*as.numeric(sigmaII[,n]))
	  rho.edge[[i]] = tmp
      rhomn[i+2, ncount] <- lfdrf(rho.edge[[i]])
	  pr.heatmapR[ncount, i] = mean(rho.edge[[i]])
   }
    tmp1.mat = sigmaij[[1]]; tmp2.mat = sigmaij[[2]]
    tmp1 = as.numeric(tmp1.mat[, ncount])
    tmp2 = as.numeric(tmp2.mat[, ncount])
    tmp = -0.5*(tmp1+tmp2)/sqrt(as.numeric(sigmaII[,m])*as.numeric(sigmaII[,n]))
    rho.edge[[3]] = tmp
    rhomn[5, ncount] <- lfdrf(rho.edge[[3]])
    pr.heatmapR[ncount, 3] = mean(rho.edge[[3]])
  }
}

########Global FDR to select edges kept for graph
for(i in 1:3){
  lfdr <- rhomn[i+2, ]
  lfdr.s <- sort(lfdr, decreasing = F)
  lfdr.o <- order(lfdr, decreasing = F)
  k = 1
  tmp.lfdr = lfdr.s[k]  
  while(tmp.lfdr/k <= typeIerr){
    k = k + 1
	  if(k > length(lfdr.s)) break
    tmp.lfdr = tmp.lfdr + lfdr.s[k]
  }
  k = k - 1
  id.sl = lfdr.o[1:k]
  pr.heatmapR[-id.sl, i] = 0
}

for(i in 1:3) print(sum(pr.heatmapR[,i] != 0))
for(i in 1:3) print(sum(pr.heatmapR[,i] != 0)/(nG*(nG-1)/2))

save(pr.heatmapR, rhomn, file = paste0('fdrm', '_',as.character(args[1]), '_', as.character(args[2]), '.rda'))








