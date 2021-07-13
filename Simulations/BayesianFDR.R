args = commandArgs(trailingOnly = T)
# 1: type 1 error; 2: threshold for edge; 
library('coda')
id = as.character(args[3])
load(paste0('post',id,'.rda'))
nG = 20 
##load the data from result
lc <<- 1000
ncount = 0
sl1 = c(0.05, 0.10, 0.20) 
typeI = as.numeric(as.character(args[1]))
typeIerr = sl1[typeI]
#threshold
sl2 = c(0.05, 0.1, 0.2)
thresh <<- sl2[as.numeric(as.character(args[2]))]
#function for local fdr
lfdrf <- function(x) return(sum(abs(x) <= thresh)/lc) 
##graph1: normal; graph2: tumor
##set up a matrix to store everything
rhomn <- matrix(0, ncol = (nG-1)*nG/2, nrow = 4) 
pr.heatmapR <- matrix(0, nrow = (nG-1)*nG/2, ncol = 2) #posterior mean of two graphs


sigmaij = list()
sigmaij[[1]] = sigmaIJ1
sigmaij[[2]] = sigmaIJ2
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
  }
}

########Global FDR to select edges kept for graph
for(i in 1:2){
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
save(pr.heatmapR, rhomn, file = paste0('fdr',id, '_', as.character(args[1]), '_', as.character(args[2]), '.rda'))







