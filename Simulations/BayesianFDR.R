args = commandArgs(trailingOnly = T)
# 1: type 1 error; 2: threshold for edge; 
library('coda')
id = as.character(args[3])
load(paste0('post',id,'.rda'))
nG = 20 ## protein number
##load the data from result
lc <<- 1000
ncount = 0
#type I error = 0.9
sl1 = c(0.05, 0.10, 0.20) # three type I errors
typeI = as.numeric(as.character(args[1]))
typeIerr = sl1[typeI]
#threshold
sl2 = c(0.05, 0.1, 0.2)
thresh <<- sl2[as.numeric(as.character(args[2]))]
#function for local fdr
lfdrf <- function(x) return(sum(abs(x) <= thresh)/lc) #just calculate the probability that below thresh

###local fdr function for differential and common graph
lfdrfc <- function(x, y) return(1 - sum( (abs(x) > thresh) & (abs(y) > thresh))/lc) #common graph
lfdrfd <- function(x, y) return(1 - sum( (abs(x) > thresh) & (abs(y) <= thresh) )/lc) #differential graph with x from y

##graph1: normal; graph2: tumor
##set up a matrix to store everything
rhomn <- matrix(0, ncol = (nG-1)*nG/2, nrow = 4) #m; n; graph1, graph2 plus two columns for index of edge; store probability of rhomn
pr.heatmapR <- matrix(0, nrow = (nG-1)*nG/2, ncol = 2) #posterior mean of two graphs
#set the differential
rhodf <- rep(0, (nG-1)*nG/2) # 1 differential graph
rhocm <- rep(0, (nG-1)*nG/2) # 1 common graph

#####first calculate the rho of edge from sigmaij
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
      #common graph
	  rhocm[ncount] <- lfdrfc(rho.edge[[1]], rho.edge[[2]]) #common

      #differential graph
	  rhodf[ncount] <- lfdrfd(rho.edge[[1]], rho.edge[[2]]) 
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

####################################################This is a global FDR strategy for analysis of differential and common edges############
#####################################################This controls a common and differential graph among four groups##################################
df.heatmapR = rep(1, nrow = nG*(nG - 1)/2)
	lfdr <- rhodf
	lfdr.s <- sort(lfdr, decreasing = F)
	lfdr.o <- order(lfdr, decreasing = F)
	k = 1
	tmp.lfdr = lfdr.s[k]  
	while(tmp.lfdr/k <= typeIerr){
		k = k + 1
		tmp.lfdr = tmp.lfdr + lfdr.s[k]
	}
	k = k - 1
	id.sl = lfdr.o[1:k]
	df.heatmapR[-id.sl] = 0
#####################################################This controls a common and differential graph among four groups##################################
cm.heatmapR = rep(1, nG*(nG - 1)/2)
	lfdr <- rhocm
	lfdr.s <- sort(lfdr, decreasing = F)
	lfdr.o <- order(lfdr, decreasing = F)
	k = 1
	tmp.lfdr = lfdr.s[k]  
	while(tmp.lfdr/k <= typeIerr){
		k = k + 1
		tmp.lfdr = tmp.lfdr + lfdr.s[k]
	}
	k = k - 1
	id.sl = lfdr.o[1:k]
	cm.heatmapR[-id.sl] = 0

#######################################################################################################################################################
save(pr.heatmapR, df.heatmapR, cm.heatmapR, rhomn, rhocm, rhodf, file = paste0('fdr',id, '_', as.character(args[1]), '_', as.character(args[2]), '.rda'))








