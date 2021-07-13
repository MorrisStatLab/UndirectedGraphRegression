args = commandArgs(trailingOnly = T)
library('coda')
library(xtable)
load('Pathway.rda')

id1 = which((pathu %in% path1) & !(pathu %in% path3) & !(pathu %in% path5))
id3 = which((pathu %in% path3) & !(pathu %in% path1) & !(pathu %in% path5))
id5 = which((pathu %in% path5) & !(pathu %in% path1) & !(pathu %in% path3))
id13 = which((pathu %in% path1) & (pathu %in% path3))
id35 = which((pathu %in% path3) & (pathu %in% path5))

load('postm.rda')
nG = ncol(sigmaII) 
##load the data from result
lc <<- 1000
typeIerr = 0.1
#threshold
sl2 = c(0.1,0.15,0.2,0.25,0.3) 
thresh <<- sl2[as.numeric(as.character(args[1]))]
#function for local fdr
pur = seq(0, 1, by = 0.05)
lfdrf <- function(x) return(sum(abs(x) <= thresh)/lc) #just calculate the probability that below thresh
LL = length(pur)
rhomn <- matrix(0, ncol = (nG-1)*nG/2, nrow = 2+LL)

pr.heatmapR <- matrix(0, nrow = (nG-1)*nG/2, ncol = LL) #posterior mean
ncount = 0
for(m in 1:(nG-1)){
  for(n in (m+1):nG){
    ncount = ncount + 1
    ###calculate the rho
    rhomn[1, ncount] = m; rhomn[2, ncount] = n 
    ###calculate the rho.edge
    for(i in 1:LL) {
        tmp1 = as.numeric(sigmaIJN[, ncount])
        tmp2 = as.numeric(sigmaIJT[, ncount])
        tmp = -(tmp1*(1-pur[i]) + tmp2*pur[i])/sqrt(as.numeric(sigmaII[,m])*as.numeric(sigmaII[,n]))
        rhomn[i+2,ncount] <- lfdrf(tmp)
        pr.heatmapR[ncount, i] = mean(tmp)
     }
  }
}
########Global FDR to select edges kept for graph
for(i in 1:LL){
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
heatmapR = pr.heatmapR; row.names(heatmapR) = 1:nrow(heatmapR)
ncount=0
for(m in 1:(nG-1)){
  for(n in (m+1):nG){
    ncount = ncount + 1
    row.names(heatmapR)[ncount] = paste(m,n,sep='_')
  }
}


res = matrix(0, nrow = 6, ncol = LL)
llc = c(11,12,13,22,23,33)
#take an iteration to count the edges within/cross subpathways for all different purity
for(l in 1:LL){
ngl = row.names(heatmapR)[which(heatmapR[,l] != 0)]
## GH, Angiogenesis, Immune
for(i in 1:length(ngl)){
  tmp = ngl[i]
  tmp = strsplit(tmp, split = '_')[[1]]
  a = as.numeric(tmp[1]); b = as.numeric(tmp[2])
  if((a == 2)||(b == 2)) next;
  a = 1*(a %in% id1) + 2*(a %in% id3) + 3*(a %in% id5) + 12*(a %in% id13) + 23*(a %in% id35)
  b = 1*(b %in% id1) + 2*(b %in% id3) + 3*(b %in% id5) + 12*(b %in% id13) + 23*(b %in% id35)
  seta = c(a%/%10, a%%10); setb = c(b%/%10, b%%10);
  seta = seta[seta != 0]; setb = setb[setb != 0]
  for(ia in 1:length(seta)){
    for(ib in 1:length(setb)){
      va = seta[ia]; vb = setb[ib]
      vc = min(va,vb); vb = max(va,vb); va = vc
      vv = va*10 + vb
      nnc = which(llc == vv); res[nnc, l] = res[nnc, l] + 1
    }
  }
   }
}

############get the overall count of edges for each pathway and cross-pathway
l01 = length(id1); l02 = length(id3); l03 = length(id5)
l1 = l01 + length(id13); l2 = l02 + length(id13) + length(id35); l3 = l03 + length(id35) 
tt = rep(0,6)

tt[which(llc == 11)] = l1*(l1-1)/2
tt[which(llc == 22)] = l2*(l2-1)/2
tt[which(llc == 33)] = l3*(l3-1)/2
tt[which(llc == 12)] = l1*l2
tt[which(llc == 13)] = l1*l3
tt[which(llc == 23)] = l2*l3
rest = res/tt
row.names(res) = c("GH","GH/Angiogenesis", "GH/Immune", "Angiogenesis", "Angiogenesis/Immune", "Immune")
row.names(rest) = c("GH","GH/Angiogenesis", "GH/Immune", "Angiogenesis", "Angiogenesis/Immune", "Immune")
colnames(res) = pur
colnames(rest) = pur
print(res)
print(rest)

res = as.matrix(res)
rest = as.matrix(rest)

library(gplots)
source('heatmap.3.R')
heatmap.3(x = rest, Rowv = FALSE, Colv = FALSE, dendrogram = "none", cellnote = res, notecol = "black",
cexRow = 2.5, cexCol = 3, col = bluered,key.title = NA, key.ylab = NA, key.xlab = 'CirBAS', keysize = 1.1, key.par = list(cex.lab=2.1, cex.axis=2.2, mar = c(3,4,3,3)),
notecex = 3.5,trace = "none", key = TRUE, margins = c(5, 12))



