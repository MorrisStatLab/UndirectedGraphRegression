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
lc <<- 1000
typeIerr = 0.1
sl2 = c(0.1,0.15,0.2,0.25,0.3) 
thresh <<- sl2[as.numeric(as.character(args[1]))]

lfdrf <- function(x) return(sum(abs(x) <= thresh)/lc) 
rhomn <- matrix(0, ncol = (nG-1)*nG/2, nrow = 5) 
pr.heatmapR <- matrix(0, nrow = (nG-1)*nG/2, ncol = 3) 
sigmaij = list()
sigmaij[[1]] = sigmaIJN; sigmaij[[2]] = sigmaIJT
ncount = 0
for(m in 1:(nG-1)){
  for(n in (m+1):nG){
    rho.edge = list()
    ncount = ncount + 1
    rhomn[1, ncount] = m; rhomn[2, ncount] = n 
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
    rhomn[5,ncount] <- lfdrf(rho.edge[[3]])
    pr.heatmapR[ncount, 3] = mean(rho.edge[[3]])
    }
}
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
heatmapR = pr.heatmapR; row.names(heatmapR) = 1:nrow(heatmapR)
ncount=0
for(m in 1:(nG-1)){
  for(n in (m+1):nG){
    ncount = ncount + 1
    row.names(heatmapR)[ncount] = paste(m,n,sep='_')
  }
}
ngn = row.names(heatmapR)[which(heatmapR[,1] != 0)] #normal
tgn = row.names(heatmapR)[which(heatmapR[,2] != 0)] #tumor
hgn = row.names(heatmapR)[which(heatmapR[,3] != 0)] #half


ccn = matrix(0, nrow = 3, ncol = 3)
for(i in 1:length(ngn)){
  tmp = ngn[i]
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
      ccn[va, vb] = ccn[va, vb] + 1
      if(va != vb) ccn[vb, va] = ccn[vb, va] + 1
    }
  }
}
#tumor
cct = matrix(0, nrow = 3, ncol = 3)
for(i in 1:length(tgn)){
  tmp = tgn[i]
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
      cct[va, vb] = cct[va, vb] + 1
      if(va != vb) cct[vb, va] = cct[vb, va] + 1
    }
  }
}
#half
cch = matrix(0, nrow = 3, ncol = 3)
for(i in 1:length(hgn)){
  tmp = hgn[i]
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
      cch[va, vb] = cch[va, vb] + 1
      if(va != vb) cch[vb, va] = cch[vb, va] + 1
    }
  }
}



############get the overall count of edges for each pathway and cross-pathway
l01 = length(id1); l02 = length(id3); l03 = length(id5)
l1 = l01 + length(id13); l2 = l02 + length(id13) + length(id35); l3 = l03 + length(id35)
tt = matrix(0, nrow = 3, ncol = 3)
tt[1,1] = l1*(l1-1)/2; tt[2,2] = l2*(l2-1)/2; tt[3,3] = l3*(l3-1)/2
tt[1,2] = l1*l2; tt[2,1] = l1*l2
tt[1,3] = l1*l3; tt[3,1] = l3*l1
tt[2,3] = l2*l3; tt[3,2] = l3*l2

row.names(ccn) = c("GH", "Angiogenesis", "Immune")
colnames(ccn) = c("GH", "Angiogenesis", "Immune")
row.names(cct) = c("GH", "Angiogenesis", "Immune")
colnames(cct) = c("GH", "Angiogenesis", "Immune")
row.names(cch) = c("GH", "Angiogenesis", "Immune")
colnames(cch) = c("GH", "Angiogenesis", "Immune")
print(thresh)
cctb = cbind(ccn, cch, cct)
print(xtable(cctb))
tt = cbind(tt,tt,tt)
print(xtable(cctb/tt))





