EAUC <- function(Sens, Spec){
  phi_grid <- seq(0,1,by=.05)  # for x axis of ROC curve
  # ROC with expected sensitivity
  ESens <- rep(NA,length(phi_grid)-1)
  for (phi in 2:length(phi_grid)){
    temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
    ESens[phi-1] <- mean(Sens[temp])
  }
  # interpolate if missing ESens value
  if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
  for (ind in which(is.na(ESens))) {
    low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
    high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
    slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
    int <- ESens[high] - slope*(phi_grid[high]+0.025)
    ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
  }
  # Calcuate "Expected AUC"
  return(0.05*sum(ESens))
}
#########################function########################
lc <<- 1000
TPR <- function(fit, truth){ 
	tp <- sum((fit == 1)&(truth == 1))
	fn <- sum((fit == 0)&(truth == 1))
	return(tp/(tp+fn))
}

FPR <- function(fit, truth){ 
	fp <- sum((fit == 1)&(truth == 0))
	tn <- sum((fit == 0)&(truth == 0))
	return(fp/(fp+tn))
}
lfdrf <- function(x, thresh) return(sum(abs(x) <= thresh)/lc) 
#########################################################
args = commandArgs(trailingOnly = T)
args = as.character(args)
load(paste0('icl',args[1],'.rda'))
load(paste0('post',args[1],'.rda'))
nG = 20 
sigmaij = list()
sigmaij[[1]] = sigmaIJ1;
sigmaij[[2]] = sigmaIJ2; 
tyerr = seq(0, 1, by = 0.01)
thval = seq(0, 1, by = 0.01)
#1
Sens1 <- matrix(NA,nrow=length(tyerr),ncol=length(thval)) 
rownames(Sens1) <- paste("typeI",tyerr,sep='')
colnames(Sens1) <- paste("thresh",thval,sep='')
Spec1 <- matrix(NA,nrow=length(tyerr),ncol=length(thval))
rownames(Spec1) <- paste("typeI",tyerr,sep='')
colnames(Spec1) <- paste("thresh",thval,sep='')
#2
Sens2 <- matrix(NA,nrow=length(tyerr),ncol=length(thval))  
rownames(Sens2) <- paste("typeI",tyerr,sep='')
colnames(Sens2) <- paste("thresh",thval,sep='')
Spec2 <- matrix(NA,nrow=length(tyerr),ncol=length(thval))
rownames(Spec2) <- paste("typeI",tyerr,sep='')
colnames(Spec2) <- paste("thresh",thval,sep='')

for (a in 1:length(tyerr)){
	for(b in 1:length(thval)){
	  typeIerr <- tyerr[a]; thresh <- thval[b]
      rhomn <- matrix(0, ncol = (nG-1)*nG/2, nrow = 4) 
	  pr.heatmapR <- matrix(1, nrow = (nG-1)*nG/2, ncol = 2) 
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
			rhomn[i+2, ncount] <- lfdrf(rho.edge[[i]], thresh)
		} 
	}
}
#############Global FDR to select edges kept for graph###########################
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
		
	pr1 = pr.heatmapR[,1]; pr2 = pr.heatmapR[,2]

Sens1[a,b] = TPR(pr1, icl[[1]]); Spec1[a,b] = 1 - FPR(pr1, icl[[1]])
Sens2[a,b] = TPR(pr2, icl[[2]]); Spec2[a,b] = 1 - FPR(pr2, icl[[2]])
	}
}
eauc1 <- EAUC(Sens1, Spec1)
eauc2 <- EAUC(Sens2, Spec2)
eauc <- c(eauc1, eauc2)
save(eauc,Sens1, Spec1, Sens2, Spec2,file = paste0('auc',args[1],'.rda'))




