library(JGL)
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
#######a function to calculate the approximate AIC of fgl objects
jaic <- function(fglo, yy){
  K = length(yy)
  aicv = 0
  for(k in 1:K){
    htheta <- fglo$theta[[k]]
    xx <- yy[[k]]
    n <- nrow(xx); p <- ncol(xx)
    S <- cov(xx)*(n-1)/n
    Ek <- (sum(htheta != 0) - p)/2
    tr <- sum(diag(S%*%htheta))
    tmp = n*tr - n*log(det(htheta)) + 2*Ek
    aicv = aicv + tmp
  }
  return(aicv)
}
#evaluation function
TPR <- function(fit, truth){ #given truth = 1, fit = 1
  tp <- sum((fit == 1)&(truth == 1))
  fn <- sum((fit == 0)&(truth == 1))
  return(tp/(tp+fn))
}
FPR <- function(fit, truth){ #given truth = 0, fit = 1
  fp <- sum((fit == 1)&(truth == 0))
  tn <- sum((fit == 0)&(truth == 0))
  return(fp/(fp+tn))
}

##########################################################################
id = commandArgs(trailingOnly = T)
id = as.character(id[1])
load(paste0('sim',id,'.rda'))
load(paste0('icl',id,'.rda'))
y1 = scale(yn, center=T, scale=T)
y2 = scale(ye, center=T, scale=T)
# lambda1v; lambda2v
nG = 20
lambda1v = seq(0.01, 2.01, length.out = 101)
lambda2v = seq(0.01, 2.01, length.out = 101)
aicf = c(); aicg = c()
lam = c()
auc1f = c(); auc2f = c();
auc1g = c(); auc2g = c();
###################define fSensitivity and fSpecificity matrix######################
#fused
#1
fSens1 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v)) 
rownames(fSens1) <- paste("l1",lambda1v,sep='')
colnames(fSens1) <- paste("l2",lambda2v,sep='')
fSpec1 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))
rownames(fSpec1) <- paste("l1",lambda1v,sep='')
colnames(fSpec1) <- paste("l2",lambda2v,sep='')
#2
fSens2 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v)) 
rownames(fSens2) <- paste("l1",lambda1v,sep='')
colnames(fSens2) <- paste("l2",lambda2v,sep='')
fSpec2 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))
rownames(fSpec2) <- paste("l1",lambda1v,sep='')
colnames(fSpec2) <- paste("l2",lambda2v,sep='')

#group
###################define gSensitivity and gSpecificity matrix######################
#
#1
gSens1 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))  
rownames(gSens1) <- paste("l1",lambda1v,sep='')
colnames(gSens1) <- paste("l2",lambda2v,sep='')
gSpec1 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))
rownames(gSpec1) <- paste("l1",lambda1v,sep='')
colnames(gSpec1) <- paste("l2",lambda2v,sep='')
#2
gSens2 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v)) 
rownames(gSens2) <- paste("l1",lambda1v,sep='')
colnames(gSens2) <- paste("l2",lambda2v,sep='')
gSpec2 <- matrix(NA,nrow=length(lambda1v),ncol=length(lambda2v))
rownames(gSpec2) <- paste("l1",lambda1v,sep='')
colnames(gSpec2) <- paste("l2",lambda2v,sep='')

###################################################################################
for(m in 1:length(lambda1v)){
for(n in 1:length(lambda2v)){
  lambda1 = lambda1v[m] 
    lambda2 = lambda2v[n] 
    yy = list(y1,y2)
    fglff = JGL(Y = yy, penalty = "fused", lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = T)
    fglfg = JGL(Y = yy, penalty = "group", lambda1 = lambda1, lambda2 = lambda2, return.whole.theta = T)
    aicvf = jaic(fglff, yy)
    aicvg = jaic(fglfg, yy)
    aicf = c(aicf, aicvf)
    aicg = c(aicg, aicvg)
    lam = rbind(lam, c(lambda1, lambda2))
  fc1 = fglff$theta[[1]]; fc2 = fglff$theta[[2]]
  gc1 = fglfg$theta[[1]]; gc2 = fglfg$theta[[2]]
  f.heatmapR = matrix(0, ncol = 2, nrow = nG*(nG - 1)/2)
  g.heatmapR = matrix(0, ncol = 2, nrow = nG*(nG - 1)/2)
  h = 0 
  for(i in 1:(nG-1)){
    for(j in (i+1):nG){
      h = h + 1
      f.heatmapR[h, 1] = sum(fc1[i,j] != 0)
      f.heatmapR[h, 2] = sum(fc2[i,j] != 0)
      g.heatmapR[h, 1] = sum(gc1[i,j] != 0)
      g.heatmapR[h, 2] = sum(gc2[i,j] != 0)
    }
  }
  #######################
  f1 = f.heatmapR[,1]; f2 = f.heatmapR[,2]
  g1 = g.heatmapR[,1]; g2 = g.heatmapR[,2]
 
  fSens1[m, n] = TPR(f1, icl[[1]]); fSpec1[m, n] = 1 - FPR(f1, icl[[1]])
  fSens2[m, n] = TPR(f2, icl[[2]]); fSpec2[m, n] = 1 - FPR(f2, icl[[2]])
  #
  gSens1[m, n] = TPR(g1, icl[[1]]); gSpec1[m, n] = 1 - FPR(g1, icl[[1]])
  gSens2[m, n] = TPR(g2, icl[[2]]); gSpec2[m, n] = 1 - FPR(g2, icl[[2]])
}
}
  
feauc1 <- EAUC(fSens1, fSpec1)
feauc2 <- EAUC(fSens2, fSpec2)
geauc1 <- EAUC(gSens1, gSpec1)
geauc2 <- EAUC(gSens2, gSpec2)



iaf = which.min(aicf); iag = which.min(aicg)
lamf = lam[iaf, ]; lamg = lam[iag, ]
#
fglff = JGL(Y = yy, penalty = "fused", lambda1 = lamf[1], lambda2 = lamf[2], return.whole.theta = T)
fglfg = JGL(Y = yy, penalty = "group", lambda1 = lamg[1], lambda2 = lamg[2], return.whole.theta = T)

fc1 = fglff$theta[[1]]; fc2 = fglff$theta[[2]]
gc1 = fglfg$theta[[1]]; gc2 = fglfg$theta[[2]]


f.heatmapR = matrix(0, ncol = 2, nrow = nG*(nG - 1)/2)
g.heatmapR = matrix(0, ncol = 2, nrow = nG*(nG - 1)/2)

h = 0 
for(i in 1:(nG-1)){
  for(j in (i+1):nG){
    h = h + 1
    f.heatmapR[h, 1] = sum(fc1[i,j] != 0)
    f.heatmapR[h, 2] = sum(fc2[i,j] != 0)
    g.heatmapR[h, 1] = sum(gc1[i,j] != 0)
    g.heatmapR[h, 2] = sum(gc2[i,j] != 0)
  }
}

#
######################
f1 = f.heatmapR[,1]; f2 = f.heatmapR[,2]
g1 = g.heatmapR[,1]; g2 = g.heatmapR[,2]

f <-rbind(
  c(TPR(f1, icl[[1]]), FPR(f1, icl[[1]]), feauc1),
  c(TPR(f2, icl[[2]]), FPR(f2, icl[[2]]), feauc2)
  )
colnames(f) <- c('TPR', 'FPR','AUC')
rownames(f) <- c('Group 1', 'Group 2')
#
g <-rbind(
  c(TPR(g1, icl[[1]]), FPR(g1, icl[[1]]), geauc1),
  c(TPR(g2, icl[[2]]), FPR(g2, icl[[2]]), geauc2)
)
colnames(g) <- c('TPR', 'FPR','AUC')
rownames(g) <- c('Group 1', 'Group 2')
save.image(file = paste0('glf', id, '.rda'))


