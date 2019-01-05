args = commandArgs(trailingOnly = T)
args = as.character(args)

AUC <- function(proba, truth)
{
  r <- rank(proba)
  n_pos <- sum(truth==1)
  n_neg <- length(truth) - n_pos
  auc <- (sum(r[truth==1]) - n_pos*(n_pos+1)/2) / (n_pos*n_neg)
  auc
}

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

load(paste0('icl',args[3],'.rda'))
load(paste0('fdr',args[3], '_', args[1],'_', args[2], '.rda'))
pr1 = apply(as.matrix(pr.heatmapR[,1] != 0), 1, sum)
pr2 = apply(as.matrix(pr.heatmapR[,2] != 0), 1, sum)
pr3 = apply(as.matrix(pr.heatmapR[,3] != 0), 1, sum)

perf <-rbind(
  c(TPR(pr1, icl[[1]]), FPR(pr1, icl[[1]]), AUC(1-rhomn[3,], icl[[1]])),
  c(TPR(pr2, icl[[2]]), FPR(pr2, icl[[2]]), AUC(1-rhomn[4,], icl[[2]])),
  c(TPR(pr3, icl[[3]]), FPR(pr3, icl[[3]]), AUC(1-rhomn[5,], icl[[3]])))
colnames(perf) <- c('TPR', 'FPR', 'AUC')
rownames(perf) <- c('Group 1', 'Group 2', 'Group 3')
save(perf, file = paste0('perf_', args[1],'_',args[2],'_',args[3], '.rda'))


