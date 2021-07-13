args = commandArgs(trailingOnly = T)
args = as.character(args)


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

load(paste0('icl',args[3],'.rda'))
load(paste0('fdr',args[3], '_', args[1],'_', args[2], '.rda'))
pr1 = apply(as.matrix(pr.heatmapR[,1] != 0), 1, sum)
pr2 = apply(as.matrix(pr.heatmapR[,2] != 0), 1, sum)

perf <-rbind(
  c(TPR(pr1, icl[[1]]), FPR(pr1, icl[[1]])),
  c(TPR(pr2, icl[[2]]), FPR(pr2, icl[[2]])))
colnames(perf) <- c('TPR', 'FPR')
rownames(perf) <- c('Normal Graph', 'Tumor Graph')
save(perf, file = paste0('perf_', args[1],'_',args[2],'_',args[3], '.rda'))


