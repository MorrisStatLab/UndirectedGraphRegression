args <- commandArgs(trailingOnly = TRUE)
id <- as.character(args[1])
load(paste0('res',id,'.rda'))
sigmaIJ1 = t(as.matrix(res[[10010]]$sigmaIJ[1,]))
sigmaIJ2 = t(as.matrix(res[[10010]]$sigmaIJ[2,]))
sigmaII = t(as.matrix(res[[10010]]$sigmaII))
for(i in 1:999){
  index = 10010 + i*10
  sigmaIJ1 = rbind(sigmaIJ1, res[[index]]$sigmaIJ[1,])
  sigmaIJ2 = rbind(sigmaIJ2, res[[index]]$sigmaIJ[2,])
  sigmaII =  rbind(sigmaII, res[[index]]$sigmaII)
}
save(sigmaIJ1, sigmaIJ2, 
     sigmaII, file = paste0('post',id,'.rda'))

