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

	  
EAUCF <- function(Sens, Spec,d0,fname){
phi_grid <- seq(0,1,by=.05)  # for x axis of ROC curve
# ROC with expected sensitivity
ESens <- rep(NA,length(phi_grid)-1)
png(paste0('ROC',fname,'.png'), height = 800, width = 1200)
plot(phi_grid[2:length(phi_grid)] - 0.025, ESens, type="p",ylim=c(0,1),xlab="1-Spec",main=paste("Exp. ROC for true delta ",d0)) # just plots the frame -- no values yet
for (phi in 2:length(phi_grid)){
temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
ESens[phi-1] <- mean(Sens[temp])
points(rep(phi_grid[phi]-0.025,length(temp)),Sens[temp],pch="*") 
}
points(phi_grid[2:length(phi_grid)] - 0.025, ESens, pch=1,type="b") 
dev.off()
}
