load('rresauc.rda')

phi_grid <- seq(0,1,by=.05)  # for x axis of ROC curve
png('RROC1_1_all.png', height = 1000, width = 1000)
par(mfrow=c(1,1), pty ='s', cex = 1.5, cex.axis=2.3,cex.lab = 2.3,cex.main =2.3, mar = c(5,5,4,3), lwd=5)
# ROC with expected sensitivity
ESens <- rep(NA,length(phi_grid)-1)
plot(c(0, phi_grid[2:length(phi_grid)] - 0.025), c(NA,ESens), type="p",ylim=c(0,1), xlim = c(0,1),xlab="Average false positive rate", 
     ylab = 'Average true positive rate', main=bquote(paste(bold("Simulated normal graph "), Omega[N], bold(" (simulation 1)"))), 
     font.lab = 2, font.axis=2, font.main=2, font.sub=2, font=2) 
abline(0,1,lty=5,lwd=4)
#begin calculation
#1
Sens = se1; Spec = sp1
for (phi in 2:length(phi_grid)){
temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
      low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
      high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
      slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
      int <- ESens[high] - slope*(phi_grid[high]+0.025)
      ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
 }
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'red', lty = 1, lwd = 5) 
#2
ESens <- rep(NA,length(phi_grid)-1)
Sens = gse1; Spec = gsp1
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col ='green', lty = 2, lwd = 5) 
#3
ESens <- rep(NA,length(phi_grid)-1)
Sens = fse1; Spec = fsp1
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'blue', lty = 3, lwd = 5) 
#4
ESens <- rep(NA,length(phi_grid)-1)
Sens = lse1; Spec = lsp1
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'orange', lty = 4, lwd = 5) 
#5
ESens <- rep(NA,length(phi_grid)-1)
Sens = nse1; Spec = nsp1
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'brown', lty = 5, lwd = 5) 
#6
ESens <- rep(NA,length(phi_grid)-1)
Sens = mse1; Spec = msp1
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'grey', lty = 1, lwd = 5) 
#7
ESens <- rep(NA,length(phi_grid)-1)
Sens = jse1; Spec = jsp1
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'yellow', lty = 6, lwd = 5) 

#
legend("bottomright", 
        c(as.expression(bquote(bold("Bayesian edge regression"))),
          as.expression(bquote(bold("Bayesian edge regression (Group)"))),
          as.expression(bquote(bold("Group graphical lasso"))),
          as.expression(bquote(bold("Fused graphical lasso"))), 
          as.expression(bquote(bold("LASICH"))), 
          as.expression(bquote(bold("NFMGGM"))), 
          as.expression(bquote(bold("BIMGGM")))), 
       lty=c(1,1,2,3,4,5,6), col=c('red', 'gray','green','blue','orange', 'brown', 'yellow'), lwd = 5,merge=FALSE, cex = 1.75,bty = "n")

dev.off()
#####################################################################################################
phi_grid <- seq(0,1,by=.05)  # for x axis of ROC curve
png('RROC1_2_all.png', height = 1000, width = 1000)
par(mfrow=c(1,1), pty ='s', cex = 1.5, cex.axis=2.3,cex.lab = 2.3,cex.main =2.3, mar = c(5,5,4,3) ,lwd=5)
# ROC with expected sensitivity
ESens <- rep(NA,length(phi_grid)-1)
plot(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(NA,ESens), type="p",ylim=c(0,1), xlim = c(0,1),xlab="Average false positive rate", 
     ylab = 'Average true positive rate', main=bquote(paste(bold("Simulated tumor graph "), Omega[T], bold(" (simulation 1)"))),
     font.lab = 2, font.axis=2, font.main=2, font.sub=2, font=2)  
abline(0,1,lty=5,lwd=4)
#begin calculation
#1
Sens = se2; Spec = sp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'red', lty = 1, lwd = 5) 
#2
ESens <- rep(NA,length(phi_grid)-1)
Sens = gse2; Spec = gsp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col ='green', lty = 2, lwd = 5) 
#3
ESens <- rep(NA,length(phi_grid)-1)
Sens = fse2; Spec = fsp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'blue', lty = 3, lwd = 5) 
#4
ESens <- rep(NA,length(phi_grid)-1)
Sens = lse2; Spec = lsp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'orange', lty = 4, lwd = 5) 
#5
ESens <- rep(NA,length(phi_grid)-1)
Sens = nse2; Spec = nsp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'brown', lty = 5, lwd = 5) 
#6
ESens <- rep(NA,length(phi_grid)-1)
Sens = mse2; Spec = msp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'gray', lty = 1, lwd = 5) 
#7
ESens <- rep(NA,length(phi_grid)-1)
Sens = jse2; Spec = jsp2
for (phi in 2:length(phi_grid)){
  temp <- which((c(1-Spec) >= phi_grid[phi-1]) & (c(1-Spec) < phi_grid[phi]))
  ESens[phi-1] <- mean(Sens[temp])
}
if (is.na(ESens[length(ESens)])) ESens[length(ESens)] <- 1  # if missing boundary
for (ind in which(is.na(ESens))) {
  low <- rev( which(!is.na(ESens[1:ind])) )[1]  # closest index that's not NaN
  high <- (which(!is.na(ESens[(ind+1):length(ESens)]))+ind)[1] # same thing
  slope <- (ESens[high] - ESens[low])/(phi_grid[high]+0.025 - (phi_grid[low]+0.025))
  int <- ESens[high] - slope*(phi_grid[high]+0.025)
  ESens[ind] <- slope*(phi_grid[ind]+0.025) + int   #replace NaN with interpolated value
}
points(c(0,phi_grid[2:length(phi_grid)] - 0.025), c(0,ESens), pch=1,type="l", col = 'yellow', lty = 6, lwd = 5) 

legend("bottomright", 
        c(as.expression(bquote(bold("Bayesian edge regression"))),
          as.expression(bquote(bold("Bayesian edge regression (Group)"))),
          as.expression(bquote(bold("Group graphical lasso"))),
          as.expression(bquote(bold("Fused graphical lasso"))), 
          as.expression(bquote(bold("LASICH"))), 
          as.expression(bquote(bold("NFMGGM"))), 
          as.expression(bquote(bold("BIMGGM")))), 
       lty=c(1,1,2,3,4,5,6), col=c('red', 'gray','green','blue','orange', 'brown', 'yellow'), lwd = 5,merge=FALSE, cex = 1.75,bty = "n")
dev.off()


