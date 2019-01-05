##This function is to calculate M.
library(MASS)
Mout <- function(yy, x){
	p = ncol(yy)
    #x is the tumor purity vector for yy
	yy.tumor <- yy[which(x >= 0.5), ]
	yy.normal <- yy[which(x < 0.5), ]
    # then we calculate the MLE for the precision matrix (inverse of sample covariance matrix)
	n.t = nrow(yy.tumor)
	n.n = nrow(yy.normal)
	prec.tumor <- (solve(cov(yy.tumor)*(n.t - 1)/n.t))^2
	prec.normal <- (solve(cov(yy.normal)*(n.n - 1)/n.n))^2
	M.t <- (sum(sum(prec.tumor)) - sum(diag(prec.tumor)))/2
	M.n <- (sum(sum(prec.normal)) - sum(diag(prec.normal)))/2
	M <- c(M.n, M.t)/(p*(p-1)/2)
	return(M)
}
