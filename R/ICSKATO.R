
# assumes that you already ran ICSKAT.R
# weightVec is a p*1 vector, gMat is an n*p matrix, xMat is an n*q matrix with intercept
ICSKATO <- function(icskatOut, liu=TRUE, rhoVec=c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1)) {
  
  # check the rhoVec
  largeRho <- which(rhoVec >= 1)
  if (length(largeRho) > 1) {
    stop("Nonsense rho")
  } else if (length(largeRho == 1)) {
    rhoVec[which(rhoVec >= 1)] <- 0.999
  }
  
  # get the Qrho, p-value of Qrho, its distribution parameters
  QrhoDF <- QrhoIC(rhoVec = rhoVec, icskatOut = icskatOut, liu=liu)
  
  # calculate the distribution of \kappa 
  # sometimes the machine precision is a little off so this matrix isn't symmetric even though it should be
  # (just numerical rounding errors)
  sig_mat <- icskatOut$sig_mat
  sig_mat[lower.tri(sig_mat)] = t(sig_mat)[lower.tri(sig_mat)]
  # sometimes sig_mat is not full rank, so you have to use pivoting
  zPrelim <- chol(sig_mat, pivot=TRUE)
  # if not full rank, have to modify it a bit
  r <- attr(zMat, 'rank')
  p <- nrow(sig_mat)
  if (r < p) {
    sig_mat[(r+1):p, (r+1):p] <- 0
    oo <- order(attr(zMat, 'pivot'))
    zMat <- zPrelim[, oo]
  } else {zPrelim <- zMat}
  
  # done with decomposing sig_mat
  zBar <- apply(zMat, 1, mean)
  #kappaSubtract <- zBar %*% t(zBar) %*% zMat / (sum(zBar^2))
  # faster, saves n^2 multiplications in getting kappaSubtract
  forDiag <- t(zBar) %*% zMat / sum(zBar^2)
  kappaSubtract <- matrix(data=rep(zBar, p), ncol=p, byrow=FALSE) %*% diag(x = forDiag[1, ]) 
  kappaHalf <- zMat - kappaSubtract
  kappaMat <- t(kappaHalf) %*% kappaHalf
  
  # keep according to SKAT package procedure
  kappaLambda <- eigen(kappaMat, symmetric = TRUE, only.values = TRUE)$values 
  idx1 <- which(kappaLambda >= 0)
  idx2 <- which(kappaLambda > mean(kappaLambda[idx1])/100000)
  if (length(idx2) < 1) {stop("Issue finding eigenvalues for kappa")}
  kappaLambda <- kappaLambda[idx2]
 
  # the moments of the kappa term 
  muK1 <- sum(kappaLambda)
  sigmaZeta <- 2 * sqrt(sum(t(kappaHalf) %*% kappaHalf * t(kappaSubtract) %*% kappaSubtract))
  sigmaK1 <- sqrt(2 * sum(kappaLambda^2) + sigmaZeta^2)
  
  # the \tau(\rho) value that varies with rho
  tauVec <- rep(NA, length(rhoVec))
  for (rho_it in 1:length(rhoVec)) {
    tempRho <- rhoVec[rho_it]
    #tauVec[rho_it] <- p^2 * tempRho * sum(zBar^2) + (1 - tempRho) * t(zBar) %*% zMat %*% t(zMat) %*% zBar / (t(zBar) %*% zBar)
    term1 <- p^2 * tempRho + sum(forDiag[1, ]^2) * (1 - tempRho)
    tauVec[rho_it] <- sum(term1) * sum(zBar^2)
  }
  
  # T statistic
  if (liu) {
    Tstat <- min(QrhoDF$liuPval)
  } else {
    Tstat <- min(QrhoDF$daviesPval)
  }
  
  # need to find qmin(rhov) for all rhov
  qMinVec <- rep(NA, length(rhoVec))
  for (rho_it in 1:length(rhoVec)) {
    # the SKAT package seems to ignore the matched value of \delta, just always set it to 0
    #tempQ <- qchisq(p = 1 - Tstat, ncp = liuDF$delta[rho_it], df = liuDF$df[rho_it])
    tempQ <- qchisq(p = 1 - Tstat, ncp = 0, df = QrhoDF$df[rho_it])
    #muX <- liuDF$delta[rho_it] + liuDF$df[rho_it]
    muX <- QrhoDF$df[rho_it]
    #sigmaX <- sqrt(2 * liuDF$df[rho_it] + 4 * liuDF$delta[rho_it]) 
    sigmaX <- sqrt(2 * QrhoDF$df[rho_it]) 
    qMinVec[rho_it] <- (tempQ - muX) * (QrhoDF$sigmaQrho[rho_it] / sigmaX) + QrhoDF$muQrho[rho_it] 
  }
  
  # append to QrhoDF
  QrhoDF <- QrhoDF %>% mutate(rhoVec = rhoVec, tauVec = tauVec, qMinVec = qMinVec)
  
  # integrate
  intOut <-  integrate(f = fIntegrate, lower=0, upper=40, subdivisions = 1000, muK1 = muK1, sigmaK1 = sigmaK1, sigmaZeta = sigmaZeta, 
                       kappaLambda = kappaLambda, QrhoDF = QrhoDF)
  # final SKATO pvalue
  skatoPval <- 1 - intOut[1]$value
  
  # return
  return(list(pval = skatoPval, QrhoDF=QrhoDF))
}

