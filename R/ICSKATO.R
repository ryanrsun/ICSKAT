#' ICSKATO.R
#'
#' Calculate SKATO test for ICSKAT.
#'
#' @param rhoVec Vector of rhos to search over.
#' @param icskatOut The output list from ICSKAT().
#' @param liu Boolean for whether to use Liu moment matching approximation for p-value of each Qrho (as opposed to Davies).
#' @param liuIntegrate Boolean for whether to use Liu moment matching approximation integration in SKATO p-value (as opposed to Davies).
#' @param kurtQvec Vector of kurtosis of Qrho, from bootstrapping.
#' @param kurtKappa Kurtosis for first part of kappa term, from bootstrapping.
#' @param alwaysCentral A boolean, if TRUE, follow SKAT package practice of always setting delta=0 in chi-square moment matching.
#'
#' @return A list with the elements:
#' \item{pval}{SKATO p-value}
#' \item{QrhoDF}{Data frame containing the distribution and p-value for each Krho.}
#' \item{r}{The rank of the cholesky decomposition of the kappa part}
#' \item{intDavies}{Boolean denoting whether integration was with Davies (true) or Liu method (false)}
#' \item{err}{Boolean, did the integration run into an error}
#'
#' @export
#'
ICSKATO <- function(rhoVec=c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1), icskatOut=icskatOut, useMixtureKurt = FALSE,
                    liu=TRUE, liuIntegrate=FALSE, bootstrapOut = NULL,  alwaysCentral=FALSE) {

  if (!is.null(bootstrapOut)) {
    kurtQvec = bootstrapOut$kurtQvec
    varQvec = bootstrapOut$varQvec
    meanQvec = bootstrapOut$meanQvec
    kurtKappa = bootstrapOut$kurtKappa
    kurtKappaAll = bootstrapOut$kurtKappaAll
    varKappaAll = bootstrapOut$varKappaAll
    sigmaKappaAll = bootstrapOut$varKappaAll
    meanKappaAll = bootstrapOut$meanKappaAll
  } else {
    kurtQvec <- NA; varQvec <- NA; meanQvec <- NA; kurtKappa <- NA;
    kurtKappaAll <- NA; varKappaAll <- NA; sigmaKappaAll <- NA; meanKappaAll <- NA
  }

  # check the rhoVec
  largeRho <- which(rhoVec >= 1)
  if (length(largeRho) > 1) {
    stop("Nonsense rho")
  } else if (length(largeRho == 1)) {
    rhoVec[which(rhoVec >= 1)] <- 0.999
  }

	# calculate the distribution of \kappa
  # sometimes the machine precision is a little off so this matrix isn't symmetric even though it should be
  # (just numerical rounding errors)
  sig_mat <- icskatOut$sig_mat
  sig_mat[lower.tri(sig_mat)] = t(sig_mat)[lower.tri(sig_mat)]
  # sometimes sig_mat is not full rank, so you have to use pivoting
  zPrelim <- chol(sig_mat, pivot=TRUE)
  # if not full rank, have to modify it a bit
  r <- attr(zPrelim, 'rank')
  p <- nrow(sig_mat)
  if (r < p) {
    sig_mat[(r+1):p, (r+1):p] <- 0
    oo <- order(attr(zPrelim, 'pivot'))
    zMat <- zPrelim[, oo]
  } else {zMat <- zPrelim}

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
  lambdaMuK1 <-sum(kappaLambda)
  sigmaZeta <- 2 * sqrt(sum(t(kappaHalf) %*% kappaHalf * t(kappaSubtract) %*% kappaSubtract))
  lambdaSigmaK1 <- sqrt(2 * sum(kappaLambda^2) + sigmaZeta^2)
  lambdaKurtK1 <- 12 * sum(kappaLambda^4) / (sum(kappaLambda^2))^2
  if (is.null(bootstrapOut)) {
    muK1 <- lambdaMuK1
    sigmaK1 <- lambdaSigmaK1
    kurtK1 <- lambdaKurtK1
  } else {  # if bootstrap available, always use those for moments
    muK1 <- meanKappaAll
    sigmaK1 <- sqrt(varKappaAll)
    kurtK1 <- kurtKappaAll
  }
  dfK1 <- 12 / kurtK1

  # the \tau(\rho) value that varies with rho - need this before kappa in case kappa
  tauVec <- rep(NA, length(rhoVec))
  for (rho_it in 1:length(rhoVec)) {
    tempRho <- rhoVec[rho_it]
    #tauVec[rho_it] <- p^2 * tempRho * sum(zBar^2) + (1 - tempRho) * t(zBar) %*% zMat %*% t(zMat) %*% zBar / (t(zBar) %*% zBar)
    term1 <- p^2 * tempRho + sum(forDiag[1, ]^2) * (1 - tempRho)
    tauVec[rho_it] <- sum(term1) * sum(zBar^2)
  }

  # get the Qrho, p-value of Qrho, its distribution parameters
  QrhoDF <- QrhoIC(rhoVec = rhoVec, icskatOut = icskatOut, liu=liu, kurtQvec=kurtQvec, sigmaZeta = sigmaZeta,
                   tauVec = tauVec, alwaysCentral=alwaysCentral)
  # sometimes numerically we just get weird things like only one eigenvalue
  if (class(QrhoDF)[1] == "numeric") { return(list(pval = NA, QrhoDF=NA, r=NA, intDavies = NA, err=1)) }

  # T statistic
  if (liu) {
    pRhoVec <- QrhoDF$liuPval
  } else {pRhoVec <- QrhoDF$daviesPval}
  Tstat <- min(pRhoVec)

  # need to find qmin(rhov) for all rhov
  qMinVec <- rep(NA, length(rhoVec))
  mixDFVec <- rep(NA, length(rhoVec))
  for (rho_it in 1:length(rhoVec)) {

    # it's not clear why when we do bootstrapping, we can't just use the kurtosis from bootstrapping each Qrho, as
    # we did to find the original p-values for each Qrho, but here the SKAT package uses the mixture method
    # for whatever reason for the bootstrap case (see else part).
    if (is.null(bootstrapOut) | !useMixtureKurt) {
      # no bootstrap, use eigenvalues
      if (is.null(bootstrapOut)) {
        tempQ <- qchisq(p = 1 - Tstat, ncp = 0, df = QrhoDF$dfLambda[rho_it])
        muX <- QrhoDF$dfLambda[rho_it]
        sigmaX <- sqrt(2 * QrhoDF$dfLambda[rho_it])
      } else { # bootstrap available
        tempQ <- qchisq(p = 1 - Tstat, ncp = 0, df = QrhoDF$dfBoot[rho_it])
        muX <- QrhoDF$dfBoot[rho_it]
        sigmaX <- sqrt(2 * QrhoDF$dfBoot[rho_it])
      }
    } else {  # mixture requires bootstrap
      tempv1 <- QrhoDF$sigmaQrho[rho_it]^2 + sigmaZeta^2
      mixKurt <- mixture_kurtosis(tempDF1 = kurtKappa, tempDF2 = 1, v1 = tempv1, a1 = 1 - rhoVec[rho_it],
                                  a2 = tauVec[rho_it])
      mixDF <- 12 / mixKurt
      mixDFVec[rho_it] <- mixDF
      tempQ <- qchisq(p = 1 - Tstat, ncp = 0, df = mixDF)
      muX <- mixDF
      sigmaX <- sqrt(2 * mixDF)
    }

    # if bootstrap mean and variance are there, use it
    if (is.null(bootstrapOut)) {
      qMinVec[rho_it] <- (tempQ - muX) * (QrhoDF$sigmaQrho[rho_it] / sigmaX) + QrhoDF$muQrho[rho_it]
    } else {
      qMinVec[rho_it] <- (tempQ - muX) * (QrhoDF$sigmaQrhoBoot[rho_it] / sigmaX) + QrhoDF$muQrhoBoot[rho_it]
    }
  }

  # append to QrhoDF
  QrhoDF <- QrhoDF %>% mutate(rhoVec = rhoVec, tauVec = tauVec, qMinVec = qMinVec)

  # integrate
  if (liuIntegrate) {
    intOut <- tryCatch(integrate(f = fIntegrateLiu, lower=0, upper=40, subdivisions = 1000,
                                 muK1 = muK1, sigmaK1 = sigmaK1, QrhoDF = QrhoDF, dfK1 = dfK1, abs.tol = 10^(-25)), error=function(e) e)
    intDavies <- FALSE
  } else {
    intOut <-  tryCatch(integrate(f = fIntegrate, lower=0, upper=40, subdivisions = 1000,
                                  muK1 = muK1, sigmaK1 = sigmaK1, sigmaZeta = sigmaZeta, kappaLambda = kappaLambda, QrhoDF = QrhoDF, abs.tol = 10^(-25)), error=function(e) e)
    intDavies <- TRUE
    # sometimes the CompQuadForm has numerical issues
    if (class(intOut)[1] == "simpleError") {
      intOut <- tryCatch(integrate(f = fIntegrateLiu, lower=0, upper=40, subdivisions = 1000,
                                   muK1 = muK1, sigmaK1 = sigmaK1, sigmaZeta = sigmaZeta, QrhoDF = QrhoDF, dfK1 = dfK1, abs.tol = 10^(-25)), error=function(e) e)
      intDavies <- FALSE
    }
  }

	# sometimes even the liu integration doesn't work
	if (class(intOut)[1] == "simpleError") {
		return(list(pval = NA, QrhoDF=QrhoDF, r=r, intDavies = intDavies, err=1))
	}

	# ICSKATO p-value
	skatoPval <- 1 - intOut[1]$value

	# SKATO package performs this check as well
	# According to their logic, since SKAT-O is between burden and SKAT,
	# SKAT-O p-value should be <= min(p-values) * 2.
	# To correct conservatively, we use min(p-values) * 3 when number(r.all) >= 3.
	# See SKAT_Optimal_Get_Pvalue_VarMatching() in SKAT_Optimal_VarMatching.R.
	multi <- ifelse(rhoVec < 3, 2, 3)
	posPval <- which(pRhoVec > 0)
	# need correction
	if (skatoPval <= 0) {
	  # same corrections as SKAT package
	  if (length(posPval) < length(rhoVec) & length(posPval) > 0) {
	    correctedP <- min(pRhoVec[posPval])[1]
	  } else if (length(posPval) == length(rhoVec)) {
	    correctedP <- multi * min(pRhoVec[posPval])[1]
	  } else {
	    correctedP <- skatoPval
	  }
	} else {correctedP <- skatoPval}

  # return
  return(list(pval = skatoPval, correctedP = correctedP, QrhoDF=QrhoDF, r=r,
              intDavies = intDavies, err=0,  lambdaKurtK1 = lambdaKurtK1, lambdaSigmaK1 = lambdaSigmaK1,
              lambdaMuK1 = lambdaMuK1, bootKurtKappaAll = kurtKappaAll, bootSigmaKappaAll = sigmaKappaAll,
              bootMuKappaAll = meanKappaAll, mixDFVec = mixDFVec))
}
