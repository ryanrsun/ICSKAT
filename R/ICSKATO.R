#' ICSKATO.R
#'
#' Calculate SKATO test for ICSKAT.
#'
#' @param rhoVec Vector of rhos to search over.
#' @param icskatOut The output list from ICSKAT().
#' @param useMixtureKurt Boolean for whether to use the mixture formula to estimate the kurtosis of Qrho when we
#' have bootstrap results. Default is false, instead we just use the bootstraped kurtosis of Qrho.
#' @param liu Boolean for whether to use Liu moment matching approximation for p-value of each Qrho (as opposed to Davies).
#' If Davies, cannot use bootstrapped moments of Qrho.
#' @param liuIntegrate Boolean for whether to use Liu moment matching approximation integration in SKATO p-value (as opposed to Davies).
#' @param bootstrapOut Output list from call to ICSKATO_bootstrap().
#' @param alwaysCentral A boolean, if TRUE, follow SKAT package practice of always setting delta=0 in chi-square moment matching.
#'
#' @return A list with the elements:
#' \item{pval}{SKATO p-value.}
#' \item{correctedP}{Corrected SKATO p-value, which will be the same as pval when not all Qrho values produce
#' a p-value between 0 and 1 (e.g. sometimes it will be 0). Correction is same as SKAT package correction..}
#' \item{QrhoDF}{Data frame containing the distribution and p-value for each Qrho.}
#' \item{r}{The rank of the cholesky decomposition of the sig_mat returned from ICSKAT(), i.e. V^-1/2 or Z.}
#' \item{intDavies}{Boolean denoting whether integration was with Davies (true) or Liu method (false).}
#' \item{err}{0 is no error, 1 is early error like possibly only one eigenvalue/issue with sigmat/issue with kappaMat/issue with QrhoDF,
#' 2 is corrected p-value (fine), 3 is integration error, 9 is no positive p-values (so SKATOp should be 0 unless burden is 1).}
#' \item{lambdaKurtK1}{Kurtosis of kappa term minus zeta using eigenvalues, we use it to approximate the kurtosis of the entire kappa.}
#' \item{lambdaSigmaK1}{Standard deviation of kappa term, including zeta, using eigenvalues.}
#' \item{lambdaMuK1}{Mean of kappa term using eigenvalues.}
#' \item{bootKurtKappaAll}{Kurtosis of entire kappa term, including zeta, using bootstrap data}
#' \item{bootSigmaKappaAll}{Standard deviation of entire kappa term using bootstrap data.}
#' \item{bootMuKappaAll}{Mean of entire kappa term using bootstrap data.}
#' \item{mixDFVec}{Degrees of freedom of Qrho if useMixtureKurt is true, only here to match SKAT package, not really used.}
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom stats qchisq
#' @importFrom stats integrate
#' @export
#' @examples
#' set.seed(1)
#' gMat <- matrix(data=rbinom(n=200, size=2, prob=0.3), nrow=100)
#' xMat <- matrix(data=rnorm(200), nrow=100)
#' bhFunInv <- function(x) {x}
#' obsTimes <- 1:5
#' etaVec <- rep(0, 100)
#' outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, windowHalf = 0.1,
#' probMiss = 0.1, etaVec = etaVec)
#' lt <- outcomeDat$leftTimes
#' rt <- outcomeDat$rightTimes
#' tpos_ind <- as.numeric(lt > 0)
#' obs_ind <- as.numeric(rt != Inf)
#' dmats <- make_IC_dmat(xMat, lt, rt, obs_ind, tpos_ind)
#' nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind = obs_ind, tpos_ind = tpos_ind,
#' lt = lt, rt = rt)
#' icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
#' lt = lt, rt = rt, obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat,
#' null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
#' ICSKATO(icskatOut = icskatOut)
ICSKATO <- function(rhoVec=c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1), icskatOut , useMixtureKurt = FALSE,
                    liu=TRUE, liuIntegrate=FALSE, bootstrapOut = NULL,  alwaysCentral=FALSE) {

  # load bootstrap information
  if (!is.null(bootstrapOut)) {
    kurtQvec = bootstrapOut$kurtQvec
    varQvec = bootstrapOut$varQvec
    meanQvec = bootstrapOut$meanQvec
    kurtKappa = bootstrapOut$kurtKappa
    kurtKappaAll = bootstrapOut$kurtKappaAll
    varKappaAll = bootstrapOut$varKappaAll
    sigmaKappaAll = bootstrapOut$varKappaAll
    meanKappaAll = bootstrapOut$meanKappaAll
    if (kurtKappa <= 0) {kurtKappa <- 1e-3}
    if (kurtKappaAll <= 0) {kurtKappaAll <- 1e-3}
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
  # just write it out or trust me that kappaHalf = (I-M)Z
  kappaHalf <- zMat - kappaSubtract
  # we need the eigenvalues of Z^T(I-M)Z which is equal to Z^T(I-M)(I-M)Z
  kappaMat <- t(kappaHalf) %*% kappaHalf
  # sometimes sig_mat is just 0
  if (length(which(is.na(kappaMat))) > 0) { return(list(pval = NA, QrhoDF=NA, r=NA, intDavies = NA, err=1)) }

  # keep according to SKAT package procedure
  # it's a little confusing because kappa lambda does not consider the zeta part of kappa
  kappaLambda <- eigen(kappaMat, symmetric = TRUE, only.values = TRUE)$values
  idx1 <- which(kappaLambda >= 0)
  idx2 <- which(kappaLambda > mean(kappaLambda[idx1])/100000)
  if (length(idx2) < 1) {
    return(list(pval = NA, QrhoDF=NA, r=r, intDavies = FALSE, err=1))
  }
  kappaLambda <- kappaLambda[idx2]

  # the moments of the kappa term minus the zeta term
  lambdaMuK1 <-sum(kappaLambda)
  sigmaZeta <- 2 * sqrt(sum(t(kappaHalf) %*% kappaHalf * t(kappaSubtract) %*% kappaSubtract))
  lambdaSigmaK1 <- sqrt(2 * sum(kappaLambda^2) + sigmaZeta^2)
  lambdaKurtK1 <- 12 * sum(kappaLambda^4) / (sum(kappaLambda^2))^2
  if (is.null(bootstrapOut)) {
    muK1 <- lambdaMuK1
    sigmaK1 <- lambdaSigmaK1
    kurtK1 <- lambdaKurtK1
  } else {  # if bootstrap available, always use those for moments
    # here we are getting the moments of the entire kappa term, with the zeta term, so likely more accruate
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
  QrhoDF <- QrhoIC(rhoVec = rhoVec, icskatOut = icskatOut, liu=liu, bootstrapOut = bootstrapOut,
                   alwaysCentral=alwaysCentral)
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
        tempQ <- stats::qchisq(p = 1 - Tstat, ncp = 0, df = QrhoDF$dfLambda[rho_it])
        muX <- QrhoDF$dfLambda[rho_it]
        sigmaX <- sqrt(2 * QrhoDF$dfLambda[rho_it])
      } else { # bootstrap available
        tempQ <- stats::qchisq(p = 1 - Tstat, ncp = 0, df = QrhoDF$dfBoot[rho_it])
        muX <- QrhoDF$dfBoot[rho_it]
        sigmaX <- sqrt(2 * QrhoDF$dfBoot[rho_it])
      }
    } else {  # mixture requires bootstrap
      tempv1 <- QrhoDF$sigmaQrho[rho_it]^2 + sigmaZeta^2
      mixKurt <- mixture_kurtosis(tempDF1 = kurtKappa, tempDF2 = 1, v1 = tempv1, a1 = 1 - rhoVec[rho_it],
                                  a2 = tauVec[rho_it])
      mixDF <- 12 / mixKurt
      mixDFVec[rho_it] <- mixDF
      tempQ <- stats::qchisq(p = 1 - Tstat, ncp = 0, df = mixDF)
      muX <- mixDF
      sigmaX <- sqrt(2 * mixDF)
    }

    # if bootstrap mean and variance are there, use it
    if (is.null(bootstrapOut)) {
      qMinVec[rho_it] <- (tempQ - muX) * (QrhoDF$sigmaQrhoLambda[rho_it] / sigmaX) + QrhoDF$muQrhoLambda[rho_it]
    } else {
      qMinVec[rho_it] <- (tempQ - muX) * (QrhoDF$sigmaQrhoBoot[rho_it] / sigmaX) + QrhoDF$muQrhoBoot[rho_it]
    }
    # the reason we are using so many different moments: we are approximating a centered and scaled Qrho
    # (where the Qrho moments come in) by another centered and scaled chi-square variable (call it X_df),
    # where X_df has a degree of freedom matched to the kurtosis of Qrho (and thus the mean of X_df is different
    # from the mean of Qrho, because the mean of Qrho is calculated from bootstrapping Qrho, and the mean of X_df
    # is calculated from bootstrapping the kurtosis of Qrho and dividing that by 12).
    # so for something like Pr(Q < q1) = a, we need to center and scale Q and q1, then replace standardized Q
    # by standardized X_df.
  }

  # append to QrhoDF
  QrhoDF <- QrhoDF %>% dplyr::mutate(rhoVec = rhoVec, tauVec = tauVec, qMinVec = qMinVec)

  # integrate
  if (liuIntegrate) {
    intOut <- tryCatch(stats::integrate(f = fIntegrateLiu, lower=0, upper=40, subdivisions = 1000,
                                 muK1 = muK1, sigmaK1 = sigmaK1, QrhoDF = QrhoDF, dfK1 = dfK1, abs.tol = 10^(-25)), error=function(e) e)
    intDavies <- FALSE
  } else {
    intOut <-  tryCatch(stats::integrate(f = fIntegrate, lower=0, upper=40, subdivisions = 1000,
                                  muK1 = muK1, sigmaK1 = sigmaK1, sigmaZeta = sigmaZeta, kappaLambda = kappaLambda, QrhoDF = QrhoDF, abs.tol = 10^(-25)), error=function(e) e)
    intDavies <- TRUE
    # sometimes the CompQuadForm has numerical issues
    if (class(intOut)[1] == "simpleError") {
      intOut <- tryCatch(stats::integrate(f = fIntegrateLiu, lower=0, upper=40, subdivisions = 1000,
                                   muK1 = muK1, sigmaK1 = sigmaK1, QrhoDF = QrhoDF, dfK1 = dfK1, abs.tol = 10^(-25)), error=function(e) e)
      intDavies <- FALSE
    }
  }

  # A check in case some of the pRhoVec values are 0 or close to zero, in which case
  # the qMinVec values can be Inf and cause errors in the integration.
  # SKATO package performs this check as well.
  # According to their logic, since SKAT-O is between burden and SKAT,
  # SKAT-O p-value should be <= min(p-values) * 2.
  # See SKAT_Optimal_Get_Pvalue_VarMatching() in SKAT_Optimal_VarMatching.R.
  multi <- ifelse(length(rhoVec) < 3, 2, 3)
  posPval <- which(pRhoVec > 0)

  # if there is a zero or negative p-value in pRhoVec, corrected should be minimum positive p-value
  if (length(posPval) < length(rhoVec) & length(posPval) > 0) {
    correctedP <- min(pRhoVec[posPval])[1]
  } else if (length(posPval) == length(rhoVec)) {
    # here all the pRhoVec values are positive
    correctedP <- multi * min(pRhoVec[posPval])[1]
  } else {
    # here there are no positive p-values in pRhoVec
    # return error 9
    return(list(pval = NA, QrhoDF=QrhoDF, r=r, intDavies = intDavies, err=9))
  }

  # sometimes the liu integration doesn't work
  if (class(intOut)[1] == "simpleError") {
    # if the reason it didn't work is because some of the pRhoVec values are too small, give it a corrected p-value
    # error 2 means corrected p-value
    if (length(which(qMinVec == Inf)) > 0) {
      return(list(pval = correctedP, QrhoDF=QrhoDF, r=r, intDavies = intDavies, err=2))
    } else {
      # error 3 is all other errors
      return(list(pval = NA, QrhoDF=QrhoDF, r=r, intDavies = intDavies, err=3))
    }
  }

  # ICSKATO p-value
  skatoPval <- 1 - intOut[1]$value

  # return
  return(list(pval = skatoPval, correctedP = correctedP, QrhoDF=QrhoDF, r=r,
              intDavies = intDavies, err=0,  lambdaKurtK1 = lambdaKurtK1, lambdaSigmaK1 = lambdaSigmaK1,
              lambdaMuK1 = lambdaMuK1, bootKurtKappaAll = kurtKappaAll, bootSigmaKappaAll = sigmaKappaAll,
              bootMuKappaAll = meanKappaAll, mixDFVec = mixDFVec))
}
