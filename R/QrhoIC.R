#' QrhoIC.R
#'
#' Calculate the test statistic, distribution, and p-value for each value of Krho in SKATO.
#'
#' @param rhoVec Numeric vector of the rho values to use in SKATO.
#' @param icskatOut The output list returned from  a call to ICSKAT().
#' @param liu Boolean for whether to use Liu (TRUE) or Davies (FALSE) method in calculating p-values for each Qrho.
#' Default is Liu, following SKAT package. If wanting to use bootstrap moments for Qrho, need to use Liu method.
#' @param bootstrapOut The output (a list) from a call the ICSKATO_bootstrap() function, holding moments for Qrho.
#' @param alwaysCentral A boolean, if TRUE, follow SKAT package practice of always setting delta=0 in chi-square moment matching.
#' @return Data frame holding the SKAT pvalue + test statistic for each fixed rho, the matched noncentrality + degree of freedom parameters
#' for each fixed rho (using both bootstrap and analytic calculation), and the mean and variance of each Qrho using both
#' bootstrap and analytic calculation.
#'
#' @importFrom stats pchisq
#'
#' @export
QrhoIC <- function(rhoVec, icskatOut, liu=TRUE, bootstrapOut=NULL, alwaysCentral=FALSE) {
  nRho <- length(rhoVec)
  p <- nrow(icskatOut$sig_mat)

  # this DF will also hold the Davies p-value and Qrho for each rho
  # in addition to the liu p-value and moment-matching parameters
  liuDF <- data.frame(liuPval = rep(NA, length(rhoVec)), delta=NA,
                      daviesPval=NA, Qrho=NA, alwaysCentral=alwaysCentral)

  # there should only be a few rho, so ok to use a for loop
  for (rho_it in 1:nRho) {
    tempRho <- rhoVec[rho_it]
    # calculate Q
    liuDF$Qrho[rho_it] <- (1 - tempRho) * icskatOut$skatQ + tempRho * icskatOut$burdenQ

    # first, the eigenvalues to get the distribution of Qrho
    Rrho <- matrix(data=tempRho, nrow=p, ncol=p)
    diag(Rrho) <- 1
    #toDecomp <- Rrho %*% icskatOut$sig_mat
    toDecomp <- eigenMapMatMultCrossTwo(Rrho, icskatOut$sig_mat)
    Aeigen <- eigen(toDecomp, symmetric = TRUE, only.values = TRUE)
    tempLambda <- Aeigen$values

    # so many weird things can happen with the eigenvalues
    idx1 <- which(tempLambda >= 0)
    idx2 <- which(tempLambda > mean(tempLambda[idx1])/100000)
    if (length(idx2) <= 1) {return(-1)}
    tempLambda <- tempLambda[idx2]

    # liu moment matching to get the distribution
    liuMatch <- chiSqMatchFast(lambdaVec = tempLambda, alwaysCentral=alwaysCentral)
    if (class(liuMatch)[1] == "numeric") {return(-1)}

    # SKATO uses liu by default
    if (liu) {
      # use kurtosis from eigenvalues or from bootstrapping
      if (is.null(bootstrapOut)) {
        tempDF <- liuMatch$l
        tempDelta <- liuMatch$delta
        muQrho <- liuMatch$muQrho
        sigmaQrho <- liuMatch$sigmaQrho

        # record
        liuDF$muQrhoBoot[rho_it] <- NA; liuDF$sigmaQrhoBoot[rho_it] <- NA; liuDF$dfBoot[rho_it] <- NA
      } else {
        tempDelta <- 0
        tempKurt <- bootstrapOut$kurtQvec[rho_it]
        if (tempKurt <= 0) {tempKurt <- 1e-3}
        tempDF <- 12 / tempKurt
        muQrho <- bootstrapOut$meanQvec[rho_it]
        sigmaQrho <- sqrt(bootstrapOut$varQvec[rho_it])

        # record
        liuDF$muQrhoBoot[rho_it] <- muQrho
        liuDF$sigmaQrhoBoot[rho_it] <- sigmaQrho
        liuDF$dfBoot[rho_it] <- tempDF
      }

      # the minimum of these p-values is the test statistic.
      # also, should I ignore the delta and set it to 0 always?
      # if so, should actually change thsi in chiSqMatchFast so that it always returns 0 for delta,
      # that will do it more cleanly.
      liuDF$liuPval[rho_it] <- stats::pchisq(q = (liuDF$Qrho[rho_it] - muQrho) / sigmaQrho * sqrt(2 * tempDF + 4 * tempDelta) + tempDF + tempDelta,
                                      df = tempDF, ncp=tempDelta, lower.tail=F)
    } else {
      # p-value of Qrho
      liuDF$daviesPval[rho_it] <- CompQuadForm::davies(q=liuDF$Qrho[rho_it], lambda=tempLambda)$Qq
      tempDF <- liuMatch$l
    }

    liuDF$muQrhoLambda[rho_it] <- liuMatch$muQrho
    liuDF$sigmaQrhoLambda[rho_it] <- liuMatch$sigmaQrho
    liuDF$dfLambda[rho_it] <- liuMatch$l
    liuDF$delta[rho_it] <- liuMatch$delta

    #cat(rho_it)
  }
  return(liuDF)
}
