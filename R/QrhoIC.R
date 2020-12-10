#' QrhoIC.R
#'
#' Calculate the test statistic, distribution, and p-value for each value of Krho in SKATO.
#'
#' @param rhoVec Numeric vector of the rho values to use in SKATO.
#' @param null_beta The output (a list) from a call to the ICSKAT() function.
#' @param liu Boolean for whether to use Liu (TRUE) or Davies (FALSE) method in calculating p-values.
#' @param kurtQvec Vector of kurtosis of Qrho, from bootstrapping.
#' @param sigmaZeta Standard deviation of zeta term, needed if using kurtosis from bootstrapping.
#' @param tauVec Vector of tau(rho) terms, needed if using kurtosis from bootstrapping.
#' @param alwaysCentral A boolean, if TRUE, follow SKAT package practice of always setting delta=0 in chi-square moment matching.
#' @return Data frame holding the pvalue + test statistic for each fixed rho, the matched noncentrality + degree of freedom parameters
#' according to the Liu method, and the mean and variance of each Qrho.
#'
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' null_fit <- skat_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#' myoutput <- ICskat(left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind=rep(1, n),
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$beta_fit, Itt=null_fit$Itt,
#' gMat=matrix(data=rbinom(n=200*10, size=2, prob=0.3), nrow=200))
#' QrhoIC(rhoVec = seq(from=0, to=1, by=0.1), icskatOut = myoutput)
#'
QrhoIC <- function(rhoVec, icskatOut, liu=TRUE, bootstrapOut=NULL, sigmaZeta=NULL, tauVec=NULL, alwaysCentral=FALSE) {
  nRho <- length(rhoVec)
	p <- nrow(icskatOut$sig_mat)

  # this DF will also hold the Davies p-value and Qrho for each rho
  # in addition to the liu p-value and moment-matching parameters
  liuDF <- data.frame(liuPval = rep(NA, length(rhoVec)), delta=NA, df=NA, muQrho = NA,
                      sigmaQrho = NA, daviesPval=NA, Qrho=NA, alwaysCentral=alwaysCentral)

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
        tempKurt <- kurtQvec[rho_it]
        if (tempKurt <= 0) {tempKurt <- 1e-3}
        tempDF <- 12 / tempKurt
        muQrho <- meanQvec[rho_it]
        sigmaQrho <- sqrt(varQvec[rho_it])

        # record
        liuDF$muQrhoBoot[rho_it] <- muQrho
        liuDF$sigmaQrhoBoot[rho_it] <- sigmaQrho
        liuDF$dfBoot[rho_it] <- tempDF
      }

      # the minimum of these p-values is the test statistic.
      # also, should I ignore the delta and set it to 0 always?
      # if so, should actually change thsi in chiSqMatchFast so that it always returns 0 for delta,
      # that will do it more cleanly.
      liuDF$liuPval[rho_it] <- pchisq(q = (liuDF$Qrho[rho_it] - muQrho) / sigmaQrho * sqrt(2 * tempDF + 4 * tempDelta) + tempDF + tempDelta,
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
