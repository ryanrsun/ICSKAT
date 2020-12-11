#' ICSKATO_bootstrap.R
#'
#' The version of ICSKATO to run when bootstrapping to match kurtosis.
#'
#' @param icskatOut The output list from ICSKAT().
#' @param B Number of bootstrap replications.
#' @param intervalProbs n*(s+1) matrix where n is number of subjects and s is number of visits, probability of falling in each interval.
#' @param allVisits n*s matrix with all the visit times for each subject.
#' @param quant_r Quantiles of time from make_IC_dmats, to keep them constant through bootstrapping.
#' @param null_fit The null fit output from ICSKAT_fit_null.
#' @param gMat Genotype matrix used in original test.
#' @param fitAgain Boolean, whether to fit the null model again in each bootstrap.
#' @param checkpoint Boolean, whether to print every time 100 bootstraps finish.
#' @param rhoVec Vector of rhos to search over in SKATO.
#'
#' @return A list with the elements:
#' \item{Qrho}{Vector of SKAT test statistics, one for each rho}
#' \item{kappaPart}{Scalar first term in kappa}
#'
#' @export
#'
ICSKATO_bootstrap <- function(icskatOut, B, intervalProbs, allVisits, quant_r, seed = NULL,
                              null_fit, gMat, fitAgain, checkpoint=FALSE, downsample=1, rhoVec) {

  # set seed
  if (!is.null(seed)) { set.seed(seed) }

  # bootstrap
  # really the only quantity we need is kappa and QrhoBoot
  bootDF <- data.frame(mu=rep(NA, B), sigma=NA, testStat=NA, pSKAT=NA, kappa=NA, kappaAll=NA)
  QrhoBoot <- matrix(data=NA, nrow=B, ncol=7)
  for (boot_it in 1:B) {

    # draw from model
    draws <- t(sapply(data.frame(t(intervalProbs)), FUN = rmultinom, n=1, size=1))
    modelTimes <- t(mapply(matchVisit, data.frame(t(draws)), data.frame(t(allVisits))))
    newLT <- modelTimes[, 1]
    newRT <- modelTimes[, 2]

    # new design matrix - feed it the old quant_r
    newDmat <- make_IC_dmat(xMat=xMat, lt=newLT, rt=newRT, quant_r = quant_r)
    newLeft <- newDmat$left_dmat
    newRight <- newDmat$right_dmat
    obs_ind <- as.numeric(newRT < 999)
    tpos_ind <- as.numeric(newLT> 0)

    # fit the null to get the new Itt
    # give it the old beta
    # or should I do the entire null fit again?
    if (fitAgain) {
      newNull <- ICSKAT_fit_null(init_beta=as.numeric(null_fit$beta_fit), lt=newLT, rt=newRT,
                                 left_dmat=newLeft, right_dmat=newRight, runOnce = FALSE,
                                 obs_ind=obs_ind, tpos_ind=tpos_ind)

      bootSKAT <- ICskat(left_dmat=newLeft, tpos_ind=tpos_ind, obs_ind=obs_ind,
                         right_dmat=newRight, gMat=gMat, lt=newLT, rt=newRT,
                         null_beta=as.numeric(newNull$beta_fit), Itt=newNull$Itt)
    } else {

      newNull <- ICSKAT_fit_null(init_beta=as.numeric(null_fit$beta_fit), lt=newLT, rt=newRT,
                                 left_dmat=newLeft, right_dmat=newRight, runOnce = TRUE,
                                 obs_ind=obs_ind, tpos_ind=tpos_ind)

      # run skat on bootstrapped data
      # give it the old beta but the null Itt
      bootSKAT <- ICskat(left_dmat=newLeft, tpos_ind=tpos_ind, obs_ind=obs_ind,
                         right_dmat=newRight, gMat=gMat, lt=newLT, rt=newRT,
                         null_beta=as.numeric(null_fit$beta_fit), Itt=newNull$Itt)
    }

    # fix eigenvalues
    newLambda <- bootSKAT$lambdaQ
    idx1 <- which(newLambda >= 0)
    idx2 <- which(newLambda > mean(newLambda[idx1])/100000)
    newLambda <- newLambda[idx2]

    # record
    bootDF$mu[boot_it] <- sum(newLambda)
    bootDF$sigma[boot_it] <- sqrt(2 * sum(newLambda^2))
    bootDF$testStat[boot_it] <- bootSKAT$skatQ
    bootDF$pSKAT[boot_it] <- bootSKAT$p_SKAT

    # it turns out the kappa part can be calculated with just sig_mat and the Ugamma
    kappaPortion <- matrix(data=rep(colSums(bootSKAT$sig_mat), ncol(bootSKAT$sig_mat)), nrow=ncol(bootSKAT$sig_mat), byrow=TRUE) /
      sum(bootSKAT$sig_mat)
    UgammaKappa <- t(as.numeric(bootSKAT$Ugamma))  %*% kappaPortion
    kappaVec <- t(as.numeric(bootSKAT$Ugamma)) -
      UgammaKappa
    bootDF$kappa[boot_it] <- sum(kappaVec^2)
    bootDF$kappaAll[boot_it] <- bootDF$kappa[boot_it] + 2 * sum(kappaVec * UgammaKappa)
    QrhoBoot[boot_it, ] <-  (1 - rhoVec) * as.numeric(bootSKAT$skatQ) +
      rhoVec * as.numeric(bootSKAT$burdenQ)

    if (boot_it%%100 == 0) {
      if (checkpoint) {cat(boot_it)}
    }
  }

  # calculate the kurtosis of Qrho and of kappa
  # now also getting the mean and variance
  kurtQvec <- rep(0, ncol(QrhoBoot))
  varQvec <- rep(0, ncol(QrhoBoot))
  meanQvec <- rep(0, ncol(QrhoBoot))
  for (i in 1:length(kurtQvec)) {
    kurtQvec[i] <- mean((QrhoBoot[, i] - mean(QrhoBoot[, i]))^4) / mean((QrhoBoot[, i] - mean(QrhoBoot[, i]))^2)^2 - 3
    varQvec[i] <- mean((QrhoBoot[, i] - mean(QrhoBoot[, i]))^2)
    meanQvec[i] <- mean(QrhoBoot[, i])
  }
  kurtKappa <- mean( (bootDF$kappa - mean(bootDF$kappa))^4 ) / mean( (bootDF$kappa - mean(bootDF$kappa))^2 )^2 - 3
  kurtKappaAll <- mean( (bootDF$kappaAll - mean(bootDF$kappaAll))^4 ) / mean( (bootDF$kappaAll - mean(bootDF$kappaAll))^2 )^2 - 3
  varKappaAll <- mean( (bootDF$kappaAll - mean(bootDF$kappaAll))^2 )
  meanKappaAll <- mean(bootDF$kappaAll)

  # if downsample < 1, give results with less bootstrapping to see if it's feasible
  if (downsample < 1) {
    dsIdx <- round(downsample * B)
    kurtQvecDS <- rep(0, ncol(QrhoBoot))
    varQvecDS <- rep(0, ncol(QrhoBoot))
    meanQvecDS <- rep(0, ncol(QrhoBoot))
    for (i in 1:length(kurtQvec)) {
      kurtQvecDS[i] <- mean((QrhoBoot[1:dsIdx, i] - mean(QrhoBoot[1:dsIdx, i]))^4) / mean((QrhoBoot[1:dsIdx, i] - mean(QrhoBoot[1:dsIdx, i]))^2)^2 - 3
      varQvecDS[i] <- mean((QrhoBoot[1:dsIdx, i] - mean(QrhoBoot[1:dsIdx, i]))^2)
      meanQvecDS[i] <- mean(QrhoBoot[1:dsIdx, i])
    }
    kurtKappaDS <- mean( (bootDF$kappa[1:dsIdx] - mean(bootDF$kappa[1:dsIdx]))^4 ) / mean( (bootDF$kappa[1:dsIdx] - mean(bootDF$kappa[1:dsIdx]))^2 )^2 - 3
    kurtKappaAllDS <- mean( (bootDF$kappaAll[1:dsIdx]  - mean(bootDF$kappaAll[1:dsIdx]))^4 ) / mean( (bootDF$kappaAll[1:dsIdx]- mean(bootDF$kappaAll[1:dsIdx]))^2 )^2 - 3
    varKappaAllDS <- mean( (bootDF$kappaAll[1:dsIdx]  - mean(bootDF$kappaAll[1:dsIdx]))^2 )
    meanKappaAllDS <- mean(bootDF$kappaAll[1:dsIdx])

    # package it all inside another list
    listDS <- list(kurtQvec = kurtQvecDS, varQvec = varQvecDS, meanQvec = meanQvecDS, kurtKappa = kurtKappaDS,
                   kurtKappaAll = kurtKappaAllDS, varKappaAll = varKappaAllDS, meanKappaAll = meanKappaAllDS,
                   bootDF = bootDF[1:dsIdx, ], QrhoBoot = QrhoBoot[1:dsIdx, ])
  }

  # return
  return(list(kurtQvec = kurtQvec, varQvec = varQvec, meanQvec = meanQvec, kurtKappa = kurtKappa,
              kurtKappaAll = kurtKappaAll, varKappaAll = varKappaAll, meanKappaAll = meanKappaAll,
              bootDF = bootDF, QrhoBoot = QrhoBoot, listDS))
}

