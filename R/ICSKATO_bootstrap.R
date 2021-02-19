#' ICSKATO_bootstrap.R
#'
#' The version of ICSKATO to run when bootstrapping to match kurtosis.
#'
#' @param icskatOut The output list from ICSKAT().
#' @param B Number of bootstrap replications.
#' @param intervalProbs n*(s+1) matrix where n is number of subjects and s is the number of visits possible, probability of falling in each interval.
#' @param allVisits n*s matrix with all the visit times for each subject.
#' @param quant_r Quantiles of time from make_IC_dmats, to keep them constant through bootstrapping.
#' @param seed Seed to start the bootstrapping.
#' @param null_fit The null fit output from ICSKAT_fit_null.
#' @param gMat Genotype matrix used in original test.
#' @param fitAgain Boolean, whether to fit the null model again in each bootstrap.
#' @param checkpoint Boolean, whether to print every time 100 bootstraps finish.
#' @param downsample A number in (0, 1], will use this fraction of the bootstrap iterations to try running the test with fewer bootstraps.
#' @param rhoVec Vector of rhos to search over in SKATO.
#'
#' @return A list with the elements:
#' \item{kurtQvec}{Vector of bootstrapped excess kurtosis of each Qrho.}
#' \item{varQvec}{Vector of bootstrapped variance of each Qrho.}
#' \item{meanQvec}{Vector of bootstrapped mean of each Qrho.}
#' \item{kurtKappa}{Bootstrapped kurtosis of kappa term without xi.}
#' \item{kurtKappaAll}{Bootstrapped kurtosis of full kappa term with xi.}
#' \item{varKappaAll}{Bootstrapped variance of full kappa term with xi.}
#' \item{meanKappaAll}Bootstrapped mean of full kappa term with xi.}
#' \item{bootDF}{Matrix with B rows containing all the bootstrapped quantities over all iterations.}
#' \item{QrhoBoot}{Matrix with B rows containing all the bootstrapped Qrho values, one column for each rho.}
#' \item{listDS}{A list containing all of the other elements in this return list, except using the downsampled iterations.}
#' \item{nonNA}{Number of bootstraps that did not result in NA (and thus were not removed).}
#'
#' @export
#' @examples
#' gMat <- matrix(data=rbinom(n=200, size=2, prob=0.3), nrow=100)
#' xMat <- matrix(data=rnorm(200), nrow=100)
#' bhFunInv <- function(x) {x}
#' obsTimes <- 1:5
#' etaVec <- rep(0, 100)
#' outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTime = obsTime, windowHalf = 0.1,
#' probMiss = 0.1, etaVec = etaVec)
#' lt <- outcomeDat$leftTimes
#' rt <- outcomeDat$rightTimes
#' tpos_ind <- as.numeric(lt > 0)
#' obs_ind <- as.numeric(rt != Inf)
#' dmats <- make_IC_dmat(xMat, lt, rt)
#' nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, 
#' obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
#' icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
#' intervalProbOutput <- construct_interval_probs(allTimes = outcomeDat$allVisits, dmats = dmats,
#' nullBeta = nullFit$beta_fit, p = ncol(xMat), nKnots=1)
#' ICSKATO_bootstrap(icskatOut = icSkatOut, B = 100, intervalProbs = intervalProbOutput$probMat,
#' allVisits = intervalProbOutput$allTimesFilled, quant_r = dmats$quant_r, seed = 0, 
#' null_fit = nullFit, gMat = gMat, fitAgain = TRUE, rhoVec=c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1))
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
      PHres <- ICSKATwrapper(left_dmat = newLeft, right_dmat = newRight, initValues = as.numeric(null_fit$beta_fit), 
        lt=newLT, rt=newRT, runOnce = FALSE, obs_ind=obs_ind, tpos_ind=tpos_ind, gMat=gMat, PH = TRUE, nKnots=1, maxIter = 5, returnNull=TRUE)
      newNull <- PHres$nullFit
      bootSKAT <- PHres$skatOutput 
    } else {
      PHres <- ICSKATwrapper(left_dmat = newLeft, right_dmat = newRight, initValues = as.numeric(null_fit$beta_fit), 
        lt=newLT, rt=newRT, runOnce = TRUE, obs_ind=obs_ind, tpos_ind=tpos_ind, gMat=gMat, PH = TRUE, nKnots=1, maxIter = 1, returnNull=TRUE)
      newNull <- PHres$nullFit
      bootSKAT <- PHres$skatOutput
    }
    # if error, go on
    if (bootSKAT$err != 0 & bootSKAT$err != 22) {next}

    # fix eigenvalues
    newLambda <- as.numeric(bootSKAT$lambdaQ)
    idx1 <- which(newLambda >= 0)
    idx2 <- which(newLambda > mean(newLambda[idx1])/100000)
    newLambda <- newLambda[idx2]

    # record
    bootDF$mu[boot_it] <- sum(newLambda)
    bootDF$sigma[boot_it] <- sqrt(2 * sum(newLambda^2))
    bootDF$testStat[boot_it] <- bootSKAT$skatQ
    bootDF$pSKAT[boot_it] <- bootSKAT$p_SKAT

    # it turns out the kappa part can be calculated with just sig_mat and the Ugamma
    # we're calculating u^T(I-M)ZZ^T(I-M)u, half of this is u^T(I-M)Z
    # remember u^TZ is just the score vector since u=V^-1/2WG^Tr and Z is the upper triangular chol decomp L^T
    # MZ is just L^TJ(J^TLL^T)^-1J^TLL^T, the first L^T multiplies with U to just result in the score vector,
    # the inverse part is just the sum of all the elements of V, JJ^T is a matrix of ones, and then LL^T = V,
    # so the matrix of ones times V is just a matrix with each row the same, the column sums of V.
    kappaPortion <- matrix(data=rep(colSums(bootSKAT$sig_mat), ncol(bootSKAT$sig_mat)), nrow=ncol(bootSKAT$sig_mat), byrow=TRUE) /
      sum(bootSKAT$sig_mat)
    UgammaKappa <- t(as.numeric(bootSKAT$Ugamma))  %*% kappaPortion
    kappaVec <- t(as.numeric(bootSKAT$Ugamma)) -
      UgammaKappa
    bootDF$kappa[boot_it] <- sum(kappaVec^2)
    # the sum(kappaVec * UgammaKappa) is just 2 * u^T(I-M)Z %*% Z^TMu which is the xi part we need,
    # remember that kappaVec is u^T(I-M)Z  and UgammaKappa is Z^TMu.
    bootDF$kappaAll[boot_it] <- bootDF$kappa[boot_it] + 2 * sum(kappaVec * UgammaKappa)
    QrhoBoot[boot_it, ] <-  (1 - rhoVec) * as.numeric(bootSKAT$skatQ) +
      rhoVec * as.numeric(bootSKAT$burdenQ)

    if (boot_it%%10 == 0) {
      if (checkpoint) {cat(boot_it)}
    }
  }
  # remove the NAs
  isNA <- which(is.na(bootDF$kappa))
  if (length(isNA) > 0) {
    bootDF <- bootDF[-isNA, ]
    QrhoBoot <- QrhoBoot[-isNA, ]
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
  } else {listDS <- NA}

  # return
  return(list(kurtQvec = kurtQvec, varQvec = varQvec, meanQvec = meanQvec, kurtKappa = kurtKappa,
              kurtKappaAll = kurtKappaAll, varKappaAll = varKappaAll, meanKappaAll = meanKappaAll,
              bootDF = bootDF, QrhoBoot = QrhoBoot, listDS = listDS, nonNA = B - length(isNA)))
}

