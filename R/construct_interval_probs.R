#' construct_interval_probs.R
#'
#' Construct the probabilities of falling into each time interval for bootstrapping of interval-censored data.
#'
#' @param allTimes n*s matrix where n is number of subjects and s is all visit times for that subjects.
#' @param dmats Output from make_IC_dmats, a list holding left_dmat and right_dmat.
#' @param nullBeta Vector of coefficients under the null model.
#' @param p Number of covariates in the null model.
#' @param nKnots Number of knots in the spline.
#' @param infVal The numeric value representing time 0 (left-censored observation).
#' @param zeroVal The numeric value representing time infinity (right-censored observation).
#'
#' @return n*(s+1) matrix where element (i,j) holds the probability that subject i will fail in interval j.
#'
#' @importFrom zoo na.locf
#'
#' @export
#' @examples
#' set.seed(2)
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
#' intervalProbOutput <- construct_interval_probs(allTimes = outcomeDat$allVisits,
#' dmats = dmats, nullBeta = nullFit$beta_fit, p = ncol(xMat), nKnots=1)
#'
construct_interval_probs <- function(allTimes, dmats, nullBeta, p, nKnots, infVal=999, zeroVal=0) {

  # number of subjects
  n <- nrow(allTimes)

  # replace the 0s with NAs, then apply na.locf() to fill the missing visits with the last visit time
  # don't have to do this for the first column, because 0 in the first column will naturally just give 0 probability
  for (col_it in 2:ncol(allTimes)) {
    zeroIdx <- which(allTimes[, col_it] == 0)
    if (length(zeroIdx) > 0) {allTimes[zeroIdx, col_it] <- NA}
  }
  filledTimes <- t(apply(allTimes, 1, zoo::na.locf))
  allTimes[, 2:ncol(allTimes)] <- filledTimes[, 2:ncol(filledTimes)]

  # holds the fitted null survival for each visit time
  fittedSurv <- matrix(data=NA, nrow=n, ncol=ncol(allTimes))
  # holds the probability of falling into each interval
  probMat <- matrix(data=NA, nrow=n, ncol=ncol(fittedSurv) + 1)
  # covariate adjustment to baseline hazard
  covariateH <- exp(dmats$left_dmat[, 1:p] %*% as.numeric(nullBeta[1:p]))

  # need to set the quantiles to be the same so we get the same design matrix if we're using the same null model!
  quant_r <- dmats$quant_r

  # loop through and fill probabilities of each interval
  for (time_it in 1:ncol(fittedSurv)) {

    # make design matrix
    temp_obs_ind <- as.numeric(allTimes[, time_it] != infVal)
    temp_tpos_ind <- as.numeric(allTimes[, time_it] != zeroVal)
    tempDmat <- make_IC_dmat(xMat=NULL, lt=allTimes[, time_it], rt=allTimes[, time_it],
                             quant_r=quant_r, obs_ind = temp_obs_ind, tpos_ind = temp_tpos_ind)

    # total baseline hazard
    tempH <- exp(tempDmat$left_dmat %*% as.numeric(nullBeta[(p+1):(p+nKnots+2)])) * covariateH

    # there can be time 0, in that case manually fix to survival = 1
    fittedSurv[, time_it] <- ifelse(temp_tpos_ind == 0, 1, exp(-tempH))
    # if there is time 999, then that survival is 0
    fittedSurv[, time_it] <- ifelse(temp_obs_ind == 0, 0, fittedSurv[, time_it])

    # prob of falling in interval
    if (time_it == 1) {
      probVec <- 1 - fittedSurv[, time_it]
    } else {
      probVec <- fittedSurv[, time_it - 1] - fittedSurv[, time_it]
    }

    # make sure it can't be zero
    if (length(which(probVec < 0)) > 0) {
      probVec[which(probVec < 0)] <- 0
    }
    probMat[, time_it] <- probVec
  }

  # fill the last interval probability
  probVec <- fittedSurv[, ncol(fittedSurv)]
  if (length(which(probVec < 0)) > 0) {
    probVec[which(probVec < 0)] <- 0
  }
  probMat[, ncol(probMat)] <- probVec

  # return
  return(list(probMat = probMat, allTimesFilled = allTimes))
}

