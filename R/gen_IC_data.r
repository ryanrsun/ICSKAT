#' gen_IC_data.R
#'
#' Generate interval-censored data under the proportional odds/PH model given a baseline hazard function and
#' some information about observation times. 
#'
#' @param bhFunInv A function, the inverse of the baseline hazard function.
#' @param obsTimes Vector of the intended observation times.
#' @param windowHalf The amount of time before or after the intended obsTimes that a visit might take place.
#' @param etaVec n*1 linear predictor in either the proportional odds or proportional hzards model.
#' @param probMiss The probability of missing any given visit.
#'
#' @return A list with the elements:
#' \item{deltaVec}{n*1 indicator of whether the event has happened for each subject}
#' \item{tVec}{n*1 vector of exact event times}
#' \item{leftTimes}{n*1 vector of left side of interval for each subject}
#' \item{rightTimes}{n*1 vector of right side of interval for each subject}
#'
#' @export
#' @examples
#' bhFunInv <- function(x) {sqrt(2*x)}
#' obsTimes <- seq(from=0.5, to=3.5, by=0.5)
#'  xMat <- cbind(rnorm(n), rbinom(n=n, size=1, prob=0.5))
#' betaVec <- c(-0.5, 0.5)
#' etaVec <- xMat %*% betaVec
#' outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, windowHalf = 0.25, etaVec = etaVec, n=n)
#'
gen_IC_data <- function(bhFunInv, obsTimes, windowHalf, etaVec, probMiss=0.1) {
  
  # sample size
  n <- length(etaVec)
  
  # generate the exact times
  uVec <- runif(n=n, min=0, max=1)
  tVec <- bhFunInv( -log(1 - uVec) / exp(etaVec) )
  
  # there is a 10% chance of missing a visit
  nVisits <- length(obsTimes)
  madeVisit <- matrix(data=rbinom(n=n*nVisits, size=1, prob=0.9), nrow=n, ncol=7)
  # your visit is uniformly distributed around the intended obsTime, windowHalf on each side
  visitTime <- sweep(matrix(data=runif(n=n*nVisits, min=-windowHalf, max=windowHalf), nrow=n, ncol=7),
                     MARGIN=2, STATS=obsTimes, FUN="+")
  
  # get all visits for each subject
  allVisits <- madeVisit * visitTime
  # make the interval for each subject
  allInts <- t(mapply(FUN=createInt, obsTimes = data.frame(t(allVisits)), eventTime=tVec))
  leftTimes <- allInts[, 1]
  rightTimes <- allInts[, 2]
  # event time indicators
  deltaVec <- ifelse(rightTimes == 999, 0, 1)
  
  # return
  return(list(deltaVec = deltaVec, tVec = tVec, leftTimes = leftTimes, rightTimes = rightTimes))
}

#' Called by gen_IC_data() to turn the actual outcome times and observation times into interval-censored
#' outcomes for each subject. Apply this with mapply over a data.frame of visit times, pass in the exact times.
#' Returns 999 instead of Inf.
#'
#' @param obsTimes A vector of all the times a subject is observed.
#' @param eventTime The exact event time for the subject.
#'
#' @return A 2*1 vector which is the interval of the event time
#'
#' @export
#' @examples
#' obsTimes <- 1:10
#' eventTime <- 7.7
#' createInt(obsTimes, eventTime)
#'
createInt <- function(obsTimes, eventTime) {
  # order the times in case the random portion causes them to go out of order
  orderedTimes <- sort(obsTimes)
  # left end of interval
  minIdx <- which(orderedTimes < eventTime)
  if (length(minIdx) == 0) {
    minTime <- 0
  } else {
    minTime <- orderedTimes[max(minIdx)]
  }
  # right end of interval
  maxIdx <- which(orderedTimes >= eventTime)
  if (length(maxIdx) == 0) {
    maxTime <- 999
  } else {
    maxTime <- orderedTimes[min(maxIdx)]
  }
  
  return(c(minTime, maxTime))
}