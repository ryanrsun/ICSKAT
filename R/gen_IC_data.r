#' gen_IC_data.R
#'
#' Generate interval-censored data under the proportional odds/PH model given a baseline hazard function and
#' some information about observation times.
#'
#' @param bhFunInv A function, the inverse of the baseline hazard function.
#' @param obsTimes Vector of the intended observation times.
#' @param windowHalf The amount of time before or after the intended obsTimes that a visit might take place.
#' @param etaVec n*1 linear predictor in either the proportional odds or proportional hzards model.
#' @param mod Either "PH" to generate under PH model or "PO" to generate under PO model.
#' @param probMiss The probability of missing any given visit.
#'
#' @return n*2 matrix with the interval times for each subject.
#'
#' @export
#' @examples
#' bhFunInv <- function(x) {x}
#' obsTimes <- 1:5
#' xMat <- matrix(data=rnorm(200), nrow=100)
#' etaVec <- rep(0, 100)
#' outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, 
#' windowHalf = 0.1, etaVec = etaVec)
#'
gen_IC_data <- function(bhFunInv, obsTimes, windowHalf, etaVec, mod = "PH", probMiss=0.1) {

  # sample size
  n <- length(etaVec)

  # generate the exact times
  uVec <- runif(n=n, min=0, max=1)
  if (mod == "PH") {
      tVec <- bhFunInv( -log(1 - uVec) / exp(etaVec) )
  } else if (mod == "PO") {
      part1 <- rje::logit(uVec) - etaVec
      tVec <- bhFunInv( -log(exp(-part1) / (1 + exp(-part1))) )
  } else {
      stop("Bad model")
  }

  # 1 - probMiss is the chance of making it to the visit
  nVisits <- length(obsTimes)
  madeVisit <- matrix(data=rbinom(n=n*nVisits, size=1, prob=(1 - probMiss)), nrow=n, ncol=nVisits)
  
  # make sure there is at least one visit for each subject
  nMadeVisits <- apply(madeVisit, 1, sum)
  zeroVisits <- which(nMadeVisits == 0)
  while (length(zeroVisits) > 0) {
    madeVisit[zeroVisits, ] <- matrix(data=rbinom(n=length(zeroVisits) * nVisits, size=1, 
                                                  prob=(1 - probMiss)), nrow=length(zeroVisits), ncol=nVisits)
    nMadeVisits <- apply(madeVisit, 1, sum)
    zeroVisits <- which(nMadeVisits == 0)
  }
		
	# actual visit time is uniformly distributed around the intended obsTime, windowHalf on each side
  visitTime <- sweep(matrix(data=runif(n=n*nVisits, min=-windowHalf, max=windowHalf), nrow=n, ncol=nVisits),
                     MARGIN=2, STATS=obsTimes, FUN="+")

  # get all visits for each subject
  allVisits <- madeVisit * visitTime
  # make the interval for each subject
  allInts <- t(mapply(FUN=createInt, obsTimes = data.frame(t(allVisits)), eventTime=tVec))
  leftTimes <- allInts[, 1]
  rightTimes <- allInts[, 2]
  # event time indicators
  obs_ind <- ifelse(rightTimes == Inf, 0, 1)

  # return
  return(list(deltaVec = deltaVec, tVec = tVec, leftTimes = leftTimes, rightTimes = rightTimes, allVisits=allVisits))
}

#' Called by gen_IC_data() to turn the actual outcome times and observation times into interval-censored
#' outcomes for each subject. Apply this with mapply over a data.frame of visit times, pass in the exact times.
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
    maxTime <- Inf
  } else {
    maxTime <- orderedTimes[min(maxIdx)]
  }

  return(c(minTime, maxTime))
}
