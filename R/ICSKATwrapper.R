#' ICSKATwrapper.R
#'
#' Wrapper to fit the null model and run ICSKAT all in one instead of separately - offers some functionality
#' for error checking or using different initial values when fit fails to converge.
#'
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param initValues (p+nknots+2)*1 vector of coefficients to initialize the Newton-Raphson.
#' @param lt Left side of interval times.
#' @param rt Right side of interval times.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param gMat n*q matrix of genotypes.
#' @param PH Boolean for whether to fit PH model (TRUE) or PO model (FALSE).
#' @param nKnots Number of knots in the spline.
#' @param maxIter Number of times to try the fit if initial values do not lead to convergence.
#' @param eps Difference in L2 norm of fitted null coefficients that stops the Newton Raphson.
#' @param runOnce Boolean, if true then just go through the algorithm once with the initial 
#' values for coefficients, updating the variance matrix, useful for bootstrapping.
#' @param returnNull Return a list with the skat output and null model, or just return the skat output (FALSE).
#'
#' @return Either a list with skatOutput and nullFit (two lists), or just skatOutput.
#'
#' @export
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
#' ICSKATwrapper(left_dmat = dmats$left_dmat, right_dmat = dmats$right_dmat, initValues = rep(0, ncol(xMat) + 3),
#' lt = lt, rt = rt, obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, returnNull = TRUE)
#'
ICSKATwrapper <- function(left_dmat, right_dmat, initValues, lt, rt, obs_ind, tpos_ind, gMat, 
                          PH=TRUE, nKnots=1, maxIter=3, eps=10^(-6), runOnce = FALSE, returnNull = FALSE) {

  xMat <- left_dmat[, 1:(ncol(left_dmat) - nKnots - 2)]	
  counter <- 0
  pass <- FALSE
  while (counter < maxIter) {

    counter <- counter + 1
	
    # initial values 
    if (counter == 1) {
      init_beta <- initValues
    } else {
      init_beta <- runif(n=ncol(left_dmat), min = -1, max = 1)
    }

    # null fit
    if (PH) {
      nullFit <- ICSKAT_fit_null(init_beta=init_beta, lt=lt, rt=rt,
                                 left_dmat=left_dmat, right_dmat=right_dmat,
                                 obs_ind=obs_ind, tpos_ind=tpos_ind, eps=eps, runOnce=runOnce)
    }	else {
      nullFit <- ICSKAT_fit_null_PO(init_beta=init_beta,
                                    lt=lt, rt=rt,
                                    ZL=left_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)],
                                    ZR=right_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)], xMat = xMat,
                                    obs_ind=obs_ind, tpos_ind=tpos_ind, eps=eps, runOnce=runOnce)
    }

    # if null fit not good, go to next loop, unless runOnce
    if ( (nullFit$err == 1 | nullFit$diff_beta > eps) & runOnce == FALSE) {
      next
    }
	
    # test
    if (PH) {
      skatOutput <- ICskat(left_dmat=left_dmat, tpos_ind=tpos_ind, obs_ind=obs_ind,
                           right_dmat=right_dmat, gMat=gMat, lt=lt, rt=rt,
                           null_beta=as.numeric(nullFit$beta_fit), Itt=nullFit$Itt)
    } else {
      skatOutput <- ICskatPO(ZL=left_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)],
                             ZR=right_dmat[, (ncol(xMat) + 1):(ncol(xMat) + nKnots + 2)], xMat = xMat,
                             lt=lt, rt=rt,
                             obs_ind=obs_ind, tpos_ind=tpos_ind, gMat=gMat, nullCoef=nullFit$nullCoef, Itt = nullFit$Itt)
    }

    # if worked or runOnce, break
    # if asking for runOnce, may return with unexpected errors
    if ( skatOutput$err == 0 | skatOutput$err == 22 | runOnce == TRUE ) {
      pass <- TRUE	
      break
    }
  } # end while loop

  # did it work
  if (!pass) {
    if (nullFit$err == 1 | nullFit$diff_beta > eps) {
      skatOutput <- list(p_SKAT=NA, p_burden=NA, complex=NA, err=1, errMsg="Failed null fit")
    } else {
      # if it failed at testing, then just let it return the testing error
      a <- 1
      #skatOutput <- list(p_SKAT=NA, p_burden=NA, complex=NA, err=1, errMsg="Failed testing")
    }	
  }

  # return
  if (returnNull) {
    return(list(skatOutput = skatOutput, nullFit = nullFit))
  } else {
    return(skatOutput)
  }

}

