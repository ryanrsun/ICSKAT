#' ICSKAT_fit_null_PO.R
#'
#' Fit the null model (cubic basis spline for baseline cumulative hazard and coefficients
#' for non-genetic coefficients) for interval-censored skat with PO model.
#'
#' @param init_beta (p+nknots+2)*1 vector of coefficients to initialize the Newton-Raphson.
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param lt n*1 vector of left side of interval times.
#' @param rt n*1 vector of right side of interval times.
#' @param eps Stop when the L2 norm of the difference in model coefficients reaches this limit.
#' @param checkpoint Boolean tells the function to print when each iteration completes.
#'
#' @return A list with the elements:
#' \item{beta_fit}{(p+nknots+2)*1 vector of fitted coefficients under null model.}
#' \item{iter}{Number of iterations needed to converge.}
#' \item{Itt}{Fisher information matrix for the fitted coefficients.}
#' \item{diff_beta}{Difference between beta_fit and previous iteration of the vector, can be checked for errors.}
#' \item{err}{err=1 if NA shows up in the calculation.}
#' \item{IterrMsg}{Describes the error.}
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
#' dmats <- make_IC_dmat(xMat = xMat, lt = lt, rt = rt, obs_ind = obs_ind,
#' tpos_ind = tpos_ind)
#' ICSKAT_fit_null_PO(init_beta = rep(0.1, 5), left_dmat = dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
ICSKAT_fit_null_PO <- function(init_beta, left_dmat, right_dmat, obs_ind, tpos_ind, lt, rt, checkpoint=FALSE, eps=10^(-6)) {

  diff_beta <- 1
  iter <- 0
  temp_beta <- init_beta

  stopSolve <- FALSE
  while( !is.nan(diff_beta) & stopSolve == FALSE & diff_beta<1000) {

    etaL <- left_dmat %*% temp_beta
    etaR <- right_dmat %*% temp_beta

    # Survival term
    SL <- ifelse(tpos_ind == 0, 1, rje::expit(-etaL))
    SR <- ifelse(obs_ind == 0, 0, rje::expit(-etaR))
    A <- SL - SR
    # sometimes A is 0
    A[which(A == 0)] <- min(A[which(A > 0)])

    # expit(-eta) = S(t) = 1 when t=0, so 1 - expit(-eta) = 1 - S(t) = 0
    # we could also use -expit(-etaL) * (1 - expit(-etaL))
    dGinvLdEta <- ifelse(tpos_ind == 0, 0, -SL * (1 - SL))
    dGinvRdEta <- ifelse(obs_ind == 0, 0, -SR * (1 - SR))

    # second derivative terms
    # again, there's an S(T) and a 1 - S(T) in these terms
    # we could also use (exp(-etaL) - exp(-2 * etaL)) / (1 + exp(-etaL))^3
    d2GinvLdEta2 <- ifelse(tpos_ind == 0, 0, 2 * SL^3 - 3 * SL^2 + SL)
    d2GinvRdEta2 <- ifelse(obs_ind == 0, 0, 2 * SR^3 - 3 * SR^2 + SR)

    # to sweep, just divide by A
    UsweepL <- ifelse(tpos_ind == 0, 0, dGinvLdEta / A)
    UsweepR <- ifelse(obs_ind == 0, 0, dGinvRdEta / A)
    # the score vector
    uVec <- rowSums( sweep(t(left_dmat), MARGIN = 2, STATS = UsweepL, FUN = "*") -
                       sweep(t(right_dmat), MARGIN = 2, STATS = UsweepR, FUN = "*") )


    # Remember that A^TB != B^TA
    IttSweepL <- ifelse(tpos_ind == 0, 0, d2GinvLdEta2 / A)
    IttSweepR <- ifelse(obs_ind == 0, 0, d2GinvRdEta2 / A)
    IttL <- sweep(t(left_dmat), 2, IttSweepL, FUN="*") %*% left_dmat
    IttR <- sweep(t(right_dmat), 2, IttSweepR, FUN="*") %*% right_dmat
    IttBothHalf <- sweep(t(left_dmat), 2, UsweepL, FUN="*")  - sweep(t(right_dmat), 2, UsweepR, FUN="*")
    IttBoth <- IttBothHalf %*% t(IttBothHalf)
    iMat <- IttL - IttR - IttBoth

    # sometimes fullImat is singular
    solvedImat <- tryCatch(solve(iMat), error=function(e) e)
    if (class(solvedImat)[1] %in% c("simpleError")) {
      return(list(beta_fit=NA, iter=iter, diff_beta=diff_beta, Itt=NA, err=1, errMsg = "iMat singular, try different initial values"))
    }
    if (length(which(is.na(solvedImat))) > 0)
    {
      return(list(beta_fit=NA, iter=iter, diff_beta=diff_beta, Itt=NA, err=1, errMsg="iMat inverse has NaN, try different initial values"))
    }

    # update
    beta_new <- temp_beta - t(uVec) %*% solvedImat
    diff_beta <- (beta_new - temp_beta) %*% t(beta_new - temp_beta)
    temp_beta <- as.numeric(beta_new)
    iter <- iter + 1

    # stop? the second clause checks if we're stuck at a local max, if so then keep going.
    stopSolve <- ifelse(diff_beta < eps & sum(uVec^2) < length(uVec), TRUE, FALSE)

    # if too many iterations, stop
    if (iter > 100) {
      return(list(beta_fit=NA, iter=iter, diff_beta=diff_beta, Itt=NA, err=1, errMsg="Too many iterations, try different initial values"))
    }

    # checkpoint prints iterations
    if(checkpoint) {cat("iter ", iter, "diff", diff_beta, "\n")}
  }
  return(list(beta_fit=as.numeric(beta_new), iter=iter, Itt=iMat, diff_beta=diff_beta, err=0, errMsg=""))
}
