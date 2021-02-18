#' ICSKAT_fit_null.R
#'
#' Fit the null model (cubic basis spline for baseline cumulative hazard and coefficients
#' for non-genetic coefficiens) for interval-censored skat.
#'
#' @param init_beta (p+nknots+2)*1 vector of coefficients to initialize the Newton-Raphson.
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param lt n*1 vector of left interval times.
#' @param rt n*1 vector of right interval times.
#' @param runOnce Boolean tells the function to just go through the loop once instead of converging (to get quantites for bootstrapping).
#' @param eps Stop when the L2 norm of the difference in model coefficients reaches this limit.
#' @param checkpoint Boolean tells the function to print when each iteration completes.
#'
#' @return A list with the elements:
#' \item{beta_fit}{(p+nknots+2)*1 vector of fitted coefficients under null model.}
#' \item{iter}{Number of iterations needed to converge.}
#' \item{Itt}{Fisher information matrix for the fitted coefficients.}
#' \item{diff_beta}{Difference between beta_fit and previous iteration of the vector, can be checked for errors}
#' \item{err}{Value is 0 if no errors and 1 if Itt is singular, can't perform fit}
#' \item{err}{Empty string if err=0, explains error if there is one}
#'
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' ICSKAT_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat, lt=lt, rt=rt,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#'
ICSKAT_fit_null <- function(init_beta, left_dmat, obs_ind, tpos_ind, right_dmat, lt, rt, runOnce=FALSE, checkpoint=FALSE, eps=10^(-6)) {

  diff_beta <- 1
  iter <- 0
  temp_beta <- init_beta
  while( !is.nan(diff_beta) & diff_beta > eps & diff_beta<1000) {
    # Cumulative hazard under null
    H_L <- exp(left_dmat %*% temp_beta)
    H_R <- exp(right_dmat %*% temp_beta)
    # Survival term
    SL <- ifelse(lt == 0, 1, exp(-H_L))
    SR <- ifelse(rt == 999, 0, exp(-H_R))
    SR[!is.finite(SR)] <- 0
    A <- SL - SR
    # sometimes A is 0
    A[which(A == 0)] <- min(A[which(A > 0)])

    # score vector
    U_term1 <- sweep(t(left_dmat), 2, ifelse(lt == 0, 0, exp(-H_L) * -H_L), FUN="*")
    U_term2 <- sweep(t(right_dmat), 2, ifelse(rt == 999, 0, exp(-H_R) * -H_R), FUN="*")
    U_term2[is.na(U_term2)] <- 0
    uVec <- rowSums( sweep(U_term1 - U_term2, 2, A, FUN="/") )

    # information matrix
    I_term1 <- crossprod(left_dmat, sweep(left_dmat, 1, tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A), FUN="*"))
    # sometimes H_R is so large that when it gets squared it goes to Inf and then multiplied
    # by exp(-H_R) it goes to NaN
    check_Iterm2 <- (H_R * exp(-H_R) - H_R^2 * exp(-H_R))
    check_Iterm2[which(is.nan(check_Iterm2))] <- 0
    I_term2 <- crossprod(right_dmat, sweep(right_dmat, 1, obs_ind * as.numeric(check_Iterm2 / A), FUN="*"))
    check_tempterm <- (H_R * exp(-H_R))
    check_tempterm[which(is.nan(check_tempterm))] <- 0
    temp_term <- sweep(t(left_dmat), 2, tpos_ind * as.numeric((H_L * exp(-H_L)) / A), FUN="*") -
    sweep(t(right_dmat), 2, obs_ind * as.numeric(check_tempterm / A), FUN="*")
    I_term3 <- temp_term %*% t(temp_term)
    iMat <- I_term1 + I_term2  - I_term3

    # sometimes iMat is singular
    solvedImat <- tryCatch(solve(iMat), error=function(e) e)
    if (class(solvedImat)[1] %in% c("simpleError")) {
      return(list(beta_fit=NA, iter=iter, diff_beta=diff_beta, Itt=NA, err=1, errMsg="iMat singular, try different initial values"))
    }
    if (length(which(is.na(solvedImat))) > 0)
    {
      return(list(beta_fit=NA, iter=iter, diff_beta=diff_beta, Itt=NA, err=1, errMsg="iMat inverse has NaN, try different initial values"))
    }
    beta_new <- temp_beta - t(uVec) %*% solve(iMat)
    diff_beta <- (beta_new - temp_beta) %*% t(beta_new - temp_beta)
    temp_beta <- as.numeric(beta_new)
    iter <- iter + 1
    if (runOnce) {break}

    # very rarely does it get stuck
    if (iter > 100) {
      return(list(beta_fit=NA, iter=iter, diff_beta=diff_beta, Itt=NA, err=1, errMsg="Too many iterations, try different initial values"))
    }

    if(checkpoint) {cat("iter ", iter, "\n")}
  }

  return(list(beta_fit=beta_new, iter=iter, diff_beta=diff_beta, Itt=iMat, err=0, errMsg=""))
}
