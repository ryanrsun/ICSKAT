#' make_IC_dmat.R
#'
#' Puts together the entire design matrix for both the left and right ends of the
#' interval, pasting together the non-genetic covariates with the cubic spline basis.
#'
#' @param xMat n*p matrix of non-genetic covariates.
#' @param lt n*1 vector with left end of intervals (min is 0).
#' @param rt n*1 vector with right end of intervals.
#' @param quant_r Quantiles of time to use in constructing the spline, pass in if doing bootstrap.
#' @param nknots Number of knots to use for cubic spline basis (default is 1).
#'
#' @return A list with the elements:
#' \item{right_dmat}{n*(p+nknots+2) design matrix for right end of interval.}
#' \item{left_dmat}{n*(p+nknots+2) design matrix for left end of interval.}
#' \item{quant_r}{Quantiles used for constructing spline.}
#'
#' @export
#' @examples
#' xMat <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' make_IC_dmat(xMat=xMat, lt=lt, rt=rt)
#'
make_IC_dmat <- function(xMat, lt, rt, quant_r=NULL, nknots=1) {
  
  # place the knots at equally spaced quantiles
  obs_ind <- as.numeric(rt < 999)
  if (is.null(quant_r)) {
    quant_r <- quantile(log(rt[obs_ind == 1]), probs=seq(from=0, to=1, length.out=nknots+2))
  }

  # a0 and a1 are always there
  right_a0 <- 1
  right_a1 <- log(rt)
  if (is.null(xMat)) {
    right_dmat <- cbind(right_a0, right_a1)
  } else {
    right_dmat <- cbind(xMat, right_a0, right_a1)
  }

  # if lt = 0, then the cumulative hazard is necessarily 0 so set it all to 0
  left_a0 <- ifelse(lt == 0, 0, 1)
  left_a1 <- ifelse(lt == 0, 0, log(lt))
  if (is.null(xMat)) {
    left_dmat <- cbind(left_a0, left_a1)
  } else {
    left_dmat <- cbind(xMat, left_a0, left_a1)
  }

  kmax <- max(quant_r)
  kmin <- min(quant_r)
  # place the knots
  for (j in 1:nknots) {
    ej <- (kmax - quant_r[j+1]) / (kmax - kmin)

    right_aj <-  ifelse(right_a1 == log(999), 999,
                        pmax(0, (right_a1-quant_r[j+1])**3) - ej*pmax(0, (right_a1-kmin)**3) -
                          (1-ej)*pmax(0, (right_a1-kmax)**3))
    right_dmat <- cbind(right_dmat, right_aj)

    left_aj <-ifelse(lt == 0, 0, pmax(0, (left_a1-quant_r[j+1])**3) - ej*pmax(0, (left_a1-kmin)**3) -
                             (1-ej)*pmax(0, (left_a1-kmax)**3))
    left_dmat <- cbind(left_dmat, left_aj)
  }

  return(list(right_dmat=right_dmat, left_dmat=left_dmat, quant_r = quant_r))
}

