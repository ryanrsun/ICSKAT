#' fIntegrateLiu.R
#'
#' The integrand in the SKATO p-value when using Liu instead of Davies method, pass it to a numerical integration function like integrate().
#'
#' @param x Variable to be integrated over, can be a vector or scalar.
#' @param muK1 Mean of the mixture of chi-squares that are the first part of the kappa variable.
#' @param sigmaK1 Standard deviation of the mixture of chi-squares that are the first part of the kappa variable.
#' @param QrhoDF The data frame output from calling QrhoIC().
#' @param dfK1 The degrees of freedom from the approximated chi-square.
#'
#' @return The value of the integrand at x.
#'
#' @importFrom stats dchisq
#' @importFrom stats pchisq
#'
#' @export
fIntegrateLiu <- function(x, muK1, sigmaK1, QrhoDF, dfK1) {

  # needs to be able to take vectorized x
  nr <- nrow(QrhoDF)
  nx <- length(x)
  tauX <- QrhoDF$tauVec %*% t(x)
  allChoiceMat <- (QrhoDF$qMinVec - tauX) / (1 - QrhoDF$rhoVec)
  minVec <- apply(allChoiceMat, 2, min)

  # no need for looping because pchisq is vectorized, unlike CompQuadForm::davies().
  # here the new chi-square RV that we are approximating kappa by has its moments determined
  # by its degrees of freedom, which are matched to its kurtosis.
  deltaX <- sqrt(2 * dfK1) * (minVec - muK1) / sigmaK1 + dfK1
  retVec <- stats::pchisq(deltaX, df = dfK1) * stats::dchisq(x, df=1)
  return(retVec)
}
