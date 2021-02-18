#' mixture_kurtosis.R
#'
#' Calculate the kurtosis of Qrho when performing SKATO with bootstrapped moments. This function is 
#' included to allow for the potential to match the SKAT package, however we generally don't call it
#' because we can just bootstrap the kurtosis of Qrho directly if we are already doing bootstrap,
#' thus avoiding this calculation. Also it's only used in calculating the qmin values, not in 
#' the final p-value calculation, which uses a kappa expression that is only the first two terms of Qrho.
#'
#' @param tempDF1 Generally the bootstrapped kurtosis of the mixture of chi-squares in kappa.
#' @param tempDF2 Generally 1 because it's for the chi-square1 RV in kappa.
#' @param v1 Generally the variance of the mixture of chi-squares plus the variance of zeta in kappa.
#' @param a1 Generally the 1-rho in front of the first part of the kappa term.
#' @param a2 Generally the tau(rho) term in front of the chi-square1 RV in kappa.
#' @return Kurtosis (excess kurtosis to be more precise), use df = 12 / kurtosis.
#'
#' @export
#' @examples
#' mixture_kurtosis(1, 1, 1, 1, 1)
#'
mixture_kurtosis <- function(tempDF1, tempDF2, v1, a1, a2) {
  # df1 is for the first quadratic form plus xi
  # df2 <- for the chisq_1
  # v1 is for the first quadratic form plus xi
  # a1 is 1-rho, a2 is tau(rho)
  v2 <- 2 * tempDF2

  # fourth central moments
  S4_1 <- (12 / tempDF1 + 3) * v1^2
  S4_2 <- (12 / tempDF2 + 3) * v2^2
  # unlike variance, which can be added for the first quadratic form plus xi + chisq1,
  # fourth central moment can't just be added
  S4 <- a1^4 * S4_1 + a2^4 * S4_2 + 6 * a1^2 * a2^2 * v1 * v2
  S2 <- a1^2 * v1 + a2^2 * v2
  tempKurt <- S4 / S2^2 - 3
  # same check as SKAT
  if (tempKurt < 0) {tempKurt <- 1e-4}

  return(tempKurt)
}

