#' chiSqMatchFast.R
#'
#' Match the moments of a mixture of scaled chi-square random variables to a single non-central chi-square,
#' assumes the quadratic form case where the mean of the multivariate normal V=RV is 0.
#'
#' @param lambdaVec Numeric vector holding the eigenvalues of the A term, where we are interested in x^TAX
#' and x is multivariate normal.
#' @param alwaysCentral Boolean determining whether to always set the noncentrality parameter to 0, as in SKAT package.
#'
#' @return A list with the elements:
#' \item{sigmaQrho}{Standard deviation of the mixture distribution}
#' \item{muQrho}{Mean of the mixture distribution}
#' \item{delta}{Noncentrality parameter of the matched distribution}
#' \item{l}{Degrees of freedom of the matched distribution}
#'
#' @export
#' @examples
#' set.seed(0)
#' gMat <- matrix(data=rbinom(n=200, size=2, prob=0.3), nrow=100)
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
#' nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
#' icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
#' kRho <- matrix(data=0.5, nrow=10, ncol=10)
#' diag(kRho) <- 1
#' toDecomp <- t(icskatOut$Ugamma) %*% kRho %*% icskatOut$Ugamma
#' tempEvals <- eigen(toDecomp, symmetric = TRUE, only.values = TRUE)$values
#' chiSqMatchFast(lambdaVec = tempEvals)
#'
chiSqMatchFast <- function(lambdaVec, alwaysCentral=FALSE) {
  c1 <- sum(lambdaVec)
  c2 <- sum(lambdaVec^2)
  c3 <- sum(lambdaVec^3)
  c4 <- sum(lambdaVec^4)
  muQrho <- c1
  sigmaQrho <- sqrt(2 * c2)
  s1 <- c3 / c2^(1.5)
  s2 <- c4 / c2^2
  # sometimes the eigenvalues are so small, c2 is too small and s2 is NaN
  if (is.na(s2)) {return(-1)}

  if (alwaysCentral | s2 >= s1^2) {
    delta <- 0
    l <- 1 / s2
  } else {
    a <- 1 / (s1 - sqrt(s1^2 - s2))
    delta <- s1 * a^3 - a^2
    l <- a^2 - 2 * delta
  }
  return(list(sigmaQrho = sigmaQrho, muQrho = muQrho, delta = delta, l = l))
}
