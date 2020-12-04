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
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' null_fit <- skat_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#' myoutput <- ICskat(left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind=rep(1, n),
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$beta_fit, Itt=null_fit$Itt,
#' gMat=matrix(data=rbinom(n=200*10, size=2, prob=0.3), nrow=200))
#' QrhoIC(rhoVec = seq(from=0, to=1, by=0.1), icskatOut = myoutput)
#' fIntegrateLiu(x=1, muK1=10, sigmaK1=2, kappaLambda = rep(1, 10), QrhoDF=QrhoDF, dfK1 = 10)
#'
fIntegrateLiu <- function(x, muK1, sigmaK1, QrhoDF, dfK1) {

  # needs to be able to take vectorized x
  nr <- nrow(QrhoDF)
  nx <- length(x)
  tauX <- QrhoDF$tauVec %*% t(x)
  allChoiceMat <- (QrhoDF$qMinVec - tauX) / (1 - QrhoDF$rhoVec)
  minVec <- apply(allChoiceMat, 2, min)

  # no need for looping because pchisq is vectorized, unlike CompQuadForm::davies()
	deltaX <- sqrt(2 * dfK1) * (minVec - muK1) / sigmaK1 + dfK1
	retVec <- pchisq(deltaX, df = dfK1) * dchisq(x, df=1)
	return(retVec)
}
