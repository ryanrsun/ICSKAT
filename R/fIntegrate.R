#' fIntegrate.R
#'
#' The integrand in the SKATO p-value, pass it to a numerical integration function like integrate(), uses
#' Davies method instead of Liu to calculate the probability in the integrand.
#'
#' @param x Variable to be integrated over, can be a vector or scalar.
#' @param muK1 Mean of the mixture of chi-squares that are the first part of the kappa variable.
#' When we do bootstrap we often pass in the mean of the entire kappa, since the mean of zeta is supposed to be 0.
#' @param sigmaK1 Standard deviation of the entire kappa term.
#' @param sigmaZeta Standard deviation of the zeta part of the kappa variable.
#' @param kappaLambda Eigenvalues that weight the mixture of chi-squares that are the first part of the kappa variable.
#' @param QrhoDF The data frame output from calling QrhoIC().
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
#' fIntegrate(x=1, muK1=10, sigmaK1=2, sigmaZeta=2, kappaLambda = rep(1, 10), QrhoDF=QrhoDF)
#'
fIntegrate <- function(x, muK1, sigmaK1, sigmaZeta, kappaLambda, QrhoDF, kurtK1=NULL, matchKurt=FALSE) {

  # needs to be able to take vectorized x
  nr <- nrow(QrhoDF)
  nx <- length(x)
  tauX <- QrhoDF$tauVec %*% t(x)
  allChoiceMat <- (QrhoDF$qMinVec - tauX) / (1 - QrhoDF$rhoVec)
  minVec <- apply(allChoiceMat, 2, min)

  # loop
  retVec <- rep(0, nx)
  for (i in 1:nx) {
    tempMin <- minVec[i]

    if(tempMin > sum(kappaLambda) * 10^4){
      tempQuantile <- 0
    } else {
      # the davies deltaX is different from the liu deltaX, this is normal.
      # the difference here is in the variance terms, we are approximating a
      # centered and scaled kappa by a centered and scaled sum of chi-square RVs that
      # are multiplied by the eigenvalues of the first part of kappa. We could arguably
      # replace the sigmaK1^2 and + muK1 terms by their values calculated based on kappaLambda. 
      deltaX <- (tempMin - muK1) * sqrt(sigmaK1^2 - sigmaZeta^2) / sigmaK1 + muK1
      daviesOutput <- CompQuadForm::davies(q=deltaX, lambda=kappaLambda, acc=10^(-6))
      Sx <- daviesOutput$Qq
      if (daviesOutput$ifault != 0) {stop("daviesOutput$ifault is not 0")}

      if (Sx > 1) {Sx <- 1}
      
      # approximate kappa by a chi-squared RV that has the kurtosis matched to kappa
      if (matchKurt) {
        deltaXkurt<- (tempMin - muK1) * sqrt(2 * dfK1) / sigmaK1 + dfK1
        daviesOutputKurt <- CompQuadForm::davies(q=deltaX, lambda=kappaLambda, acc=10^(-6))
        Sx <- daviesOutput$Qq
        SxKurt <- 
        if (Sx > 1) {Sx <- 1}
      }
    }
    retVec[i]<-(1-Sx) * dchisq(x[i],df=1, ncp=0)
  }
  return(retVec)
}
