#' ICSKATPO.R
#'
#' Calculate the test statistic and p-value for interval-censored skat with PO model.
#'
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param gMat n*q genotype matrix.
#' @param null_beta (p+nknots+2)*1 vector of coefficients for null model.
#' @param Itt (p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.
#'
#' @return A list with the elements:
#' \item{p_SKAT}{ICSKAT p-value}
#' \item{p_burden}{IC burden test p-value}
#' \item{complex}{Indicator of whether the SKAT variance matrix was positive definite}
#' \item{sig_mat}{The covariance matrix of the score equations for genetic effects when treated as fixed effects}
#' \item{skatQ}{SKAT test statistic}
#' \item{burdenQ}{Burden test statistic}
#'
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' null_fit <- skat_fit_null_PO(init_beta=rep(0, 5), left_dmat=dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#' ICskatPO(left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind=rep(1, n),
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$nullCoef, Itt=null_fit$Itt,
#' gMat=matrix(data=rbinom(n=200*10, size=2, prob=0.3), nrow=200))
#'
ICskatPO <- function(ZL, ZR, xMat, lt, rt, obs_ind, tpos_ind, gMat, nullCoef, Itt) {

	# assume that you put the spline coefficients last  
  betaVec <- nullCoef[1:ncol(xMat)]
  alphaVec <- nullCoef[(ncol(xMat)+1):length(nullCoef)]

  # corrected
  etaL <- ZL %*% alphaVec + xMat %*% betaVec
  etaR <- ZR %*% alphaVec + xMat %*% betaVec
  
  # Survival term
  SL <- ifelse(lt == 0, 1, expit(-etaL))
  SR <- ifelse(rt == 999, 0, expit(-etaR))
  SLSR <- SL - SR
  
  # if lt=0 or rt=0, we swap 0 because -expit(etaL) is the survival and should be 0 or 1. 
  dGinvLdEta <- ifelse(lt == 0, 0, -expit(-etaL) * (1 - expit(-etaL)))
  dGinvRdEta <- ifelse(rt == 999, 0, -expit(-etaR) * (1 - expit(-etaR))) 
  # put a 1 in here to remind that we've taken the ZL and ZR out of these two terms
  # these are all terms that will be put in the final matrix as sweep(t(Z), 2, toSweep) %*% Z
  dEtaLdAlpha <- 1 
  dEtaRdAlpha <- 1 
  
  # second derivative terms
  # if lt=0 or rt=0, we swap 0 because -expit(etaL) is the survival and should be 0 or 1. 
  d2GinvLdEta2 <- ifelse(lt == 0, 0, (exp(-etaL) - exp(-2 * etaL)) / (1 + exp(-etaL))^3)  
  d2GinvRdEta2 <- ifelse(rt == 999, 0, (exp(-etaR) - exp(-2 * etaR)) / (1 + exp(-etaR))^3)
  
  # score vector
  UgammaSweep <- (dGinvLdEta - dGinvRdEta) / SLSR
  Ugamma <- rowSums( sweep(t(gMat), MARGIN = 2, STATS = UgammaSweep, FUN = "*") )

  # Igg
  IggSweep <- (d2GinvLdEta2 - d2GinvRdEta2) / SLSR - UgammaSweep^2
  Igg <- sweep(t(gMat), 2, IggSweep, FUN = "*") %*% gMat
  
  # Iga
  IgaSweepL <- ifelse(lt == 0, 0, d2GinvLdEta2 * dEtaLdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvLdEta * dEtaLdAlpha / SLSR^2)
  IgaSweepR <- ifelse(rt == 999, 0, d2GinvRdEta2 * dEtaRdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvRdEta * dEtaRdAlpha / SLSR^2)
  IgaL <- sweep(t(gMat), 2, IgaSweepL, FUN = "*") %*% ZL
  IgaR <- sweep(t(gMat), 2, IgaSweepR, FUN = "*") %*% ZR
  Iga <- IgaL - IgaR
  
  # Igb
  Igb <- sweep(t(gMat), 2, IggSweep, FUN = "*") %*% xMat
  
  # off diagonal
  offDiag <- cbind(Igb, Iga)
  sig_mat <- (-Igg) - (-offDiag) %*% solve(-Itt) %*% t(-offDiag)
  skatQ <- t(Ugamma) %*% Ugamma
  lambdaQ <- eigen(sig_mat)$values
  p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
  
  burdenQ <- (sum(Ugamma))^2
  B_burden= burdenQ / sum(sig_mat);
  p_burden=1-pchisq(B_burden,df=1)
  
  return(list(p_SKAT=p_SKAT, p_burden=p_burden, skatQ=skatQ, burdenQ=burdenQ, sig_mat = sig_mat, complex=is.complex(lambdaQ)))
}


