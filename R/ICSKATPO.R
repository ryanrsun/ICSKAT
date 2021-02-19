#' ICSKATPO.R
#'
#' Calculate the test statistic and p-value for interval-censored skat with PO model.
#'
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param lt n*1 vector of left side of interval times.
#' @param rt n*1 vector of right side of interval times.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param gMat n*q genotype matrix.
#' @param null_beta (p+nknots+2)*1 vector of coefficients for null model.
#' @param Itt (p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.
#'
#' @return A list with the elements:
#' \item{p_SKAT}{ICSKAT p-value for PO model.}
#' \item{p_burden}{IC burden test p-value for PO model.}
#' \item{complex}{Indicator of whether the SKAT variance matrix was positive definite}
#' \item{sig_mat}{The covariance matrix of the score equations for genetic effects when treated as fixed effects}
#' \item{skatQ}{SKAT test statistic.}
#' \item{burdenQ}{Burden test statistic.}
#' \item{err}{err=1 for a bad null fit.}
#' \item{errMsg}{Describes the error.}
3''
#' @export
#' @examples
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
#' nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, 
#' obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
#' ICskatPO(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
#'
ICskatPO <- function(left_dmat, right_dmat, lt, rt, obs_ind, tpos_ind, gMat, null_beta, Itt) {

  # corrected
  etaL <- left_dmat %*% null_beta 
  etaR <- right_dmat %*% null_beta
  
  # Survival term
  SL <- ifelse(tpos_ind == 0, 1, expit(-etaL))
  SR <- ifelse(obs_ind == 0, 0, expit(-etaR))
  SR[!is.finite(SR)] <- 0
  SLSR <- SL - SR
  # sometimes A is 0
  # set it at 10^(-100) because you can square that and R will not round it to 0
  #SLSR[which(SLSR <= 0)] <- min(SLSR[which(SLSR > 0)])
  #SLSR[which(SLSR < 10^(-100))] <- 10^(-100)
  SLSR[which(SLSR == 0)] <- min(SLSR[which(SLSR > 0)])
 
  # if lt=0 or rt=999, we swap 0 because -expit(etaL) is the survival and should be 0 or 1. 
  dGinvLdEta <- ifelse(tpos_ind == 0, 0, -expit(-etaL) * (1 - expit(-etaL)))
  dGinvRdEta <- ifelse(obs_ind == 0, 0, -expit(-etaR) * (1 - expit(-etaR))) 
  # put a 1 in here to remind that we've taken the ZL and ZR out of these two terms
  # these are all terms that will be put in the final matrix as sweep(t(Z), 2, toSweep) %*% Z
  dEtaLdAlpha <- 1 
  dEtaRdAlpha <- 1 
  
  # second derivative terms
  # if lt=0 or rt=999, we swap 0 because -expit(etaL) is the survival and should be 0 or 1. 
  # also, if checkNumR is -Inf, that means the denominator is like -Inf^2, so that's why we set it to 0 so the 
  # quotient will be correctly 0.	
  checkNumR <- exp(-etaR) - exp(-2 * etaR)
  numR <- ifelse(exp(-2 * etaR) == Inf, 0, checkNumR)
  checkNumL <- exp(-etaL) - exp(-2 * etaL)
  numL <- ifelse(exp(-2 * etaL) == Inf, 0, checkNumL)
  d2GinvLdEta2 <- ifelse(tpos_ind == 0, 0, numL / (1 + exp(-etaL))^3)  
  d2GinvRdEta2 <- ifelse(obs_ind == 0, 0, numR / (1 + exp(-etaR))^3)
 
  # score vector
  # this works fine with minimal checking because when SL or SR is very small, both the numerator and denominator are on the
  # same scale as long as SLSR is not allowed to be 0.
  UgammaSweep <- (dGinvLdEta - dGinvRdEta) / SLSR
  Ugamma <- rowSums( sweep(t(gMat), MARGIN = 2, STATS = UgammaSweep, FUN = "*") )
  if (length(which(is.na(UgammaSweep))) > 0 | length(which(is.na(Ugamma))) > 0) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  }

  # Igg
  #IggSweepCheck <- (d2GinvLdEta2 - d2GinvRdEta2) / SLSR - UgammaSweep^2
  #IggMax <- min( max(abs(IggSweepCheck[which(abs(IggSweepCheck) != Inf)])), 10^20 )
  #IggSweep <- ifelse(abs(IggSweepCheck) > IggMax, IggMax * sign(IggSweepCheck), IggSweepCheck)
  IggSweep <- (d2GinvLdEta2 - d2GinvRdEta2) / SLSR - UgammaSweep^2	
  Igg <- sweep(t(gMat), 2, IggSweep, FUN = "*") %*% gMat
  if (length(which(is.na(Igg))) > 0) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  }
  
  # Iga
  # we can get in trouble when SL - SR is very small, in that case the second term is kind of like dGinvLdEta / SLSR
  #checkIgaSweepL <- ifelse(lt == 0, 0, d2GinvLdEta2 * dEtaLdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvLdEta * dEtaLdAlpha / SLSR^2)
  #IgaSweepL1 <- ifelse(checkIgaSweepL > 10^100 | is.na(checkIgaSweepL), d2GinvLdEta2 * dEtaLdAlpha / SLSR -  dGinvLdEta * dEtaLdAlpha / SLSR, checkIgaSweepL)
  #IgaSweepL <- ifelse(abs(IgaSweepL1) > 10^30, sign(IgaSweepL1) * 10^30, IgaSweepL1) 	
  #checkIgaSweepR <- ifelse(rt == 999, 0, d2GinvRdEta2 * dEtaRdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvRdEta * dEtaRdAlpha / SLSR^2)
  #IgaSweepR1 <- ifelse(checkIgaSweepR > 10^100 | is.na(checkIgaSweepR), d2GinvRdEta2 * dEtaRdAlpha / SLSR - dGinvRdEta * dEtaRdAlpha / SLSR, checkIgaSweepR)
  #IgaSweepR	<- ifelse(abs(IgaSweepR1) > 10^30, sign(IgaSweepR1) * 10^30, IgaSweepR1)
  IgaSweepL <- ifelse(tpos_ind == 0, 0, d2GinvLdEta2 * dEtaLdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvLdEta * dEtaLdAlpha / SLSR^2)
  IgaSweepR <- ifelse(obs_ind == 0, 0, d2GinvRdEta2 * dEtaRdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvRdEta * dEtaRdAlpha / SLSR^2)
  IgaL <- sweep(t(gMat), 2, IgaSweepL, FUN = "*") %*% ZL
  IgaR <- sweep(t(gMat), 2, IgaSweepR, FUN = "*") %*% ZR
  Iga <- IgaL - IgaR
  
  # Igb
  Igb <- sweep(t(gMat), 2, IggSweep, FUN = "*") %*% xMat
 
  # remove to save RAM?
  rm(UgammaSweep)
  rm(IggSweep)
  rm(IgaSweepL)
  rm(IgaSweepR)
  rm(IgaL)
  rm(IgaR) 

  # off diagonal
  offDiag <- cbind(Igb, Iga)
  if (length(which(is.na(offDiag))) > 0) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  }
  sig_mat <- (-Igg) - (-offDiag) %*% solve(-Itt) %*% t(-offDiag)
  skatQ <- t(Ugamma) %*% Ugamma
  # final check for bad fit
  if (length(which(is.na(sig_mat))) > 0 | is.na(skatQ)) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  } 
  lambdaQ <- eigen(sig_mat)$values
  p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
  
  burdenQ <- (sum(Ugamma))^2
  B_burden= burdenQ / sum(sig_mat);
  p_burden=1-pchisq(B_burden,df=1)
  
  return(list(p_SKAT=p_SKAT, p_burden=p_burden, skatQ=skatQ, burdenQ=burdenQ, sig_mat = sig_mat, complex=is.complex(lambdaQ), err=0, errMsg=""))
}


