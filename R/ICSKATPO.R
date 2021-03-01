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
#'
#' @importFrom rje expit
#' @importFrom stats pchisq
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
#' nullFit <- ICSKAT_fit_null(init_beta = rep(0.1, 5), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
#' ICskatPO(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
#'
ICskatPO <- function(left_dmat, right_dmat, lt, rt, obs_ind, tpos_ind, gMat, null_beta, Itt) {

  # linear predictor
  etaL <- left_dmat %*% null_beta
  etaR <- right_dmat %*% null_beta

  # survival term
  SL <- ifelse(tpos_ind == 0, 1, rje::expit(-etaL))
  SR <- ifelse(obs_ind == 0, 0, rje::expit(-etaR))
  SR[!is.finite(SR)] <- 0
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

  # score vector
  UgammaSweep <- (dGinvLdEta - dGinvRdEta) / A
  Ugamma <- rowSums( sweep(t(gMat), MARGIN = 2, STATS = UgammaSweep, FUN = "*") )
  if (length(which(is.na(UgammaSweep))) > 0 | length(which(is.na(Ugamma))) > 0) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  }

  # Igg
  IggSweep <- (d2GinvLdEta2 - d2GinvRdEta2) / A - UgammaSweep^2
  Igg <- sweep(t(gMat), 2, IggSweep, FUN = "*") %*% gMat
  if (length(which(is.na(Igg))) > 0) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  }

  # Iga
  IgtSweepL <- ifelse(tpos_ind == 0, 0, d2GinvLdEta2  / A - (dGinvLdEta - dGinvRdEta) * dGinvLdEta  / A^2)
  IgtSweepR <- ifelse(obs_ind == 0, 0, d2GinvRdEta2  / A - (dGinvLdEta - dGinvRdEta) * dGinvRdEta  / A^2)
  IgtL <- sweep(t(gMat), 2, IgtSweepL, FUN = "*") %*% left_dmat
  IgtR <- sweep(t(gMat), 2, IgtSweepR, FUN = "*") %*% right_dmat
  Igt <- IgtL - IgtR

  # remove to save RAM?
  rm(UgammaSweep)
  rm(IggSweep)
  rm(IgtSweepL)
  rm(IgtSweepR)
  rm(IgtL)
  rm(IgtR)

  sig_mat <- (-Igg) - (-Igt) %*% solve(-Itt) %*% t(-Igt)
  skatQ <- t(Ugamma) %*% Ugamma
  # final check for bad fit
  if (length(which(is.na(sig_mat))) > 0 | is.na(skatQ)) {
    return(list(p_SKAT=NA, p_burden=NA, skatQ=NA, burdenQ=NA, sig_mat = NA, complex=NA, err=1, errMsg="bad null fit, try different starting values"))
  }
  lambdaQ <- eigen(sig_mat)$values
  p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq

  # errors?
  errCode <- 0
  errMsg <- ""
  if (p_SKAT > 1) {
    paramDF <- data.frame(expand.grid(lim = c(10000, 20000, 50000), acc=c(1e-7, 1e-6, 1e-5, 1e-4)))
    paramCounter <- 1
    # normal adjustment for davies function
    while(p_SKAT > 1) {
      tempLim <- paramDF$lim[paramCounter]
      tempAcc <- paramDF$acc[paramCounter]
      p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=tempAcc, lim=tempLim)$Qq
      paramCounter <- paramCounter + 1
      if (paramCounter > nrow(paramDF)) {break}
    }
    errCode <- 22
    errMsg <- "Had to adjust parameters on CompQuadForm."
  }

  # burden
  burdenQ <- (sum(Ugamma))^2
  B_burden = burdenQ / sum(sig_mat);
  p_burden = 1 - stats::pchisq(B_burden, df = 1)

  return(list(p_SKAT=p_SKAT, p_burden=p_burden, skatQ=skatQ, burdenQ=burdenQ, sig_mat = sig_mat, complex=is.complex(lambdaQ), err=errCode, errMsg=errMsg))
}


