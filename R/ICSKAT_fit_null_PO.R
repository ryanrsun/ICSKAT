#' ICSKAT_fit_null_PO.R
#'
#' Fit the null model (cubic basis spline for baseline cumulative hazard and coefficients
#' for non-genetic coefficients) for interval-censored skat with PO model.
#'
#' @param init_beta (p+nknots+2)*1 vector of coefficients to initialize the Newton-Raphson.
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param lt n*1 vector of left side of interval times.
#' @param rt n*1 vector of right side of interval times.
#' @param eps Stop when the L2 norm of the difference in model coefficients reaches this limit.
#' @param checkpoint Boolean tells the function to print when each iteration completes.
#'
#' @return A list with the elements:
#' \item{beta_fit}{(p+nknots+2)*1 vector of fitted coefficients under null model.}
#' \item{iter}{Number of iterations needed to converge.}
#' \item{Itt}{Fisher information matrix for the fitted coefficients.}
#' \item{diff_beta}{Difference between beta_fit and previous iteration of the vector, can be checked for errors.}
#' \item{err}{err=1 if NA shows up in the calculation.}
#' \item{IterrMsg}{Describes the error.}
#' @export
#' @examples
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
#' ICSKAT_fit_null_PO(init_beta = rep(0, 5), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
#' obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt
#'
ICSKAT_fit_null_PO <- function(init_beta, left_dmat, right_dmat, obs_ind, tpos_ind, lt, rt, ZL, ZR, xMat, checkpoint=FALSE, eps=10^(-6)) {

  diff_beta <- 1
  iter <- 0
  temp_beta <- init_beta

  stopSolve <- FALSE
  while( !is.nan(diff_beta) & stopSolve == FALSE & diff_beta<1000) {

    etaL <- left_dmat %*% temp_beta
    etaR <- right_dmat %*% temp_beta

    # Survival term
    SL <- ifelse(tpos_ind == 0, 1, rje::expit(-etaL))
    SR <- ifelse(obs_ind == 0, 0, rje::expit(-etaR))
    SLSR <- SL - SR

    # expit(-eta) = S(t) = 1 when t=0, so 1 - expit(-eta) = 1 - S(t) = 0
    dGinvLdEta <- ifelse(tpos_ind == 0, 0, -expit(-etaL) * (1 - expit(-etaL)))
    dGinvRdEta <- ifelse(obs_ind == 0, 0, -expit(-etaR) * (1 - expit(-etaR)))
    # put a 1 in here to remind that we've taken the ZL and ZR out of these two terms
    # these are all terms that will be put in the final matrix as sweep(t(Z), 2, toSweep) %*% Z
    dEtaLdAlpha <- 1
    dEtaRdAlpha <- 1

    # second derivative terms
    # again, there's an S(T) and a 1 - S(T) in these terms
    d2GinvLdEta2 <- ifelse(tpos_ind == 0, 0, (exp(-etaL) - exp(-2 * etaL)) / (1 + exp(-etaL))^3)
    d2GinvRdEta2 <- ifelse(obs_ind == 0, 0, (exp(-etaR) - exp(-2 * etaR)) / (1 + exp(-etaR))^3)

    # dEtaLdAlpha is 1
    UalphaSweepL <- ifelse(tpos_ind == 0, 0, (dGinvLdEta * dEtaLdAlpha) / SLSR)
    UalphaSweepR <- ifelse(obs_ind == 0, 0, (dGinvRdEta * dEtaRdAlpha) / SLSR)
    # second part of score vector
    Ualpha <- rowSums( sweep(t(ZL), MARGIN = 2, STATS = UalphaSweepL, FUN = "*") -
                         sweep(t(ZR), MARGIN = 2, STATS = UalphaSweepR, FUN = "*") )

    # Remember that A^TB != B^TA
    IaaSweepL <- ifelse(tpos_ind == 0, 0, (d2GinvLdEta2 * dEtaLdAlpha^2)/ SLSR) - UalphaSweepL^2
    IaaSweepR <- ifelse(obs_ind == 0, 0, -(d2GinvRdEta2 * dEtaRdAlpha^2) / SLSR) - UalphaSweepR^2
    IaaSweepBoth <- UalphaSweepL * UalphaSweepR

    IaaL <- sweep(t(ZL), 2, IaaSweepL, FUN="*") %*% ZL
    IaaR <- sweep(t(ZR), 2, IaaSweepR, FUN="*") %*% ZR
    IaaBoth <- sweep(t(ZL), 2, IaaSweepBoth, FUN="*") %*% ZR + sweep(t(ZR), 2, IaaSweepBoth, FUN="*") %*% ZL
    # works, compared with iMat from ICSKAT_fit_null
    Iaa <- IaaL + IaaR + IaaBoth

    # first part of score vector
    UbetaSweep <- (dGinvLdEta - dGinvRdEta) / SLSR
    Ubeta <- rowSums( sweep(t(xMat), MARGIN = 2, STATS = UbetaSweep, FUN = "*") )

    # Ibb
    IbbSweep <- (d2GinvLdEta2 - d2GinvRdEta2) / SLSR - UbetaSweep^2
    Ibb <- sweep(t(xMat), 2, IbbSweep, FUN = "*") %*% xMat

    # Iba
    IbaSweepL <- ifelse(tpos_ind == 0, 0, d2GinvLdEta2 * dEtaLdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvLdEta * dEtaLdAlpha / SLSR^2)
    IbaSweepR <- ifelse(obs_ind == 0, 0, d2GinvRdEta2 * dEtaRdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvRdEta * dEtaRdAlpha / SLSR^2)
    IbaL <- sweep(t(xMat), 2, IbaSweepL, FUN = "*") %*% ZL
    IbaR <- sweep(t(xMat), 2, IbaSweepR, FUN = "*") %*% ZR
    Iba <- IbaL - IbaR

    # full iMat for "fixed effects"
    fullImat <- rbind(cbind(Ibb, Iba), cbind(t(Iba), Iaa))
    # full score vector
    fullU <- c(Ubeta, Ualpha)

    # update
    # sometimes fullImat is singular
    solvedImat <- tryCatch(solve(fullImat), error=function(e) e)
    if (class(solvedImat)[1] %in% c("simpleError")) {
      return(list(beta_fit=NA, iter=iter, Itt=NA, err=1, errMsg = "fullImat singular, try different initial values"))
    }

    beta_new <- temp_beta - t(fullU) %*% solve(fullImat)
    diff_beta <- (beta_new - temp_beta) %*% t(beta_new - temp_beta)
    temp_beta <- as.numeric(beta_new)
    iter <- iter + 1
    if (checkpoint) {cat("iter ", iter, "diff ", diff_beta, "\n")}

    # stop?
    stopSolve <- ifelse(diff_beta < eps & sum(fullU^2) < length(fullU), TRUE, FALSE)
  }

  return(list(beta_fit=beta_new, iter=iter, Itt=fullImat, err=0, errMsg=""))
}
