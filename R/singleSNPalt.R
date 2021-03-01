#' singleSNPalt.R
#'
#' Take a matrix of SNPs and get the interval-censored regression p-value for each one separately using either
#' survreg() or coxph() with midpoint approximation.
#'
#' @param lt n*1 vector of left side of time interval.
#' @param rt n*1 vector of right side of time interval.
#' @param tpos_ind n*1 binary vector of whether the event was observed after follow-up started (takes value 1 if t>0, 0 otherwise).
#' @param obs_ind n*1 vector of whether the event was observed or right-censored (takes value 1 if observed or 0 if right-censored).
#' @param xMat non-SNP covariates matrix.
#' @param gMat n*q genotype matrix.
#' @param coxph Boolean, whether to fit Cox PH model.
#' @param survreg Boolean, whether to fit survreg() Wiibull model.
#'
#' @return A list with the elements:
#' \item{pvalCox}{q*1 vector of marginal SNP p-values with Cox model}
#' \item{pvalSurv}{q*1 vector of marginal SNP p-values with survreg Weibull model}
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
#' singleSNPalt(lt = lt, rt = rt,tpos_ind = tpos_ind, obs_ind = obs_ind, xMat = xMat, gMat = gMat
#'
singleSNPalt <- function(lt, rt, tpos_ind, obs_ind, xMat, gMat, coxph=TRUE, survreg=TRUE) {

  # make the midpoint approximation
  midTime <- ifelse(obs_ind == 0, lt, (lt + rt) / 2)
  midEvent <- obs_ind

  # make the times for Surv interval2 type
  leftTime2 <- ifelse(tpos_ind == 0, NA, lt)
  rightTime2 <- ifelse(obs_ind == 0, NA, rt)

  # now apply
  p <- ncol(xMat)
  if (coxph) {
    pvalCox <- apply(gMat, 2, coxphFn, xMat = xMat, midTime = midTime, midEvent = midEvent, p=p)
  } else {pvalCox <- NA}

  if (survreg) {
    pvalSurv <- apply(gMat, 2, survregFn, xMat = xMat, leftTime2 = leftTime2, rightTime2 = rightTime2, p=p)
  } else {pvalSurv <- NA}

  # return
  return(list(pvalCox = pvalCox, pvalSurv = pvalSurv))

}



#' survregFn
#'
#' Function to be applied over gMat to get p-values from survreg().
#'
#' @param x n*1 genotype vector.
#' @param xMat n*p matrix of non-genotype covariates.
#' @param leftTime2 n*1 vector of left interval times in the format of Surv() interval2 type, i.e NA for left or right censored observations.
#' @param rightTime2 n*1 vector of right interval times in the format of Surv() interval2 type, i.e NA for left or right censored observations.
#' @param p scalar, number of columns in xMat.
#'
#' @return A scalar p-value for testing the effect of the genotype in survreg() Weibull model.
#'
#' @export
#'
survregFn <- function(x, xMat, leftTime2, rightTime2, p) {
  tempMod <- survival::survreg(Surv(time = leftTime2, time2 = rightTime2, type="interval2") ~
                                xMat + x, dist="weibull")
  pval <- summary(tempMod)$table[p+2, 4]
  return(pval)
}

#' coxphFn
#'
#' Function to be applied over gMat to get p-values from coxPH().
#'
#' @param x n*1 genotype vector.
#' @param xMat n*p matrix of non-genotype covariates.
#' @param midTime n*1 vector of event times imputed to be right-censored times using the midpoint imputation method.
#' @param midEvent n*1 vector event indicators (0 for censoring, 1 for event) after times have been transformed to right-censored observations.
#' @param p scalar, number of columns in xMat.
#'
#' @return A scalar p-value for testing the effect of the genotype in survreg() Weibull model.
#'
#' @export
#'
coxphFn <- function(x, xMat, midTime, midEvent, p) {
  tempMod <- survival::coxph(Surv(time = midTime,  event=midEvent) ~ xMat + x)
  pval <- summary(tempMod)$coefficients[p+1, 5]
  return(pval)
}
