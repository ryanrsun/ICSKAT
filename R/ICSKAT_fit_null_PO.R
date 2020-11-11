#' ICSKAT_fit_null_PO.R
#'
#' Fit the null model (cubic basis spline for baseline cumulative hazard and coefficients
#' for non-genetic coefficiens) for interval-censored skat with PO model.
#'
#' @param init_beta (p+nknots+2)*1 vector of coefficients to initialize the Newton-Raphson.
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param eps Stop when the L2 norm of the difference in model coefficients reaches this limit.
#' @param checkpoint Boolean tells the function to print when each iteration completes.
#'
#' @return A list with the elements:
#' \item{beta_fit}{(p+nknots+2)*1 vector of fitted coefficients under null model.}
#' \item{iter}{Number of iterations needed to converge.}
#' \item{Itt}{Fisher information matrix for the fitted coefficients.}
#'
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' ICSKAT_fit_null_PO(init_beta=rep(0, 5), left_dmat=dmats$left_dmat, lt=lt, rt=rt,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#'
ICSKAT_fit_null_PO <- function(init_beta, ZL, ZR, xMat, obs_ind, tpos_ind, lt, rt, checkpoint=FALSE, eps=10^(-6)) {
  
  diffCoef <- 1
  iter <- 0
  tempCoef <- init_beta
  
  stopSolve <- FALSE
  while( !is.nan(diffCoef) & stopSolve == FALSE & diffCoef<1000) {
   
		# assume that you put the spline cofficients last 
    betaVec <- tempCoef[1:ncol(xMat)]
    alphaVec <- tempCoef[(ncol(xMat)+1):length(tempCoef)]
    
    etaL <- ZL %*% alphaVec + xMat %*% betaVec
    etaR <- ZR %*% alphaVec + xMat %*% betaVec
    
    # Survival term
    SL <- ifelse(lt == 0, 1, expit(-etaL))
    SR <- ifelse(rt == 999, 0, expit(-etaR))
    SLSR <- SL - SR
    
    # understand why it's a 0 if lt == 0! (because expit(-eta) = S0 = 1 when t=0, so 1 - expit(-eta) = 0)
    dGinvLdEta <- ifelse(lt == 0, 0, -expit(-etaL) * (1 - expit(-etaL)))
    dGinvRdEta <- ifelse(rt == 999, 0, -expit(-etaR) * (1 - expit(-etaR))) 
    # put a 1 in here to remind that we've taken the ZL and ZR out of these two terms
    # these are all terms that will be put in the final matrix as sweep(t(Z), 2, toSweep) %*% Z
    dEtaLdAlpha <- 1 
    dEtaRdAlpha <- 1 
    
    # second derivative terms
    # understand why it's a 0 if lt == 0! (because again there's a 1- expit(-eta) term in there)
    d2GinvLdEta2 <- ifelse(lt == 0, 0, (exp(-etaL) - exp(-2 * etaL)) / (1 + exp(-etaL))^3)  
    d2GinvRdEta2 <- ifelse(rt == 999, 0, (exp(-etaR) - exp(-2 * etaR)) / (1 + exp(-etaR))^3)
    
    # dEtaLdAlpha is 1
    UalphaSweepL <- ifelse(lt == 0, 0, (dGinvLdEta * dEtaLdAlpha) / SLSR)
    UalphaSweepR <- ifelse(rt == 999, 0, (dGinvRdEta * dEtaRdAlpha) / SLSR)
    # second part of score vector
    Ualpha <- rowSums( sweep(t(ZL), MARGIN = 2, STATS = UalphaSweepL, FUN = "*") - 
                         sweep(t(ZR), MARGIN = 2, STATS = UalphaSweepR, FUN = "*") )
    
    # Remember that A^TB != B^TA 
    IaaSweepL <- ifelse(lt == 0, 0, (d2GinvLdEta2 * dEtaLdAlpha^2)/ SLSR) - UalphaSweepL^2
    IaaSweepR <- ifelse(rt == 999, 0, -(d2GinvRdEta2 * dEtaRdAlpha^2) / SLSR) - UalphaSweepR^2
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
    IbaSweepL <- ifelse(lt == 0, 0, d2GinvLdEta2 * dEtaLdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvLdEta * dEtaLdAlpha / SLSR^2)
    IbaSweepR <- ifelse(rt == 999, 0, d2GinvRdEta2 * dEtaRdAlpha / SLSR - (dGinvLdEta - dGinvRdEta) * dGinvRdEta * dEtaRdAlpha / SLSR^2)
    IbaL <- sweep(t(xMat), 2, IbaSweepL, FUN = "*") %*% ZL
    IbaR <- sweep(t(xMat), 2, IbaSweepR, FUN = "*") %*% ZR
    Iba <- IbaL - IbaR
    
    # full iMat for "fixed effects"
    fullImat <- rbind(cbind(Ibb, Iba), cbind(t(Iba), Iaa))
    # full score vector
    fullU <- c(Ubeta, Ualpha)
    
    # update
    newCoef <- tempCoef - t(fullU) %*% solve(fullImat)
    diffCoef <- (newCoef - tempCoef) %*% t(newCoef - tempCoef)
    tempCoef <- as.numeric(newCoef)
    iter <- iter + 1
    if (checkpoint) {cat("iter ", iter, "diff ", diffCoef, "\n")}
    
    # stop?
    stopSolve <- ifelse(diffCoef < eps & sum(fullU^2) < length(fullU), TRUE, FALSE)
  }
  
  return(list(beta_fit=tempCoef, iter=iter, Itt=fullImat))
}
