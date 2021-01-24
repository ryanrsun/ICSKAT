#' ICSKAT.R
#'
#' Calculate the test statistic and p-value for interval-censored skat.
#'
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param gMat n*q genotype matrix.
#' @param null_beta (p+nknots+2)*1 vector of coefficients for null model.
#' @param Itt (p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.
#' @param Itt (p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.
#' @param pvalue Boolean, if TRUE then find the p-value (maybe don't need it if bootstrapping, saves eigendecomposition)
#' @return A list with the elements:
#' \item{p_SKAT}{ICSKAT p-value}
#' \item{p_burden}{IC burden test p-value}
#' \item{complex}{Indicator of whether the SKAT variance matrix was positive definite}
#' \item{sig_mat}{The covariance matrix of the score equations for genetic effects when treated as fixed effects}
#' \item{skatQ}{SKAT test statistic}
#' \item{burdenQ}{Burden test statistic}
#' \item{lambdaQ}{Vector of eigenvalues of variance matrix}
#' \item{null_beta}{The fitted null parameters}
#' \item{err}{Will be 0 for no error, 22 if had to adjust parameters on CompQuadForm (totally normal), or 99 if NA in variance matrix. ICSKATwrapper will return 1 here if the null fit has an error}
#' \item{errMsg}{Explains error code, blank string if no error}
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' null_fit <- skat_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#' ICskat(left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind=rep(1, n),
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$beta_fit, Itt=null_fit$Itt,
#' gMat=matrix(data=rbinom(n=200*10, size=2, prob=0.3), nrow=200))
#'
ICskat <- function(left_dmat, right_dmat, lt, rt, obs_ind, tpos_ind, gMat, null_beta, Itt, pvalue=TRUE) {

  # Cumulative hazard under null
  H_L <- exp(left_dmat %*% null_beta)
  H_R <- exp(right_dmat %*% null_beta)
  # Sometimes H_R goes to infinity
  infHR <- which(H_R == Inf)
  if (length(infHR) > 0) {H_R[infHR] <- max( max(H_R[which(H_R < Inf)]), 10 )}
  # Survival term
  SL <- ifelse(lt == 0, 1, exp(-H_L))
  SR <- ifelse(rt == 999, 0, exp(-H_R))
  SR[!is.finite(SR)] <- 0
  A <- SL - SR
  # sometimes A is 0
  A[which(A == 0)] <- min(A[which(A > 0)])

  # just sweep once
  Ug_sweep1 <- tpos_ind * exp(-H_L) * -H_L
  Ug_sweep2 <- obs_ind * exp(-H_R) * -H_R
  Ug_sweep2[which(is.na(Ug_sweep2))] <- 0
  Ug_sweepTerm <- (Ug_sweep1 - Ug_sweep2) / A
  Ugamma <- rowSums(sweep(t(gMat), 2, Ug_sweepTerm, FUN="*"))
	# rm to save RAM
	rm(Ug_sweep1)
	rm(Ug_sweep2)
	rm(Ug_sweepTerm)


	# The Igg term
	ggTerm1 <- tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A)
	ggTerm2 <- obs_ind * as.numeric(H_R * exp(-H_R) - H_R^2 * exp(-H_R)) / A
	ggTerm2[which(is.nan(ggTerm2))] <- 0
	ggTerm3 <- ( (tpos_ind * as.numeric((H_L * exp(-H_L))) - obs_ind * as.numeric((H_R * exp(-H_R)))) / A )^2
	ggTerm <- ggTerm1 + ggTerm2 - ggTerm3
	IggHalf <- sweep(gMat, 1, ggTerm, FUN="*")
	# have to do this to make it a double matrix
	gMat[1, 1] <- 1.0 * gMat[1, 1]
	Igg <- eigenMapMatMultCrossTwo(gMat, IggHalf)
	# rm to save ram
	rm(ggTerm1, ggTerm3, ggTerm, IggHalf)

  # off-diagonals for Igtheta
	gtTermL <- tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A) +
	  tpos_ind * as.numeric((H_L * exp(-H_L)) / A) * (tpos_ind * as.numeric((-H_L * exp(-H_L)) / A) + obs_ind * as.numeric((H_R * exp(-H_R)) / A))
	gtTermR <- obs_ind * ggTerm2 - obs_ind * as.numeric((H_R * exp(-H_R)) / A) *
	  (tpos_ind * as.numeric((-H_L * exp(-H_L)) / A) + obs_ind * as.numeric((H_R * exp(-H_R)) / A))
	gtHalfL <- sweep(left_dmat, 1, gtTermL, FUN="*")
	gtHalfR <- sweep(right_dmat, 1, gtTermR, FUN="*")
	Igt <- eigenMapMatMultCrossTwo(gMat, gtHalfL) + eigenMapMatMultCrossTwo(gMat, gtHalfR)
	# rm to save ram
	rm(gtTermL, gtTermR, ggTerm2, gtHalfL, gtHalfR)

  # we just need the Igg portion of the inverse
 	# make sure there are no NA
	if (length(which(is.na(Igg))) > 0 | length(which(is.na(Igt))) > 0 | length(which(is.na(Itt))) > 0) {
		return(list(lambdaQ=NA, p_SKAT=NA, p_burden=NA, skatQ=NA, Ugamma=Ugamma, null_beta = null_beta,
              burdenQ=NA, sig_mat = NA, complex=NA, err=99, errMsg="NA in variance matrix"))
	} 
	sig_mat <- (-Igg) - (-Igt) %*% solve(-Itt) %*% t(-Igt)
	if (length(which(is.na(sig_mat))) > 0) {
		return(list(lambdaQ=NA, p_SKAT=NA, p_burden=NA, skatQ=NA, Ugamma=Ugamma, null_beta = null_beta,
              burdenQ=NA, sig_mat = NA, complex=NA, err=99, errMsg="NA in variance matrix"))
  } 
  skatQ <- t(Ugamma) %*% Ugamma
  burdenQ <- (sum(Ugamma))^2

	errCode <- 0
	errMsg <- ""
  if (pvalue) {
    lambdaQ <- eigen(sig_mat)$values
    p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
 		# as noted in the CompQuadForm documentation, sometimes you need to play with acc or lim parameters
		if (!is.na(p_SKAT)) {
			if (p_SKAT > 1) {
				paramDF <- data.frame(expand.grid(lim = c(10000, 20000, 50000), acc=c(1e-7, 1e-6, 1e-5, 1e-4)))
				paramCounter <- 1
				while(p_SKAT > 1) {
					tempLim <- paramDF$lim[paramCounter]
					tempAcc <- paramDF$acc[paramCounter]
					p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=tempAcc, lim=tempLim)$Qq
					paramCounter <- paramCounter + 1
					if (paramCounter > nrow(paramDF)) {break}
				}
				errCode <- 22
				errMsg <- "Had to adjust parameters on CompQuadForm"
			}
		} 
		B_burden= burdenQ / sum(sig_mat);
    p_burden=1-pchisq(B_burden,df=1) 
	} else {
    pSKAT <- NA
    lambdaQ <- 1
    p_burden <- NA
  }

  return(list(lambdaQ=lambdaQ, p_SKAT=p_SKAT, p_burden=p_burden, skatQ=skatQ, Ugamma=Ugamma, null_beta = null_beta,
              burdenQ=burdenQ, sig_mat = sig_mat, complex=is.complex(lambdaQ), err=errCode, errMsg=errMsg))
}
