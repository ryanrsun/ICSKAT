#' ICsingleSNP.R
#'
#' Same as ICSKAT() except do a separate score test for each SNP in gMat, one at a time.
#'
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param gMat n*p genotype matrix.
#' @param null_beta (p+nknots+2)*1 vector of coefficients for null model.
#' @param Itt (p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.
#'
#' @return A list with the elements:
#' \item{testStatsVec}{p*1 vector of score test statistics}
#' \item{pVec}{p*1 vector of score test p-values}
#'
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' null_fit <- skat_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#' ICsingleSNP(left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind=rep(1, n),
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$beta_fit, Itt=null_fit$Itt, 
#' gMat=matrix(data=rbinom(n=200*10, size=2, prob=0.3), nrow=200))
#'
ICsingleSNP <- function(left_dmat, right_dmat, lt, rt, obs_ind, tpos_ind, gMat, null_beta, solveItt, p) {

    # You need these terms for all the SNPs
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

    # multiply each SNP vector by this term
    Ug_term1 <- tpos_ind * exp(-H_L) * -H_L 
    Ug_term2 <- obs_ind * exp(-H_R) * -H_R
    Ug_term2[which(is.na(Ug_term2))] <- 0
    UgTerm <- (Ug_term1 - Ug_term2) / A
		
	  # The Igg term
	  ggTerm1 <- tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A)
	  ggTerm2 <- obs_ind * as.numeric(H_R * exp(-H_R) - H_R^2 * exp(-H_R)) / A
	  ggTerm2[which(is.nan(ggTerm2))] <- 0
	  ggTerm3 <- ( (tpos_ind * as.numeric((H_L * exp(-H_L))) - obs_ind * as.numeric((H_R * exp(-H_R)))) / A )^2
	  ggTerm <- ggTerm1 + ggTerm2 - ggTerm3
	
    # off-diagonals for Igtheta
	  gtTermL <- tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A) + 
	    tpos_ind * as.numeric((H_L * exp(-H_L)) / A) * (tpos_ind * as.numeric((-H_L * exp(-H_L)) / A) + obs_ind * as.numeric((H_R * exp(-H_R)) / A))
	  gtTermR <- obs_ind * ggTerm2 - obs_ind * as.numeric((H_R * exp(-H_R)) / A) * 
	    (tpos_ind * as.numeric((-H_L * exp(-H_L)) / A) + obs_ind * as.numeric((H_R * exp(-H_R)) / A))
	  gtTermCombined <- gtTermL + gtTermR
	  nKnots <- ncol(left_dmat) - p
	  gtTermCommon <- sweep(left_dmat[, 1:p], 1, gtTermCombined, FUN="*")
	  gtHalfL <- sweep(left_dmat[, (p+1):(p+nKnots)], 1, gtTermL, FUN="*")
	  gtHalfR <- sweep(right_dmat[, (p+1):(p+nKnots)], 1, gtTermR, FUN="*")
	
    # apply across all columns of gMat
	  scoreStatsOutput <- apply(gMat, 2, calcScoreStats, UgTerm=UgTerm, ggTerm=ggTerm, 
	                            gtTermCommon=gtTermCommon, gtHalfL=gtHalfL, gtHalfR=gtHalfR, 
	                            solveItt=solveItt)

    return(list(testStats = scoreStatsOutput[1, ], pvals = scoreStatsOutput[2, ]))
}
