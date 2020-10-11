#' ICSKAT.R
#'
#' Calculate the test statistic and p-value for interval-censored skat.
#'
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param G n*q genotype matrix.
#' @param null_beta (p+nknots+2)*1 vector of coefficients for null model.
#' @param Itt (p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.
#'
#' @return A list with the elements:
#' \item{p_SKAT}{ICSKAT p-value.}
#' \item{p_burden}{IC burden test p-value.}
#' \item{complex}{Indicator of whether the SKAT variance matrix was positive definite}
#'
#' @export
#' @examples
#' X <- matrix(data=rnorm(200), nrow=100)
#' lt <- runif(n=100, min=0, max=5)
#' rt <- lt + runif(n=100, min=0, max=5)
#' dmats <- make_IC_dmat(X=X, lt=lt, rt=rt)
#' null_fit <- skat_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#' ICskat(left_dmat=dmats$left_dmat, right_dmat=dmats$right_dmat, obs_ind=rep(1, n),
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$beta_fit, Itt=null_fit$Itt)
#'
ICskat <- function(left_dmat, right_dmat, lt, rt, obs_ind, tpos_ind, G, null_beta, Itt) {

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

    # score vector for gamma
    Ug_term1 <- sweep(t(G), 2, ifelse(lt == 0, 0, exp(-H_L) * -H_L), FUN="*")
    Ug_term2 <- sweep(t(G), 2, ifelse(rt == 999, 0, exp(-H_R) * -H_R), FUN="*")
    Ug_term2[which(is.na(Ug_term2))] <- 0
    Ugamma <- rowSums( sweep(Ug_term1 - Ug_term2, 2, A, FUN="/") )
		# rm to save RAM
		rm(Ug_term1)
		rm(Ug_term2)

    # bottom right corner for information matrix corresponding to Igg
    Igg_term1 <- crossprod(G, sweep(G, 1, tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A), FUN="*"))
    # sometimes H_R is so large that when it gets squared it goes to Inf and then multiplied
    # by exp(-H_R) it goes to NaN
    check_Iggterm2 <- (H_R * exp(-H_R) - H_R^2 * exp(-H_R))
    check_Iggterm2[which(is.nan(check_Iggterm2))] <- 0
    Igg_term2 <- crossprod(G, sweep(G, 1, obs_ind * as.numeric(check_Iggterm2 / A), FUN="*"))
    temp_termgg <- sweep(t(G), 2, tpos_ind * as.numeric((H_L * exp(-H_L)) / A), FUN="*") -
        sweep(t(G), 2, obs_ind * as.numeric((H_R * exp(-H_R)) / A), FUN="*")
    Igg_term3 <- temp_termgg %*% t(temp_termgg)
    Igg <- Igg_term1 + Igg_term2  - Igg_term3
		# rm to save RAM
		rm(Igg_term1)
		rm(Igg_term2)
		rm(Igg_term3)

    # off-diagonals for Igtheta
    Igt_term1 <- crossprod(G, sweep(left_dmat, 1, tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A), FUN="*"))
    Igt_term2 <- crossprod(G, sweep(right_dmat, 1, obs_ind * as.numeric(check_Iggterm2 / A), FUN="*"))
    temp_term_gt1 <- sweep(t(left_dmat), 2, tpos_ind * as.numeric((H_L * exp(-H_L)) / A), FUN="*") -
        sweep(t(right_dmat), 2, obs_ind * as.numeric((H_R * exp(-H_R)) / A), FUN="*")
    temp_term_gt2 <- sweep(t(G), 2, tpos_ind * as.numeric((-H_L * exp(-H_L)) / A), FUN="*") +
        sweep(t(G), 2, obs_ind * as.numeric((H_R * exp(-H_R)) / A), FUN="*")
   	Igt_term3 <- temp_term_gt2 %*% t(temp_term_gt1) 
    Igt <- Igt_term1 + Igt_term2 + Igt_term3
		# rm to save RAM
		rm(Igt_term1)
		rm(Igt_term2)
		rm(temp_term_gt1)
		rm(temp_term_gt2)
		rm(Igt_term3)

    # we just need the Igg portion of the inverse
    sig_mat <- (-Igg) - (-Igt) %*% solve(-Itt) %*% t(-Igt)
    skatQ <- t(Ugamma) %*% Ugamma
    lambdaQ <- eigen(sig_mat)$values
    p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq

    burdenQ <- (sum(Ugamma))^2
    B_burden= burdenQ / sum(sig_mat);
    p_burden=1-pchisq(B_burden,df=1)

    return(list(p_SKAT=p_SKAT, p_burden=p_burden, skatQ=skatQ, burdenQ=burdenQ, sig_mat = sig_mat, complex=is.complex(lambdaQ)))
}
