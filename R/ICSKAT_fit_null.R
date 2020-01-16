#' ICSKAT_fit_null.R
#'
#' Fit the null model (cubic basis spline for baseline cumulative hazard and coefficients
#' for non-genetic coefficiens) for interval-censored skat.
#'
#'
#' @param init_beta (p+nknots+2)*1 vector of coefficients to initialize the Newton-Raphson.
#' @param left_dmat n*(p+nknots+2) design matrix for left end of interval.
#' @param obs_ind n*1 vector of whether the event was observed before last follow-up.
#' @param tpos_ind n*1 vector of whether the event was observed after follow-up started (t>0).
#' @param right_dmat n*(p+nknots+2) design matrix for right end of interval.
#' @param eps Stop when the L2 norm of the difference in model coefficients reaches this limit.
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
#' ICSKAT_fit_null(init_beta=rep(0, 5), left_dmat=dmats$left_dmat, lt=lt, rt=rt,
#' right_dmat=dmats$right_dmat, obs_ind=rep(1, n), tpos_ind = as.numeric(lt > 0))
#'
ICSKAT_fit_null <- function(init_beta, left_dmat, obs_ind, tpos_ind, right_dmat, lt, rt, eps=10^(-6)) {

    diff_beta <- 1
    iter <- 0
    temp_beta <- init_beta
    while( !is.nan(diff_beta) & diff_beta > eps & diff_beta<1000) {
        # Cumulative hazard under null
        H_L <- exp(left_dmat %*% temp_beta)
        H_R <- exp(right_dmat %*% temp_beta)
        # Survival term
        SL <- ifelse(lt == 0, 1, exp(-H_L))
        SR <- ifelse(rt == 999, 0, exp(-H_R))
        SR[!is.finite(SR)] <- 0
        A <- SL - SR
        # sometimes A is 0
        A[which(A == 0)] <- min(A[which(A > 0)])

        # score vector
        U_term1 <- t(left_dmat) %*% diag(ifelse(lt == 0, 0, exp(-H_L) * -H_L))
        U_term2 <- t(right_dmat) %*% diag(ifelse(rt == 999, 0, exp(-H_R) * -H_R))
        U_term2[is.na(U_term2)] <- 0
        U <- apply( (U_term1 - U_term2) %*% diag(1/A), 1, sum)

        # information matrix
        I_term1 <- t(left_dmat) %*% diag( tpos_ind * as.numeric((-H_L * exp(-H_L) + H_L^2 * exp(-H_L)) / A) ) %*% left_dmat
        # sometimes H_R is so large that when it gets squared it goes to Inf and then multiplied
        # by exp(-H_R) it goes to NaN
        check_Iterm2 <- (H_R * exp(-H_R) - H_R^2 * exp(-H_R))
        check_Iterm2[which(is.nan(check_Iterm2))] <- 0
        I_term2 <- t(right_dmat) %*% diag( obs_ind * as.numeric(check_Iterm2 / A) ) %*% right_dmat
        check_tempterm <- (H_R * exp(-H_R))
        check_tempterm[which(is.nan(check_tempterm))] <- 0
        temp_term <- t(left_dmat) %*% diag( tpos_ind * as.numeric((H_L * exp(-H_L)) / A) ) -
            t(right_dmat) %*% diag( obs_ind * as.numeric( check_tempterm / A) )
        I_term3 <- temp_term %*% t(temp_term)
        I <- I_term1 + I_term2  - I_term3

        beta_new <- temp_beta - t(U) %*% solve(I)
        diff_beta <- (beta_new - temp_beta) %*% t(beta_new - temp_beta)
        temp_beta <- as.numeric(beta_new)
        iter <- iter + 1
    }

    return(list(beta_fit=beta_new, iter=iter, Itt=I))
}
