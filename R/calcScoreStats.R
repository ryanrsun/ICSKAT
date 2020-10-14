#' calcScoreStats.R
#'
#' Function that is applied in ICsingleSNP() to calculate a score statistic and p-value for each
#' column of an n*p genotype matrix.
#'
#' @param x n*1 vector of genotypes.
#' @param UgTerm n*1 vector multiplier for the score statistic.
#' @param ggTerm n*1 vector multiplier for the Igg term of the variance.
#' @param gtTermCommon n*p matrix multiplier for the common part of the Igt term of the variance.
#' @param gtHalfL n*(nknots+1) matrix multiplier for one half of the unique part of the Igt term of the variance.
#' @param gtHalfR n*(nknots+1) matrix multiplier for one half of the unique part of the Igt term of the variance.
#' @param solveItt p*p inverse of the Itt matrix from ICSKAT_fit_null().
#' @return A 2*1 vector with the test statistic and then p-value.
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
#' tpos_ind = as.numeric(lt > 0), null_beta=null_fit$beta_fit, Itt=null_fit$Itt, 
#' gMat=matrix(data=rbinom(n=200*10, size=2, prob=0.3), nrow=200))
#'
calcScoreStats <- function(x, UgTerm, ggTerm, gtTermCommon, gtHalfL, gtHalfR, solveItt) {

    # numerator of test statistic
    numerator <- sum(UgTerm * x)
    
    # variance terms
    Igg <- sum(x^2 * ggTerm)
    IgtCommon <- t(x) %*% gtTermCommon
    IgtLeft <- t(x) %*% gtHalfL
    IgtRight <- t(x) %*% gtHalfR
    Igt <- c(IgtCommon, IgtLeft + IgtRight)
    # variance
    varTerm <- -Igg + t(Igt) %*% solveItt %*% Igt
		
    # test statistic
    testStat <- numerator^2 / varTerm
    pval <- 1 - pchisq(testStat, df=1, lower.tail = TRUE)
    return(c(testStat, pval))
}