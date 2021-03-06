% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICSKAT.R
\name{ICskat}
\alias{ICskat}
\title{ICSKAT.R}
\usage{
ICskat(
  left_dmat,
  right_dmat,
  lt,
  rt,
  obs_ind,
  tpos_ind,
  gMat,
  null_beta,
  Itt,
  pvalue = TRUE
)
}
\arguments{
\item{left_dmat}{n*(p+nknots+2) design matrix for left end of interval.}

\item{right_dmat}{n*(p+nknots+2) design matrix for right end of interval.}

\item{lt}{n*1 vector of left side of interval times.}

\item{rt}{n*1 vector of right side of interval times.}

\item{obs_ind}{n*1 vector of whether the event was observed before last follow-up.}

\item{tpos_ind}{n*1 vector of whether the event was observed after follow-up started (t>0).}

\item{gMat}{n*q genotype matrix.}

\item{null_beta}{(p+nknots+2)*1 vector of coefficients for null model.}

\item{Itt}{(p+nknots+2)*(p+nknots+2) Fisher information matrix for null model coefficients.}

\item{pvalue}{Boolean, if TRUE then find the p-value (maybe don't need it if bootstrapping, saves eigendecomposition)}
}
\value{
A list with the elements:
\item{p_SKAT}{ICSKAT p-value}
\item{p_burden}{IC burden test p-value}
\item{complex}{Indicator of whether the SKAT variance matrix was positive definite}
\item{sig_mat}{The covariance matrix of the score equations for genetic effects when treated as fixed effects}
\item{skatQ}{SKAT test statistic}
\item{burdenQ}{Burden test statistic}
\item{Ugamma}{Score vector}
\item{lambdaQ}{Vector of eigenvalues of variance matrix}
\item{null_beta}{The fitted null parameters}
\item{err}{Will be 0 for no error, 22 if had to adjust parameters on CompQuadForm (totally normal), or 99 if NA in variance matrix. ICSKATwrapper will return 1 here if the null fit has an error}
\item{errMsg}{Explains error code, blank string if no error}
}
\description{
Calculate the test statistic and p-value for interval-censored SKAT.
}
\examples{
set.seed(2)
gMat <- matrix(data=rbinom(n=2000, size=2, prob=0.3), nrow=100)
xMat <- matrix(data=rnorm(200), nrow=100)
bhFunInv <- function(x) {x}
obsTimes <- 1:5
etaVec <- rep(0, 100)
outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, windowHalf = 0.1,
probMiss = 0.1, etaVec = etaVec)
lt <- outcomeDat$leftTimes
rt <- outcomeDat$rightTimes
tpos_ind <- as.numeric(lt > 0)
obs_ind <- as.numeric(rt != Inf)
dmats <- make_IC_dmat(xMat, lt, rt, obs_ind, tpos_ind)
nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat,
right_dmat=dmats$right_dmat, obs_ind = obs_ind, tpos_ind = tpos_ind,
lt = lt, rt = rt)
ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
lt = lt, rt = rt, obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat,
null_beta = nullFit$beta_fit, Itt = nullFit$Itt)

}
