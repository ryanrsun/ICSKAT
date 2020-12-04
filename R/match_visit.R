#' match_visit.R
#'
#' Match visit to a time for model-based bootstrap with interval-censored data.
#'
#' @param draw (s+1)-length vector of all 0s except for one 1, which is the failure interval.
#' @param visitTimes s-length vector where is the number of inspection times for a subject.
#'
#' @return n*(s+1) matrix where element (i,j) holds the probability that subject i will fail in interval j.
#'
#' @export
#'
matchVisit <- function(draw, visitTimes) {
  # draw should be all 0s and one 1
  if (sum(abs(draw)) > 1) {stop("Only 1 draw per patient")}
  # number of visits should be one less than length of draw vector
  nVisits <- length(visitTimes)
  if (length(draw) != nVisits + 1) {stop("nDraws = nVisits + 1")}
  idx <- which(draw == 1)

  # draw on the ends or in the middle
  if (idx == 1) {
    return(c(0, visitTimes[1]))
  } else if (idx == nVisits + 1) {
    return(c(visitTimes[nVisits], 999))
  } else {
    return(c(visitTimes[idx-1], c(visitTimes[idx])))
  }
}
