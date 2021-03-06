# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' matrixMultCpp.cpp
#'
#' Faster matrix multiplication for ICSKAT tests.
#'
#' @param A First matrix to be multiplied, will be transposed for t(A) %*% B
#' @param B Second matrix to be multiplied
#' @return C Product of A and B
#' @export
NULL

eigenMapMatMultCrossTwo <- function(A, B) {
    .Call('_ICSKAT_eigenMapMatMultCrossTwo', PACKAGE = 'ICSKAT', A, B)
}

