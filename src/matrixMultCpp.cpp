// [[Rcpp::depends(RcppEigen)]]
//' matrixMultCpp.cpp
//'
//' Faster matrix multiplication for ICSKAT tests.
//'
//' @param A First matrix to be multiplied, will be transposed for t(A) %*% B
//' @param B Second matrix to be multiplied
//' @return C Product of A and B
//' @export


#include <RcppEigen.h>
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
SEXP eigenMapMatMultCrossTwo(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A.adjoint() * B;
  return Rcpp::wrap(C);
}

