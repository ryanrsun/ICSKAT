// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
using namespace Eigen;
using namespace std;

// [[Rcpp::export]]
SEXP eigenMapMatMultC(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMultCrossprod(const Eigen::Map<Eigen::MatrixXd> A){
  const int n(A.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
  return Rcpp::wrap(AtA);
}

// [[Rcpp::export]]
SEXP eigenMapMatMultCrossTwo(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A.adjoint() * B;
  return Rcpp::wrap(C);
}

