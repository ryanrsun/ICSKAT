// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace std;

using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
  arma::mat C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

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


// Below here are tests

// [[Rcpp::export]]
VectorXd getEigenValues(Map<MatrixXd> M) {
  SelfAdjointEigenSolver<MatrixXd> es(M);
  return es.eigenvalues();
}

// [[Rcpp::export]]
SEXP myfunc(const Eigen::Map<Eigen::MatrixXi> A) {
  const int n(A.cols());
  const MatrixXi At(A.transpose());
  return Rcpp::wrap(At);
}

// [[Rcpp::export]]
SEXP myfunc2(const Eigen::Map<Eigen::MatrixXd> A) {
  const MatrixXd At(A.adjoint());
  return Rcpp::wrap(At);
}


// [[Rcpp::export]]
SEXP myfunc3(const Eigen::Map<Eigen::MatrixXd> A) {
  return Rcpp::wrap(A.adjoint());
}


// [[Rcpp::export]]
SEXP myfunc4(const Eigen::Map<Eigen::MatrixXd> A) {
  return Rcpp::wrap(A.adjoint() * A);
}


// [[Rcpp::export]]
SEXP myfunc5(const Eigen::Map<Eigen::MatrixXi> A) {
  const int n(A.cols());
  MatrixXi AtA(MatrixXi(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
  return Rcpp::wrap(AtA);
}

// [[Rcpp::export]]
SEXP myfunc6(const Eigen::Map<Eigen::MatrixXd> A) {
  const int n(A.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
  return Rcpp::wrap(AtA);
}


