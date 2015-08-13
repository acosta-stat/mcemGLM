/**
 * \file mcemGLM-util.cpp
 * \author Felipe Acosta
 * \date 2014-12-08
 * \brief This file contains various functions that are used in the package.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double min0(double a, double b) {
  if (a < b)
    return a;
  return b;
}

// [[Rcpp::export]]
/** Evaluate the log density of a multivariate normal distribution with mean vector 0 */
double ldmn(const arma::vec& x, const arma::mat& sigma) {
  int kDIM = sigma.n_cols;
  for (int i = 0; i < kDIM; i++) {
    for (int j = 0; j < kDIM; j++) {
      if(sigma(i, j) < 0) {
        return -INFINITY;
      }
    }
  }
  
  double VALUE = -0.5 * kDIM * log(2 * M_PI) - 0.5 * log(arma::det(sigma));
  
  arma::mat sigmainv;
  sigmainv = inv(sigma);
  //std::cout<<sigmainv(1,1);
  
  double tmp0 = 0;
  NumericVector tmpVector(kDIM); /** stores the product of x^t and sigma^-1 */
  for (int i = 0; i < kDIM; i++) {
    for (int j = 0; j < kDIM; j++) {
      tmpVector(i) += x(j) * sigmainv(j, i);
    }
  }
  for (int i = 0; i < kDIM; i++) {
    tmp0 += tmpVector(i) * x(i);
  }
  VALUE += - 0.5 * tmp0;
  return VALUE;
  
}

// [[Rcpp::export]]
/** Evaluate the log density of a multivariate t distribution with mean vector 0*/
double ldmt(arma::vec x, double df, arma::mat sigma, int sigmaType) {
  if (df <= 0) {
    return -INFINITY;
  }
  int n = sigma.n_cols;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (sigma(i, j) < 0) {
        return -INFINITY;
      }
    }
  }
  double value = 0; /** Density value */
  int k = x.size();
  //arma::mat sigma0 = as<arma::mat>(sigma);
  value += lgamma(0.5 * (df + k)) - lgamma(0.5 * df) - 0.5 * k * log(df) - 0.5 * k * log(M_PI) - 0.5 * log(arma::det(sigma));
  
  /** If sigma is diagonal it's easier to do the multiplication */
  double tmp0 = 0;
  if (sigmaType == 0) {
    for (int i = 0; i < k; i++) {
      tmp0 += x(i) * x(i) / sigma(i,i);
    }
    value += - 0.5 * (df + k) * log(1 + tmp0/df);
    return value;
  }
  
  arma::mat sigmainv;
  sigmainv = inv_sympd(sigma);
  //std::cout<<sigmainv(1,1);
  NumericVector tmpVector(k); /** stores the product of x^t and sigma^-1 */
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      tmpVector(i) += x(j) * sigmainv(j, i);
    }
  }
  for (int i = 0; i < k; i++) {
    tmp0 += tmpVector(i) * x(i);
  }
  value += - 0.5 * (df + k) * log(1 + tmp0/df);
  return value;
}

