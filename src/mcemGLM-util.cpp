/**
 * \file mcemGLM-util.cpp
 * \author Felipe Acosta
 * \date 2014-12-08
 * \brief This function contains various functions that are used in the package.
 */


#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

// #include <Rcpp.h>
// [[Rcpp::export]]
/** Evaluate the log density of a multivariate t distribution with mean vector 0*/
double ldmt(arma::vec x, double df, arma::mat sigma, int sigmaType) {
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
