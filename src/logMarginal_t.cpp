/**
 * \file logMarginal_t.cpp
 * \author Felipe Acosta
 * \date 2014-12-02
 * \brief This function evaluates the LogLikelihood function for the logistic regression case with t 
 * random effects.
 * Arguments:
 * sigma:     Matrix with r rows. The covariance matrices for the random effects. There are 'r' 
 *            different covariance 
 *            matrices, one matrix per row. The first number of each row is the dimension of each 
 *            covariance matrix. The matrix is reconstructed with the function 'getSigma'.
 * sigmaType: Covariance matrix types, in case of a diagonal matrix the determinant and inverse have 
 *            closed forms. The types are:
 *            0 - diagonal
 * u:         Vector of random effects.
 * df:        Degrees of freedom for the different groups.
 * kKi        Number of random effects in each variance component. Its length is equal to the number 
 *            of variance components. Its sum is equal to the length of 'u'.
 * kLh:       Number of sub-variance components in each variance component. These have a common 
 *            covariance structure but different degrees of freedom.
 * kLhi:      Number of random effects in each subvariance component.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]


// [[Rcpp::export]]
double logMarginalCpp_t(const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi) {
  double value = 0; /** The value to be returned */
  int kR = kKi.n_elem; /** Number of variance components */
  int from = 0;
  int to = - 1;
  int counter = 0;
  for (int i = 0; i < kR; i++) {
    for (int j = 0; j < kLh(i); j++) {
      // std::cout<<i<<"\n";
      to += kLhi(counter);
      value += ldmt(u.subvec(from, to), df(counter), sigma.submat(from, from, to, to), sigmaType(i));
      from = to + 1;
      counter += 1;
    }
  }
  return value;
}
