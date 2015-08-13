/**
 * \file MCMCloglikelihoodNegBinom_t.cpp
 * \author Felipe Acosta
 * \date 2015-8-02
 * \brief This function evaluates the LogLikelihood function for the negative binomial
 * regression case with t random effects up to a constant (the factorial term is ommited.)
 * Arguments:
 * beta:      The fixed effects coefficients.
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
 * kY:        Observations, 0 for failure and 1 for success.
 * kX:        Design matrix for fixed effects.
 * kZ:        Design matrix for random effects.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]


// [[Rcpp::export]]
arma::vec MCMCloglikelihoodNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  int kM = u.n_rows;  /** Number of MCMC samples */
  arma::vec loglike(kM);
  loglike.fill(0);
  
  for (int i = 0; i < kM; i++) {
    loglike(i) = loglikelihoodNegBinomCpp_t(beta, sigma, alpha, sigmaType, u.row(i).t(), df, kKi, kLh, kLhi, kY, kX, kZ);
  }
  
  return loglike;
}
