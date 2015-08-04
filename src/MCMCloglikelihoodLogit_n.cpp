/**
 * \file MCMCloglikelihoodLogit_n.cpp
 * \author Felipe Acosta
 * \date 2015-08-03
 * \brief This function evaluates the LogLikelihood function for the logistic regression case with normal 
 * random effects.
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
 * kY:        Observations, 0 for failure and 1 for success.
 * kX:        Design matrix for fixed effects.
 * kZ:        Design matrix for random effects.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]


// [[Rcpp::export]]
arma::vec MCMCloglikelihoodLogitCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::mat& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  int kM = u.n_rows;
  arma::vec loglike(kM);
  loglike.fill(0);
  
  for (int i = 0; i < kM; i++) {
    loglike(i) = loglikelihoodLogitCpp_n(beta, sigma, u.row(i).t(), kY, kX, kZ);
  }
  
  return loglike;
}
