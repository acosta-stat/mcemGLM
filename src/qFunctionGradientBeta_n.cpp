/**
 * \file qFunctionGradientBeta_n.cpp
 * \author Felipe Acosta
 * \date 2014-12-30
 * \brief This function evaluates the Q function gradient of the algorithm. It performs a loop on the gradient of the loglikelihood function 
 * evaluated on different values of the random effects.
 * Arguments:
 * beta:      The fixed effects coefficients.
 * sigma:     Matrix with r rows. The covariance matrices for the random effects. There are 'r' 
 *            different covariance 
 *            matrices, one matrix per row. The first number of each row is the dimension of each 
 *            covariance matrix. The matrix is reconstructed with the function 'getSigma'.
 * sigmaType: Covariance matrix types, in case of a diagonal matrix the determinant and inverse have 
 *            closed forms. The types are:
 *            0 - diagonal
 * u:         Matrix of MCMC iterations for the random effects. Each row corresponds to one vector of observations.
 * df:        Degrees of freedom for the different groups.
 * kKi        Number of random effects in each variance component. Its length is equal to the number 
 *            of variance components. Its sum is equal to the length of 'u'.
 * kLh:       Number of sub-variance components in each variance component. These have a common 
 *            covariance matrix but different degrees of freedom.
 * kLhi:      Number of random effects in each subvariance componenet.
 * kY:        Observations, 0 for failure and 1 for success.
 * kX:        Design matrix for fixed effects.
 * kZ:        Design matrix for random effects.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::vec qFunctionGradientBetaCpp_n(const arma::vec& beta, const arma::mat& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  int kM = u.n_rows;
  int kP = kX.n_cols;  /** Dimension of Beta */
  arma::vec gradient(kP); /** The value to be returned */
  gradient.fill(0);
  for (int i = 0; i < kM; i++) {
    gradient += loglikelihoodLogitGradientBetaCpp_n(beta, u.row(i).t(), kY, kX, kZ);
  }
  return gradient / kM;
}
