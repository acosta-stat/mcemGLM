/**
 * \file loglikelihoodGamma_n.cpp
 * \author Felipe Acosta
 * \date 2015-11-20
 * \brief This function evaluates the LogLikelihood function for the gamma 
 * regression case with normal random effects up to a constant (the factorial term is ommited.)
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
double loglikelihoodGammaCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  double value = 0; /** The value to be returned */
  
  int nObs = kY.n_elem;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kK = kZ.n_cols;  /** Dimension of U */
  
  /** sum of yij * (wij - log(1 + ...))
   *  This corresponds to the 
  */
  for (int i = 0; i < nObs; i++) {
    double wij = 0;
    for (int j = 0; j < kP; j++) {
      wij += kX(i, j) * beta(j);
    }
    
    for (int j = 0; j < kK; j++) {
      wij += kZ(i, j) * u(j);
    }
    // value += exp(wij) * alpha * log(alpha * kY(i)) - lgamma(alpha * exp(wij)) - alpha * kY(i);
    value += alpha * log(alpha) - alpha * wij - lgamma(alpha) + alpha * log(kY(i)) - alpha * kY(i) * exp(-wij);
  }
  
  value += ldmn(u, sigma);
  return value;
}
