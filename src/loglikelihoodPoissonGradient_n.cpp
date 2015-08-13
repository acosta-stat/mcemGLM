/**
 * \file loglikelihoodPoissonGradient_n.cpp
 * \author Felipe Acosta
 * \date 2014-12-02
 * \brief This function evaluates the LogLikelihood gradient function for the Poisson regression 
 * case with normal random effects in the diagonal case.
 * Arguments:
 * beta:      Fixed effects parameters.
 * sigma:     Covariance matrix.
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
arma::vec loglikelihoodPoissonGradientCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  int nObs = kY.n_elem;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kK = kZ.n_cols;  /** Dimension of U */
  int kR = kKi.n_elem; /** Number of variance components */
  
  arma::vec gradient(kP + kR); /** The value to be returned */
  gradient.fill(0);
  
  for (int i = 0; i < nObs; i++) {
    double wij = 0;
    for (int j = 0; j < kP; j++) {
      wij += kX(i, j) * beta(j);
    }
    
    for (int j = 0; j < kK; j++) {
      wij += kZ(i, j) * u(j);
    }
    for (int j = 0; j < kP; j++) {
      gradient(j) += - kX(i, j) * exp(wij) + kY(i) * kX(i, j);
    }
  }
  
  int counter = 0;
  for (int i = 0; i < kR; i++) {
    double sumU = 0;
    double lambda_i = sigma(counter, counter);
    for (int j = 0; j < kKi(i); j++) {
      sumU += u(counter) * u(counter);
      counter += 1;
    }
    gradient(kP + i) = -0.5 * kKi(i) / lambda_i + 0.5 / (lambda_i * lambda_i) * sumU;
  }
  
  return gradient;
}
