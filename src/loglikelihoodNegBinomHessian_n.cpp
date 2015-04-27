/**
 * \file loglikelihoodNegBinomHessian_n.cpp
 * \author Felipe Acosta
 * \date 2015-04-27
 * \brief This function evaluates the LogLikelihood Hessian matrix for the negative binomial
 * regression case with normal random effects for the diagonal case.
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
arma::mat loglikelihoodNegBinomHessianCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  int nObs = kY.n_elem;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kK = kZ.n_cols;  /** Dimension of U */
  int kR = kKi.n_elem; /** Number of random effects */
  
  arma::mat hessian(kP + 1 + kR, kP + 1 + kR); /** The value to be returned */
  hessian.fill(0);
  
  for (int i = 0; i < nObs; i++) {
    double wij = 0;
    for (int j = 0; j < kP; j++) {
      wij += kX(i, j) * beta(j);
    }
    
    for (int j = 0; j < kK; j++) {
      wij += kZ(i, j) * u(j);
    }
    // Hessian for Beta
    for (int j = 0; j < kP; j++) {
      for (int k = 0; k <= j; k++) {
        hessian(j, k) += -alpha * kX(i, j) * kX(i, k) * (alpha + kY(i)) * exp(wij) / pow(alpha + exp(wij), 2);
        if (k < j) {
          hessian(k, j) = hessian(j, k);
        }
      }
    }
    // Hessian for alpha+Beta
    for (int j = 0; j < kP; j++) {
      hessian(j, kP) += kX(i, j) * exp(wij) * (kY(i) - exp(wij)) / pow(alpha + exp(wij), 2.0);
      hessian(kP, j) = hessian(j, kP);
    }
    // Hessian for alpha
    hessian(kP, kP) += R::trigamma(alpha + kY(i)) - R::trigamma(alpha) + 1.0/alpha - 2.0/(alpha + exp(wij)) + (alpha + kY(i)) / pow(alpha + exp(wij), 2.0);
  }
  
  int counter = 0;
  for (int i = 0; i < kR; i++) {
    double sumU = 0;
    double lambda_i = sigma(counter, counter);
    for (int j = 0; j < kKi(i); j++) {
      sumU += u(counter) * u(counter);
      counter += 1;
    }
    hessian(kP + 1 + i, kP + 1 + i) = 0.5 * kKi(i) / (lambda_i * lambda_i) - sumU / (lambda_i * lambda_i * lambda_i);
  }
  
  return hessian;
}
