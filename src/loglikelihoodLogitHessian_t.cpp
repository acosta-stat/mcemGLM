/**
 * \file loglikelihoodLogitHessian_t.cpp
 * \author Felipe Acosta
 * \date 2015-04-07
 * \brief This function evaluates the LogLikelihood gradient function for the logistic regression case with normal 
 * random effects for the diagonal case.
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
arma::mat loglikelihoodLogitHessianCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ) {
  int nObs = kY.n_elem;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kK = kZ.n_cols;  /** Dimension of U */
  int kR = kKi.n_elem; /** Number of random effects */
  
  arma::mat hessian(kP + kR, kP + kR); /** The value to be returned */
  hessian.fill(0);
  
  for (int i = 0; i < nObs; i++) {
    double wij = 0;
    for (int j = 0; j < kP; j++) {
      wij += kX(i, j) * beta(j);
    }
    
    for (int j = 0; j < kK; j++) {
      wij += kZ(i, j) * u(j);
    }
    for (int j = 0; j < kP; j++) {
      for (int k = 0; k <= j; k++) {
        hessian(j, k) += -kX(i, j) * kX(i, k) * exp(wij) / pow((1 + exp(wij)), 2);
        if (k < j) {
          hessian(k, j) = hessian(j, k);
        }
      }
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
    hessian(kP + i, kP + i) = 0.5 * kKi(i) / (pow(lambda_i, 2.0)) - 0.5 * (df(i) + kKi(i)) *
                              (2.0 * pow(lambda_i, -3.0) * df(i) * sumU + pow(lambda_i, -4.0) * pow(sumU, 2.0)) / pow(df(i) + sumU / lambda_i, 2);
  }
  return hessian;
}
