/**
 * \file iJacobDiagNegBinom_n.cpp
 * \author Felipe Acosta
 * \date 2015-04-07
 * \brief This function calculates the inverse of I - Jacobian for the logit with normal random effects.
 * Arguments:
 * beta:      The fixed effects coefficients.
 * sigma:     Matrix with r rows. The covariance matrices for the random effects. There are 'r' 
 *            different covariance 
 *            matrices, one matrix per row. The first number of each row is the dimension of each 
 *            covariance matrix. The matrix is reconstructed with the function 'getSigma'.
 * sigmaType: Covariance matrix types, in case of a diagonal matrix the determinant and inverse have 
 *            closed forms. The types are:
 *            0 - diagonal
 * u:         Initial value for the vector of random effects.
 * kY:        Observations, 0 for failure and 1 for success.
 * kX:        Design matrix for fixed effects.
 * kZ:        Design matrix for random effects.
 */

#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::mat iJacobDiagNegBinomCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, 
const arma::mat& uSample, const arma::vec& kKi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ, int B, double sd0) {
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kR = kKi.n_elem; /** Number of random effects */
  
  //arma::mat uSample(B, kR); /** MCMC sample from U */
  //uSample = uSamplerCpp_n(beta, sigma, u, kY, kX, kZ, B, sd0);
  
  arma::vec g0(kP + 1 + kR); /** Gradient vector */
  g0.fill(0);
  arma::vec g1(kP + kR); /** Gradient vector */
  g1.fill(0);
  arma::mat h0(kP + 1 + kR, kP + 1 + kR); /** Hessian Matrix */
  h0.fill(0);
  arma::mat iMatrix(kP + 1 + kR, kP + 1 + kR); /** Information Matrix */
  iMatrix.fill(0);
  arma::mat Ix(kP + 1 + kR, kP + 1 + kR);
  Ix.fill(0);
  
  for (int i = 0; i < B; i++) {
    g0 = loglikelihoodNegBinomGradientCpp_n(beta, sigma, alpha, kKi, uSample.row(i).t(), kY, kX, kZ);
    h0 = loglikelihoodNegBinomHessianCpp_n(beta, sigma, alpha, kKi, uSample.row(i).t(), kY, kX, kZ);
    iMatrix += (-h0 - g0 * g0.t()) / (double) B;
    g1 += g0 / (double) B;
    Ix += -h0 / (double) B;
  }
  iMatrix =+ g1 * g1.t();
  return(Ix * inv(iMatrix));
}