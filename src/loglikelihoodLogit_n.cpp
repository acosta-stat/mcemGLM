/**
 * \file loglikelihoodLogit_t.cpp
 * \author Felipe Acosta
 * \date 2014-12-02
 * \brief This function evaluates the LogLikelihood function for the logistic regression case with t 
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
double loglikelihoodLogitCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  double value = 0; /** The value to be returned */
  
  int nObs = kY.n_elem;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kK = kZ.n_cols;  /** Dimension of U */
  int kR = kKi.n_elem; /** Number of variance components */
  int kL = sum(kLh);   /** Number of subvariance components */
  
  
  
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
    value += kY(i) * wij - log(1 + exp(wij));
  }
  
  value += ldmn(u, sigma);
  
  return value;
}