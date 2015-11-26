/**
 * \file qFunctionDiagGamma_t.cpp
 * \author Felipe Acosta
 * \date 2015-11-22
 * \brief This function evaluates the Q function of the algorithm for the diagonal case.
 * It performs a loop on the loglikelihood function evaluated on different values of the random effects.
 * The function returns the value, gradient, and Hessian of the Q function.
 * Arguments:
 * beta:      The fixed effects coefficients.
 * sigma:     Matrix with r rows. The covariance matrices for the random effects. There are 'r' 
 *            different covariance 
 *            matrices, one matrix per row. The first number of each row is the dimension of each 
 *            covariance matrix. The matrix is reconstructed with the function 'getSigma'.
 * u:         Matrix of MCMC iterations for the random effects. Each row corresponds to one vector of observations.
 * df:        Degrees of freedom for the different groups.
 * kKi        Number of random effects in each variance component. Its length is equal to the number 
 *            of variance components. Its sum is equal to the length of 'u'.
 * kY:        Observations, 0 for failure and 1 for success.
 * kX:        Design matrix for fixed effects.
 * kZ:        Design matrix for random effects.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
List qFunctionDiagGammaCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::mat& u, 
                                const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
                                const arma::mat& kX, const arma::mat& kZ) {
  int kM = u.n_rows;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kR = kKi.n_elem; /** Number of random effects */

  double value = 0;
  arma::vec gradient(kP + 1 + kR);
  gradient.fill(0);
  arma::mat hessian(kP + 1 + kR, kP + 1 + kR);
  hessian.fill(0);

  for (int i = 0; i < kM; i++) {
    value += loglikelihoodGammaCpp_t(beta, sigma, alpha, sigmaType, u.row(i).t(), df, kKi, kLh, kLhi, kY, kX, kZ) / kM;
    gradient += loglikelihoodGammaGradientCpp_t(beta, sigma, alpha, u.row(i).t(), df, kKi, kLh, kLhi, kY, kX, kZ) / kM;
    hessian += loglikelihoodGammaHessianCpp_t(beta, sigma, alpha, u.row(i).t(), df, kKi, kLh, kLhi, kY, kX, kZ) / kM;
  }
  return List::create(Named("value") = value, Named("gradient") = gradient, Named("hessian") = hessian);
}
