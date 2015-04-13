/** 
 * \file uSamplerPoisson_n.cpp
 * \author Felipe Acosta
 * \date 2014-12-30
 * \brief This function performs an MCMC run on the random effects. The arguments are the same arguments used in
 * loglikelihood with an extra argument 'B' which indicates the MCMC sample size and the argument 'u' now 
 * indicates the intial value for the chain.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

double logAcceptPoisson_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& ucurrent, 
const arma::vec& uproposed, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  return min0(0.0, loglikelihoodPoissonCpp_n(beta, sigma, uproposed, kY, kX, kZ) - loglikelihoodPoissonCpp_n(beta, sigma, ucurrent, kY, kX, kZ));
}

// [[Rcpp::export]]
arma::mat uSamplerPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0) {
  RNGScope scope;
  int kK = u.n_rows;
  
  arma::mat usample(B, kK);
  arma::vec ucurrent(kK);
  arma::vec uproposed(kK);
  ucurrent = u;
  usample.row(0) = ucurrent.t();
  
  for (int i = 1; i < B; i++){
    uproposed = rnorm(kK, 0, sd0);
    uproposed += ucurrent;
    if (log(R::runif(0, 1)) < logAcceptPoisson_n(beta, sigma, ucurrent, uproposed, kY, kX, kZ)) {
      ucurrent = uproposed;
    }
    usample.row(i) = ucurrent.t();
  }
  
  return usample;
}
