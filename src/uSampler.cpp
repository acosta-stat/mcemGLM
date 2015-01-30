/** 
 * \file uSampler.cpp
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

double min0(double a, double b) {
  if (a < b)
    return a;
  return b;
}

double logAccept(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& ucurrent, 
const arma::vec& uproposed, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ) {
  return min0(0.0, loglikehoodLogitCpp_t(beta, sigma, sigmaType, uproposed, df, kKi, kLh, kLhi, kY, kX, kZ)
  - loglikehoodLogitCpp_t(beta, sigma, sigmaType, ucurrent, df, kKi, kLh, kLhi, kY, kX, kZ));
}

// [[Rcpp::export]]
arma::mat uSamplerCpp(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, 
int B, double sd0) {
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
    if (log(R::runif(0, 1)) < logAccept(beta, sigma, sigmaType, ucurrent, uproposed, df, kKi, kLh, kLhi, kY, kX, kZ)) {
      ucurrent = uproposed;
    }
    usample.row(i) = ucurrent.t();
  }
  
  return usample;
  
  /*
  std::cout<<u<<"\n";
  std::cout<<ucurrent<<"\n";
  ucurrent(0) = -10;
  std::cout<<u<<"\n";
  std::cout<<ucurrent<<"\n";
  
  arma::vec xx(5);
  xx = rnorm(5,0,1);
  xx(0) = -10;
  std::cout<<xx.t()<<"\n\n";
  
  arma::vec x(7);
  for (int i = 0; i < 7; i++)
    x(i) = 1;
  std::cout<<usample.row(0);
  usample.row(0) = x.t();
  std::cout<<usample.row(0);
  */
  
  
  
}
