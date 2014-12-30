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

double logAccept(arma::vec beta, arma::mat sigma, arma::vec sigmaType, arma::vec ucurrent,arma::vec uproposed,
arma::vec df, arma::vec kKi, arma::vec kLh, arma::vec kY, arma::mat kX, arma::mat kZ) {
  return min0(0.0, loglikehoodLogitCpp_t(beta, sigma, sigmaType, uproposed, df, kKi, kLh, kY, kX, kZ)
  - loglikehoodLogitCpp_t(beta, sigma, sigmaType, ucurrent, df, kKi, kLh, kY, kX, kZ));
}

// [[Rcpp::export]]
arma::mat uSamplerCpp(arma::vec beta, arma::mat sigma, arma::vec sigmaType, const arma::vec& u, 
arma::vec df, arma::vec kKi, arma::vec kLh, arma::vec kY, const arma::mat& kX, const arma::mat& kZ, int B) {
  RNGScope scope;
  int kK = u.n_rows;
  
  arma::mat usample(B, kK);
  arma::vec ucurrent(kK);
  arma::vec uproposed(kK);
  ucurrent = u;
  
  std::cout<<u<<"\n";
  std::cout<<ucurrent<<"\n";
  ucurrent(0) = -10;
  std::cout<<u<<"\n";
  std::cout<<ucurrent<<"\n";
  
  arma::vec xx(5);
  xx = rnorm(5,0,1);
  xx(0) = -100;
  std::cout<<xx.t()<<"\n\n";
  
  arma::vec x(7);
  for (int i = 0; i < 7; i++)
    x(i) = 1;
  std::cout<<usample.row(0);
  usample.row(0) = x.t();
  std::cout<<usample.row(0);
  return usample;
}
