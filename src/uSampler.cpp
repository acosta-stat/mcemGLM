/** 
 * \file uSampler.cpp
 * \author Felipe Acosta
 * \date 2014-12-30
 * \brief This function performs an MCMC run on the random effects. The arguments are the same arguments used in
 * loglikelihood with an extra argument 'B' which indicates the MCMC sample size.
 */


#include "RcppArmadillo.h"
#include "mcemGLM.h"

using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]


// [[Rcpp::export]]
arma::mat uSamplerCpp(arma::vec beta, arma::mat sigma, arma::vec sigmaType, arma::vec u, 
arma::vec df, arma::vec kKi, arma::vec kLh, arma::vec kY, arma::mat kX, arma::mat kZ, int B) {
  RNGScope scope;
  arma::mat u_sam(B, u.n_rows);
  
  return u_sam;
}