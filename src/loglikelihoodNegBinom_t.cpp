/**
 * \file loglikelihoodNegBinom_t.cpp
 * \author Felipe Acosta
 * \date 2014-12-02
 * \brief This function evaluates the LogLikelihood function for the negative binomial
 * regression case with t random effects up to a constant (the factorial term is ommited.)
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
double loglikelihoodNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ) {
  double value = 0; /** The value to be returned */
  
  int nObs = kY.n_elem;
  int kP = kX.n_cols;  /** Dimension of Beta */
  int kK = kZ.n_cols;  /** Dimension of U */
  int kR = kKi.n_elem; /** Number of variance components */
  // int kL = sum(kLh);   /** Number of subvariance components */
  
  
  
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
    value += lgamma(kY(i) + alpha) - lgamma(alpha) + alpha * log(alpha) + kY(i) * wij - (kY(i) + alpha) * log(alpha + exp(wij));
  }
  
  int from = 0;
  int to = - 1;
  int counter = 0;
  for (int i = 0; i < kR; i++) {
    for (int j = 0; j < kLh(i); j++) {
      // std::cout<<i<<"\n";
      to += kLhi(counter);
      // std::cout<<"from:"<<from<<'\n';
      // std::cout<<"to:"<<to<<'\n';
      // std::cout<<sigmaType(i)<<"\n";
      // std::cout<<kron(arma::mat(kLhi(counter), kLhi(counter), arma::fill::eye), getSigma(sigma.row(i).t()))<<"\n";
      value += ldmt(u.subvec(from, to), df(counter), sigma.submat(from, from, to, to), sigmaType(i));
      from = to + 1;
      counter += 1;
    }
  }
  
  
  return value;
}
