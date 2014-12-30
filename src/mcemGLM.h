#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

double ldmt(arma::vec x, double df, arma::mat sigma, int sigmaType);

arma::mat getSigma(arma::vec x);

double loglikehoodLogitCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double qFunctionCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat uSamplerCpp(arma::vec beta, arma::mat sigma, arma::vec sigmaType, arma::vec u, 
arma::vec df, arma::vec kKi, arma::vec kLh, arma::vec kY, arma::mat kX, arma::mat kZ, int B);
