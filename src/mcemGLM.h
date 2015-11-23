#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]

double min0(double a, double b);

double ldmt(arma::vec x, double df, arma::mat sigma, int sigmaType);

double ldmn(const arma::vec& x, const arma::mat& sigma);

arma::mat getSigma(arma::vec x);

double loglikelihoodLogitCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double qFunctionCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat uSamplerCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

double loglikelihoodLogitCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, 
const arma::vec& u, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double qFunctionCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat uSamplerCpp(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::vec loglikelihoodLogitGradientCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodLogitHessianCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

List qFunctionDiagCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& kKi, const arma::mat& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::vec loglikelihoodLogitGradientCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodLogitHessianCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

List qFunctionDiagCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat iMatrixDiagCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::mat iMatrixDiagCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kKi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, 
int B, double sd0);

double loglikelihoodPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

double loglikelihoodPoissonCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ);

arma::vec loglikelihoodPoissonGradientCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::vec loglikelihoodPoissonGradientCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodPoissonHessianCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodPoissonHessianCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat uSamplerPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::mat uSamplerPoissonCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::mat iMatrixDiagPoissonCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::mat iMatrixDiagPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, 
const arma::vec& kKi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, 
int B, double sd0);

double loglikelihoodNegBinomCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double loglikelihoodNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ);

arma::vec loglikelihoodNegBinomGradientCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::vec loglikelihoodNegBinomGradientCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodNegBinomHessianCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodNegBinomHessianCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

List qFunctionDiagNegBinomCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, const arma::mat& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

List qFunctionDiagNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat uSamplerNegBinomCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::mat uSamplerNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::mat iMatrixDiagNegBinomCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& u, 
const arma::vec& kKi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, 
int B, double sd0);

arma::mat iMatrixDiagNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::vec MCMCloglikelihoodPoissonCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ);

arma::vec MCMCloglikelihoodPoissonCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::mat& u, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::vec MCMCloglikelihoodNegBinomCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, 
const arma::mat& u, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::vec MCMCloglikelihoodNegBinomCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, 
const arma::mat& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::vec MCMCloglikelihoodLogitCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ);

arma::vec MCMCloglikelihoodLogitCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::mat& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double loglikelihoodGammaCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::vec loglikelihoodGammaGradientCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat loglikelihoodGammaHessianCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, 
const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

arma::mat iMatrixDiagGammaCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::mat& uSample, 
const arma::vec& kKi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::vec MCMCloglikelihoodGammaCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, 
const arma::mat& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

List qFunctionDiagGammaCpp_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& kKi, const arma::mat& u, 
const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double logAcceptGamma_n(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& ucurrent, 
const arma::vec& uproposed, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);

double loglikelihoodGammaCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ);

arma::vec loglikelihoodGammaGradientCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, 
const arma::mat& kZ);

arma::mat loglikelihoodGammaHessianCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat iMatrixDiagGammaCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, 
const arma::mat& uSample, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);

arma::vec MCMCloglikelihoodGammaCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, 
const arma::mat& u, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, \
const arma::mat& kX, const arma::mat& kZ);

List qFunctionDiagGammaCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::mat& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ);

arma::mat uSamplerGammaCpp_t(const arma::vec& beta, const arma::mat& sigma, double alpha, const arma::vec& sigmaType, const arma::vec& u, 
const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, 
const arma::mat& kX, const arma::mat& kZ, int B, double sd0);
