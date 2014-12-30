/**
 * \file loglikelihoodLogit_t.cpp
 * \author Felipe Acosta
 * \date 2014-12-02
 * \brief This function evaluates the LogLikelihood function for the logistic regression case with t random effects.
 * Arguments:
 * beta:      The fixed coefficients.
 * sigma:     The covariance matrix for the random effects.
 * sigmaType: Covariance matrix types, in case of a diagonal matrix the determinant and inverse have closed forms. The types are:
 *            0 - diagonal
 * u:         Vector of random effects.
 * df:        Degrees of freedom for the different groups.
 * kY:        Observations, 0 for failure and 1 for success.
 * kYgroup:   Grouping vector. Groups numbering start at 0. Maximum group number is the length of 'df'.
 * kX:        Design matrix for fixed effects.
 * kZ:        Design matrix for random effects.
 */


#include <Rcpp.h>
#include "mcemGLM.h"

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]

double loglikehoodLogitCpp_t(NumericVector beta, NumericMatrix sigma, IntegerVector sigmaType, NumericVector u, NumericVector df,
IntegerVector kKi, IntegerMatrix kKih, NumericVector kY, IntegerMatrix kYgroup, NumericMatrix kX, NumericMatrix kZ) {
  double value = 0; /** The value to be returned */
  
  /** Diagonal matrix case */
  if (sigma == 0) {
    int nObs = kY.size();
    int kP = kX.ncol();
    int kR = kZ.ncol();
    int kK = sigma.ncol();
    
    /** sum of yij * (wij - log(1 + ...)) */
    for (int i = 0; i < nObs; i++) {
      double wij = 0;
      for (int j = 0; j < kP; j++) {
        wij += kX(i, j) * beta(j);
      }
      for (int j = 0; j < kR; j++) {
        wij += kZ(i, j) * u(j);
      }
      value += kY(i) * wij - log(1 + exp(wij));
    }
  }
  return value;
}