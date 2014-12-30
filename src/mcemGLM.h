#include <Rcpp.h>
using namespace Rcpp;

double ldmt(NumericVector x, double df, NumericMatrix sigma, int sigmaType);