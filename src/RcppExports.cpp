// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// loglikelihoodLogitCpp_n
double loglikelihoodLogitCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);
RcppExport SEXP mcemGLM_loglikelihoodLogitCpp_n(SEXP betaSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP kYSEXP, SEXP kXSEXP, SEXP kZSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kY(kYSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kX(kXSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kZ(kZSEXP );
        double __result = loglikelihoodLogitCpp_n(beta, sigma, u, kY, kX, kZ);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// loglikelihoodLogitCpp_t
double loglikelihoodLogitCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);
RcppExport SEXP mcemGLM_loglikelihoodLogitCpp_t(SEXP betaSEXP, SEXP sigmaSEXP, SEXP sigmaTypeSEXP, SEXP uSEXP, SEXP dfSEXP, SEXP kKiSEXP, SEXP kLhSEXP, SEXP kLhiSEXP, SEXP kYSEXP, SEXP kXSEXP, SEXP kZSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type sigmaType(sigmaTypeSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type df(dfSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kKi(kKiSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kLh(kLhSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kLhi(kLhiSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kY(kYSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kX(kXSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kZ(kZSEXP );
        double __result = loglikelihoodLogitCpp_t(beta, sigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ldmn
/** Evaluate the log density of a multivariate normal distribution with mean vector 0 */ double ldmn(const arma::vec& x, const arma::mat& sigma);
RcppExport SEXP mcemGLM_ldmn(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        /** Evaluate the log density of a multivariate normal distribution with mean vector 0 */ double __result = ldmn(x, sigma);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// ldmt
/** Evaluate the log density of a multivariate t distribution with mean vector 0*/ double ldmt(arma::vec x, double df, arma::mat sigma, int sigmaType);
RcppExport SEXP mcemGLM_ldmt(SEXP xSEXP, SEXP dfSEXP, SEXP sigmaSEXP, SEXP sigmaTypeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        Rcpp::traits::input_parameter< double >::type df(dfSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< int >::type sigmaType(sigmaTypeSEXP );
        /** Evaluate the log density of a multivariate t distribution with mean vector 0*/ double __result = ldmt(x, df, sigma, sigmaType);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// getSigma
/** Reconstruct a nxn covariance matrix from a vector */ arma::mat getSigma(arma::vec x);
RcppExport SEXP mcemGLM_getSigma(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP );
        /** Reconstruct a nxn covariance matrix from a vector */ arma::mat __result = getSigma(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// qFunctionCpp_t
double qFunctionCpp_t(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);
RcppExport SEXP mcemGLM_qFunctionCpp_t(SEXP betaSEXP, SEXP sigmaSEXP, SEXP sigmaTypeSEXP, SEXP uSEXP, SEXP dfSEXP, SEXP kKiSEXP, SEXP kLhSEXP, SEXP kLhiSEXP, SEXP kYSEXP, SEXP kXSEXP, SEXP kZSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type sigmaType(sigmaTypeSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type df(dfSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kKi(kKiSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kLh(kLhSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kLhi(kLhiSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kY(kYSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kX(kXSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kZ(kZSEXP );
        double __result = qFunctionCpp_t(beta, sigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// qFunctionCpp_n
double qFunctionCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::mat& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ);
RcppExport SEXP mcemGLM_qFunctionCpp_n(SEXP betaSEXP, SEXP sigmaSEXP, SEXP sigmaTypeSEXP, SEXP uSEXP, SEXP kYSEXP, SEXP kXSEXP, SEXP kZSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type sigmaType(sigmaTypeSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type u(uSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kY(kYSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kX(kXSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kZ(kZSEXP );
        double __result = qFunctionCpp_n(beta, sigma, sigmaType, u, kY, kX, kZ);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP mcemGLM_rcpp_hello_world() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = rcpp_hello_world();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// uSamplerCpp
arma::mat uSamplerCpp(const arma::vec& beta, const arma::mat& sigma, const arma::vec& sigmaType, const arma::vec& u, const arma::vec& df, const arma::vec& kKi, const arma::vec& kLh, const arma::vec& kLhi, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0);
RcppExport SEXP mcemGLM_uSamplerCpp(SEXP betaSEXP, SEXP sigmaSEXP, SEXP sigmaTypeSEXP, SEXP uSEXP, SEXP dfSEXP, SEXP kKiSEXP, SEXP kLhSEXP, SEXP kLhiSEXP, SEXP kYSEXP, SEXP kXSEXP, SEXP kZSEXP, SEXP BSEXP, SEXP sd0SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type sigmaType(sigmaTypeSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type df(dfSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kKi(kKiSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kLh(kLhSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kLhi(kLhiSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kY(kYSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kX(kXSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kZ(kZSEXP );
        Rcpp::traits::input_parameter< int >::type B(BSEXP );
        Rcpp::traits::input_parameter< double >::type sd0(sd0SEXP );
        arma::mat __result = uSamplerCpp(beta, sigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ, B, sd0);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// uSamplerCpp_n
arma::mat uSamplerCpp_n(const arma::vec& beta, const arma::mat& sigma, const arma::vec& u, const arma::vec& kY, const arma::mat& kX, const arma::mat& kZ, int B, double sd0);
RcppExport SEXP mcemGLM_uSamplerCpp_n(SEXP betaSEXP, SEXP sigmaSEXP, SEXP uSEXP, SEXP kYSEXP, SEXP kXSEXP, SEXP kZSEXP, SEXP BSEXP, SEXP sd0SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::vec& >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type kY(kYSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kX(kXSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type kZ(kZSEXP );
        Rcpp::traits::input_parameter< int >::type B(BSEXP );
        Rcpp::traits::input_parameter< double >::type sd0(sd0SEXP );
        arma::mat __result = uSamplerCpp_n(beta, sigma, u, kY, kX, kZ, B, sd0);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
