# File : gradientLogit_t.R
# Calculates the gradient of the loglikelihood function at a given point.


gradientLogit_t <- function(pars, u, sigmaType, kKi, kLh, kLhi, kY, kX, kZ) {
  require(numDeriv)
  kP <- ncol(kX)    # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kK <- ncol(kZ)    # Number of random effects
  kL <- sum(kLh)    # Number of subvariance components
  
  beta <- pars[1:kP]
  # The gradient can be separated into Beta and Sigma. The gradient of Beta can be computed analytically.
  # The marginal using normal and t distribution are the same so we can use the 'n' function.
  gradientBeta <- loglikelihoodLogitGradientBetaCpp_n(beta = beta, u = u, kY = kY, kX = kX, kZ = kZ)
  
  # The gradient of sigma will be approximated by Richardson's extrapolation using the numDeriv package.
  loglikelihoodSigma <- function(pars, u) {
    ovSigma <- constructSigma_n(pars = pars[-(1:kL)], sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    return(logMarginalCpp_t(sigma = ovSigma, df = pars[1:kL], sigmaType = sigmaType, u = u, kKi = kKi, kLh = kLh, kLhi = kLhi))
  }
  gradientSigma <- grad(func = loglikelihoodSigma, x = pars[-(1:kP)], u = u)
  return(c(gradientBeta, gradientSigma))
}
