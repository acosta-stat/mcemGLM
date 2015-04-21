# File : hessianLogit_t_fixed_df.R
# Calculates the hessian of the loglikelihood function at a given point.


hessianLogit_t_fixed_df <- function(pars, df = df, u, sigmaType, kKi, kLh, kLhi, kY, kX, kZ) {
  require(numDeriv)
  kP <- ncol(kX)    # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kK <- ncol(kZ)    # Number of random effects
  kL <- sum(kLh)    # Number of subvariance components
  
  beta <- pars[1:kP]
  s0 <- length(pars[-(1:kP)])
  # The gradient can be separated into Beta and Sigma. The gradient of Beta can be computed analytically.
  # Again we can use the same function for the t and normal distribution for the fixed effects.
  hessianBeta <- loglikelihoodLogitHessianBetaCpp_n(beta = beta, u = u, kY = kY, kX = kX, kZ = kZ)
  # The gradient of sigma will be approximated by Richardson's extrapolation using the numDeriv package.
  # The gradient of sigma will be approximated by Richardson's extrapolation using the numDeriv package.
  loglikelihoodSigma <- function(pars, df, u) {
    ovSigma <- constructSigma(pars = pars, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    return(logMarginalCpp_t(sigma = ovSigma, df = df, sigmaType = sigmaType, u = u, kKi = kKi, kLh = kLh, kLhi = kLhi))
  }
  hessianSigma <- hessian(func = loglikelihoodSigma, x = pars[-(1:kP)], df = df, u = u)
  
  hessian0 <- matrix(0, length(pars), length(pars))
  hessian0[1:kP, 1:kP] <- hessianBeta
  hessian0[(kP + 1):(kP + s0), (kP + 1):(kP + s0)] <- hessianSigma
  return(hessian0)
}
