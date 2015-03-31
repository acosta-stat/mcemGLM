# This function is used with optim to maximize the Q function.
# The paraters in 'pars' are: beta and sigma.
# df:         kL dimensional vector of degrees of freedom.
# u:          Matrix of MCMC ouput for the random effects.
# sigmaType:  Type of each covariance matrix:
#             0 - Diagonal
#             1 - Exchangeable
#             2 - AR(1)
# kKi:        Dimension of each variance component vector.
# kLh:        Number of subvariance components within each variance components. The
#             subvariance components share a covariance matrix but have different
#             degrees of freedom.
# KLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
toMax_t_fixed_df <- function(pars, df, u, sigmaType, sigmaDim, kKi, kLh, kLhi, kY, kX, kZ, utrust = TRUE) {
  kP <- ncol(kX)  # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kK <- ncol(kZ)    # Number of random effects in the model
  kL <- sum(kLh)    # Number of subvariance components
  
  beta <- pars[1:kP]
  s0 <- length(pars[-(1:kP)]) # Number of variance parameters
  ovSigma <- constructSigma_t(pars[-(1:kP)], sigmaType, kK, kR, kLh, kLhi)
  if (min(eigen(ovSigma)$values) <= 0) {
    return(list(value = -Inf, gradient = rep(0, length(pars)), hessian = matrix(0, length(pars), length(pars))))
  }
  
  if (utrust == FALSE) {
    return(-qFunctionCpp_t(beta, sigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ))
  }
  
  value <- qFunctionCpp_t(beta, ovSigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ)
  
  
  loglikelihoodSigma <- function(pars, df, u) {
    ovSigma <- constructSigma_n(pars = pars, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    return(logMarginalCpp_t(sigma = ovSigma, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi))
  }
  
  
  
  gradientBeta <- rep(0, kP)
  gradientSigma <- rep(0, s0)
  hessianBeta <- matrix(0, kP, kP)
  hessianSigma <- matrix(0, s0, s0)
  for (i in 1:nrow(u)) {
    gradientBeta <- gradientBeta + loglikelihoodLogitGradientBetaCpp_n(beta = beta, u = u[i, ], kY = kY, kX = kX, kZ = kZ)
    gradientSigma <- gradientSigma + grad(func = loglikelihoodSigma, x = pars[-(1:kP)],df = df, u = u[i, ])
    hessianBeta <- hessianBeta + loglikelihoodLogitHessianBetaCpp_n(beta = beta, u = u[i, ], kY = kY, kX = kX, kZ = kZ)
    hessianSigma <- hessianSigma + hessian(func = loglikelihoodSigma, x = pars[-(1:kP)], df = df, u = u[i, ])
  }
  
  gradient0 <- c(gradientBeta, gradientSigma) / nrow(u)
  hessian0 <- matrix(0, length(pars), length(pars))
  hessian0[1:kP, 1:kP] <- hessianBeta / nrow(u)
  hessian0[(kP + 1):(kP + s0), (kP + 1):(kP + s0)] <- hessianSigma / nrow(u)

  return(list(value = value, gradient = gradient0, hessian = hessian0))
}
