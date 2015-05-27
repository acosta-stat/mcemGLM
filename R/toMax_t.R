# This function is used with optim to maximize the Q function.
# The paraters in 'pars' are: beta, df, and sigma.
# u:          Matrix of MCMC ouput for the random effects.
# sigmaType:  Type of each covariance matrix:
#             0 - Diagonal
#             1 - Exchangeable
#             2 - AR(1)
# kKi:        Dimension of each variance component vector
# kLh:        Number of subvariance components within each variance components. The
#             subvariance components share a covariance matrix but have different
#             degrees of freedom.
# KLhi:       Number of random effects in each subvariance component
# kY, kX, kZ: Data and design matrices
toMax_t <- function(pars, u, sigmaType, kKi, kLh, kLhi, kY, kX, kZ) {
  kP <- dim(kX)[2]  # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kK <- ncol(kZ)    # Number of random effects in the model
  kL <- sum(kLh)    # Number of subvariance components
  
  beta <- pars[1:kP]
  df <- pars[(kP + 1):(kP + kL)]
  s0 <- length(pars[-(1:(kP + kL))]) # Number of variance parameters
  ovSigma <- constructSigma(pars[-(1:(kP + kL))], sigmaType, kK, kR, kLh, kLhi)
    
  value <- qFunctionCpp_t(beta, ovSigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ)
  
  
  loglikelihoodSigma <- function(pars, u) {
    ovSigma <- constructSigma(pars = pars[-(1:kL)], sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    return(logMarginalCpp_t(sigma = ovSigma, df = pars[1:kL], sigmaType = sigmaType, u = u, kKi = kKi, kLh = kLh, kLhi = kLhi))
  }
  # The gradient can be separated into Beta and Sigma. The gradient of Beta can be computed analytically.
  # The marginal using normal and t distribution are the same so we can use the 'n' function.
  # The functions for "Sigma" actually contains the df
  gradientBeta <- rep(0, kP)
  gradientSigma <- rep(0, kL + s0)
  hessianBeta <- matrix(0, kP, kP)
  hessianSigma <- matrix(0, kL + s0, kL + s0)
  for (i in 1:nrow(u)) {
    gradientBeta <- gradientBeta + loglikelihoodLogitGradientBetaCpp_n(beta = beta, u = u[i, ], kY = kY, kX = kX, kZ = kZ)
    gradientSigma <- gradientSigma + grad(func = loglikelihoodSigma, x = pars[-(1:kP)], u = u[i, ])
    hessianBeta <- hessianBeta + loglikelihoodLogitHessianBetaCpp_n(beta = beta, u = u[i, ], kY = kY, kX = kX, kZ = kZ)
    hessianSigma <- hessianSigma + hessian(func = loglikelihoodSigma, x = pars[-(1:kP)], u = u[i, ])
  }
  hessian0 <- matrix(0, length(pars), length(pars))
  gradient0 <- c(gradientBeta, gradientSigma) / nrow(u)
  hessian0[1:kP, 1:kP] <- hessianBeta / nrow(u)
  hessian0[(kP + 1):(kP + kL + s0), (kP + 1):(kP + kL + s0)] <- hessianSigma / nrow(u)
  
  return(list(value = value, gradient = gradient0, hessian = hessian0))  
}
