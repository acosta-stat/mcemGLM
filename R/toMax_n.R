# This function is used with optim to maximize the Q function.
# The paraters in 'pars' are: beta and sigma.
# df:         kL dimensional vector of degrees of freedom.
# u:          Matrix of MCMC ouput for the random effects.
# sigmaType:  Type of each covariance matrix:
#             0 - Diagonal
#             1 - Exchangeable
#             2 - AR(1)
# sigmaDim:   Dimensions of the sigma matrices.
# kKi:        Dimension of each variance component vector. Length equal to kR.
# kLh:        Number of subvariance components within each variance components. The
#             subvariance components share a covariance matrix. Length equal to kR.
# KLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
toMax_n <- function(pars, u, sigmaType, kKi, kLh, kLhi, kY, kX, kZ) {
  kP <- ncol(kX)    # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kK <- ncol(kZ)    # Number of random effects
  kL <- sum(kLh)    # Number of subvariance components
  
  beta <- pars[1:kP]
  s0 <- length(pars[-(1:kP)]) # Number of variance parameters
  
  # We call ovSigma the overall covariance matrix.
  ovSigma <- constructSigma(pars = pars[-(1:kP)], sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
  
  if (min(eigen(ovSigma)$values) <= 0) {
    return(list(value = -Inf, gradient = rep(0, length(pars)), hessian = matrix(0, length(pars), length(pars))))
  }
  
  # Value of the function
  value <- qFunctionCpp_n(beta, ovSigma, sigmaType, u, kY, kX, kZ)
  
  
  loglikelihoodSigma <- function(pars, u) {
    ovSigma <- constructSigma(pars = pars, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    return(ldmn(x = u, sigma = ovSigma))
  }

  # The gradient can be separated into Beta and Sigma. The gradient of Beta can be computed analytically.
  # The gradient of sigma will be approximated by Richardson's extrapolation using the numDeriv package.
  gradientBeta <- rep(0, kP)
  gradientSigma <- rep(0, s0)
  hessianBeta <- matrix(0, kP, kP)
  hessianSigma <- matrix(0, s0, s0)
  
  for (i in 1:nrow(u)) {
    gradientBeta <- gradientBeta + loglikelihoodLogitGradientBetaCpp_n(beta = beta, u = u[i, ], kY = kY, kX = kX, kZ = kZ)
    gradientSigma <- gradientSigma + grad(func = loglikelihoodSigma, x = pars[-(1:kP)], u = u[i, ])
    hessianBeta <- hessianBeta + loglikelihoodLogitHessianBetaCpp_n(beta = beta, u = u[i, ], kY = kY, kX = kX, kZ = kZ)
    hessianSigma <- hessianSigma + hessian(func = loglikelihoodSigma, x = pars[-(1:kP)], u = u[i, ])
  }
  
  gradient <- c(gradientBeta, gradientSigma) / nrow(u)
  hessian0 <- matrix(0, length(pars), length(pars))
  hessian0[1:kP, 1:kP] <- hessianBeta  / nrow(u)
  hessian0[(kP + 1):(kP + s0), (kP + 1):(kP + s0)] <- hessianSigma  / nrow(u)

  return(list(value = value, gradient = gradient, hessian = hessian0))
  
  # return(list(beta=beta, df=df, sigma=sigma))
}