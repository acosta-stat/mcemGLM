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
  
  sigma <- constructSigma_t(pars[-(1:(kP + kL))], sigmaType, kK, kR, kLh, kLhi)
  return(-qFunctionCpp_t(beta, sigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ))
  # return(list(beta=beta, df=df, sigma=sigma))
}