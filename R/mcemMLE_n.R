# This function starts the estimation of the parameters. It runs for a fixed number of EM iterations.
# sigmaType:  Structure of the sigma matrices in the model.
# kKi:        Number of random effects per variance component.
# kLh:        Number of subvariance components in each variance component.
# kLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
# EMit:       Number of EM iterations.
# MCit:       Number of intial MCMC iterations
# MCf:        Factor to increase the number of MCMC iterations.

mcemMLE_n <- function (sigmaType, kKi, kLh, kLhi, kY, kX, kZ, EMit, MCit, MCf) {
  # Number of fixed effects, random effects, variance and subvariance components.
  kP <- ncol(kX)
  beta <- rep(0, kP)
  
  kK <- ncol(kZ)
  kR <- length(kKi)
  u <- rep(0, kK)
  kL <- sum(kLh)
  
  # Parameters needed in sigma, one for diagonal, two for exchangeable and AR(1).
  sigma <- NULL
  for (i in sigmaType) {
    if (i == 0) {
      sigma <- c(sigma, 1)
    } else {
      sigma <- c(sigma, 1, 1)
    }
  }
  
  theta <- c(beta, sigma)
  outMLE <- matrix(0, EMit, length(theta))
  outMLE[1, ] <- theta
  

  for (j in 2:EMit) {
    # Obtain MCMC sample for u with the current parameter estimates. We need to give it the sigma matrix in the 'compact form'.
    ovSigma <- constructSigma_n(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    uSample <- uSamplerCpp_n(beta = beta, sigma = ovSigma, u = u, kY = kY, kX = kX, kZ = kZ, B = MCit, sd0 = 1)
    
    # Now we optimize.
    outOptim <- optim(par = theta, fn = toMax_n, control = list(fnscale = -1, maxit = 10000), u = uSample, sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    if(outOptim$convergence != 0) {
      print("Convergence issues")
      print(outOptim)
    }
    outMLE[j, ] <- outOptim$par
    
    # The current estimates are updated now
    beta <- outMLE[j, 1:kP]
    sigma <- outMLE[j, -c(1:kP)]
    theta <- c(beta, sigma)
    
    # The starting value for the next MCMC run is the mean of the previous iteration.
    u <- apply(uSample, 2, mean)
    
    # We modify the number of MCMC iterations
    MCit <- MCit * MCf
  }
  # Get a final sample from U using the last MLE estimates
  uSample <- uSamplerCpp_n(beta = beta, sigma = ovSigma, u = u, kY = kY, kX = kX, kZ = kZ, B = MCit, sd0 = 1)
  
  return(list(mcemEST=outMLE, randeff = uSample))
}