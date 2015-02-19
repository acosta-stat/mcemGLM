# This function starts the estimation of the parameters. It runs for a fixed number of EM iterations.
# sigmaType:  Structure of the sigma matrices in the model.
# sigmaDim:   Dimension of the sigma matrices in the model.
# kKi:        Number of random effects per variance component.
# kLh:        Number of subvariance components in each variance component.
# kLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
# EMit:       Number of EM iterations.
# MCit:       Number of intial MCMC iterations
# MCf:        Factor to increase the number of MCMC iterations.

mcemMLE_t_fixed_df <- function (sigmaType, sigmaDim, df, kKi, kLh, kLhi, kY, kX, kZ, EMit, MCit, MCf) {
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
  
  # (pars, u, sigmaType, sigmaDim, kKi, kLh, kLhi, kY, kX, kZ)
  for (j in 2:EMit) {
    # Obtain MCMC sample for u with the current parameter estimates. We need to give it the sigma matrix in the 'compact form'.
    sigmaMat <- matrix(-1, kR, 1 + (max(sigmaDim))^2)
    counter <- 1
    for (i in 1:kR) {
      sigmaMat[i, 1] <- sigmaDim[i]
      if (sigmaType[i] == 0) {
        # Diagonal matrix. One parameter
        tmp_mat <- sigma[counter] * diag(sigmaDim[i])
        sigmaMat[i, 2:(sigmaDim[i]^2 + 1)] <- as.vector(tmp_mat)
        counter <- counter + 1
      }
      if (sigmaType[i] == 1) {
        # Exchangeable matrix. Two parameters, diagonal and off-diagonal.
        tmp_mat <- sigma[counter] * diag(sigmaDim[i])
        counter <- counter + 1
        tmp_mat[lower.tri(tmp_mat)] <- sigma[counter]
        tmp_mat[upper.tri(tmp_mat)] <- sigma[counter]
        counter <- counter + 1
        sigmaMat[i, 2:(sigmaDim[i]^2 + 1)] <- as.vector(tmp_mat)
      }
      if (sigmaType[i] == 2) {
        # AR(1) matrix. Two parameters, sigma^2 and pho.
        sigma2 <- sigma[counter]
        counter <- counter + 1
        pho <- sigma[counter]
        counter <- counter + 1
        d0 <- abs(outer(1:sigmaDim[i], 1:sigmaDim[i], "-"))
        tmp_mat <- sigma2 * pho^d0
        sigmaMat[i, 2:(sigmaDim[i]^2 + 1)] <- as.vector(tmp_mat)
      }
    }
    # The matrix input is:
    # print("Matrix:")
    # print(sigmaMat)
    uSample <- uSamplerCpp(beta = beta, sigma = sigmaMat, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = MCit, sd0 = 1)
    
    # Now we optimize.
    outOptim <- optim(par = theta, fn = toMax_t_fixed_df, control = list(fnscale = -1, maxit = 10000), u = uSample, sigmaType = sigmaType, sigmaDim = sigmaDim, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    if(outOptim$convergence != 0) {
      print("Convergence issues")
      print(outOptim)
    }
    outMLE[j, ] <- outOptim$par
    
    # The current estimates are updated now
    beta <- outMLE[j, 1:kP]
    # df <- outMLE[j, (kP + 1):(kP + kL)]
    sigma <- outMLE[j, -c(1:kP)]
    theta <- c(beta, sigma)
    
    # The starting value for the next MCMC run is the mean of the previous iteration.
    u <- apply(uSample, 2, mean)
    
    # We modify the number of MCMC iterations
    MCit <- MCit * MCf
  }
  # Get a final sample from U using the last MLE estimates
  uSample <- uSamplerCpp(beta = beta, sigma = sigmaMat, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = MCit, sd0 = 1)
  
  return(list(mcemEST=outMLE, randeff = uSample))
}