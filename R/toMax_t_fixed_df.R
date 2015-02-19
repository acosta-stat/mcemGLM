# This function is used with optim to maximize the Q function.
# The paraters in 'pars' are: beta and sigma.
# df:         kL dimensional vector of degrees of freedom.
# u:          Matrix of MCMC ouput for the random effects.
# sigmaType:  Type of each covariance matrix:
#             0 - Diagonal
#             1 - Exchangeable
#             2 - AR(1)
# sigmaDim:   Dimensions of the sigma matrices.
# kKi:        Dimension of each variance component vector.
# kLh:        Number of subvariance components within each variance components. The
#             subvariance components share a covariance matrix but have different
#             degrees of freedom.
# KLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
toMax_t_fixed_df <- function(pars, df, u, sigmaType, sigmaDim, kKi, kLh, kLhi, kY, kX, kZ) {
  kP <- dim(kX)[2]  # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kL <- sum(kLh)    # Number of subvariance components
  
  beta <- pars[1:kP]
  # df <- pars[(kP + 1):(kP + kL)]
  sigma <- matrix(-1, kR, 1 + (max(sigmaDim))^2)
  counter <- kP + 1
  for (i in 1:kR) {
    sigma[i, 1] <- sigmaDim[i]
    if (sigmaType[i] == 0) {
      # Diagonal matrix. One parameter
      tmp_mat <- pars[counter] * diag(sigmaDim[i])
      sigma[i, 2:(sigmaDim[i]^2 + 1)] <- as.vector(tmp_mat)
      counter <- counter + 1
    }
    if (sigmaType[i] == 1) {
      # Exchangeable matrix. Two parameters, diagonal and off-diagonal.
      tmp_mat <- pars[counter] * diag(sigmaDim[i])
      counter <- counter + 1
      tmp_mat[lower.tri(tmp_mat)] <- pars[counter]
      tmp_mat[upper.tri(tmp_mat)] <- pars[counter]
      counter <- counter + 1
      sigma[i, 2:(sigmaDim[i]^2 + 1)] <- as.vector(tmp_mat)
    }
    if (sigmaType[i] == 2) {
      # AR(1) matrix. Two parameters, sigma^2 and pho.
      sigma2 <- pars[counter]
      counter <- counter + 1
      pho <- pars[counter]
      counter <- counter + 1
      d0 <- abs(outer(1:sigmaDim[i], 1:sigmaDim[i], "-"))
      tmp_mat <- sigma2 * pho^d0
      sigma[i, 2:(sigmaDim[i]^2 + 1)] <- as.vector(tmp_mat)
    }
  }
  return(qFunctionCpp_t(beta, sigma, sigmaType, u, df, kKi, kLh, kLhi, kY, kX, kZ))
  # return(list(beta=beta, df=df, sigma=sigma))
}