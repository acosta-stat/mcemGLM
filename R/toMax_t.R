# This function is used in optim to maximize the Q function.
# The paraters in 'pars' are: beta, df, and sigma.
# u: Vector composed of kR variance components.
# sigmaType: Type of each covariance matrix:
#            0 - Diagonal
#            1 - Exchangeable
#            2 - AR(1)
# kKi: Dimension of each variance component vector
# kLh: Number of subvariance components within each variance components. The
#      subvariance components share a covariance matrix but have different
#      degrees of freedom.
# kY, kX, kZ: Data and design matrices
toMax_t <- function(pars, u, sigmaType, kKi, kLh, kY, kX, kZ) {
  kP <- dim(kX)[2]  # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kL <- sum(kLh)    # Number of subvariance components
  kS <- kKi/kLh     # Dimensions of the sigma matrices
  
  beta <- pars[1:kP]
  df <- pars[(kP + 1):(kP + kL)]
  sigma <- matrix(-1, kR, 1 + (max(kS))^2)
  counter <- kP + kL + 1
  for (i in 1:kR) {
    sigma[i, 1] <- kS[i]
    if (sigmaType[i] == 0) {
      tmp_mat <- pars[counter] * diag(kS[i])
      sigma[i, 2:(kS[i]^2 + 1)] <- as.vector(tmp_mat)
      counter <- counter + 1
    }
    if (sigmaType[i] == 1) {
      tmp_mat <- pars[counter] * diag(kS[i])
      counter <- counter + 1
      tmp_mat[lower.tri(tmp_mat)] <- pars[counter]
      tmp_mat[upper.tri(tmp_mat)] <- pars[counter]
      counter <- counter + 1
      sigma[i, 2:(kS[i]^2 + 1)] <- as.vector(tmp_mat)
    }
    if (sigmaType[i] == 2) {
      sigma2 <- pars[counter]
      counter <- counter + 1
      pho <- pars[counter]
      counter <- counter + 1
      d0 <- abs(outer(1:kS[i], 1:kS[i], "-"))
      tmp_mat <- sigma2 * pho^d0
      sigma[i, 2:(kS[i]^2 + 1)] <- as.vector(tmp_mat)
    }
  }
  return(qFunctionCpp_t(beta, sigma, sigmaType, u, df, kKi, kLh, kY, kX, kZ))
  # return(list(beta=beta, df=df, sigma=sigma))
}