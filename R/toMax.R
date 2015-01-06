# This function is used in optim to maximize the Q function.
# The paraters in 'pars' are: beta, df, and sigma.
toMax <- function(pars, u, sigmaType, kKi, kLh, kY, kX, kZ) {
  kP <- dim(kX)[2]  # Number of fixed coefficients
  kR <- length(kKi) # Number of variance components, this is the number of sigma matrices
  kL <- sum(kLh)    # Number of subvariance components
  kS <- kKi/kLh     # Dimensions of the sigma matrices
  
  beta <- pars[1:kP]
  df <- pars[(kP + 1):(kP + kL + 1)]
  sigma <- matrix(0, kR, max(kS) + 1)
  
}