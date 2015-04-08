# User interface function.
# Input information:
# The random effects have to be stated individually when using more than one. For nested random effects the interaction has to be specified. Right now a data set is required (that will change next.) When using a t distribution, the parameter df must have length equal to the number of variance components in the model.

# From the input we need to extract:
# sigmaType:  Structure of the sigma matrices in the model.
# kKi:        Number of random effects per variance component.
# kLh:        Number of subvariance components in each variance component.
# kLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
# controlEM:
#   EMit:     Number of EM iterations.
#   MCit:     Number of intial MCMC iterations
#   MCf:      Factor to increase the number of MCMC iterations.
#   MCsd:     Standard deviation for the proposal step.

mcemGLMM <- function(fixed, random, data, family, vcDist = c("t", "normal"), df, corType, controlEM = list(), controlTrust = list(), controlOptim = list(), methodOptim = "Nelder-Mead", initial = NULL) {
  kY <- data[, all.vars(fixed)[1]]
  kX <- model.matrix(fixed, data = data)
  if (missing(corType)) {
    if (!is.list(random)) {
      # Case 1: One random effect with a diagonal matrix.
      kZ <- model.matrix(random, data = data)
      sigmaType <- 0
      kKi <- ncol(kZ)
      kLh <- 1
      kLhi <- kKi
      
      if (vcDist == "normal") {
        fit0 <- mcemMLE_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM, controlTrust, methodOptim, controlOptim)
      } else {
        fit0 <- mcemMLE_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM, controlTrust, methodOptim, controlOptim)
      }
      return(fit0)
    } else {
      # Case 2: Multiple random effects with diagonal matrices.
      sigmaType <- rep(0, length(random))
      kZ <- NULL
      kKi <- NULL
      kLh <- NULL
      for (i in 1:length(random)) {
        tmpZ <- model.matrix(random[[i]], data = data)
        tmpZ <- tmpZ[, colSums(tmpZ^2) !=0] # Check for columns with only zeros. This might happen with nested random effects.
        kZ <- cbind(kZ, tmpZ)
        kKi <- c(kKi, ncol(tmpZ))
        kLh <- c(kLh, 1)        
      }
      kLhi <- kKi
      # return(list(sigmaType, kY, kX, kZ, kKi, kLh, kLhi))
      if (vcDist == "normal") {
        fit0 <- mcemMLE_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM, controlTrust, methodOptim, controlOptim)
      } else {
        fit0 <- mcemMLE_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM, controlTrust, methodOptim, controlOptim)
      }
      return(fit0)
    }
  } else {
    # Different correlation types.
    return("Not yet.")
  }
  
}