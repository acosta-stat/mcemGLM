# This function starts the estimation of the parameters. It runs for a fixed number of EM iterations.
# sigmaType:  Structure of the sigma matrices in the model.
# kKi:        Number of random effects per variance component.
# kLh:        Number of subvariance components in each variance component.
# kLhi:       Number of random effects in each subvariance component.
# kY, kX, kZ: Data and design matrices.
# EMit:       Number of EM iterations.
# MCit:       Number of intial MCMC iterations
# MCf:        Factor to increase the number of MCMC iterations.
# MCsd:       Standard deviation for the proposal step.

mcemMLE_t <- function (sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial = NULL, controlEM = list(), methodOptim = "Nelder-Mead", controlOptim = list()) {
  ctrl <- list(EMit = 10, MCit = 1000, MCf = 1.04, verb = TRUE, MCsd = 0.2)
  ctrlN <- names(ctrl)
  ctrl[(controlN <- names(controlEM))] <- controlEM
  if(length(unkwn <- controlN[!controlN %in% ctrlN])){
    warning("Unknown names in control: ", paste(unkwn, collapse = ", "))
  }
  
  # Number of fixed effects, random effects, variance and subvariance components.
  kP <- ncol(kX)
  kK <- ncol(kZ)
  kR <- length(kKi)
  u <- rep(0, kK)
  kL <- sum(kLh)
  
  # Parameters needed in sigma, one for diagonal, two for exchangeable and AR(1).
  if (!is.null(initial)) {
    beta <- initial[1:kP]
    sigma <- initial[-(1:kP)]
  } else {
    beta <- rep(0, kP)
    df <- rep(20, kL)
    sigma <- NULL
    for (i in sigmaType) {
      if (i == 0) {
        sigma <- c(sigma, 1)
      } else {
        sigma <- c(sigma, 1, 1)
      }
    }
  }
  
  loglikeVal <- NULL
  theta <- c(beta, df, sigma)
  outMLE <- matrix(0, ctrl$EMit, length(theta))
  outMLE[1, ] <- theta
  
  # (pars, u, sigmaType, kKi, kLh, kLhi, kY, kX, kZ)
  for (j in 2:ctrl$EMit) {
    # Obtain MCMC sample for u with the current parameter estimates. We need to give it the sigma matrix in the 'compact form'.
    
    sigmaMat <- constructSigma(sigma, sigmaType, kK, kR, kLh, kLhi)
    
    # The matrix input is:
    # print("Matrix:")
    # print(sigmaMat)
    uSample <- uSamplerCpp(beta = beta, sigma = sigmaMat, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = ctrl$MCit, sd0 = ctrl$MCsd)
    
    # Now we optimize.
    outTrust <- trust(toMax_t, parinit = theta, rinit = 10, rmax = 20, minimize = FALSE, u = uSample, sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    print(outTrust)
    outMLE[j, ] <- outTrust$argument
    loglikeVal <- c(loglikeVal, outTrust$value)
    
    # The current estimates are updated now
    beta <- outMLE[j, 1:kP]
    df <- outMLE[j, (kP + 1):(kP + kL)]
    sigma <- outMLE[j, -c(1:(kP + kL))]
    theta <- c(beta, df, sigma)
    if (ctrl$verb == TRUE) {
      print(theta)
      print(ts.plot(uSample[, sample(1:kK, 1)]))
    }
    
    # The starting value for the next MCMC run is the mean of the previous iteration.
    u <- colMeans(uSample)
    
    # We modify the number of MCMC iterations
    ctrl$MCit <- ctrl$MCit * ctrl$MCf
  }
  
  # Get a final sample from U using the last MLE estimates to estimate the information matrix.
  uSample <- uSamplerCpp(beta = beta, sigma = sigmaMat, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = ctrl$MCit, sd0 = ctrl$MCsd)
  iMatrix <- matrix(0, length(theta), length(theta))
  for (i in 1:ctrl$MCit) {
    h0 <- hessianLogit_t(pars = theta, u = uSample[i, ], sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    g0 <- gradientLogit_t(pars = theta, u = uSample[i, ], sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    iMatrix <-  iMatrix + (h0 - g0 %*% t(g0)) / ctrl$MCit
  }
  
  return(list(mcemEST = outMLE, iMatrix = -iMatrix, randeff = uSample, loglikeVal = loglikeVal))
}