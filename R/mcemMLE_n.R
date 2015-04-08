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

mcemMLE_n <- function (sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial = NULL, controlEM = list(), controlTrust = list(), methodOptim = "Nelder-Mead", controlOptim = list()) {
  ctrl <- list(EMit = 3, MCit = 1000, MCf = 1.04, verb = TRUE, MCsd = NULL, utrust = TRUE)
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
    sigma <- NULL
    for (i in sigmaType) {
      if (i == 0) {
        sigma <- c(sigma, 1)
      } else {
        sigma <- c(sigma, 1, 0.1)
      }
    }
  }
  theta <- c(beta, sigma)
  
  outMLE <- matrix(0, ctrl$EMit, length(theta))
  outMLE[1, ] <- theta
  
  if (is.null(ctrl$MCsd)) {
    if (ctrl$verb == TRUE)
      print("Tuning acceptance rate.")
    ar <- 1
    sdtune <- 1
    ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    u <- rmvnorm(1, u, ovSigma)
    while (ar > 0.4 | ar < 0.1) {
      uSample <- uSamplerCpp_n(beta = beta, sigma = ovSigma, u = u, kY = kY, kX = kX, kZ = kZ, B = 1000, sd0 = sdtune)
      ar <- length(unique(uSample[, 1])) / 1000
      if (ar < 0.1)
        sdtune <- 0.8 * sdtune
      if (ar > 0.4)
        sdtune <- 1.2 * sdtune
    }
    if (ctrl$verb == TRUE)
      print(ar)
    ctrl$MCsd <- sdtune
  }

  for (j in 2:ctrl$EMit) {
    # Obtain MCMC sample for u with the current parameter estimates. We need to give it the sigma matrix in the 'compact form'.
    ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    if (j == 2){
      u <- rmvnorm(1, u, ovSigma)
    }
    uSample <- uSamplerCpp_n(beta = beta, sigma = ovSigma, u = u, kY = kY, kX = kX, kZ = kZ, B = ctrl$MCit, sd0 = ctrl$MCsd)
    
    
    
    
    
    
    # Now we optimize.
    if (ctrl$utrust == FALSE) {
      outOptim <- optim(par = theta, fn = toMax_n, method = methodOptim, control = controlOptim, u = uSample, sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, utrust = FALSE)    
      if(outOptim$convergence != 0) {
        print("Convergence issues")
        print(outOptim)
      }
      outMLE[j, ] <- outOptim$par
    } else {
      if (sum(sigmaType == 0)) {
        outTrust <- trust(toMaxDiag_n, parinit = theta, rinit = 10, rmax = 20, minimize = FALSE, u = uSample, sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
      } else {
        outTrust <- trust(toMax_n, parinit = theta, rinit = 10, rmax = 20, minimize = FALSE, u = uSample, sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
      }
      if (ctrl$verb == TRUE)
        print(outTrust)
      outMLE[j, ] <- outTrust$argument
    }
    
    # The current estimates are updated now
    beta <- outMLE[j, 1:kP]
    sigma <- outMLE[j, -c(1:kP)]
    theta <- c(beta, sigma)
    if (ctrl$verb == TRUE) {
      print(theta)
      print(ts.plot(uSample[, sample(1:kK, 1)]))
    }
    
    # Retuning the accepatance rate.
    ar <- length(unique(uSample[, 1]))/ctrl$MCit
    if (ar < 0.1 | ar > 0.4) {
      if (ctrl$verb == TRUE)
        print("Tuning acceptance rate.")
      ar <- 1
      sdtune <- ctrl$MCsd
      ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
      u <- rmvnorm(1, u, ovSigma)
      while (ar > 0.4 | ar < 0.1) {
        uSample <- uSamplerCpp_n(beta = beta, sigma = ovSigma, u = u, kY = kY, kX = kX, kZ = kZ, B = 1000, sd0 = sdtune)
        ar <- length(unique(uSample[, 1])) / 1000
        if (ar < 0.1)
          sdtune <- 0.9 * sdtune
        if (ar > 0.4)
          sdtune <- 1.1 * sdtune
      }
      if (ctrl$verb == TRUE)
        print(ar)
      ctrl$MCsd <- sdtune
    }
    
    # The starting value for the next MCMC run is the mean of the previous iteration.
    u <- colMeans(uSample)
    # We modify the number of MCMC iterations
    ctrl$MCit <- ctrl$MCit * ctrl$MCf
  }
  #Estimation of the information matrix.
  ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
  if (sum(sigmaType) == 0) {
    iMatrix <- iMatrixDiagCpp_n(beta = beta, sigma = ovSigma, u = u, kKi = kKi, kY = kY, kX = kX, kZ = kZ, B = ctrl$MCit, sd0 = ctrl$MCsd)
  } else {
    uSample <- uSamplerCpp_n(beta = beta, sigma = ovSigma, u = u, kY = kY, kX = kX, kZ = kZ, B = ctrl$MCit, sd0 = ctrl$MCsd)
    iMatrix <- matrix(0, length(theta), length(theta))
    for (i in 1:ctrl$MCit) {
      h0 <- hessianLogit_n(pars = theta, u = uSample[i, ], sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
      g0 <- gradientLogit_n(pars = theta, u = uSample[i, ], sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
      iMatrix <-  iMatrix + (h0 - g0 %*% t(g0)) / ctrl$MCit
    }
  }
  colnames(uSample) <- colnames(kZ)
  return(list(mcemEST = outMLE, iMatrix = -iMatrix, randeff = uSample, y = kY, x = kX, z = kZ))
}