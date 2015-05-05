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

mcemMLENegBinom_t_fixed_df <- function (sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM, controlTrust) {  
  # Number of fixed effects, random effects, variance and subvariance components.
  kP <- ncol(kX)
  kK <- ncol(kZ)
  kR <- length(kKi)
  u <- rep(0, kK)
  kL <- sum(kLh)
  
  # Parameters needed in sigma, one for diagonal, two for exchangeable and AR(1).
  if (!missing(initial)) {
    beta <- initial[1:kP]
    alpha <- initial[kP + 1]
    sigma <- initial[-(1:(kP + 1))]
  } else {
    beta <- rep(0, kP)
    alpha <- 1
    sigma <- NULL
    for (i in sigmaType) {
      if (i == 0) {
        sigma <- c(sigma, 5)
      } else {
        sigma <- c(sigma, 5, 0.1)
      }
    }
  }
  theta <- c(beta, alpha, sigma)
  ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
  
  outMLE <- matrix(0, controlEM$EMit, length(theta))
  outMLE[1, ] <- theta
  
  # MCMC step size tuning
  if (is.null(controlEM$MCsd)) {
    if (controlEM$verb == TRUE)
      print("Tuning acceptance rate.")
    ar <- 1
    sdtune <- 1
    u <- rmvnorm(1, rep(0, kK), ovSigma)
    while (ar > 0.4 | ar < 0.1) {
      uSample <- uSamplerNegBinomCpp_t(beta = beta, sigma = ovSigma, alpha = alpha, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = 1000, sd0 = sdtune)
      ar <- length(unique(uSample[, 1])) / 1000
      if (ar < 0.1)
        sdtune <- 0.8 * sdtune
      if (ar > 0.4)
        sdtune <- 1.2 * sdtune
    }
    if (controlEM$verb == TRUE)
      print(ar)
    controlEM$MCsd <- sdtune
  }
  
  # EM iterations
  j <- 2
  errorCounter <- 0
  while (j <= controlEM$EMit & sum(tail(errorCounter, 3)) < 3) {
    # Obtain MCMC sample for u with the current parameter estimates.
    u <- rmvnorm(1, rep(0, kK), ovSigma) # Initial value for u
    uSample <- uSamplerNegBinomCpp_t(beta = beta, sigma = ovSigma, alpha = alpha, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = controlEM$MCit, sd0 = controlEM$MCsd)
    
    # Now we optimize.
    if(sum(sigmaType) == 0) {
      outTrust <- trust(toMaxDiagNegBinom_t, parinit = theta, rinit = controlTrust$rinit, rmax = controlTrust$rmax, iterlim = controlTrust$iterlim, minimize = FALSE, u = uSample, sigmaType = sigmaType, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    } else {
      stop("Not implemented yet.")
      # outTrust <- trust(toMax_t_fixed_df, parinit = theta, rinit = 10, rmax = 20, minimize = FALSE, u = uSample, sigmaType = sigmaType, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    }
    if (controlEM$verb == TRUE)
      print(outTrust)
    outMLE[j, ] <- outTrust$argument
    
    # The current estimates are updated now
    beta <- outMLE[j, 1:kP]
    alpha <- outMLE[j, kP + 1]
    sigma <- outMLE[j, -c(1:(kP + 1))]
    theta <- c(beta, alpha, sigma)
    ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
    if (controlEM$verb == TRUE) {
      print(theta)
      print(ts.plot(uSample[, sample(1:kK, 1)]))
    }
    
    # Retuning the accepatance rate.
    ar <- length(unique(uSample[, 1]))/controlEM$MCit
    if (ar < 0.1 | ar > 0.4) {
      if (controlEM$verb == TRUE)
        print("Tuning acceptance rate.")
      ar <- 1
      sdtune <- controlEM$MCsd
      u <- rmvnorm(1, rep(0, kK), ovSigma)
      while (ar > 0.4 | ar < 0.1) {
        uSample <- uSamplerNegBinomCpp_t(beta = beta, sigma = ovSigma, alpha = alpha, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = 1000, sd0 = sdtune)
        ar <- length(unique(uSample[, 1])) / 1000
        if (ar < 0.1)
          sdtune <- 0.8 * sdtune
        if (ar > 0.4)
          sdtune <- 1.2 * sdtune
      }
      if (controlEM$verb == TRUE)
        print(ar)
      controlEM$MCsd <- sdtune
    }
    
    # We modify the number of MCMC iterations
    controlEM$MCit <- controlEM$MCit * controlEM$MCf
    
    # Error checking
    error <- max(abs(outMLE[j, ] - outMLE[j - 1, ])/(abs(outMLE[j, ]) + controlEM$EMdelta))
    if (error < controlEM$EMepsilon) {
      errorCounter <- c(errorCounter, 1)
    } else {
      errorCounter <- c(errorCounter, 0)
    }
    j <- j + 1
  }
  # Estimation of the information matrix.
  ovSigma <- constructSigma(pars = sigma, sigmaType = sigmaType, kK = kK, kR = kR, kLh = kLh, kLhi = kLhi)
  if (sum(sigmaType) == 0) {
    iMatrix <- iMatrixDiagNegBinomCpp_t(beta = beta, sigma = ovSigma, alpha = alpha, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = controlEM$MCit, sd0 = controlEM$MCsd)
  } else {
    stop("Not implemented yet.")
    # uSample <- uSamplerCpp(beta = beta, sigma = ovSigma, sigmaType = sigmaType, u = u, df = df, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ, B = controlEM$MCit, sd0 = controlEM$MCsd)
    # iMatrix <- matrix(0, length(theta), length(theta))
    # for (i in 1:controlEM$MCit) {
    #   h0 <- hessianLogit_t_fixed_df(pars = theta,df = df, u = uSample[i, ], sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    #   g0 <- gradientLogit_t_fixed_df(pars = theta, df= df, u = uSample[i, ], sigmaType = sigmaType, kKi = kKi, kLh = kLh, kLhi = kLhi, kY = kY, kX = kX, kZ = kZ)
    #   iMatrix <-  iMatrix + (h0 - g0 %*% t(g0)) / controlEM$MCit
    # }
  }
  
  colnames(uSample) <- colnames(kZ)
  return(list(mcemEST = outMLE, iMatrix = -iMatrix, randeff = uSample, y = kY, x = kX, z = kZ))
}