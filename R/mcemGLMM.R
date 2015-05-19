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

mcemGLMM <- function(fixed, random, data, family = c("bernoulli", "poisson"), vcDist = c("t", "normal"), df, corType, controlEM = list(), controlTrust = list(), initial) {
  # Reading Y and X.
  if (missing(data)) {
    kY <- get_all_vars(fixed)[, 1]
    kX <- model.matrix(fixed)
  } else {
    kY <- data[, all.vars(fixed)[1]]
    kX <- model.matrix(fixed, data = data)
  }
  if (attr(terms(fixed), "intercept") == 0) {
    tnames <- attr(terms(fixed), "term.labels")
  } else {
    tnames <- c("(Intercept)", attr(terms(fixed), "term.labels"))
  }
  
  # Drop columns with only zeros. Keep attributes.
  kX0 <- kX[, colSums(kX^2) !=0]
  attr(kX0, "assign") <- attr(kX, "assign")[colSums(kX^2) !=0]
  kX <- kX0
  xlabs <- colnames(kX)
  
  # Options
  ctrl <- list(EMit = 50, MCit = 5000, MCf = 1.03, verb = TRUE, MCsd = NULL, EMdelta = 0.02, EMepsilon = 0.01, utrust = TRUE)
  ctrlN <- names(ctrl)
  ctrl[(controlN <- names(controlEM))] <- controlEM
  if(length(unkwn <- controlN[!controlN %in% ctrlN])){
    warning("Unknown names in control: ", paste(unkwn, collapse = ", "))
  }
  
  # Options for trust
  cTrust <- list(rinit = 20, rmax = 200, iterlim = 100)
  cTrustNames <- names(cTrust)
  cTrust[(controlN <- names(controlTrust))] <- controlTrust
  if(length(unkwn <- controlN[!controlN %in% cTrustNames])){
    warning("Unknown names in control: ", paste(unkwn, collapse = ", "))
  }
  
  if (missing(initial)) {
    if (family == "bernoulli") {
      if(!missing(data)) {
        initial0 <- glm(fixed, family = binomial, data = data)$coefficients
      } else {
        initial0 <- glm(fixed, family = binomial)$coefficients
      }
    }
    if (family == "poisson") {
      if(!missing(data)) {
        initial0 <- glm(fixed, family = poisson, data = data)$coefficients
      } else {
        initial0 <- glm(fixed, family = poisson)$coefficients
      }
    }
    if (family == "negbinom") {
      if(!missing(data)) {
        initial0 <- c(glm(fixed, family = poisson, data = data)$coefficients, 100)
      } else {
        initial0 <-c(glm(fixed, family = poisson)$coefficients, 100)
      }
    }
    
    
    if(!is.list(random)) {
      initial <- c(initial0, 5)
    } else {
      initial <- c(initial0, rep(5, length(random)))
    }
  }
  
  # Missing corType corresponds to the diagonal covariance matrix cases.
  if (missing(corType)) {
    if (!is.list(random)) {
      # Case 1: One random effect with a diagonal matrix.
      if (missing(data)) {
        kZ <- model.matrix(random)
      } else {
        kZ <- model.matrix(random, data = data)
      }
      zlabs <- attr(terms(random), "term.labels")
      sigmaType <- 0
      kKi <- ncol(kZ)
      kLh <- 1
      kLhi <- kKi
      
      if (vcDist == "normal") {
        if (family == "bernoulli") {
          fit0 <- mcemMLE_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
        if (family == "poisson") {
          fit0 <- mcemMLEPoisson_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
        if (family == "negbinom") {
          fit0 <- mcemMLENegBinom_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
      } else {
        if (length(df) > 1) {
          stop("The number of variance components and the length of df must me equal.")
        }
        if (family == "bernoulli") {
          fit0 <- mcemMLE_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        } 
        if (family == "poisson") {
          fit0 <- mcemMLEPoisson_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
        if (family == "negbinom") {
          fit0 <- mcemMLENegBinom_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
      }
      # return(fit0)
    } else {
      # Case 2: Possible multiple random effects with diagonal matrices.
      sigmaType <- rep(0, length(random))
      kZ <- NULL
      kKi <- NULL
      kLh <- NULL
      for (i in 1:length(random)) {
        if (missing(data)) {
          tmpZ <- model.matrix(random[[i]])
        } else {
          tmpZ <- model.matrix(random[[i]], data = data)
        }
        tmpZ <- tmpZ[, colSums(tmpZ^2) !=0] # Check for columns with only zeros. This might happen with nested random effects.
        kZ <- cbind(kZ, tmpZ)
        kKi <- c(kKi, ncol(tmpZ))
        kLh <- c(kLh, 1)        
      }
      kLhi <- kKi
      zlabs <- as.character(lapply(lapply(random, terms), attr, which="term.labels"))
      # Normal random effects
      if (vcDist == "normal") {
        if (family == "bernoulli") {
          fit0 <- mcemMLE_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        } 
        if (family == "poisson") {
          fit0 <- mcemMLEPoisson_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
        if (family == "negbinom") {
          fit0 <- mcemMLENegBinom_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
      } 
      # t random effects
      if (vcDist == "t") {
        if (length(df) != length(random)) {
          stop("The number of variance components and the length of df must me equal.")
        }
        if (family == "bernoulli") {
          fit0 <- mcemMLE_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        } 
        if (family == "poisson") {
          fit0 <- mcemMLEPoisson_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
        if (family == "negbinom") {
          fit0 <- mcemMLENegBinom_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, initial, controlEM = ctrl, controlTrust = cTrust)
        }
      }
    }
    fit0$tnames <- tnames
    class(fit0) <- "mcemGLMM"
    if (family != "negbinom") {
      colnames(fit0$mcemEST) <- c(xlabs, zlabs)
    } else {
      colnames(fit0$mcemEST) <- c(xlabs, "alpha", zlabs)
    }
    fit0$mcemEST <- fit0$mcemEST[rowSums(fit0$mcemEST^2) !=0, ]
    fit0$family <- family
    colnames(fit0$x) <- colnames(kX)
    if (det(fit0$iMatrix) < .Machine$double.eps) {
      warning("Information matrix is not invertible.")
    }
    return(fit0)
  } else {
    # Different correlation types.
    return("Not yet.")
  }
  
}