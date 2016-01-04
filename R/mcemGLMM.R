#' Generalized linear mixed models estimation
#' 
#' Maximum likelihood estimation for logistic, Poisson, negative binomial,
#' and gamma models with random effects using a Monte Carlo EM algorithm.
#' 
#' The function \code{mcemGLMM} allows the fitting of generalized linear mixed
#' models when the random effects are normal or a t distributed.  The supported
#' models are the logistic, Poisson and negative binomial.  The degrees of
#' freedom for the t case must be supplied by the user with a vector in the
#' \code{df} argument. The length of the vector must be equal to the number of
#' variance components. For normal random effects the argument \code{df} does
#' not need to be included.
#' 
#' To fit a model with one random effect a formula must be supplied in the
#' \code{random} argument. Note that it is necessary that the variable is a
#' factor and to specify that there is no intercept for the random part of the
#' model. To use more than one random effects a list of formulas must be
#' supplied. Each member must be formula with in which the variables involved
#' must be factors and it also it is necessary to specify that there is no
#' intercept. To fit crossed random effects each variable must no appear in its
#' own formula. To fit nested random effects a formula with the highest level
#' variable must be specified and each subsequent variable must be specified
#' with an interaction of the variables above it. See examples below.
#' 
#' A note on the negative binomial overdispersion:
#' 
#' The variance for the negative binomial model is set equal to \eqn{(1 +
#' \mu/\alpha)} so we have that there is no overdispersion as \eqn{\alpha} goes
#' to infinity.
#' 
#' The variance for the gamma distributionis set to \eqn{\mu^2/\alpha}. The
#' value of \eqn{\alpha = 1} corresponds to an exponential regression.
#' 
#' Stopping rules and convergence criteria:
#' 
#' The algorithm runs for a maximum of \code{EMit} iterations or until the
#' criteria \deqn{\max_i \left\{ \frac{ |\theta_i^{(t)} -
#' \theta_i^{(t-1)}|}{|\theta_i^{(t)}| + \delta}\right\} <
#' \epsilon}{\max_i|\theta_i^(t)-\theta_i^(t-1)|/(|\theta_i^(t)| + \delta) <
#' \epsilon} is satisfied three times in a row for pre-defined values of
#' \eqn{\epsilon} and \eqn{\delta}. Once this criterion has been achieved two
#' times we increase the Monte Carlo sample size more rapidly to have a better
#' estimation of the model's information matrix. For a detailed discussion on
#' convergence diagnostics see Neath R.C. (2012).  After fitting a model it is
#' recommended to plot the EM estimates at each step to assess convergence.
#' 
#' Control options for EM: \describe{ \item{EMit}{maximum number of EM
#' iterations.} \item{MCit}{initial number of Monte Carlo iterations for the
#' MCMC step.} \item{ranefsam}{logical. If true the full random effect sample 
#' is returned.} \item{MCf}{factor in which the MC iterations increase in each EM
#' iteration.} \item{verb}{logical. If TRUE at each EM iteration the function
#' will print convergence information and a trace plot for on of the random
#' effects. This is useful to assess the performance of the algorithm but it
#' can impact the actual running time.} \item{MCsd}{initial standard deviation
#' for the proposal density of the MCMC step. If zero (default) an auto-tuning
#' step will be performed.} \item{EMdelta}{constant for the EM error
#' assessment.} \item{EMepsilon}{constant for the EM error assessment.} }
#' 
#' Control options for trust, see \code{help(trust)} for more details:
#' \describe{ \item{rinit}{starting trust region radius. Default value set to
#' 20.} \item{rmax}{maximum allowed trust region radius. Default value set to
#' 200.} \item{iterlim}{maximum number of iterations. Default value set to
#' 100.} }
#' 
#' @param fixed the fixed effects model. This is specified by a \code{formula}
#' object.
#' @param random the random effects models. This is specified by a
#' \code{formula} object or a list of \code{formula} objects.  See details
#' below.
#' @param data an optional data frame containing the variables in the model. 
#' If missing, the variables are taken from the current environment.
#' @param family a string indicating the type of model to be fitted.  The
#' options are "bernoulli" for logistic regression, "Poisson" for Poisson count
#' regression, and "negbinom" for negative binomial count regression.
#' @param vcDist a string indicating the distribution of the marginal variance
#' components. The options are "normal" and "t" for normal and t distributed
#' random effects respectively.
#' @param df a vector of degrees of freedom of the random effects when these
#' are t distributed. The length of the vector must be equal to the number of
#' variance components in the model.
#' @param controlEM a list of options for the algorithm. See Details below.
#' @param controlTrust a list of options to be passed to the \code{trust}
#' optimizer. See details below.
#' @param initial optional initial values for the parameters. If missing the
#' initial values for the fixed effects are taken from a generalized linear
#' model fitted without random effects and the initial values for the variance
#' components are set to 5.
#' @return A list of class "mcemGLMM" with the following items: \describe{
#' \item{mcemEST}{a matrix with the value of the maximum likelihood estimators
#' at the end of each EM step.} \item{iMatrix}{Fisher's information matrix.}
#' \item{QfunVal}{ value (up to a constant) of the Q function.  Used to perform
#' likelihood ratio tests.} \item{QfunMCMC}{Q function MCMC sample.}
#' \item{ranefsam}{a sample from the conditional distribution of the random
#' effects given the data and the maximum likelihood estimators.}
#' \item{ranef}{the predicted random effects.}
#' \item{y}{vector of observations.} \item{x}{design matrix for the fixed
#' effects.} \item{z}{design matrix for the random effects.}
#' \item{EMerror}{relative error at the last iteration. See details.}
#' \item{MCsd}{last value for MCMC step size.} \item{call}{original call.}
#' \item{MCit}{last value for MCMC sample size.} \item{call}{original call.} }
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @references Neath, R. C. (2012) On Convergence Properties of the Monte Carlo
#' EM Algorithm In Advances in Modern Statistical Theory and Applications: A
#' Festschrift in Honor of Morris L. Eaton. \emph{Institute of Mathematical
#' Statistics} 43--62
#' @keywords glmm
#' @examples
#' 
#' \donttest{
#' # Data set for a logistic model with one binary fixed effects and two 
#' # possible random effects.
#' # Initial values and MC iterations are given to speed up the examples 
#' # but these are not necessary in general.
#' set.seed(0123210)
#' data(exdata)
#' 
#' # To fit a model with one random effect
#' fit.1 <- mcemGLMM(obs ~ x, random = ~ 0 + z1, data = exdata, 
#'                 family = "bernoulli", vcDist = "normal", 
#'                 controlEM = list(MCit = 10000), 
#'                 initial = c(0.27, -0.13, 0.003))
#' summary(fit.1)
#' 
#' # We can assess convergence by looking at a trace plot of the EM estimates
#' # and the loglikelihood values
#' ts.plot(fit.1$mcemEST)
#' ts.plot(fit.1$QfunVal)
#' 
#' # To fit a model with crossed random effects
#' fit.crossed <- mcemGLMM(obs ~ x, random = list(~ 0 + z1, ~ 0 + z2), 
#'                 data = exdata, 
#'                 family = "bernoulli", vcDist = "normal", 
#'                 controlEM = list(EMit = 10, MCit = 10000), 
#'                 initial = c(0.28, -0.15, 0.001, 0.4))
#' summary(fit.crossed)
#' 
#' 
#' # To fit a model with crossed random effects
#' fit.nested <- mcemGLMM(obs ~ x, random = list(~ 0 + z2, ~ 0 + z2:z1), 
#'                 data = exdata, 
#'                 family = "bernoulli", vcDist = "normal", 
#'                 controlEM = list(EMit = 10, MCit = 10000), 
#'                 initial = c(0.31, -0.15, 0.29, 0.27))
#' summary(fit.nested)
#' 
#' # Fit a Poisson model
#' fit.pois <- mcemGLMM(obs2 ~ x, random = ~ 0 + z1, data = exdata, 
#'                 family = "poisson", vcDist = "normal", 
#'                 controlEM = list(EMit = 10, MCit = 10000), 
#'                 initial = c(1.95, 0.03, 0.003))
#' summary(fit.pois)
#' }
#' 
#' @export
#' 

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

mcemGLMM <- function(fixed, random, data, 
                     family = c("bernoulli", "poisson", "negbinom", "gamma"), 
                     vcDist = c("normal", "t"), df, controlEM = list(), 
                     controlTrust = list(), initial) {
  # Reading Y and X.
  call0 <- match.call()
  if (missing(data)) {
    kY <- get_all_vars(fixed)[, 1]
    kX <- model.matrix(fixed)
  } else {
    kY <- data[, all.vars(fixed)[1]]
    kX <- model.matrix(fixed, data = data)
  }
  
  # Drop columns with only zeros. Keep attributes.
  kX0 <- kX[, colSums(kX^2) !=0]
  attr(kX0, "assign") <- attr(kX, "assign")[colSums(kX^2) !=0]
  kX <- kX0
  xlabs <- colnames(kX)
  
  # Options
  ctrl <- list(EMit = 90, MCit = 2500, ranefsam = FALSE, MCf = 1.05, verb = 0, 
               MCsd = 0, EMdelta = 0.05, EMepsilon = 0.015)
  ctrlN <- names(ctrl)
  ctrl[(controlN <- names(controlEM))] <- controlEM
  if(length(unkwn <- controlN[!controlN %in% ctrlN])){
    warning("Unknown names in control: ", paste(unkwn, collapse = ", "))
  }
  
  # Options for trust
  cTrust <- list(rinit = 10, rmax = 200, iterlim = 100)
  cTrustNames <- names(cTrust)
  cTrust[(controlN <- names(controlTrust))] <- controlTrust
  if(length(unkwn <- controlN[!controlN %in% cTrustNames])){
    warning("Unknown names in control: ", paste(unkwn, collapse = ", "))
  }
  
  # Initial values
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
    if (family == "gamma") {
      if(!missing(data)) {
        initial0 <- c(lm(fixed, data = data)$coefficients, 1)
      } else {
        initial0 <-c(lm(fixed)$coefficients, 1)
      }
    }
    if(!is.list(random)) {
      initial <- c(initial0, 4)
    } else {
      initial <- c(initial0, rep(4, length(random)))
    }
  }
  
  # Fitting starts. Case 1: One random effect not in a list.
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
        fit0 <- mcemMLE_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, 
                          controlEM = ctrl, controlTrust = cTrust)
      }
      if (family == "poisson") {
        fit0 <- mcemMLEPoisson_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial,
                                 controlEM = ctrl, controlTrust = cTrust)
      }
      if (family == "negbinom") {
        fit0 <- mcemMLENegBinom_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial,
                                  controlEM = ctrl, controlTrust = cTrust)
      }
      if (family == "gamma") {
        fit0 <- mcemMLEGamma_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, 
                               controlEM = ctrl, controlTrust = cTrust)
      }
    } else {
      if (length(df) > 1) {
        stop("The number of variance components and the length of df must me equal.")
      }
      if (family == "bernoulli") {
        fit0 <- mcemMLE_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, kZ, 
                                   initial, controlEM = ctrl, controlTrust = cTrust)
      } 
      if (family == "poisson") {
        fit0 <- mcemMLEPoisson_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, 
                                          kX, kZ, initial, controlEM = ctrl, 
                                          controlTrust = cTrust)
      }
      if (family == "negbinom") {
        fit0 <- mcemMLENegBinom_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, 
                                           kX, kZ, initial, controlEM = ctrl, 
                                           controlTrust = cTrust)
      }
      if (family == "gamma") {
        fit0 <- mcemMLEGamma_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, kX, 
                                        kZ, initial, controlEM = ctrl, 
                                        controlTrust = cTrust)
      }
    }
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
      # Check for columns with only zeros. This might happen with nested random effects.
      tmpZ <- tmpZ[, colSums(tmpZ^2) !=0] 
      kZ <- cbind(kZ, tmpZ)
      kKi <- c(kKi, ncol(tmpZ))
      kLh <- c(kLh, 1)        
    }
    kLhi <- kKi
    zlabs <- as.character(lapply(lapply(random, terms), attr, which="term.labels"))
    
    # Normal random effects
    if (vcDist == "normal") {
      if (family == "bernoulli") {
        fit0 <- mcemMLE_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, 
                          controlEM = ctrl, controlTrust = cTrust)
      } 
      if (family == "poisson") {
        fit0 <- mcemMLEPoisson_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, 
                                 initial, controlEM = ctrl, controlTrust = cTrust)
      }
      if (family == "negbinom") {
        fit0 <- mcemMLENegBinom_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, 
                                  initial, controlEM = ctrl, controlTrust = cTrust)
      }
      if (family == "gamma") {
        fit0 <- mcemMLEGamma_n(sigmaType, kKi, kLh, kLhi, kY, kX, kZ, initial, 
                               controlEM = ctrl, controlTrust = cTrust)
      }
    }
    
    # t random effects
    if (vcDist == "t") {
      if (length(df) != length(random)) {
        stop("The number of variance components and the length of df must me equal.")
      }
      if (family == "bernoulli") {
        fit0 <- mcemMLE_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, kY, 
                                   kX, kZ, initial, controlEM = ctrl, 
                                   controlTrust = cTrust)
      } 
      if (family == "poisson") {
        fit0 <- mcemMLEPoisson_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, 
                                          kY, kX, kZ, initial, controlEM = ctrl,
                                          controlTrust = cTrust)
      }
      if (family == "negbinom") {
        fit0 <- mcemMLENegBinom_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, 
                                           kY, kX, kZ, initial, controlEM = ctrl,
                                           controlTrust = cTrust)
      }
      if (family == "gamma") {
        fit0 <- mcemMLEGamma_t_fixed_df(sigmaType, df, kKi, kLh, kLhi, 
                                        kY, kX, kZ, initial, controlEM = ctrl, 
                                        controlTrust = cTrust)
      }
    }
  }
  
  # Some cleaning up.
  class(fit0) <- "mcemGLMM"
  if (family == "bernoulli" | family == "poisson") {
    colnames(fit0$mcemEST) <- c(xlabs, zlabs)
  } else {
    colnames(fit0$mcemEST) <- c(xlabs, "alpha", zlabs)
  }
  
  # Remove rows of zeros (stopping before max iterations)
  fit0$mcemEST <- fit0$mcemEST[rowSums(fit0$mcemEST^2) != 0, ]
  
  # Save call
  fit0$call <- call0
  
  # Names for fixed effect's design matrix
  colnames(fit0$x) <- colnames(kX)
  
  # Write warnings or other messages
  if (abs(det(fit0$iMatrix)) < .Machine$double.eps) {
    warning("Information matrix is not invertible.")
  }
  
  if (sum(diag(solve(fit0$iMatrix)) < rep(0, length(diag(fit0$iMatrix))))) {
    warning("Negative standard error estimate. \n This is possible due to Monte Carlo error. Extending the model with mcemGLMMext is recommended. See help(mcemGLMMext) for details.")
  }
  
  return(fit0)
}
