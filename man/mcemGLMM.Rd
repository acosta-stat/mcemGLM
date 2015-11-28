\name{mcemGLMM}
\alias{mcemGLMM}
\title{Generalized Linear Mixed Models Estimation}
\description{
  Maximum likelihood estimation for logistic, Poisson, and negative binomial models with random effects using a Monte Carlo EM algorithm.
}
\usage{
mcemGLMM(fixed, random, data, family = c("bernoulli", "poisson", 
    "negbinom", "gamma"), vcDist = c("normal", "t"), df, 
    controlEM = list(), controlTrust = list(), initial)
}
\arguments{
  \item{fixed}{the fixed effects model. This is specified 
    by a \code{formula} object.}
  \item{random}{the random effects models. This is specified
    by a \code{formula} object or a list of \code{formula} objects.
    See details below.}
  \item{data}{an optional data frame containing the variables in 
    the model. If missing the variables are taken from the current 
    environment.}
  \item{family}{a string indicating the type of model to be fitted. 
    The options are "bernoulli" for logistic regression, "Poisson"
    for Poisson count regression, and "negbinom" for negative 
    binomial count regression.}
  \item{vcDist}{a string indicating the distribution of the marginal 
    variance components. The options are "normal" and "t" for normal 
    and t distributed random effects respectively.}
  \item{df}{a vector of degrees of freedom of the random effects 
    when these are t distributed. The length of the vector must be
    equal to the number of variance components in the model.}
  \item{controlEM}{a list of options for the algorithm. See Details below.}
  \item{controlTrust}{a list of options to be passed to the \code{trust}
    optimizer. See details below.}
  \item{initial}{optional initial values for the parameters. If missing 
    the initial values for the fixed effects are taken from a generalized 
    linear model fitted without random effects and the initial values for 
    the variance components are set to 5.}
}
\value{
  A list of class "mcemGLMM" with the following items:
  \describe{
    \item{mcemEST}{a matrix with the value of the maximum likelihood 
    estimators at the end of each EM step.}
    \item{iMatrix}{Fisher's information matrix.}
    \item{QfunVal}{ value (up to a constant) of the Q function. 
      Used to perform likelihood ratio tests.}
    \item{QfunMCMC}{Q function MCMC sample.}
    \item{randeff}{a sample from the conditional distribution of the 
      random effects given the data and the maximum likelihood estimators.}
    \item{y}{vector of observations.}
    \item{x}{design matrix for the fixed effects.}
    \item{z}{design matrix for the random effects.}
    \item{EMerror}{relative error at the last iteration. See details.}
    \item{MCsd}{last value for MCMC step size.}
    \item{call}{original call.}
  }
}
\details{
  The function \code{mcemGLMM} allows the fitting of generalized linear 
  mixed models when the random effects are normal or a t distributed.
  The supported models are the logistic, Poisson and negative binomial. 
  The degrees of freedom for the t case must be supplied by the user 
  with a vector in the \code{df} argument. The length of the vector 
  must be equal to the number of variance components. For normal random 
  effects the argument \code{df} does not need to be included.
  
  To fit a model with one random effect a formula must be supplied in 
  the \code{random} argument. Note that it is necessary that the 
  variable is a factor and to specify that there is no intercept for 
  the random part of the model. To use more than one random effects a 
  list of formulas must be supplied. Each member must be formula with 
  in which the variables involved must be factors and it also it is 
  necessary to specify that there is no intercept. To fit crossed 
  random effects each variable must no appear in its own formula. To 
  fit nested random effects a formula with the highest level variable 
  must be specified and each subsequent variable must be specified with 
  an interaction of the variables above it. See examples below.
  
  
  A note on the negative binomial overdispersion:
  
  The variance for the negative binomial model is set equal to
  \eqn{(1 + \mu/\alpha)} so we have that there is no overdispersion
  as \eqn{\alpha} goes to infinity.
  
  The variance for the gamma distributionis set to \eqn{\mu^2/\alpha}. The value
  of \eqn{\alpha = 1} corresponds to an exponential regression.
  
  Stopping rules and convergence criteria:
  
  The algorithm runs for a maximum of \code{EMit} iterations or until
  the criteria 
  \deqn{\max_i \left\{ \frac{ |\theta_i^{(t)} - \theta_i^{(t-1)}|}{|\theta_i^{(t)}| + \delta}\right\} < \epsilon}{\max_i|\theta_i^(t)-\theta_i^(t-1)|/(|\theta_i^(t)| + \delta) < \epsilon}
  is satisfied three times in a row for pre-defined values of 
  \eqn{\epsilon} and \eqn{\delta}. Once this criterion has been achieved 
  two times we increase the Monte Carlo sample size more rapidly to 
  have a better estimation of the model's information matrix. For 
  a detailed discussion on convergence diagnostics see Neath R.C. (2012). 
  After fitting a model it is recommended to plot the EM estimates 
  at each step to assess convergence.
 
  Control options for EM:
  \describe{
    \item{EMit}{maximum number of EM iterations.}
    \item{MCit}{initial number of Monte Carlo iterations for the MCMC 
    step.}
    \item{MCf}{factor in which the MC iterations increase in each EM 
    iteration.}
    \item{verb}{logical. If TRUE at each EM iteration the function will 
    print convergence information and a trace plot for on of the random 
    effects. This is useful to assess the performance of the algorithm 
    but it can impact the actual running time.}
    \item{MCsd}{initial standard deviation for the proposal density of 
    the MCMC step. If zero (default) an auto-tuning step will be 
    performed.}
    \item{EMdelta}{constant for the EM error assessment.}
    \item{EMepsilon}{constant for the EM error assessment.}
  }
  
  Control options for trust, see \code{help(trust)} for more details:
  \describe{
    \item{rinit}{starting trust region radius. Default value set to 20.}
    \item{rmax}{maximum allowed trust region radius. Default value set to 200.}
    \item{iterlim}{maximum number of iterations. Default value set to 100.}
  }
  
}
\examples{
\donttest{
# Data set for a logistic model with one binary fixed effects and two 
# possible random effects.
# Initial values and MC iterations are given to speed up the examples 
# but these are not necessary in general.
set.seed(0123210)
data(exdata)

# To fit a model with one random effect
fit.1 <- mcemGLMM(obs ~ x, random = ~ 0 + z1, data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(MCit = 10000), 
                initial = c(0.27, -0.13, 0.003))
summary(fit.1)

# We can assess convergence by looking at a trace plot of the EM estimates
# and the loglikelihood values
ts.plot(fit.1$mcemEST)
ts.plot(fit.1$QfunVal)

# To fit a model with crossed random effects
fit.crossed <- mcemGLMM(obs ~ x, random = list(~ 0 + z1, ~ 0 + z2), 
                data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(EMit = 10, MCit = 10000), 
                initial = c(0.28, -0.15, 0.001, 0.4))
summary(fit.crossed)


# To fit a model with crossed random effects
fit.nested <- mcemGLMM(obs ~ x, random = list(~ 0 + z2, ~ 0 + z2:z1), 
                data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(EMit = 10, MCit = 10000), 
                initial = c(0.31, -0.15, 0.29, 0.27))
summary(fit.nested)

# Fit a Poisson model
fit.pois <- mcemGLMM(obs2 ~ x, random = ~ 0 + z1, data = exdata, 
                family = "poisson", vcDist = "normal", 
                controlEM = list(EMit = 10, MCit = 10000), 
                initial = c(1.95, 0.03, 0.003))
summary(fit.pois)
}
}
\references{
Neath, R. C. (2012) On Convergence Properties of the 
Monte Carlo EM Algorithm In Advances in Modern Statistical 
Theory and Applications: A Festschrift in Honor of Morris 
L. Eaton. \emph{Institute of Mathematical Statistics} 43--62  
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}

\keyword{glmm}
