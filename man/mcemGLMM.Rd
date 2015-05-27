\name{mcemGLMM}
\alias{mcemGLMM}
\title{Generalized Linear Mixed Models Estimation}
\description{
  Maximum likelihood estimation for logistic, Poisson, and negative binomial models with random effects using a Monte Carlo EM algorithm.
}
\usage{
mcemGLMM(fixed, random, data, family = c("bernoulli", "poisson", 
    "negbinom"), vcDist = c("normal", "t"), df, corType, 
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
    for Poisson count regression, and "negbinomial" for negative 
    binomial count regression.}
  \item{vcDist}{a string indicating the distribution of the marginal 
    variance components.}
  \item{df}{a vector of degrees of freedom of the random effects 
    when these are t distributed. The length of the vector must be
    equal to the number of variance components in the model.}
  \item{corType}{correlation type. This does nothing at the moment.}
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
    \item{mcemEST}{a matrix with the value of the estimators at the end 
      of each EM step.}
    \item{iMatrix}{Fisher's information matrix.}
    \item{randeff}{a sample from the conditional distribution of the 
      random effects given the data and the maximum likelihood estimators.}
    \item{y}{vector of observations.}
    \item{x}{design matrix for the fixed effects.}
    \item{z}{design matrix for the random effects.}
    \item{EMerror}{relative error at the last iteration. See details.}
    \item{call}{original call.}
  }
}
\details{
  The function \code{mcemGLMM} allows the fitting of generalized linear 
  mixed models when the random effects are normal or a t distributed. 
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
}
\examples{
# Data set for a logistic model with one binary fixed effects and two 
# possible random effects.
data(exdata)

# To fit a model with one random effect
fit.1 <- mcemGLMM(obs ~ x, random = ~ 0 + z1, data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(verb = FALSE, MCit = 10000), 
                initial = c(0.27, -0.13, 0.003))
summary(fit.1)

# To fit a model with crossed random effects
fit.crossed <- mcemGLMM(obs ~ x, random = list(~ 0 + z1, ~ 0 + z2), 
                data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(verb = FALSE, MCit = 10000), 
                initial = c(0.28, -0.15, 0.001, 0.4))
summary(fit.crossed)


# To fit a model with crossed random effects
fit.nested <- mcemGLMM(obs ~ x, random = list(~ 0 + z2, ~ 0 + z2:z1), 
                data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(verb = FALSE, MCit = 10000), 
                initial = c(0.31, -0.15, 0.29, 0.27))
summary(fit.nested)

}
\keyword{glmm}
