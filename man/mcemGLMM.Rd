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
    linear model fitted without ranom effects and the initial values for 
    the variance components are set to 5.}
}
\value{
  A list of class "mcemGLMM".
}
\details{
This will talk about the random effects options, the em options list, and the trust list.
}
\examples{
# Same as before.
set.seed(123)
x <- rnorm(30, 10, 1)
z <- factor(rep(1:6, each = 5))
obs <- sample(0:1, 30, TRUE)
fit <- mcemGLMM(obs ~ x, random = ~ 0 + z, family = "bernoulli",
vcDist = "normal", controlEM = list(verb = FALSE, MCit = 10000), 
initial = c(3.30, -0.35, 0.005))
summary(fit)
anova(fit)
}
\keyword{glmm}
