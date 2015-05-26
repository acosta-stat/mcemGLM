\name{summary.mcemGLMM}
\alias{summary}
\alias{summary.mcemGLMM}
\title{Summary method for mcemGLMM objects}
\description{
  Maximum likelihood estimation for logistic, Poisson, and negative binomial models with random effects using a Monte Carlo EM algorithm.
}
\usage{
\method{summary}{mcemGLMM}(object, ...)
}
\arguments{
  \item{object}{a model fitted with the mcemGLMM function.}
  \item{...}{aditional arguments.}
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
}
\keyword{glmm}
