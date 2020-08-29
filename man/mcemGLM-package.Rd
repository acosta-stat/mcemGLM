\name{mcemGLM-package}
\alias{mcemGLM-package}
\alias{mcemGLM}
\docType{package}
\title{
Generalized Linear Mixed Model Estimation via Monte Carlo EM
}
\description{
\code{mcemGLM} performs maximum likelihood estimation for logistic, 
Poisson, and negative binomial regression when random effects are 
present. The package uses an MCEM algorithm to estimate the model's 
fixed parameters and variance components with their respective 
standard errors.

A Wald test based \code{anova} is available to test significance of 
multi-leveled variables and for multiple contrast testing.
}
\details{
\tabular{ll}{
Package: \tab mcemGLM\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2015-05-22\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Felipe Acosta Archila

Maintainer: Felipe Acosta Archila <acosta@umn.edu>
}
\keyword{ glmm }
\examples{
\donttest{
set.seed(123)
x <- rnorm(30, 10, 1)
z <- factor(rep(1:6, each = 5))
obs <- sample(0:1, 30, TRUE)
fit <- mcemGLMM(obs ~ x, random = ~ 0 + z, family = "bernoulli",
vcDist = "normal", controlEM = list(EMit = 15, MCit = 10000), 
initial = c(3.30, -0.35, 0.005))
summary(fit)
anova(fit)
}
}
