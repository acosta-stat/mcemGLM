\name{contrasts.mcemGLMM}
\alias{contrasts.mcemGLMM}
\title{Contrast estimation for mcemGLMM objects}
\description{
  Contrast testing for a model fitted with \code{mcemGLMM}.
}
\usage{
contrasts.mcemGLMM(object, ctr.mat)
}
\arguments{
  \item{object}{a model fitted with the mcemGLMM function.}
  \item{ctr.mat}{contrast matrix. Each row corresponds to a 
    linear combination of the fixed effect coefficients.}
}
\value{
  A matrix with each row corresponding to a contrast and with the following columns:
  \describe{
    \item{Estimate}{contrasts's point estimator.}
    \item{Std. Err.}{standard error for each fitted contrast.}
    \item{Wald}{Wald statistic for each fitted contrast.}
    \item{Adj. p-value}{p-value adjusted for multiple comparison (Bonferroni.)}
  }
}
\examples{
\donttest{
set.seed(72327)
data(exdata)
fit <- mcemGLMM(obs ~ z2 + x, random = ~ 0 + z1, 
                data = exdata, 
                family = "bernoulli", vcDist = "normal", 
                controlEM = list(verb = FALSE, EMit = 15, MCit = 10000), 
                initial = c(-0.13, -0.15, -0.21, 1.59, 0.002))
mat <- rbind("1 - 2" = c(0, -1, 0, 0), "1 - 3" = c(0, 0, -1, 0), "2 - 3" = c(0, 1, -1, 0))
contrasts.mcemGLMM(fit, mat)
}
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}
\keyword{glmm}
