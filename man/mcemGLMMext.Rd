\name{mcemGLMMext}
\alias{mcemGLMMext}

\title{Extending the Iterations of a Model Fitted with mcemGLMM}
\description{
Given a model fitted with the function \code{mcemGLMM} this function will add iterations and update the model estimates for more accurate results.

This is recommended if the initial fitting seems to have a large Monte Carlo error. This function will use the previous maximum likelihood estimate as its initial point and will also start with a Monte Carlo sample size equal to the sample size used in the last iteration of the previous fitting.
}
\usage{
mcemGLMMext(object, it = 20, controlEM)
}
\arguments{
  \item{object}{an model fitted with \code{mcemGLMM}}
  \item{it}{the maximum number of iterations to be performed.}
  \item{controlEM}{a list. New set of options for the EM algorithm. Can be missing}
}
\value{
An updated object of class \code{mcemGLMM}.
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}
\note{
If \code{controlEM} is supplied it is important that the value for \code{MCit} is at least equal to number of Monte Carlo iterations used in the last EM step to fit \code{object} since providing a lower number will increase the Monte Carlo error.
}

\seealso{
\code{\link{mcemGLMM}}
}
\examples{
\donttest{
set.seed(72327)
data(exdata)
fit1 <- mcemGLMM(obs ~ z2 + x, random = ~ 0 + z1, 
                 data = exdata, 
                 family = "bernoulli", vcDist = "normal", 
                 controlEM = list(verb = FALSE, EMit = 5, MCit = 8000), 
                 initial = c(-0.13, -0.15, -0.21, 1.59, 0.002))
                 
# Now we extend the algorithm to do at least another 10 iterations
fit2 <- mcemGLMMext(fit1, it = 10)
}
}
\keyword{ glmm }
\keyword{ mcemGLMM }
