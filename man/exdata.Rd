\name{exdata}
\alias{exdata}
\docType{data}
\title{
Example data
}
\description{
Simulated binary response dataset to use in examples.
}
\usage{data("exdata")}
\format{
  A data frame with 60 observations on the following 4 variables.
  \describe{
    \item{\code{obs}}{a numeric binary vector}
    \item{\code{obs2}}{a numeric count vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{z1}}{a factor with levels \code{1}, and \code{2}}
    \item{\code{z2}}{a factor with levels \code{1}, \code{2}, and \code{3}}
  }
  The observations were generated independently with the code shown 
  in the examples section.
}
\examples{
set.seed(123)
obs <- c(sample(0:1, 30, TRUE, prob = c(0.5, 0.5)), sample(0:1, 30, TRUE, prob = c(0.3, 0.7)))
obs2 <- c(rpois(30, 5), rpois(30, 10))
}
\keyword{datasets}
