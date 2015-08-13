\name{simData}
\alias{simData}
\docType{data}
\title{
Data used for fitting examples 
}
\description{
Example data for logistic, Poisson and negative binomial models.
}
\usage{data("simData")}
\format{
  A data frame with 200 observations on the following 8 variables.
  \describe{
    \item{\code{obs}}{a binary vector. Used as a response for a logistic model.}
    \item{\code{x1}}{a numeric vector.}
    \item{\code{x2}}{a numeric vector.}
    \item{\code{x3}}{a categorical vector with levels \code{blue}, \code{red}, and \code{yellow}.}
    \item{\code{z1}}{a categorical vector with levels \code{D1}, \code{D2}, \code{D3}, \code{D4}, and \code{D5}.}
    \item{\code{z2}}{a categorical vector with levels \code{1}, \code{2}, \code{3}, \code{4}, and \code{5}.}
    \item{\code{z3}}{a categorical vectir with levels \code{A}, \code{B}, and \code{D}.}
    \item{\code{count}}{a numeric vector. Used as a response for a Poisson model.}
  }
}
\details{
The levels of \code{z2} can be nested within $\code{z1}$. The observations were generated with the code shown in the examples section.
}
\examples{
set.seed(47819)
x1 <- rnorm(200, 10, 1)
x2 <- rnorm(200, 5, 1)
x3 <- sample(c("red", "blue", "yellow"), size = 200, replace = TRUE)
z1 <- factor(rep(c("D1", "D2", "D3", "D4", "D5"), each = 40))
z2 <- factor(rep(rep(1:4, each = 5), 10))
z3 <- factor(c(rep("A", 100), rep("B", 60), rep("D", 40)))
kX <- model.matrix(~x1 + x2 + x3)
kZ <- cbind(model.matrix(~ 0+z1), model.matrix(~ 0+z1:z2), model.matrix(~ 0+z3))
kBeta <- c(5, -4, 5, 0, 8)
kU <- 3 * rt(28, 5)
linf0 <- kX %*% kBeta + kZ %*% kU
prob0 <- exp(linf0)/(1+exp(linf0))
obs <- as.numeric(runif(100) < prob0)
simData <- data.frame(obs, x1, x2, x3, z1, z2, z3)
}
\keyword{datasets}
