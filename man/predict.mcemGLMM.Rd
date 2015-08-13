\name{predict.mcemGLMM}
\alias{predict}
\alias{predict.mcemGLMM}
\title{Predict method for mcemGLMM objects}
\description{
  This functions returns predicted link function of observations for a model fitted with \code{mcemGLMM}.
}
\usage{
  \method{predict}{mcemGLMM}(object, newdata, type = c("link", "response"), ...)
}
\arguments{
  \item{object}{a model fitted with the mcemGLMM function.}
  \item{newdata}{optional data frame with new data. The variable names must match the original variables. If missing, the function will return predicted values at each observation.}
  \item{type}{character string. The type of predictions to be returned. Either "link" or "response" predictions are available.}
  \item{...}{additional arguments.}
}
\value{
  A vector with the predictions from the observed data or by using the supplied new data.
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}
\keyword{glmm, prediction}
