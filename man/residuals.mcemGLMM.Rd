\name{residuals.mcemGLMM}
\alias{residuals}
\alias{residuals.mcemGLMM}
\title{Residual extraction method for mcemGLMM objects}
\description{
  This functions returns the residuals of a model fitted with \code{mcemGLMM}.
}
\usage{
\method{residuals}{mcemGLMM}(object, type = c("deviance", "pearson"), ...)
}
\arguments{
  \item{object}{a model fitted with the mcemGLMM function.}
  \item{type}{character string. The type of residuals to be returned.}
  \item{...}{additional arguments.}
}
\value{
  A vector with the residuals of the model.
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}
\concept{residuals}
