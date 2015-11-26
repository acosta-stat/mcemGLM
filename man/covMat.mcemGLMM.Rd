\name{covMat.mcemGLMM}
\alias{covMat}
\alias{covMat.mcemGLMM}
\title{Anova method for mcemGLMM objects}
\description{
  Extract the fixed effect's covariance matrix from a model fitted with \code{mcemGLMM}.
}
\usage{
covMat.mcemGLMM(object, ...)
}
\arguments{
  \item{object}{a model fitted with the mcemGLMM function.}
  \item{...}{additional arguments.}
}
\value{
  A matrix corresponding to the covariance matrix of the fixed effects.
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}
\keyword{glmm}
