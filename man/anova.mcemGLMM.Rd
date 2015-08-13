\name{anova.mcemGLMM}
\alias{anova}
\alias{anova.mcemGLMM}
\title{Anova method for mcemGLMM objects}
\description{
  ANOVA table based on Wald tests for a model fitted with \code{mcemGLMM}.
}
\usage{
\method{anova}{mcemGLMM}(object, opt, ...)
}
\arguments{
  \item{object}{a model fitted with the mcemGLMM function.}
  \item{opt}{optional model fitted with the mcemGLMM function. If supplied a likelihood ratio test will be performed.}
  \item{...}{additional arguments.}
}
\value{
  A matrix with the rows corresponding to a test for the different 
  terms of the model and the following columns:
  \describe{
    \item{Df}{degrees of freedom for the term.}
    \item{Wald}{Wald's chi squared statistic.}
    \item{Pr}{p value for the test statistic.}
  }
  
  In the case a likelihood ratio test is performed a list with the 
  following items will be returned:
  \describe{
    \item{Chi.Sq}{test statistc value.}
    \item{Df}{degrees of freedom for the test.}
    \item{p.value}{p value corresponding to the test statistic.}
  }
}
\author{
Felipe Acosta Archila <acosta@umn.edu>
}
\keyword{glmm}
