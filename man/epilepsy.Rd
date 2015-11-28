\name{epilepsy}
\alias{epilepsy}
\docType{data}
\title{
Epilesy Data
}
\description{
Data from an experiment with $i=1,\dots,59$ epilepsy patients. Each of the 
patients was assigned to a control group or a treatment group. The experiment
recorded the number of seizures experienced by each patient over four 
two-week periods. The experiment also recorded a baseline count of the 
number of seizures the patients had experienced during the previous eight weeks.
}
\usage{data("epilepsy")}
\format{
  A data frame with 236 observations on the following 6 variables.
  \describe{
    \item{\code{count}}{Number of seizures experienced.}
    \item{\code{id}}{Patient ID.}
    \item{\code{group}}{Treatment or control groups.}
    \item{\code{age}}{The logarithm of the patient's age.}
    \item{\code{base}}{The logarithm of baseline/4.}
    \item{\code{visit}}{The corresponding time period.}
  }
}

\references{
Thall, P. F. and Vail, S. C. (1990) Some covariance models for
longitudinal count data with overdispersion In \emph{Biometrics} 46, 657--671
}
\keyword{datasets}
