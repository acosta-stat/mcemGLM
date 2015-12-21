#' Generalized Linear Mixed Model Estimation via Monte Carlo EM
#' 
#' \code{mcemGLM} performs maximum likelihood estimation for logistic, Poisson,
#' and negative binomial regression when random effects are present. The
#' package uses an MCEM algorithm to estimate the model's fixed paramters and
#' variance components with their respective standard errors.
#' 
#' A Wald test based \code{anova} is available to test significance of
#' multi-leveled variables and for multiple contrast testing.
#' 
#' \tabular{ll}{ Package: \tab mcemGLM\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2015-05-22\cr License: \tab GPL (>= 2)\cr }
#' 
#' @name mcemGLM-package
#' @aliases mcemGLM-package mcemGLM
#' @docType package
#' @author Felipe Acosta Archila
#' 
#' Maintainer: Felipe Acosta Archila <acosta@@umn.edu>
#' @keywords glmm
#' @examples
#' 
#' \donttest{
#' set.seed(123)
#' x <- rnorm(30, 10, 1)
#' z <- factor(rep(1:6, each = 5))
#' obs <- sample(0:1, 30, TRUE)
#' fit <- mcemGLMM(obs ~ x, random = ~ 0 + z, family = "bernoulli",
#' vcDist = "normal", controlEM = list(EMit = 15, MCit = 10000), 
#' initial = c(3.30, -0.35, 0.005))
#' summary(fit)
#' anova(fit)
#' }
#' 
#' @importFrom trust trust
#' @importFrom stats rnorm as.formula binomial get_all_vars glm model.matrix pchisq pnorm poisson terms ts.plot lm
#' @importFrom utils tail
#' @importFrom Rcpp evalCpp
#' @useDynLib mcemGLM
"_PACKAGE"
#> [1] "_PACKAGE"