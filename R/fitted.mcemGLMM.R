#' Extract fitted values for an mcemGLMM object
#' 
#' This function returns the fitted values of a model fitted with \code{mcemGLMM}.
#' 
#' @aliases fitted fitted.values
#' @param object a model fitted with the mcemGLMM function.
#' @param ... additional arguments.
#' @return A vector of fitted values.
#' 
#' @export
#' 
fitted.mcemGLMM <- function(object, ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  u0 <- colMeans(object$randeff)
  lin0 <- object$x %*% coef0 + object$z %*% u0
  if (object$call$family == "bernoulli") {
    mu0 <- exp(lin0) / (1 + exp(lin0))
    res0 <- sqrt(mu0 * (1 - mu0))
  }
  
  if (object$call$family == "poisson") {
    mu0 <- exp(lin0)
    res0 <- sqrt(mu0)
  }
  
  if (object$call$family == "negbinom") {
    mu0 <- exp(lin0)
    a0 <- tail(object$mcemEST, 1)[kP + 1]
    res0 <- sqrt(mu0 * (1 + 1/a0))
  }
  
  if (object$call$family == "gamma") {
    mu0 <- exp(lin0)
    a0 <- tail(object$mcemEST, 1)[kP + 1]
    res0 <- mu0 / sqrt(a0)
  }
  return(as.vector(mu0))
}