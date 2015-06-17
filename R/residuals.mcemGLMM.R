residuals.mcemGLMM <- function(object, type = c("deviance", "pearson"), ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  u0 <- colMeans(object$randeff)
  lin0 <- object$x %*% coef0 + object$z %*% u0
  if (object$call$family == "bernoulli") {
    prob0 <- exp(lin0) / (1 + exp(lin0))
    if (type == "pearson") {
      res0 <- (object$y - prob0) / sqrt(prob0 * (1 - prob0))
    }
    if (type == "deviance") {
      res0 <- ifelse(object$y, 1, -1) * sqrt(-2 * (object$y * log(prob0) + (1 - object$y) * log(1 - prob0)))
    }
  }
  
  if (object$call$family == "poisson") {
    mu0 <- exp(lin0)
    if (type == "pearson") {
      res0 <- (object$y - mu0) / sqrt(mu0)
    }
    if (type == "deviance") {
      print(mu0)
      res0 <- sign(object$y - mu0) * sqrt(2 * object$y * log(object$y / mu0) - 2 * (object$y - mu0))
    }
    
  }
  
  return(res0)
}