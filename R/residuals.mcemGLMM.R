residuals.mcemGLMM <- function(object, ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  u0 <- colMeans(object$randeff)
  lin0 <- object$x %*% coef0 + object$z %*% u0
  prob0 <- exp(lin0) / (1 + exp(lin0))
  res0 <- object$y - prob0
  return(res0)
}