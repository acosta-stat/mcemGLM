predict.mcemGLMM <- function(object, newdata, type = c("link", "response"), ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  if (missing(newdata)) {
    lin0 <- as.vector(object$x %*% coef0)
    if (type[1] == "link") {
      return(lin0)
    } 
    if (type == "response") {
      if (object$call$family == "bernoulli") {
        return(exp(lin0) / (1 + exp(lin0)))
      }
      if (object$call$family == "poisson") {
        return(exp(lin0))
      }
      if (object$call$family == "negbinom") {
        return(exp(lin0))
      }
    }
  } else {
    newdata
    tmp.x <- model.matrix(as.formula(object$call$fixed)[-2], data = newdata)
    if (!prod(colnames(tmp.x) == colnames(object$x))) {
      stop("Incorrect new data.")
    }
    lin0 <- as.vector(tmp.x %*% coef0)
    if (type[1] == "link") {
      return(lin0)
    } 
    if (type == "response") {
      if (object$call$family == "bernoulli") {
        return(exp(lin0) / (1 + exp(lin0)))
      }
      if (object$call$family == "poisson") {
        return(exp(lin0))
      }
      if (object$call$family == "negbinom") {
        return(exp(lin0))
      }
    }
  }
}