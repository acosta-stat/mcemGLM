#' Predict method for mcemGLMM objects
#' 
#' This functions returns predicted link function of observations for a model
#' fitted with \code{mcemGLMM}.
#' 
#' @aliases predict predict.mcemGLMM
#' @param object a model fitted with the mcemGLMM function.
#' @param newdata optional data frame with new data. The variable names must
#' match the original variables. If missing, the function will return predicted
#' values at each observation.
#' @param type character string. The type of predictions to be returned. Either
#' "link" or "response" predictions are available.
#' @param se.fit logical. If true, standard erros will be returned.
#' @param ... additional arguments.
#' @return A vector with the predictions from the observed data or by using the
#' supplied new data.
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @keywords glmm, prediction
#' 
#' @export
#' 
predict.mcemGLMM <- function(object, newdata, type = c("link", "response"), se.fit = FALSE, ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  if (missing(newdata)) {
    lin0 <- as.vector(object$x %*% coef0)
    if (type[1] == "link") {
      if (se.fit == FALSE)
        return(lin0)
      else {
        kN <- length(object$y)
        tmp <- matrix(0, kN, 2)
        colnames(tmp) <- c("Estimate", "SE")
        cmat <- vcov(object)
        for (i in 1:kN) {
          tmp[i, 2] <- sqrt(object$x[i, ] %*% cmat %*% (object$x[i, ]))
        }
        return(tmp)
      }
    } 
    if (type == "response") {
      if (object$call$family == "bernoulli") {
        if (se.fit == FALSE)
          return(exp(lin0) / (1 + exp(lin0)))
        else {
          kN <- length(object$y)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0) / (1 + exp(lin0))
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- vcov(object)
          for (i in 1:kN) {
            tmp[i, 2] <- sqrt(exp(2 * lin0[i])/(1 + exp(lin0[i]))^4 * object$x[i, ] %*% cmat %*% object$x[i, ])
          }
          return(tmp)
        }
      }
      if (object$call$family %in% c("poisson", "negbinom", "gamma")) {
        if (se.fit == FALSE)
          return(exp(lin0))
        else {
          kN <- length(object$y)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0)
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- vcov(object)
          for (i in 1:kN) {
            tmp[i, 2] <- sqrt(exp(2 * lin0[i]) * object$x[i, ] %*% cmat %*% object$x[i, ])
          }
          return(tmp)
        }
      }
    }
  } else {
    # newdata
    tmp.x <- model.matrix(as.formula(object$call$fixed)[-2], data = newdata)
    if (!prod(colnames(tmp.x) == colnames(object$x))) {
      stop("Incorrect new data.")
    }
    lin0 <- as.vector(tmp.x %*% coef0)
    if (type[1] == "link") {
      if (se.fit == FALSE)
        return(lin0)
      else {
        kN <- nrow(tmp.x)
        tmp <- matrix(0, kN, 2)
        tmp[, 1] <- lin0
        colnames(tmp) <- c("Estimate", "SE")
        cmat <- vcov(object)
        for (i in 1:kN) {
          tmp[i, 2] <- sqrt(tmp.x[i, ] %*% cmat %*% (tmp.x[i, ]))
        }
        return(tmp)
      }
    } 
    if (type == "response") {
      if (object$call$family == "bernoulli") {
        if (se.fit == FALSE)
          return(exp(lin0) / (1 + exp(lin0)))
        else {
          kN <- nrow(tmp.x)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0) / (1 + exp(lin0))
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- vcov(object)
          for (i in 1:kN) {
            tmp[i, 2] <- sqrt(exp(2 * lin0[i])/(1 + exp(lin0[i]))^4 * tmp.x[i, ] %*% cmat %*% tmp.x[i, ])
          }
          return(tmp)
        }
      }
      if (object$call$family %in% c("poisson", "negbinom", "gamma")) {
        if (se.fit == FALSE)
          return(exp(lin0))
        else {
          kN <- nrow(tmp.x)
          tmp <- matrix(0, kN, 2)
          tmp[, 1] <- exp(lin0)
          colnames(tmp) <- c("Estimate", "SE")
          cmat <- vcov(object)
          for (i in 1:kN) {
            tmp[i, 2] <- sqrt(exp(2 * lin0[i]) * tmp.x[i, ] %*% cmat %*% tmp.x[i, ])
          }
          return(tmp)
        }
      }
    }
  }
}