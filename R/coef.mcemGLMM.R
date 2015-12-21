#' Fixed effect coefficients extraction for mcemGLMM objects
#' 
#' Extract the fixed effect coeffcients from a model fitted with
#' \code{mcemGLMM}.
#' 
#' 
#' @aliases coef coef.mcemGLMM
#' @param object a model fitted with the mcemGLMM function.
#' @param ... additional arguments.
#' @return A vector with the fixed effect coefficients.
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @keywords glmm
#' 
#' @export
#' 
coef.mcemGLMM <- function(object,...) {
  coef0 <- tail(object$mcemEST, n = 1)[1:ncol(object$x)]
  names(coef0) <- colnames(object$mcemEST)[1:ncol(object$x)]
  return(coef0)
}
