#' Random effect prediction for mcemGLMM objects
#' 
#' Extract the random effect predictions from a model fitted with
#' \code{mcemGLMM}.
#' 
#' @aliases ranef ranef.mcemGLMM
#' @param object a model fitted with the mcemGLMM function.
#' @param ... additional arguments.
#' @return A vector with the random effect predictions.
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @keywords glmm
#' 
#' @export
#' 
ranef.mcemGLMM <- function(object, ...) {
  return(colMeans(object$randeff))
}