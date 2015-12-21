#' Fixed effect covariance matrix for an mcemGLMM object
#' 
#' Extract the fixed effect's covariance matrix from a model fitted with
#' \code{mcemGLMM}.
#' 
#' @aliases vcov vcov.mcemGLMM
#' @param object a model fitted with the mcemGLMM function.
#' @param ... additional arguments.
#' @return A matrix corresponding to the covariance matrix of the fixed
#' effects.
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @keywords glmm
#' 
#' @export
#' 
vcov.mcemGLMM <- function(object, ...) {
  kP <- ncol(object$x)
  fMat <- solve(object$iMatrix)
  fMat <- fMat[1:kP, 1:kP]
  colnames(fMat) <- colnames(object$x)
  rownames(fMat) <- colnames(object$x)
  return(fMat)
}