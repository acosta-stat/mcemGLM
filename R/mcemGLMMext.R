#' Extending the iterations of a model fitted with mcemGLMM
#' 
#' Given a model fitted with the function \code{mcemGLMM} this function will
#' add iterations and update the model estimates for more accurate results.
#' 
#' This is recommended if the initial fitting seems to have a large Monte Carlo
#' error. This function will use the previous maximum likelihood estimate as
#' its initial point and will also start with a Monte Carlo sample size equal
#' to the sample size used in the last iteration of the previous fitting.
#' 
#' 
#' @param object an model fitted with \code{mcemGLMM}
#' @param it the maximum number of iterations to be performed.
#' @param controlEM a list. New set of options for the EM algorithm. Can be
#' missing
#' @return An updated object of class \code{mcemGLMM}.
#' @note If \code{controlEM} is supplied it is important that the value for
#' \code{MCit} is at least equal to number of Monte Carlo iterations used in
#' the last EM step to fit \code{object} since providing a lower number will
#' increase the Monte Carlo error.
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @seealso \code{\link{mcemGLMM}}
#' @keywords glmm mcemGLMM
#' @examples
#' 
#' \donttest{
#' set.seed(72327)
#' data(exdata)
#' fit1 <- mcemGLMM(obs ~ z2 + x, random = ~ 0 + z1, 
#'                  data = exdata, 
#'                  family = "bernoulli", vcDist = "normal", 
#'                  controlEM = list(verb = FALSE, EMit = 5, MCit = 8000), 
#'                  initial = c(-0.13, -0.15, -0.21, 1.59, 0.002))
#'                  
#' # Now we extend the algorithm to do at least another 10 iterations
#' fit2 <- mcemGLMMext(fit1, it = 10)
#' }
#' 
#' @export
#' 
mcemGLMMext <- function(object, it = 20, controlEM) {
  if (class(object) != "mcemGLMM") {
    stop("Wrong object type.")
  }
  
  if (missing(controlEM)) {
    ctrl <- list(     EMit = it,
                      MCit = eval(object$MCit),
                       MCf = ifelse(      is.null(object$call$controlEM$MCf),  1.25,       object$call$controlEM$MCf),
                      verb = ifelse(     is.null(object$call$controlEM$verb), FALSE,      object$call$controlEM$verb),
                      MCsd = object$MCsd, 
                   EMdelta = ifelse(  is.null(object$call$controlEM$EMdelta),  0.025,   object$call$controlEM$EMdelta),
                 EMepsilon = ifelse(is.null(object$call$controlEM$EMepsilon), 0.001, object$call$controlEM$EMepsilon))
    
    fit0 <- mcemGLMM(       fixed = eval(object$call$fixed),
                           random = eval(object$call$random), 
                             data = eval(object$call$data), 
                           family = eval(object$call$family), 
                           vcDist = eval(object$call$vcDist), 
                               df = eval(object$call$df), 
                        controlEM = ctrl, 
                     controlTrust = eval(object$call$controlTrust), 
                          initial = tail(eval(object$mcemEST), 1))
  } else {
    fit0 <- mcemGLMM(       fixed = eval(object$call$fixed),
                           random = eval(object$call$random),
                             data = eval(object$call$data),
                           family = eval(object$call$family),
                           vcDist = eval(object$call$vcDist),
                               df = eval(object$call$df),
                        controlEM = controlEM,
                     controlTrust = eval(object$call$controlTrust),
                          initial = tail(eval(object$mcemEST), 1))
  }
  fit0$call <- object$call
  return(fit0)
}