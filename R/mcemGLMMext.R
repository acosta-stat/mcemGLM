mcemGLMMext <- function(object, minIt = 25, controlEM) {
  if (class(object) != "mcemGLMM") {
    stop("Wrong object type.")
  }
  
  if (missing(controlEM)) {
    ctrl <- list(     EMit = max(nrow(object$mcemEST), minIt),
                      MCit = nrow(object$randeff),
                       MCf = ifelse(      is.null(object$call$controlEM$MCf),  1.04,       object$call$controlEM$MCf),
                      verb = ifelse(     is.null(object$call$controlEM$verb), FALSE,      object$call$controlEM$verb),
                      MCsd = ifelse(     is.null(object$call$controlEM$MCsd),     0,      object$call$controlEM$MCsd), 
                   EMdelta = ifelse(  is.null(object$call$controlEM$EMdelta),  0.01,   object$call$controlEM$EMdelta),
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
                           family = object$call$family,
                           vcDist = object$call$vcDist,
                               df = object$call$df,
                        controlEM = controlEM,
                     controlTrust = object$call$controlTrust,
                          initial = tail(object$mcemEST, 1))
  }
  fit0$call <- object$call
  return(fit0)
}