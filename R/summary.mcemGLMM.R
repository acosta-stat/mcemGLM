summary.mcemGLMM <- function(obj, covm = FALSE) {
  # Fixed effects
  coef0 <- tail(obj$mcemEST, n = 1)[1:ncol(obj$x)]
  names(coef0) <- colnames(obj$mcemEST)[1:ncol(obj$x)]
  
  ran.eff0 <- colMeans(obj$randeff)
  
  # Variance estimates
  var.est0 <-tail(obj$mcemEST, n = 1)[-(1:ncol(obj$x))]
  names(var.est0) <- colnames(obj$mcemEST)[-(1:ncol(obj$x))]
  
  # Covariance matrix and standard errors
  cmat0 <- solve(obj$iMatrix)
  std.err0 <- sqrt(diag(cmat0))
  std.err1 <- std.err0[-(1:ncol(obj$x))]
  std.err0 <- std.err0[1:ncol(obj$x)]
  
  # z values
  zval0 <- coef0/std.err0[1:ncol(obj$x)]
  zval1 <- var.est0/std.err1
  
  # p values
  pval0 <- 2 * pnorm(-abs(zval0))
  pval1 <- pnorm(-abs(zval1))
  
  resultsFixed <- matrix(0, length(coef0), 4)
  resultsFixed[, 1] <- coef0
  resultsFixed[, 2] <- std.err0[1:ncol(obj$x)]
  resultsFixed[, 3] <- zval0
  resultsFixed[, 4] <- pval0
  rownames(resultsFixed) <- names(coef0)
  colnames(resultsFixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  resultsVar <- matrix(0, length(var.est0), 4)
  resultsVar[, 1] <- var.est0
  resultsVar[, 2] <- std.err1
  resultsVar[, 3] <- zval1
  resultsVar[, 4] <- pval1
  rownames(resultsVar) <- names(var.est0)
  colnames(resultsVar)  <- c("Estimate", "Std. Error", "z value", "Pr(>z)")
  
  cat("   Two sided Wald tests for fixed effects coefficients:\n\n")
  print(resultsFixed)
  
  cat("\n\n   One sided Wald tests for variance components:\n\n")
  print(resultsVar)
  
  tbr <- list(coefficients = list(fixed = coef0, random = ran.eff0), var.est = var.est0, std.err = c(std.err0, std.err1), z.val = c(zval0, zval1))
  invisible(tbr)
}