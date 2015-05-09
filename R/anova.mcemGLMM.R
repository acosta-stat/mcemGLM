anova.mcemGLMM <- function(obj) {
  # Fixed effects
  coef0 <- tail(obj$mcemEST, n = 1)[1:ncol(obj$x)]
  names(coef0) <- colnames(obj$mcemEST)[1:ncol(obj$x)]
    
  # Covariance matrix and standard errors
  cmat0 <- solve(obj$iMatrix)
  # std.err0 <- sqrt(diag(cmat0))
  
  # z2 values
  
  # Coefficients positions
  pred0 <- attr(obj$x, "assign")
  
  imat0 <- solve(obj$iMatrix)
  
  names0 <- obj$tnames
  # Drop intercept
  if(pred0[1] == 0) {
    coef0 <- coef0[-1]
    pred0 <- pred0[-1]
    imat0 <- imat0[-1, -1]
    names0 <- names0[-1]
  }
  
  # Number of predictors 
  npred <- unique(pred0)
  wald0 <- rep(0, length(npred))
  for (i in npred) {
    ind0 <- which(pred0 == i)
    wald0[i] <- t(coef0[ind0]) %*% solve(imat0[ind0, ind0]) %*% coef0[ind0]
  }
  
  df0 <- table(pred0)
  pval0 <- pchisq(wald0, df0, lower.tail = FALSE)
  tbr <- cbind(df0, wald0, pval0)
  cat("   Wald's Chi-squared Type II Anova table\n\n")
  colnames(tbr) <- c("Df", "Wald Stat.", "Pr(>W)")
  rownames(tbr) <- names0
  return(tbr)
}