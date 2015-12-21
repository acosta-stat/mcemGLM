#' Contrast estimation for mcemGLMM objects
#' 
#' Contrast testing for a model fitted with \code{mcemGLMM}.
#' 
#' 
#' @param object a model fitted with the mcemGLMM function.
#' @param ctr.mat contrast matrix. Each row corresponds to a linear combination
#' of the fixed effect coefficients.
#' @param ... additional arguments.
#' @return A matrix with each row corresponding to a contrast and with the
#' following columns: \describe{ \item{Estimate}{contrasts's point estimator.}
#' \item{Std. Err.}{standard error for each fitted contrast.} \item{Wald}{Wald
#' statistic for each fitted contrast.} \item{Adj. p-value}{p-value adjusted
#' for multiple comparison (Bonferroni.)} }
#' @author Felipe Acosta Archila <acosta@@umn.edu>
#' @keywords glmm
#' @examples
#' 
#' \donttest{
#' set.seed(72327)
#' data(exdata)
#' fit <- mcemGLMM(obs ~ z2 + x, random = ~ 0 + z1, 
#'                 data = exdata, 
#'                 family = "bernoulli", vcDist = "normal", 
#'                 controlEM = list(verb = FALSE, EMit = 15, MCit = 10000), 
#'                 initial = c(-0.13, -0.15, -0.21, 1.59, 0.002))
#' mat <- rbind("1 - 2" = c(0, -1, 0, 0), "1 - 3" = c(0, 0, -1, 0), "2 - 3" = c(0, 1, -1, 0))
#' contrasts.mcemGLMM(fit, mat)
#' }
#' 
#' @export
#' 
contrasts.mcemGLMM <- function(object, ctr.mat) {
  names0 <- rownames(ctr.mat)
  coef0 <- as.vector(tail(object$mcemEST, n = 1))
  kP <- ncol(object$x)
  kR <- length(coef0) - kP
  for (i in 1:kR) {
    ctr.mat <- cbind(ctr.mat, 0)
  }
  wald0 <- rep(0, nrow(ctr.mat))
  ctr.est0 <- rep(0, nrow(ctr.mat))
  ctr.err0 <- rep(0, nrow(ctr.mat))
  for (i in 1:nrow(ctr.mat)) {
    cm <- t(ctr.mat[i, ])
    ctr.est0[i] <- cm %*% coef0
    ctr.err0[i] <- solve(cm %*% solve(object$iMatrix) %*% t(cm))
    wald0[i] <- ctr.est0[i]^2 * ctr.err0[i]
  }
  std.err <- 1/sqrt(ctr.err0)
  pval0 <- pchisq(wald0, 1, lower.tail = FALSE) * nrow(ctr.mat)
  pval0 <- ifelse(pval0 < 1, pval0, 1)
  tbr <- cbind(ctr.est0, std.err, wald0, pval0)
  colnames(tbr) <- c("Estimate", "Std. Err.", "Wald", "Adj. p-value")
  rownames(tbr) <- rownames(ctr.mat)
  return(tbr)
}