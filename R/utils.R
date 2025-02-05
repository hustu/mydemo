
betaRA <- function(M, X, Y, COV, pscore.c=NULL) {
  b <- summary(glm(Y ~ M + X + COV))[["coefficients"]][2, 4]
  return(b)
}

betaPS <- function(M, X, Y, COV, pscore.c) {
  b <- summary(glm(Y ~ M + X + pscore.c))[["coefficients"]][2, 4]
  return(b)
}

betaIPW <- function(M, X, Y, COV, pscore.c) {
  wt4x <- ifelse(X==1, 1/pscore.c, 1/(1-pscore.c))    # weight for X
  coef_x <-  glm(Y ~ X + M, weights = wt4x)[["coefficients"]][2]
  b <- summary(glm((Y - coef_x*X) ~ M))[["coefficients"]][2,4]
  return(b)
}

betaOW <- function(M, X, Y, COV, pscore.c) {
  wt4x <- ifelse(X==1, 1-pscore.c, pscore.c)    # weight for X
  coef_x <- glm(Y ~ X + M, weights = wt4x)[["coefficients"]][2]
  b <- summary(glm((Y - coef_x*X) ~ M))[["coefficients"]][2,4]
  return(b)
}
