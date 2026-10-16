perm.reg <- function(y, x, rads = TRUE, type = "vm", tol = 1e-6, maxiters = 500, R = 999) {

  if ( !rads )  y <- y * pi / 180
  n <- length(y)
  
  m0 <- circda::circ.mle(y, rads = TRUE, type = type, tol = tol, maxiters = maxiters)
  if ( type == "vm" ) {
    m <- m0$param[1]
  } else if ( type == "pn" ) {
    m <- m0$mumu
  } else if ( type == "cipc" ) {
    m <- m0$circmu 
  }  

  stat <- circda::circ.reg(y = y, x = x, rads = TRUE, type = type, tol = tol, maxiters = maxiters)$loglik - m0$loglik
  bstat <- numeric(R)
  for ( i in 1:R ) {
    yb <- y[sample.int(n, n)]
    bstat[i] <- circda::circ.reg(y = yb, x = x, rads = TRUE, type = type, tol = tol, maxiters = maxiters)$loglik - m0$loglik
  }
  pv <- ( sum(bstat >= stat) + 1 ) / (R + 1)
  
  res <- c(stat, pv)
  names(res) <- c("statistic", "p-value")
  res
}