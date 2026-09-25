perm.reg <- function(y, x, rads = TRUE, type = "vm", tol = 1e-6, maxiters = 500, R = 999) {

  if ( !rads )  y <- y * pi / 180
  n <- length(y)
  bstat <- numeric(R)
  
  m0 <- circda::circ.mle(y, rads = TRUE, type = type, tol = tol, maxiters = maxiters)
  if ( type == "vm" ) {
    m <- m0$param[1]
  } else if ( type == "pn" ) {
    m <- m0$mumu
  } else if ( type == "cipc" ) {
    m <- m0$circmu 
  }  

  stat <- circ.reg(y = y, x = x, rads = rads, type = type, tol = tol, maxiters = maxiters)$loglik -
          circda::circ.mle(y, rads = TRUE, type = type, tol = tol, maxiters = maxiters)$loglik
  ei <- y - m
  ei <- cbind( cos(ei), sin(ei) )
  ei <- ( atan(ei[, 2]/ei[, 1]) + pi * I(ei[, 1] < 0) ) %% (2 * pi)  
  for ( i in 1:R ) {
    eib <- ei[sample.int(n, n)]     
    yb <- (m + eib) %% (2 * pi)
    bstat[i] <- circda::circ.reg(y = yb, x = x, rads = TRUE, type = type, tol = tol, maxiters = maxiters)$loglik - 
                circda::circ.mle(yb, rads = TRUE, type = type, tol = tol, maxiters = maxiters)$loglik
  }
  pv <- ( sum(bstat >= stat) + 1 ) / (R + 1) 
  
  res <- c(stat, pv)
  names(res) <- c("statistic", "p-value")
  res
}