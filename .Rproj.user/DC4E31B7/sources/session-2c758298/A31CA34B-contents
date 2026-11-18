vmf.kernel.da <- function(ynew, y, ina, h = "rot", ncores = 1) {
  ynew <- as.matrix(ynew)
  if ( dim(ynew)[2] == 1 )  ynew <- t(ynew)
  y <- as.matrix(y)
  p <- dim(y)[2]   ;   nu <- p/2 - 1
  ina <- as.numeric(ina)  ;  k <- max(ina)  ;  ni <- tabulate(ina)
  prior <- ni / length(ina)
  if ( is.character(h) ) {
    hrule <- h  ;  h <- numeric(k)
    for (j in 1:k)  h[j] <- Directionalkernel::vmf.kernel(y[ina == j, , drop = FALSE], h = hrule, ncores = ncores)$h
  } else if ( length(h) == 1 )  h <- rep(h, k)
  kappa <- 1 / h^2
  ## log C_p(kappa)
  logcon <- nu * log(kappa) - 0.5 * p * log(2 * pi) -
    log( besselI(kappa, nu, expon.scaled = TRUE) ) - kappa
  cons <- log(prior) - log(ni) + logcon
  res <- kernel_da_cpp( t(ynew), t(y), as.integer(ina), kappa, cons, as.integer(ncores) )
  list( h = h, class = res$class, score = res$score )
}
