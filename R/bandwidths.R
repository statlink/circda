## Rule of thumb bandwidths. They need only the data summaries, not the kernel
## density estimate, so the discriminant analysis functions call them directly.

## u: angles in radians. h: "tay" (Taylor, 2008) or "rot" (Garcia-Portugues, 2013).
## The value returned is the concentration of the von Mises kernel.
.vm.bw <- function(u, h) {
  if ( !h %in% c("tay", "rot") )  stop( 'h must be "tay", "rot" or numeric' )
  n <- length(u)
  est <- Rfast::rowmeans( rbind( cos(u), sin(u) ) )
  Rbar <- sqrt( sum(est^2) )  ## mean resultant length
  lik <- function(k)  k * Rbar - log( besselI(k, 0, expon.scaled = TRUE) ) - k
  k <- optimize(lik, c(0, 50000), maximum = TRUE, tol = 1e-6)$maximum
  i0 <- besselI(k, 0, expon.scaled = TRUE)^2
  i2 <- besselI(2 * k, 2, expon.scaled = TRUE)
  if ( h == "tay" ) {
    res <- 3 * n * k^2 * i2 / ( 4 * pi^0.5 * i0 )
  } else {
    i1 <- besselI(2 * k, 1, expon.scaled = TRUE)
    res <- k * n * ( 2 * i1 + 3 * k * i2 ) / ( 4 * pi^0.5 * i0 )
  }
  res^0.4
}

## y: matrix of unit vectors, one per row. h: "rot" (Garcia-Portugues, 2013).
## The value returned is h, where the concentration of the kernel is 1 / h^2.
.vmf.bw <- function(y, h) {
  if ( h != "rot" )  stop( 'h must be "rot" or numeric' )
  n <- dim(y)[1]
  q <- dim(y)[2] - 1
  k <- Rfast::vmf.mle(y)$kappa
  up <- 4 * pi^0.5 * besselI(k, (q - 1)/2, expon.scaled = TRUE)^2
  down <- k^( (q + 1)/2 ) * n * ( 2 * q * besselI(2 * k, (q + 1)/2, expon.scaled = TRUE) +
          (2 + q) * k * besselI(2 * k, (q + 3)/2, expon.scaled = TRUE) )
  ( up / down )^( 1/(4 + q) )
}
