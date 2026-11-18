vm.kernel <- function(y, h = "tay", rads = TRUE, ncores = 1) {
  ## u is the data
  ## h is the bandwidth you want
  ## thumb is either 'none' (defualt), or 'tay' (Taylor, 2008) or
  ## 'rot' (Garcia-Portugues, 2013)
  ## if the data are in degrees we transform them into radians
  if ( !rads )  y <- y / 180 * pi
  n <- length(y)  ## sample size
  u <- y
  y <- rbind( cos(y), sin(y) )

  if ( h == "tay"  |  h == "rot" ) {
    est <- Rfast::rowmeans(y)
    mesos <- atan2(est[2], est[1]) %% (2 * pi)
    lik <- function(k)   k * sum( cos(u - mesos) ) - n * ( log(besselI( k, 0, expon.scaled = TRUE) ) + k )
    k <- optimize(lik, c(0, 50000), maximum = TRUE, tol = 1e-6)$maximum
    if ( h == "tay" ) {
      h <- ( (3 * n * k^2 * besselI(2 * k, 2)) / ( 4 * pi^0.5 * besselI(k, 0)^2 ) )^0.4
    } else if ( h == "rot" ) {
      h <- ( ( k * n * ( 2 * besselI(2 * k, 1) + 3 * k * besselI(2 * k, 2)) ) / ( 4 * pi^0.5 * besselI(k, 0)^2 ) )^0.4
    }
  } else h <- h
  com <- (n - 1) * 2 * pi * besselI(h, 0)

  f <- vm_kernel_cpp(y, h, com, as.integer(ncores))

  list( h = h, f = f )
}
