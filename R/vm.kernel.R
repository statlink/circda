vm.kernel <- function(y, h = "tay", rads = TRUE, ncores = 1) {
  ## y is the data
  ## h is "tay" (Taylor, 2008), "rot" (Garcia-Portugues, 2013) or a number
  ## if the data are in degrees we transform them into radians
  if ( !rads )  y <- y / 180 * pi
  n <- length(y)  ## sample size
  if ( is.character(h) )  h <- .vm.bw(y, h)
  ## exponentially scaled, vm_kernel_cpp() works with exp( h * (cos - 1) )
  com <- (n - 1) * 2 * pi * besselI(h, 0, expon.scaled = TRUE)

  f <- vm_kernel_cpp( rbind( cos(y), sin(y) ), h, com, as.integer(ncores) )

  list( h = h, f = f )
}
