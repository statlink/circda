vmf.kernel <- function(y, h = "rot", ncores = 1) {
  ## y is the data
  ## h is "rot" (Garcia-Portugues, 2013) or a number
  p <- dim(y)[2]  ## dimensionality of the data
  n <- dim(y)[1]  ## sample size of the data
  if ( is.character(h) )  h <- .vmf.bw(y, h)
  ## exponentially scaled, vmf_kernel_cpp() works with exp( (y'y - 1) / h^2 )
  cpk <- (1/h^2)^(p/2 - 1) / ( (2 * pi)^(p/2) * besselI(1/h^2, p/2 - 1, expon.scaled = TRUE) )
  com <- cpk / (n - 1)

  f <- vmf_kernel_cpp(t(y), h, com, as.integer(ncores))
  list( h = h, f = f )
}
