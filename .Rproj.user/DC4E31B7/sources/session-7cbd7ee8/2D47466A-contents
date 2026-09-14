colcirc.mle <- function(u, rads = TRUE, type = "vm", tol = 1e-07, maxiters = 100, parallel = FALSE) {

  if ( !rads )  u <- u * pi / 180

  if ( type == "vm" ) {
    res <- Rfast::colvm.mle(u, tol = tol)
  } else if ( type == "spml" ) {
    res <- Rfast2::colspml.mle(u, tol = tol, maxiters = maxiters, parallel = parallel)
  } else if ( type == "cp" ) {
    res <- .colpurka.mle(u, tol = tol)
  } else if ( type == "cipc" ) {
    res <- .colcipc.mle(u, tol = tol)
  } else if ( type == "gcpc" ) {
    res <- .colgcpc.mle(u)
  }

  res
}



.colpurka.mle <- function(x, tol) {
  res <- matrix(NA, dim(x)[2], 5)
  for ( i in 1:dim(x)[2] )  res[i, ] <- unlist( Rfast2::purka.mle(x[, i]) )
  colnames(res) <- c("theta1", "theta2", "alpha", "loglik", "alpha.sd")
  res
}


.colgcpc.mle <- function(x) {
  res <- matrix(NA, dim(x)[2], 6)
  for ( i in 1:dim(x)[2] )  res[i, ] <- unlist( Directional::gcpc.mle(x[, i], rads = TRUE) )
  colnames(res) <- c("mu1", "mu2", "circmu", "gamma", "rho", "loglik")
  res
}


.colcipc.mle <- function(x, tol) {
  res <- matrix(NA, dim(x)[2], 5)
  for ( i in 1:dim(x)[2] )  res[i, ] <- unlist( Directional::cipc.mle(x[, i], rads = TRUE, tol = tol) )
  colnames(res) <- c("mu1", "mu2", "circmu", "gamma", "loglik")
  res
}


