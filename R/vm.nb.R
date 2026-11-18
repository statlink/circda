vm.nb <- function(ynew, y, ina, h = "tay", rads = TRUE, ncores = 1) {
  d <- dim(y)[2]
  ina <- as.numeric(ina)
  k <- max(ina)
  ynew <- as.matrix(ynew)
  if ( dim(ynew)[2] == 1 )  ynew <- t(ynew)
  score <- matrix(nrow = dim(ynew)[1], ncol = k)
  for ( j in 1:d ) {
    score <- score + Directionalkernel::vm.kernel.da(ynew[, j], y[, j], ina,  h, rads, ncores)$score
  }
  list(class = Rfast::rowMaxs( score ), score = score )
}

