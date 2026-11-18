vm.nb <- function(ynew, y, ina, h = "tay", rads = TRUE, ncores = 1) {
  d <- dim(y)[2]
  ynew <- as.matrix(ynew)
  if ( dim(ynew)[2] == 1 )  ynew <- t(ynew)
  score <- numeric( dim(ynew)[1] )
  for ( j in 1:d ) {
    score <- score + Directionalkernel::vm.kernel.da(ynew[, j], y[, j], ina,  h, rads, ncores)$score
  }
  list(class = Rfast::rowMaxs( score ), score = score )
}

