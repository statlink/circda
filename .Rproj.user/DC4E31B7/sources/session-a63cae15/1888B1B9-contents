circ.da <- function(unew, u, ina, rads = TRUE, type = c("vm", "cp", "pn", "gcpc", "cipc") ) {

  ina <- as.numeric(ina)
  pj <- tabulate(ina) / length(ina)
  logpj <- log(pj)
  g <- max(ina)
  if ( !rads )  u <- u * pi/180
  mat <- matrix(0, length(unew), g)
  est <- matrix(NA, nrow = length(unew), 5)
  x <- cbind( cos(u), sin(u) )
  xnew <- cbind( cos(unew), sin(unew) )

  if ( sum( type == "vm") == 1 ) {
    for ( j in 1:g ) {
      mod <- Rfast::vm.mle( u[ina == j], tol = 1e-6 )
      mat[, j] <- Directional::dvm(unew, mod$param[1], mod$param[2], rads = TRUE, logden = TRUE ) + logpj[j]
    }
    est[, 1] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "cp") == 1 ) {
    for ( j in 1:g ) {
      mod <- Directional::purka.mle( x[ina == j, ], tol = 1e-6 )
      mat[, j] <- Directional::dcircpurka(unew, mod$circtheta, mod$alpha, rads = TRUE, logden = TRUE ) + logpj[j]
    }
    est[, 2] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "pn") == 1 ) {
    for ( j in 1:g ) {
      mod <- Directional::spml.mle( u[ina == j], rads = TRUE, tol = 1e-6)
      mat[, j] <- Directional::dspml(unew, mod$mu, rads = TRUE, logden = TRUE ) + logpj[j]
    }
    est[, 3] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "gcpc") == 1 ) {
    for ( j in 1:g ) {
      mod <- Directional::gcpc.mle2( u[ina == j], rads = TRUE)
      mat[, j] <- Directional::dgcpc(unew, mod$circmu, mod$gamma, mod$rho, rads = TRUE, logden = TRUE) + logpj[j]
    }
    est[, 4] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "cipc") == 1 ) {
    for ( j in 1:g ) {
      mod <- Directional::cipc.mle( u[ina == j], rads = TRUE)
      mat[, j] <- Directional::dcipc(unew, mod$circmu, mod$gamma, rads = TRUE, logden = TRUE) + logpj[j]
    }
    est[, 5] <- Rfast::rowMaxs(mat)
  }

  colnames(est) <- c("vM", "PN", "Purka", "GCPC", "CIPC")
  est
}




