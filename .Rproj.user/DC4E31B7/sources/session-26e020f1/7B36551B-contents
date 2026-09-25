circ.gomp <- function(y, x, rads = TRUE, type = "vm", xstand = TRUE, thresh = qchisq(0.95, 1), tol = 1e-6, maxiters = 100) {

  tic <- proc.time()

  if ( xstand )  x <- Rfast::standardise(x)

  d <- dim(x)[2]
  ind <- 1:d
  can <- which( is.na( Rfast::colsums(x) ) )
  ind[can] <- 0

  if ( !is.matrix(y) ) {
    if ( !rads )  y <- y * pi / 180
    y <- cbind( cos(y), sin(y) )
  }
  u <- ( atan(y[, 2]/y[, 1]) + pi * I(y[, 1] < 0) ) %% (2 * pi)

  rho <- circda::circ.mle(u, type = type, tol = tol, maxiters = maxiters)$loglik
  ela <- as.vector( cov(u, x) )
  sel <- which.max( abs(ela) )
  sela <- sel
  names(sela) <- NULL
  mod <- circda::circ.reg( y, x[, sela], type = type, xnew = x[, sela], tol = tol, maxiters = maxiters )
  res <- u - mod$est
  rho[2] <-  mod$loglik
  ind[sel] <- 0
  i <- 2
  while ( 2 * (rho[i] - rho[i - 1]) > thresh ) {
    r <- rep(NA, d)
    i <- i + 1
    r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
    sel <- which.max( abs(r) )
    sela <- c(sela, sel)
    mod <- circda::circ.reg( y, x[, sela], type = type, xnew = x[, sela], tol = tol, maxiters = maxiters )
    res <- u - mod$est
    rho[i] <-  mod$loglik
    ind[sela] <- 0
  } ## end while ( 2 * (rho[i - 1] - rho[i]) > tol )
  runtime <- proc.time() - tic
  len <- length(sela)
  result <- cbind(c(0, sela[-len]), rho[1:len])
  colnames(result) <- c("Selected Vars", "Log-likelihood")
  list(runtime = runtime, result = result)
}




