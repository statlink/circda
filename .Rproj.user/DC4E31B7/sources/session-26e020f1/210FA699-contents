circ.ridge <- function(y, x, rads = TRUE, type = "vm", lambda = NULL, nlambda = 100, xnew = NULL, tol = 1e-6, maxiters = 100) {
  if ( type == "vm" )  {
    res <- .vm.ridge(y, x, rads, lambda, nlambda, xnew, tol, maxiters)
  } else if ( type == "cipc" ) {
    res <- .cipc.ridge(y, x, rads, lambda, nlambda, xnew, tol, maxiters)
  } else if ( type == "pn" )  {
    res <- .pn.ridge(y, x, rads, lambda, nlambda, xnew, tol, maxiters)
  }
  res
}


.vm.ridge <- function(y, x, rads = TRUE, lambda = NULL, nlambda = 100, xnew = NULL, 
                     tol = 1e-6, maxiters = 100) {

  runtime <- proc.time()

  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }
  n <- dim(y)[1]    
  x <- as.matrix(x)
  p <- dim(x)[2] + 1
  nam <- colnames(x)
  if ( !is.null(nam) ) {
    nam <- c( "Intercept", nam )
  } else   nam <- c( "Intercept", paste("X", 1:(p - 1), sep = "" ) )
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  x <- t( ( t(x) - m ) / s )
  x <- cbind(1, x)
  pen <- rep(1, p)   ;   pen[1] <- 0    ;    pen <- rep(pen, 2)
  ridge <- glmnet::glmnet(x[, -1], y, alpha = 0, nlambda = nlambda, lambda = lambda, 
                          family = "mgaussian", standardize = FALSE)
  lambda <- ridge$lambda
  nlambda <- length(lambda)
  beta <- coef(ridge)
  ini <- cbind( beta[[ 1 ]][, 1], beta[[ 2 ]][, 1] )
  beta <- ini
  be <- list()   ;   loglik_pen <- numeric(nlambda)

  for ( vim in 1:nlambda ) {
    mod <- .vm(y, x, beta, lambda[vim], n, pen, tol, maxiters)
    beta <- mod$be
    colnames(beta) <- paste( c("cos(y)", "sin(y)") )
    rownames(beta) <- nam
    be[[ vim ]] <- beta
    loglik_pen[vim] <- mod$loglik_pen
  }

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- list()
    xnew <- as.matrix(xnew)
    if ( nrow(xnew) == 1 )  xnew <- t(xnew)
    xnew <- t( ( t(xnew) - m ) / s )
    xnew <- cbind(1, xnew)
    for ( i in 1:nlambda ) {
      mu <- xnew %*% be[[ i ]]
      est[[ i ]] <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
      if ( !rads )  est[[ i ]] <- est[[ i ]] * 180 / pi
    }
    names(est) <- paste("lambda=", lambda, sep = "")
  }

  info <- cbind(lambda, loglik_pen)
  colnames(info) <- c("lambda", "loglik")
  runtime <- proc.time() - runtime

  list( runtime = runtime, info = info, be = be, est = est )
}




.vm <- function(y, x, be, lambda, n, pen, tol, maxiters) {

  lambda_vec <- lambda * pen
  lamI <- diag( lambda_vec, nrow = length(lambda_vec) )

  mu <- x %*% be
  ki <- sqrt( Rfast::rowsums(mu^2) )
  lik1 <- sum(mu * y) - sum( log(besselI(ki, 0, expon.scaled = TRUE) ) + ki ) - 
          0.5 * sum( lambda_vec * as.vector(be)^2 )

  A1  <- besselI(ki, 1) / besselI(ki, 0)                 
  A1d <- 1 - A1 / ki - A1^2                               
  resid <- y - (A1 / ki) * mu                             
  grad <- crossprod(x, resid)                            
  mu_ki  <- mu / ki                                       
  d11 <- A1d * mu_ki[, 1]^2 + (A1 / ki) * (1 - mu_ki[, 1]^2)
  d22 <- A1d * mu_ki[, 2]^2 + (A1 / ki) * (1 - mu_ki[, 2]^2)
  d12 <- (A1d - A1/ki) * mu_ki[, 1] * mu_ki[, 2]
  H11 <-  -crossprod(x * d11, x)
  H22 <-  -crossprod(x * d22, x)
  H12 <-  -crossprod(x * d12, x)
  H <- rbind( cbind(H11, H12), cbind(H12, H22) ) - lamI 
  be <- be - solve(H, as.vector(grad) - lambda_vec * as.vector(be) )
  be <- matrix(be, ncol = 2)
  mu <- x %*% be
  ki <- sqrt( Rfast::rowsums(mu^2) )
  lik2 <- sum(mu * y) - sum( log(besselI(ki, 0, expon.scaled = TRUE) ) + ki ) - 
          0.5 * sum( lambda_vec * as.vector(be)^2 )

  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {
    i <- i + 1
    lik1 <- lik2
    A1  <- besselI(ki, 1) / besselI(ki, 0)                  
    A1d <- 1 - A1 / ki - A1^2                               
    resid <- y - (A1 / ki) * mu                             
    grad <- crossprod(x, resid)                            
    mu_ki  <- mu / ki                                      
    d11 <- A1d * mu_ki[, 1]^2 + (A1 / ki) * (1 - mu_ki[, 1]^2)
    d22 <- A1d * mu_ki[, 2]^2 + (A1 / ki) * (1 - mu_ki[, 2]^2)
    d12 <- (A1d - A1/ki) * mu_ki[, 1] * mu_ki[, 2]
    H11 <-  -crossprod(x * d11, x)
    H22 <-  -crossprod(x * d22, x)
    H12 <-  -crossprod(x * d12, x)
    H <- rbind( cbind(H11, H12), cbind(H12, H22) ) - lamI
    be <- be - solve(H, as.vector(grad) - lambda_vec * as.vector(be) )
    be <- matrix(be, ncol = 2)
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    lik2 <- sum(mu * y) - sum( log(besselI(ki, 0, expon.scaled = TRUE) ) + ki ) - 
            0.5 * sum( lambda_vec * as.vector(be)^2 )
  }

  list(be = be, loglik_pen = lik2 - n * log(2 * pi) )
}


.cipc.ridge <- function(y, x, rads = TRUE, lambda = NULL, nlambda = 100, xnew = NULL, 
                     tol = 1e-6, maxiters = 100) {

  runtime <- proc.time()

  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }
  n <- dim(y)[1]    
  x <- as.matrix(x)
  p <- dim(x)[2] + 1
  nam <- colnames(x)
  if ( !is.null(nam) ) {
    nam <- c( "Intercept", nam )
  } else   nam <- c( "Intercept", paste("X", 1:(p - 1), sep = "" ) )
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  x <- t( ( t(x) - m ) / s )
  x <- cbind(1, x)
  pen <- rep(1, p)   ;   pen[1] <- 0    ;    pen <- rep(pen, 2)
  ridge <- glmnet::glmnet(x[, -1], y, alpha = 0, nlambda = nlambda, lambda = lambda, 
                          family = "mgaussian", standardize = FALSE)
  lambda <- ridge$lambda
  nlambda <- length(lambda)
  beta <- coef(ridge)
  ini <- cbind( beta[[ 1 ]][, 1], beta[[ 2 ]][, 1] )
  beta <- ini
  be <- list()   ;   loglik_pen <- numeric(nlambda)

  for ( vim in 1:nlambda ) {
    mod <- .cipc(y, x, beta, lambda[vim], n, p, pen, tol, maxiters)
    beta <- mod$be
    colnames(beta) <- paste( c("cos(y)", "sin(y)") )
    rownames(beta) <- nam
    be[[ vim ]] <- beta
    loglik_pen[vim] <- mod$loglik_pen
  }

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- list()
    xnew <- as.matrix(xnew)
    if ( nrow(xnew) == 1 )  xnew <- t(xnew)
    xnew <- t( ( t(xnew) - m ) / s )
    xnew <- cbind(1, xnew)
    for ( i in 1:nlambda ) {
      mu <- xnew %*% be[[ i ]]
      est[[ i ]] <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
      if ( !rads )  est[[ i ]] <- est[[ i ]] * 180 / pi
    }
    names(est) <- paste("lambda=", lambda, sep = "")
  }

  info <- cbind(lambda, loglik_pen)
  colnames(info) <- c("lambda", "loglik")
  runtime <- proc.time() - runtime

  list( runtime = runtime, info = info, be = be, est = est )
}




.cipc <- function(y, x, be, lambda, n, p, pen, tol, maxiters) {

  lambda_vec <- lambda * pen
  lamI <- diag( lambda_vec, nrow = length(lambda_vec) )
  H <- matrix(0, 2 * p, 2 * p)

  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  a <- Rfast::rowsums(y * mu)
  com <- sqrt(g2 + 1)
  com2 <- com - a
  lik <-  - sum( log( com2 ) ) - 0.5 * sum( lambda_vec * as.vector(be)^2 ) 

  muc_y <- mu / com - y
  der1 <- Rfast::eachcol.apply(x, muc_y[, 1] / com2 )
  der2 <- Rfast::eachcol.apply(x, muc_y[, 2] / com2 )
  ### Jacobian of b1
  a1 <- ( com - mu[, 1]^2 / com ) / ( com^2 * com2 )
  up1 <- crossprod(x, x * a1)
  up2 <- crossprod(x * muc_y[, 1]/com2)
  H[1:p, 1:p] <- up2 - up1
  ### Jacobian of b2
  a1 <- ( com - mu[, 2]^2 / com ) / ( com^2 * com2 )
  up1 <- crossprod(x, x * a1)
  up2 <- crossprod(x * muc_y[, 2]/com2)
  H[(p + 1):(2*p), (p + 1):(2*p)] <- up2 - up1
  ### Jacobian of b12
  a1 <- mu[, 1] * mu[, 2] / ( com^3 * com2)
  up1 <- crossprod(x, x * a1)
  up2 <- crossprod(x * muc_y[, 1]/com2, x * muc_y[, 2]/com2)
  H[1:p, (p + 1):(2*p)] <- H[(p + 1):(2*p), 1:p] <- up2 + up1

  be <- be + solve(H - lamI , c(der1, der2))
  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  a <- Rfast::rowsums(y * mu)
  com <- sqrt(g2 + 1)
  com2 <- com - a
  lik[2] <-  - sum( log( com2 ) ) - 0.5 * sum( lambda_vec * as.vector(be)^2 ) 

  i <- 2
  while ( lik[i] - lik[i-1] > tol  &  i < maxiters ) {
    i <- i + 1
    muc_y <- mu / com - y
    der1 <- Rfast::eachcol.apply(x, muc_y[, 1] / com2 )
    der2 <- Rfast::eachcol.apply(x, muc_y[, 2] / com2 )
    ### Jacobian of b1
    a1 <- ( com - mu[, 1]^2 / com ) / ( com^2 * com2 )
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(x * muc_y[, 1]/com2)
    H[1:p, 1:p] <- up2 - up1
    ### Jacobian of b2
    a1 <- ( com - mu[, 2]^2 / com ) / ( com^2 * com2 )
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(x * muc_y[, 2]/com2)
    H[(p + 1):(2*p), (p + 1):(2*p)] <- up2 - up1
    ### Jacobian of b12
    a1 <- mu[, 1] * mu[, 2] / ( com^3 * com2)
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(x * muc_y[, 1]/com2, x * muc_y[, 2]/com2)
    H[1:p, (p + 1):(2*p)] <- H[(p + 1):(2*p), 1:p] <- up2 + up1

    be <- be + solve(H - lamI, c(der1, der2))
    mu <- x %*% be
    g2 <- Rfast::rowsums(mu^2)
    a <- Rfast::rowsums(y * mu)
    com <- sqrt(g2 + 1)
    com2 <- com - a
    lik[i] <-  - sum( log( com2 ) ) - 0.5 * sum( lambda_vec * as.vector(be)^2 ) 
  }

  list(be = be, loglik_pen = lik[i] - n * log(2 * pi) )
}


.pn.ridge <- function(y, x, rads = TRUE, lambda = NULL, nlambda = 100, xnew = NULL, 
                     tol = 1e-6, maxiters = 100) {

  runtime <- proc.time()

  if ( is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    u <- y
    ci <- u[, 1]   ;   si <- u[, 2]
  } else {
    ci <- cos(y)   ;   si <- sin(y)
    u <- cbind(ci, si)
  }

  n <- dim(y)[1]    
  x <- as.matrix(x)
  p <- dim(x)[2] + 1
  nam <- colnames(x)
  if ( !is.null(nam) ) {
    nam <- c( "Intercept", nam )
  } else   nam <- c( "Intercept", paste("X", 1:(p - 1), sep = "" ) )
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  x <- t( ( t(x) - m ) / s )
  x <- cbind(1, x)
  pen <- rep(1, p)   ;   pen[1] <- 0    ;    pen <- rep(pen, 2)
  ridge <- glmnet::glmnet(x[, -1], y, alpha = 0, nlambda = nlambda, lambda = lambda, 
                          family = "mgaussian", standardize = FALSE)
  lambda <- ridge$lambda
  nlambda <- length(lambda)
  beta <- coef(ridge)
  ini <- cbind( beta[[ 1 ]][, 1], beta[[ 2 ]][, 1] )
  beta <- ini
  be <- list()   ;   loglik_pen <- numeric(nlambda)

  for ( vim in 1:nlambda ) {
    mod <- .pn(y, ci, si, x, beta, lambda[vim], n, p, pen, tol, maxiters)
    beta <- mod$be
    colnames(beta) <- paste( c("cos(y)", "sin(y)") )
    rownames(beta) <- nam
    be[[ vim ]] <- beta
    loglik_pen[vim] <- mod$loglik_pen
  }

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- list()
    xnew <- as.matrix(xnew)
    if ( nrow(xnew) == 1 )  xnew <- t(xnew)
    xnew <- t( ( t(xnew) - m ) / s )
    xnew <- cbind(1, xnew)
    for ( i in 1:nlambda ) {
      mu <- xnew %*% be[[ i ]]
      est[[ i ]] <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
      if ( !rads )  est[[ i ]] <- est[[ i ]] * 180 / pi
    }
    names(est) <- paste("lambda=", lambda, sep = "")
  }

  info <- cbind(lambda, loglik_pen)
  colnames(info) <- c("lambda", "loglik")
  runtime <- proc.time() - runtime

  list( runtime = runtime, info = info, be = be, est = est )
}




.pn <- function(u, ci, si, x, be, lambda, n, p, pen, tol, maxiters) {

  lambda_vec <- lambda * pen
  lamI <- diag( lambda_vec, nrow = length(lambda_vec) )
  f <-  - 0.5   ;   con <- sqrt(2 * pi)
 
  mu <- x %*% be
  tau <- rowsums(u * mu)
  ptau <- pnorm(tau)
  lik1 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) ) - 
               0.5 * sum( lambda_vec * as.vector(be)^2 ) 
  rat <- ptau / ( exp(f * tau^2)/con + tau * ptau )
  psit <- tau + rat
  psit2 <- 2 - tau * rat - rat^2
  der <- as.vector( crossprod(x, - mu + psit * u) )
  a11 <- crossprod(x, x * (psit2 * ci^2 - 1) )
  a12 <- crossprod(x, x * (psit2 * ci * si ) )
  a22 <- crossprod(x, x * (psit2 * si^2 - 1 ) )
  der2 <- cbind( rbind(a11, a12), rbind(a12, a22) )
  be <- be - solve(der2 - lamI, der)
  mu <- x %*% be
  tau <- rowsums(u * mu)
  ptau <- pnorm(tau)
  lik2 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) ) - 
             0.5 * sum( lambda_vec * as.vector(be)^2 ) 
  i <- 2
  ## mono th while
  while ( abs(lik2 - lik1) > tol  & i < maxiters ) {
    lik1 <- lik2
    i <- i + 1
    rat <- ptau / ( exp(f * tau^2)/con + tau * ptau )
    psit <- tau + rat
    psit2 <- 2 - tau * rat - rat^2
    der <- as.vector( crossprod(x, - mu + psit * u) )
    a11 <- crossprod(x, x * (psit2 * ci^2 - 1) )
    a12 <- crossprod(x, x * (psit2 * ci * si ) )
    a22 <- crossprod(x, x * (psit2 * si^2 - 1 ) )
    der2 <- cbind( rbind(a11, a12), rbind(a12, a22) )
    be <- be - solve(der2 - lamI, der)
    mu <- x %*% be
    tau <- rowsums(u * mu)
    ptau <- pnorm(tau)
    lik2 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) ) - 
               0.5 * sum( lambda_vec * as.vector(be)^2 ) 

  }

  list(be = be, loglik_pen = lik2 - n * log(2 * pi) )
}