mixcirc.mle <- function(u, type = "vm", rads = TRUE, g = 2, tol = 1e-4, maxiters = 500) {

  if ( type == "vm" ) {
    res <- .mixvm.mle(u, rads, g, tol, maxiters)
    mu <- res$param[, 3:4]
    mu <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
    res$param <- cbind( res$param[, 1], mu, res$param[, 2] )
    colnames(res$param) <- c("probs", "mu", "kappa")
  } else if ( type == "cp" ) {
    res <- .mixpurka.mle(u, rads, g, tol, maxiters)
    mu <- res$param[, 2:3]
    mu <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
    res$param <- cbind( res$param[, 1], mu, res$param[, 4] )
    colnames(res$param) <- c("probs", "theta", "alpha")
  } else if ( type == "pn" ) {
    res <- .mixpn.mle(u, rads, g, tol, maxiters)
    mu <- res$param[, 2:3]
    mu <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
    res$param <- cbind( res$param[, 1], mu, res$param[, 4] )
    colnames(res$param) <- c("probs", "mu", "gamma")
  } else if ( type == "gcpc" ) {
    res <- .mixgcpc.mle(u, rads, g, tol, maxiters)
  } else if ( type == "cipc" ) {
    res <- .mixcipc.mle(u, rads, g, tol, maxiters)
    mu <- res$param[, 3:4]
    mu <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
    res$param <- cbind( res$param[, 1], mu, res$param[, 2] )
    colnames(res$param) <- c( "probs", "omega", "gamma")
  }
  res
}


.mixvm.mle <- function(u, rads = TRUE, g = 2, tol = 1e-4, maxiters = 100) {

  if ( ! rads )  u <- u * pi / 180
  x <- cbind( cos(u), sin(u) )
  res <- Directional::mixvmf.mle(x, g, n.start = 10, tol = tol, maxiters = maxiters)
  names(res)[[ 4 ]] <- "probs"
  res
}


.mixcipc.mle <- function(u, rads = TRUE, g = 2, tol = 1e-4, maxiters = 100) {

  if ( ! rads )  u <- u * pi / 180
  x <- cbind( cos(u), sin(u) )
  res <- Directional::mixspcauchy.mle(x, g, n.start = 10, tol = tol, maxiters = maxiters)
  names(res)[[ 5 ]] <- "probs"
  res
}


.mixpurka.mle <- function(u, rads = TRUE, g = 2, tol = 1e-4, maxiters = 100) {

  fun2 <- function(wlika, rswlika, x, g, param, lika) {

    wij <- wlika / rswlika  ## weights
    pj <- Rfast::colmeans(wij)   # PANOS

    for ( j in 1:g ) {
      mod <- .wpurka.mle(x, w = wij[, j])
      theta <- mod$theta
      a <- mod$alpha
      param[j, ] <- c(theta, a)
      v <- drop(x %*% theta)
      v <- pmin(pmax(v, -1), 1)
      lika[, j] <- log(a) - log( 1 - exp(-a * pi) ) - a * acos(v) + log( pj[j] )
    }

    wlika <- exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika / rswlika
    list(wij = wij, param = param, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  if ( !rads )  u <- u * pi / 180
  x <- cbind( cos(u), sin(u) )
  n <- dim(x)[1]

  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  param <- matrix(nrow = g, ncol = 3)

  runtime <- proc.time()
  ## Step 1

  if ( g > 1 ) {
    cl <- kmeans(x, g, nstart = 10)$cl

  } else  cl <- rep(1, n)
  wij <- tabulate(cl)

  while ( min(wij) <= 6 ) {
    g <- g - 1
    lika <- matrix(nrow = n, ncol = g)
    param <- matrix(nrow = g, ncol = 3)
    cl <- kmeans(x, g, nstart = 10)$cl
    wij <- tabulate(cl)
  }

  for ( j in 1:g ) {
    mod <- Directional::purka.mle(x[cl == j, ])
    theta <- mod$theta
    a <- mod$alpha
    param[j, ] <- c(theta, a)
    v <- drop(x %*% theta)
    v <- pmin(pmax(v, -1), 1) 
    lika[, j] <- log(a) - log(1 - exp(-a * pi)) - a * acos(v)
  }

  wlika <- exp(lika)
  rswlika <- Rfast::rowsums(wlika)

  ep <- fun2(wlika, rswlika, x, g, param, lika)
  lik[1] <- ep$lik
  ep2 <- fun2(ep$wlika, ep$rswlika, x, g, ep$param, ep$lika)
  lik[2] <- ep2$lik

  i <- 2
  while ( abs(lik[i] - lik[i - 1] ) > tol & i < maxiters ) {
    i <- i + 1
    ep <- ep2
    ep2 <- fun2(ep$wlika, ep$rswlika, x, g, ep$param, ep$lika)
    lik[i] <- ep2$lik
  }
  res <- ep2
  if ( ep$lik > ep2$lik )  res <- ep

  pj <- Rfast::colmeans(res$wij)
  probs <- res$wij
  loglik <- res$lik
  ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
  param <- cbind(pj, res$param)
  colnames(param) <- c( "probs", paste("theta", 1:2, sep = ""), "alpha" )
  rownames(param) <- paste("Cluster", 1:g, sep = " ")

  runtime <- proc.time() - runtime
  list( param = param, loglik = loglik - n * log(2),
        pred = ta, probs = probs, iters = i, runtime = runtime )
}



.wpurka.mle <- function(x, w = NULL, tol = 1e-07) {
  n <- dim(x)[1]
  sw <- sum(w)
  eps <- 1e-6

  theta <- Rfast::mediandir(x)  # unweighted start
  for ( iter in 1:500 ) {
    ct <- drop(x %*% theta)
    ct <- pmin(pmax(ct, -1 + eps), 1 - eps)   # keep away from ±1
    d <- pmax(acos(ct), eps)
    wts <- w / d
    g <- Rfast::eachcol.apply(x, wts)
    g_norm  <- sqrt( sum(g^2) )
    if ( g_norm < eps ) break          # already at a Weiszfeld fixed point
    theta_new <- g / g_norm
    if ( sqrt( sum( (theta_new - theta)^2 ) ) < tol ) {
      theta <- theta_new
      break
    }
    theta <- theta_new
  }

  a <- drop(x %*% theta)
  a[ is.na(a) ] <- 0
  a <- pmin(pmax(a, -1), 1)
  A <- sum( w * acos(a) )           # WEIGHTED sum of geodesic distances
  circle <- function(a, A, sw)  sw * log(a) - sw * log( 1 - exp(-a * pi) ) - a * A

  lika <- optimize(circle, c(0.001, 1000), maximum = TRUE, A = A, sw = sw, tol = tol)
  list(theta = theta, alpha = lika$maximum, loglik = lika$objective - sw * log(2) )
}


.mixpn.mle <- function(u, rads = TRUE, g = 2, tol = 1e-4, maxiters = 100) {

  fun2 <- function(wlika, rswlika, x, g, param, lika) {

    wij <- wlika / rswlika  ## weights
    pj <- Rfast::colmeans(wij)   # PANOS
    f <-  - 0.5   ;   con <- sqrt(2 * pi)

    for ( j in 1:g ) {
      mod <- .wpn.mle(u, w = wij[, j])
      mu <- drop(mod$mu)
      gam <- sum(mu^2)
      tau <- drop( x %*% mu )
      ptau <- pnorm(tau)
      param[j, ] <- c(mu, sqrt(gam) )
      lika[, j] <-  - 0.5 * gam + log1p( tau * ptau * con / exp(f * tau^2) ) + log( pj[j] )
    }

    wlika <- exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika / rswlika
    list(wij = wij, param = param, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  if ( ! rads )  u <- u * pi / 180
  x <- cbind( cos(u), sin(u) )
  n <- dim(x)[1]
  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  param <- matrix(nrow = g, ncol = 3)

  runtime <- proc.time()
  ## Step 1

  if ( g > 1 ) {
    cl <- kmeans(x, g, nstart = 10)$cl

  } else  cl <- rep(1, n)
  wij <- tabulate(cl)

  while ( min(wij) <= 6 ) {
    g <- g - 1
    lika <- matrix(nrow = n, ncol = g)
    param <- matrix(nrow = g, ncol = 3)
    cl <- kmeans(x, g, nstart = 10)$cl
    wij <- tabulate(cl)
  }

  f <-  - 0.5   ;   con <- sqrt(2 * pi)
  for ( j in 1:g ) {
    mod <- Directional::spml.mle(u[cl == j], rads = TRUE)
    mu <- drop(mod$mu)
    gam <- sum(mu^2)
    tau <- drop( x %*% mu )
    ptau <- pnorm(tau)
    param[j, ] <- c( mu, sqrt(gam) )
    lika[, j] <- - 0.5 * gam + log1p( tau * ptau * con / exp(f * tau^2) )
  }

  wlika <- exp(lika)
  rswlika <- Rfast::rowsums(wlika)

  RES <- try({
    ep <- fun2(wlika, rswlika, x, g, param, lika)
    lik[1] <- ep$lik
    ep2 <- fun2(ep$wlika, ep$rswlika, x, g, ep$param, ep$lika)
    lik[2] <- ep2$lik

    i <- 2
    while ( lik[i] - lik[i - 1] > tol & i < maxiters ) {
      i <- i + 1
      ep <- ep2
      ep2 <- fun2(ep$wlika, ep$rswlika, x, g, ep$param, ep$lika)
      lik[i] <- ep2$lik
    }
    res <- ep2
    if ( ep$lik > ep2$lik )  res <- ep
  }, silent = TRUE)

  if ( !identical(class(RES), "try-error") ) {
    probs <- res$wij
    pj <- Rfast::colmeans(probs)
    loglik <- res$lik
    ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
    param <- cbind(pj, res$param)
    colnames(param) <- c( "probs", paste("mu", 1:2, sep = ""), "gamma" )
    rownames(param) <- paste("Cluster", 1:g, sep = " ")
  } else  pj <- loglik <- ta <- param <- probs <- iters <- NA

  runtime <- proc.time() - runtime
  list( param = param, loglik = loglik - n * log(2 * pi),
        pred = ta, probs = probs, iters = i, runtime = runtime )
}


#### Projected multivariate normal MLE (weighted version)
#### Presnell, Morrison and Littell (1998), JASA
################################
.wpn.mle <- function(u, w = NULL, tol = 1e-09) {
  ci <- cos(u)  ;   si <- sin(u)
  x <- cbind(ci, si)  ## bring the data onto the circle
  n <- dim(x)[1]

  if ( is.null(w) )   w <- rep(1, n)
  sw <- sum(w)

  ini <- Rfast::vmf.mle(x)  ## weighted vMF starting values
  mu1 <- ini$mu * ini$kappa
  f <-  - 0.5   ;   con <- sqrt(2 * pi)
  tau <- drop( x %*% mu1 )
  ptau <- pnorm(tau)
  rat <- ptau / ( exp(f * tau^2)/con + tau * ptau )
  psit <- tau + rat
  psit2 <- 2 - tau * rat - rat^2
  der <- Rfast::eachcol.apply(x, w * psit) - sw * mu1
  dera <- der[1]   ;  derb <- der[2]
  dera2 <- sum( w * (psit2 * ci^2 - 1) )
  derab <- sum( w * psit2 * ci * si )
  derb2 <- sum( w * (psit2 * si^2 - 1) )
  mu2 <- mu1 - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )

  i <- 2
  while ( sum( abs(mu2 - mu1) ) > tol ) {
    i <- i + 1
    mu1 <- mu2
    tau <- drop( x %*% mu1 )
    ptau <- pnorm(tau)
    rat <- ptau / ( exp(f * tau^2)/con + tau * ptau )
    psit <- tau + rat
    psit2 <- 2 - tau * rat - rat^2
    der <- Rfast::eachcol.apply(x, w * psit) - sw * mu1
    dera <- der[1]   ;  derb <- der[2]
    dera2 <- sum( w * (psit2 * ci^2 - 1) )
    derab <- sum( w * psit2 * ci * si )
    derb2 <- sum( w * (psit2 * si^2 - 1) )
    mu2 <- mu1 - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  }
  gam <- sum(mu2^2)
  loglik <-  - 0.5 * sw * gam + sum( w * log1p( tau * ptau * con / exp(f * tau^2) ) ) - sw * log(2 * pi)
  list(iters = i, loglik = loglik, gamma = gam, mu = mu2)
}


.mixgcpc.mle <- function(u, rads = TRUE, g = 2, tol = 1e-4, maxiters = 100) {

  fun2 <- function(wlika, rswlika, u, g, param, lika) {

    wij <- wlika / rswlika  ## weights
    pj <- Rfast::colmeans(wij)   # PANOS

    for ( j in 1:g ) {
      mod <- .wgcpc.mle2(u, w = wij[, j])
      mod1 <- .wgcpc.mle(u, w = wij[, j])
      if ( mod1$loglik > mod$loglik )  mod <- mod1
      omega <- mod$circmu
      gam <- mod$gamma
      rho <- mod$rho
      param[j, ] <- c(omega, gam, rho)
      phi <- u - omega
      a <- gam * cos(phi)
      b <- cos(phi)^2 + sin(phi)^2/rho
      lika[, j] <- 0.5 * log(rho) + log( b * sqrt(gam^2 + 1) - a * sqrt(b) ) + log( pj[j] )
    }

    wlika <- exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika / rswlika
    list(wij = wij, param = param, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  if ( !rads )  u <- u * pi / 180
  n <- length(u)

  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  param <- matrix(nrow = g, ncol = 3)

  runtime <- proc.time()
  ## Step 1

  if ( g > 1 ) {
    x <- cbind( cos(u), sin(u) )
    cl <- kmeans(x, g, nstart = 10)$cl
  } else  cl <- rep(1, n)
  wij <- tabulate(cl)

  while ( min(wij) <= 6 ) {
    g <- g - 1
    lika <- matrix(nrow = n, ncol = g)
    param <- matrix(nrow = g, ncol = 3)
    cl <- kmeans(x, g, nstart = 10)$cl
    wij <- tabulate(cl)
  }

  for ( j in 1:g ) {
    mod <- Directional::gcpc.mle2(u[cl == j], rads = TRUE)
    mod1 <- Directional::gcpc.mle(u[cl == j], rads = TRUE)
    if ( mod1$loglik > mod$loglik )  mod <- mod1
    omega <- mod$circmu
    gam <- mod$gamma
    rho <- mod$rho
    param[j, ] <- c(omega, gam, rho)
    phi <- u - omega
    a <- gam * cos(phi)
    b <- cos(phi)^2 + sin(phi)^2/rho
    lika[, j] <- 0.5 * log(rho) + log( b * sqrt(gam^2 + 1) - a * sqrt(b) )
  }

  wlika <- exp(lika)
  rswlika <- Rfast::rowsums(wlika)

  ep <- fun2(wlika, rswlika, u, g, param, lika)
  lik[1] <- ep$lik
  ep2 <- fun2(ep$wlika, ep$rswlika, u, g, ep$param, ep$lika)
  lik[2] <- ep2$lik

  i <- 2
  while ( lik[i] - lik[i - 1] > tol & i < maxiters ) {
    i <- i + 1
    ep <- ep2
    ep2 <- fun2(ep$wlika, ep$rswlika, u, g, ep$param, ep$lika)
    lik[i] <- ep2$lik
  }
  res <- ep2
  if ( ep$lik > ep2$lik )  res <- ep

  pj <- Rfast::colmeans(res$wij)
  loglik <- res$lik
  ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
  param <- cbind(pj, res$param)
  runtime <- proc.time() - runtime
  colnames(param) <- c( "probs", "omega", "gamma", "rho" )
  rownames(param) <- paste("Cluster", 1:g, sep = " ")

  list( param = param, loglik = loglik - n * log(2 * pi),
        pred = ta, probs = res$wij, iters = i, runtime = runtime )
}



.wgcpc.mle2 <- function(u, w) {
  lik <- function(para, u, w, sw) {
    omega <- para[1]
    r <- max(exp(para[2]), 1000)  ;  g <- max(exp(para[3]), 1000)
    phi <- u - omega
    a <- g * cos(phi)
    b <- cos(phi)^2 + sin(phi)^2/r
    0.5 * sw * log(r) + sum( w * log( b * sqrt(g^2 + 1) - a * sqrt(b) ) )
  }

  n <- length(u)
  sw <- sum(w)
  mod <- Directional::cipc.mle(u, TRUE)
  suppressWarnings({
    ini <- c(mod$mu, log( mod$gamma^2), rnorm(1) )
    mod <- optim(ini, lik, u = u, w = w, sw = sw, control = list(maxit = 5000))
    mod <- optim(mod$par, lik, u = u, w = w, sw = sw, control = list(maxit = 5000))
    mod <- optim(mod$par, lik, u = u, w = w, sw = sw, control = list(maxit = 5000))
  })
  circmu <- mod$par[1] %% (2 *pi)
  rho <- exp(mod$par[2])  ;  gama <- exp(mod$par[3])
  mu <- c( cos(circmu), sin(circmu) ) * gama
  list(mu = mu, circmu = circmu, gamma = gama, rho = rho, loglik = -mod$value - n * log(2 * pi) )
}


.wgcpc.mle <- function(x, w = NULL) {
  likint <- function(mu, rho, x, w, sw) {
    g2 <- sum(mu^2)
    ksi <- mu / sqrt(g2)
    sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
                       ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
    a <- as.vector(x %*% mu)
    b <- Rfast::rowsums( x %*% sinv * x )
    sw * 0.5 * log(rho) + sum( w * log( b * sqrt(g2 + 1) - a * sqrt(b) ) )
  }
  lik0 <- function(rho, x, w, sw, ma) {
    m1 <- optim(ma, likint, rho = rho, x = x, w = w, sw = sw, control = list(maxit = 1000) )
    - optim( m1$par, likint, rho = rho, x = x, w = w, sw = sw, control = list(maxit = 1000) )$value
  }
  lik <- function(param, x, w, sw) {
    mu <- param[1:2]
    rho <- param[3]
    g2 <- sum(mu^2)
    ksi <- mu / sqrt(g2)
    sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
                       ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
    a <- drop(x %*% mu)
    b <- Rfast::rowsums( x %*% sinv * x )
    sw * 0.5 * log(rho) + sum( w * log( b * sqrt(g2 + 1) - a * sqrt(b) ) )
  }

  x <- cbind( cos(x), sin(x) )
  n <- dim(x)[1]
  sw <- sum(w)

  ma <- Rfast::eachcol.apply(x, w) / sw      # weighted resultant vector as starting mu

  rho <- optimize(lik0, c(0.001, 1000), x = x, w = w, sw = sw, ma = ma, maximum = TRUE)$maximum
  mod <- optim(ma, likint, rho = rho, x = x, w = w, sw = sw, control = list(maxit = 5000) )
  suppressWarnings({
    mod <- optim( c(mod$par, rho), lik, x = x, w = w, sw = sw, control = list(maxit = 5000) )
  })
  mu <- mod$par[1:2]  ;  rho = mod$par[3]
  gama <- sqrt( sum(mu^2) )
  circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  list(mu = mu, circmu = circmu, gamma = gama, rho = rho,
       loglik = -mod$value - sw * log(2 * pi) )
}
