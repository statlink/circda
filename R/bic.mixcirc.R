bic.mixcirc <- function(u, rads = TRUE, type = "vm", G = 5, tol = 1e-4, maxiters = 500) {

  if ( !rads )  u <- u * pi / 180

  if ( type == "vm" ) {
    x <- cbind( cos(u), sin(u) )
    res <- Directional::bic.mixvmf(x, G, n.start = 10, tol = tol, maxiters = maxiters)
  } else if ( type == "cp" ) {
    res <- .bic.mixpurka(u,  G, tol, maxiters)
  } else if ( type == "pn" ) {
    res <- .bic.mixpn(u, G, tol, maxiters)
  } else if ( type == "gcpc" ) {
    res <- .bic.mixgcpc(u,  G, tol, maxiters)
  } else if ( type == "cipc" ) {
    x <- cbind( cos(u), sin(u) )
    res <- Directional::bic.mixspcauchy(x, G = G, n.start = 10, tol = tol, maxiters = maxiters)
  }

  res
}


.bic.mixpurka <- function(u, G = 5, tol = 1e-4, maxiters = 100) {
  runtime <- proc.time()
  logn <- log( length(u) )  ## sample size of the data
  bic <- icl <- 1:G
  mod <- Directional::purka.mle(u)
  bic[1] <- icl[1] <-  - 2 * mod$loglik + 2 * logn  ## BIC assuming one cluster

  for ( vim in 2:G ) {
    a <- .mixpurka.mle(u, rads = TRUE, vim, tol = tol, maxiters = maxiters)
    d <- dim(a$param)[1]
    tab <- table(a$pred)
    nm <- min( tab )
    if ( d == vim  &  d == length(tab)  &  nm >= 6 ) {
      bic[vim] <-  - 2 * a$loglik + ( d - 1 + 2 * d ) * logn
      icl[vim] <- bic[vim] - 2 * sum( a$probs * log(a$probs), na.rm = TRUE )
    } else  {
      bic[vim] <- bic[vim - 1]
      icl[vim] <- icl[vim - 1]
    }
  }  ## BIC for a range of different clusters

  runtime <- proc.time() - runtime
  names(bic) <- names(icl) <- paste("g=", 1:G, sep = "")
  ina <- rep(1, G)
  char <- rep(16, G)
  ina[ which.min(icl) ] <- 3  ## chosen number of clusters will appear with red on the plot
  char[ which.min(icl) ] <- 17
  plot(1:G, icl, col = ina, xlab = "Number of components", ylab = "ICL values", cex.lab = 1.3, cex.axis = 1.3)
  abline(v = 1:G, lty = 2, col = "lightgrey")
  abline(h = seq(min(icl, na.rm = FALSE), max(icl, na.rm = FALSE), length = 10), lty = 2, col = "lightgrey" )
  lines(1:G, icl, lwd = 2)
  abline(v = which.min(icl), col = 3, lwd = 2)
  points(1:G, icl, pch = char, col = ina)
  list(bic = bic, icl = icl, runtime = runtime)
}


.bic.mixpn <- function(u, G = 5, tol = 1e-4, maxiters = 100) {
  runtime <- proc.time()
  logn <- log( length(u) )  ## sample size of the data
  bic <- icl <- 1:G
  mod <- Directional::spml.mle(u, rads = TRUE)
  bic[1] <- icl[1] <-  - 2 * mod$loglik + 2 * logn  ## BIC assuming one cluster

  for ( vim in 2:G ) {
    a <- .mixpn.mle(u, rads = TRUE, vim, tol = tol)  ## model based clustering for some possible clusters
    if ( sum( is.na(a$param) ) == 0 ) {
      d <- dim(a$param)[1]
      tab <- table(a$pred)
      nm <- min(tab)
      if ( d == vim  &  d == length(tab)  &  nm >= 6 ) {
        bic[vim] <-  - 2 * a$loglik + ( d - 1 + 2 * d ) * logn
        icl[vim] <- bic[vim] - 2 * sum( a$probs * log(a$probs), na.rm = TRUE )
      } else  {
        bic[vim] <- bic[vim - 1]
        icl[vim] <- icl[vim - 1]
      }
    } else {
      bic[vim] <- bic[vim - 1]
      icl[vim] <- icl[vim - 1]
    }
  }  ## BIC for a range of different clusters

  runtime <- proc.time() - runtime
  names(bic) <- names(icl) <- paste("g=", 1:G, sep = "")
  ina <- rep(1, G)
  char <- rep(16, G)
  ina[ which.min(icl) ] <- 3  ## chosen number of clusters will appear with red on the plot
  char[ which.min(icl) ] <- 17
  plot(1:G, icl, col = ina, xlab = "Number of components", ylab = "ICL values", cex.lab = 1.3, cex.axis = 1.3)
  abline(v = 1:G, lty = 2, col = "lightgrey")
  abline(h = seq(min(icl, na.rm = FALSE), max(icl, na.rm = FALSE), length = 10), lty = 2, col = "lightgrey" )
  lines(1:G, icl, lwd = 2)
  abline(v = which.min(icl), col = 3, lwd = 2)
  points(1:G, icl, pch = char, col = ina)
  list(bic = bic, icl = icl, runtime = runtime)
}



.bic.mixgcpc <- function(u, G = 5, tol = 1e-4, maxiters = 100) {
  runtime <- proc.time()
  logn <- log( length(u) )
  bic <- icl <- 1:G
  mod <- Directional::gcpc.mle2(u, rads = TRUE)
  bic[1] <- icl[1] <-  - 2 * mod$loglik + 3 * logn  ## BIC assuming one cluster

  for ( vim in 2:G ) {
    a <- .mixgcpc.mle(u, rads = TRUE, vim, tol = tol, maxiters = 100)
    d <- dim(a$param)[1]
    tab <- table(a$pred)
    nm <- min( tab )
    if ( d == vim  &  d == length(tab)  &  nm >= 6 ) {
      bic[vim] <-  - 2 * a$loglik + ( d - 1 + 3 * d ) * logn
      icl[vim] <- bic[vim] - 2 * sum( a$probs * log(a$probs), na.rm = TRUE )
    } else  {
      bic[vim] <- bic[vim - 1]
      icl[vim] <- icl[vim - 1]
    }
  }  ## BIC for a range of different clusters

  runtime <- proc.time() - runtime
  names(bic) <- names(icl) <- paste("g=", 1:G, sep = "")
  ina <- rep(1, G)
  char <- rep(16, G)
  ina[ which.min(icl) ] <- 3  ## chosen number of clusters will appear with red on the plot
  char[ which.min(icl) ] <- 17
  plot(1:G, icl, col = ina, xlab = "Number of components", ylab = "ICL values", cex.lab = 1.3, cex.axis = 1.3)
  abline(v = 1:G, lty = 2, col = "lightgrey")
  abline(h = seq(min(icl, na.rm = FALSE), max(icl, na.rm = FALSE), length = 10), lty = 2, col = "lightgrey" )
  lines(1:G, icl, lwd = 2)
  abline(v = which.min(icl), col = 3, lwd = 2)
  points(1:G, icl, pch = char, col = ina)
  list(bic = bic, icl = icl, runtime = runtime)
}
