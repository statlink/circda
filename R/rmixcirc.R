rmixcirc <- function(n, rads = TRUE, probs, mu, kappa, rho, type = "vm") {

  if ( type == "vm" ) {
    mu <- cbind( cos(mu), sin(mu) )
    u <- Directional::rmixvmf(n, probs, mu, kappa)
    u$x <- ( atan(u$x[, 2]/u$x[, 1]) + pi * I(u$x[, 1] < 0) ) %% (2 * pi)
    names(u)[[ 2 ]] <- "u"
  } else if ( type == "cp" ) {
    u <- .rmixpurka(n, probs, mu, kappa)
  } else if ( type == "pn" ) {
    mu <- kappa * cbind( cos(mu), sin(mu) )
    u <- .rmixpn(n, probs, mu)
  } else if ( type == "gcpc" ) {
    u <- .rmixgcpc(n, probs, mu, kappa, rho)
  } else if ( type == "cipc" ) {
    u <- .rmixcipc(n, probs, mu, kappa)
  }
  if ( !rads )  u$u <- u$u / pi * 180
  u
}

.rmixpurka <- function(n, probs, mu, kappa) {
  p <- c( 0, cumsum(probs) )
  u <- rangen::Runif(n)
  g <- length(probs)  ## how many clusters are there
  ina <- as.numeric( cut(u, breaks = p) )  ## the cluster of each observation
  ina <- sort(ina)
  nu <- tabulate(ina)  ## frequency table of each cluster
  u <- list()
  for ( j in 1:g )  u[[ j ]] <- Directional::rcircpurka(nu[j], mu[j], kappa[j])
  u <- unlist(u)
  list(id = ina, u = u)
}


.rmixpn <- function(n, probs, mu) {
  p <- c( 0, cumsum(probs) )
  u <- rangen::Runif(n)
  g <- length(probs)  ## how many clusters are there
  ina <- as.numeric( cut(u, breaks = p) )  ## the cluster of each observation
  ina <- sort(ina)
  nu <- tabulate(ina)  ## frequency table of each cluster
  u <- list()
  for ( j in 1:g )  u[[ j ]] <- Directional::rspml(nu[j], mu[j, ])
  u <- unlist(u)
  list(id = ina, u = u)
}

.rmixgcpc <- function(n, probs, mu, kappa, rho) {
  p <- c( 0, cumsum(probs) )
  u <- rangen::Runif(n)
  g <- length(probs)  ## how many clusters are there
  ina <- as.numeric( cut(u, breaks = p) )  ## the cluster of each observation
  ina <- sort(ina)
  nu <- tabulate(ina)  ## frequency table of each cluster
  u <- list()
  for ( j in 1:g )  u[[ j ]] <- Directional::rgcpc(nu[j], omega = mu[j], g = kappa[j], rho = rho[j])
  u <- unlist(u)
  list(id = ina, u = u)
}

.rmixcipc <- function(n, probs, mu, kappa) {
  p <- c( 0, cumsum(probs) )
  u <- rangen::Runif(n)
  g <- length(probs)  ## how many clusters are there
  ina <- as.numeric( cut(u, breaks = p) )  ## the cluster of each observation
  ina <- sort(ina)
  nu <- tabulate(ina)  ## frequency table of each cluster
  u <- list()
  for ( j in 1:g )  u[[ j ]] <- Directional::rcipc(nu[j], omega = mu[j], g = kappa[j])
  u <- unlist(u)
  list(id = ina, u = u)
}
