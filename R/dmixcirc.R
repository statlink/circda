dmixcirc <- function(u, rads = TRUE, probs, mu, kappa, rho, type = "vm", logden = FALSE) {

  g <- length(probs)
  den <- matrix(nrow = length(u), ncol = g)

  if ( type == "vm" ) {
    for ( j in 1:g )  den[, j] <- Directional::dvm(u, mu[j], kappa[j], rads = rads)
  } else if ( type == "cp" ) {
    for ( j in 1:g )  den[, j] <- Directional::dcircpurka(u, mu[j], kappa[j], rads = rads)
  } else if ( type == "pn" ) {
    for ( j in 1:g )  den[, j] <- Directional::dspml(u, kappa[j] * mu[j], rads = rads)
  } else if ( type == "gcpc" ) {
    for ( j in 1:g )  den[, j] <- Directional::dgcpc(u, omega = mu[j], g = kappa[j],
                                                     rho = rho[j], rads = rads)
  } else if ( type == "cipc" ) {
    for ( j in 1:g )  den[, j] <- Directional::dcipc(u, omega = mu[j], g = kappa[j], rads = rads)
  }

  den <- den %*% probs
  if ( logden )  den <- log(den)
  den
}
