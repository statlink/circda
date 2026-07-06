dcirc <- function(u, rads = TRUE, mu, kappa, rho, type = "vm", logden = FALSE) {

  if ( type == "vm" ) {
    den <- Directional::dvm(u, mu, kappa, rads = rads, logden = logden)
  } else if ( type == "cp" ) {
    den <- Directional::dcircpurka(u, mu, kappa, rads = rads, logden = logden)
  } else if ( type == "pn" ) {
    den <- Directional::dspml(u, kappa * mu, rads = rads, logden = logden)
  } else if ( type == "gcpc" ) {
    den <- Directional::dgcpc(u, omega = mu, g = kappa, rho = rho, rads = rads, logden = logden)
  } else if ( type == "cipc" ) {
    den <- Directional::dcipc(u, omega = mu, g = kappa, rads = rads, logden = logden)
  }
  den
}
