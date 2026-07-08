rcirc <- function(n, rads = TRUE, mu, kappa, rho, type = "vm") {

  if ( type == "vm" ) {
    mu <- cbind( cos(mu), sin(mu) )
    u <- Directional::rvonmises(n, mu, kappa, rads = TRUE)
    u <- ( atan(u[, 2]/u[, 1]) + pi * I(u[, 1] < 0) ) %% (2 * pi)
  } else if ( type == "cp" ) {
    u <- Directional::rcircpurka(n, mu, kappa, rads = TRUE)
  } else if ( type == "pn" ) {
    mu <- kappa * cbind( cos(mu), sin(mu) )
    u <- Directional::rspml(n, mu, rads = TRUE)
  } else if ( type == "gcpc" ) {
    u <- Directional::rgcpc(n, omega = mu, g = kappa, rho = rho, rads = TRUE)
  } else if ( type == "cipc" ) {
    u <- Directional::rcipc(n, omega = mu, g = kappa, rads = TRUE)
  }
  if ( !rads )  u <- u / pi * 180
  u
}


