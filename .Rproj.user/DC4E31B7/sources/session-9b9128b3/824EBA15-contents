circ.mle <- function(u, rads = TRUE, type = "vm", tol = 1e-6, maxiters = 100) {

  if ( !rads )  u <- u * pi / 180
  if ( type == "vm" ) {
    res <- Rfast::vm.mle(u, tol)
  } else if ( type == "cp" ) {
    res <- Directional::purka.mle(u, tol)
  } else if ( type == "pn" ) {
    res <- Directional::spml.mle(u, rads = TRUE, tol)
  } else if ( type == "gcpc" ) {
    res <- Directional::gcpc.mle2(u, rads = TRUE)
  } else if ( type == "cipc" ) {
    res <- Directional::cipc.mle(u, rads = TRUE, tol)
  }
  res
}
