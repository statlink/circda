vm.kernel.da <- function(ynew, y, ina, h = "tay", rads = TRUE, ncores = 1) {
  if ( !rads ) {  y <- y / 180 * pi  ;  ynew <- ynew / 180 * pi  }
  ina <- as.numeric(ina)  ;  k <- max(ina)  ;  ni <- tabulate(ina)
  prior <- ni / length(ina)
  if ( is.character(h) ) {
    hrule <- h  ;  h <- rep(NA_real_, k)
    for ( j in which(ni > 0) )  h[j] <- .vm.bw(y[ina == j], hrule)
  } else if ( length(h) == 1 )  h <- rep(h, k)
  cons <- log(prior) - log(ni) - log(2 * pi) - log( besselI(h, 0, expon.scaled = TRUE) ) - h
  res <- kernel_da_cpp( rbind( cos(ynew), sin(ynew) ), rbind( cos(y), sin(y) ),
                        as.integer(ina), h, cons, as.integer(ncores) )
  list( h = h, class = res$class, score = res$score )
}
