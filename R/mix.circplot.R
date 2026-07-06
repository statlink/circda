mix.circplot <- function(u, ina, mu, rads = TRUE, lwd = 2) {
  if ( !rads )  u <- u * pi / 180
  dat <- circular::circular(u)
  cx <- cos(mu)  ;  cy <- sin(mu)
  k <- length(mu)
  plot(dat[ina == 1], stack = TRUE)
  arrows(x0 = 0, y0 = 0, x1 = cx[1], y1 = cy[1], length = 0.10, lwd = lwd, col = 1)
  for ( i in 2:k ) {
    points(dat[ina == i], stack = TRUE, col = i)
    arrows(x0 = 0, y0 = 0, x1 = cx[i], y1 = cy[i], length = 0.10, lwd = lwd, col = i)
  }

}
