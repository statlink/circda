circplot <- function(u, mu  = NULL, rads = TRUE, col = 3, lwd = 2) {
  if ( !rads )  u <- u * pi / 180
  dat <- circular::circular(u)
  plot(dat, stack = TRUE)
  if ( is.null(mu) ) {
    cx <- mean( cos(u) )  ;  cy <- mean( sin(u) )
    mu <- ( atan(cy/cx) + pi * I(cx < 0) ) %% (2 * pi)
    cx <- cos(mu)  ;  cy <- sin(mu)
  } else {
    cx <- cos(mu)
    cy <- sin(mu)
  }
  arrows(x0 = 0, y0 = 0, x1 = cx, y1 = cy, length = 0.10, lwd = lwd, col = col)
}
