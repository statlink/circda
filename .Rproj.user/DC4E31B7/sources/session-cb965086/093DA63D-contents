circridge.cv <- function(y, x, rads = TRUE, type = "vm", lambda = NULL, nlambda = 100, tol = 1e-06, maxiters = 100,
                        folds = NULL, nfolds = 10, seed = NULL) {
  runtime <- proc.time()
  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }
  n <- dim(x)[1]
  if ( is.null(folds) )  folds <- Directional::makefolds(1:n, nfolds = nfolds, seed = seed, stratified = FALSE)
  nfolds <- length(folds)
  lamkld <- list()

  mod <- circda::circ.ridge(y = y, x = x, type = type, lambda = lambda, nlambda = nlambda, tol = tol, maxiters = maxiters)
  lambda <- mod$info[, 1]
  nlambda <- length(lambda)
  fit <- matrix(nrow = nfolds, ncol = nlambda )

  for ( k in 1:nfolds ) {
    xtest <- x[ folds[[ k ]], ]
    xtrain <- x[ -folds[[ k ]], ]
    ytest <- y[ folds[[ k ]], ]
    ytrain <- y[ -folds[[ k ]], ]
    est <- circda::circ.ridge(y = ytrain, x = xtrain, type = type, lambda = lambda, nlambda = nlambda, xnew = xtest, tol = tol, maxiters = maxiters)$est
    for ( j in 1:length(est) ) {
      est2 <- cbind( cos( est[[ j ]] ), sin( est[[ j ]] ) )
      fit[k, j] <- sum( ytest * est2 ) / dim(est2)[1]
    }
  }
  fit <- cbind( lambda, colMeans(fit, na.rm = TRUE) )
  colnames(fit) <- c("lambda", "Fit")
  lopt <- which.max(fit[, 2])
  runtime <- proc.time() - runtime

  list(runtime = runtime, fit = fit, lambda.opt = lambda[lopt], be = mod$be[[ lopt ]] )
}


