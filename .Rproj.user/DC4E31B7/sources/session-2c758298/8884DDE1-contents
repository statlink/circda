cv.vmkernelda <- function(y, ina, h = "rot", rads = TRUE, folds = NULL, nfolds = 10,
                          stratified = TRUE, ncores = 1) {

  if ( !rads )  y <- y / 180 * pi
  ina <- as.numeric(ina)
  if ( is.null(folds) )  folds <- .makefolds(ina, nfolds = nfolds, stratified = stratified)
  nfolds <- length(folds)
  per <- numeric(nfolds)

  for ( i in 1:nfolds ) {
    ytrain <- y[ -folds[[ i]] ]
    ytest <- y[ folds[[ i]] ]
    inatrain <- ina[ -folds[[ i ]] ]
    inatest <- ina[ folds[[ i ]] ]
    est <- Directionalkernel::vm.kernel.da(ytest, ytrain, inatrain, h = h, rads = TRUE, ncores = ncores)
    per[i] <- mean( est$class == inatest )
  }

  list(per = per, perf = mean(per) )
}
