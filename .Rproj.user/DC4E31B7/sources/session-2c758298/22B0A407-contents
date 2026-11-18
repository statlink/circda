cv.vmnb <- function(y, ina, h = "rot", rads = TRUE, folds = NULL, nfolds = 10,
                    stratified = TRUE, ncores = 1) {

  ina <- as.numeric(ina)
  if ( is.null(folds) )  folds <- .makefolds(ina, nfolds = nfolds, stratified = stratified)
  nfolds <- length(folds)
  per <- numeric(nfolds)

  for ( i in 1:nfolds ) {
    ytrain <- y[ -folds[[ i]], ]
    ytest <- y[ folds[[ i]], ]
    inatrain <- ina[ -folds[[ i ]] ]
    inatest <- ina[ folds[[ i ]] ]
    est <- Directionalkernel::vm.nb(ytest, ytrain, inatrain, h, rads, ncores)
    per[i] <- mean( est$class == inatest )
  }

  list(per = per, perf = mean(per) )
}
