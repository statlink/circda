circda.cv <- function(u, ina, rads = TRUE, folds = NULL, nfolds = 10, stratified = FALSE,
                      type = c("vm", "cp", "pn", "gcpc", "cipc"), seed = NULL) {

  if ( is.null(folds) )  folds <- Directional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)

  if ( !rads )  u <- u * pi/180
  per <- matrix(0, nfolds, 5)
  colnames(per) <- c("vM", "cp", "PN", "GCPC", "CIPC")
  for ( i in 1:nfolds ) {
    utrain <- u[ - folds[[ i ]] ]
    utest <- u[ folds[[ i ]] ]
    inatrain <- ina[ - folds[[ i ]] ]
    inatest <- ina[ folds[[ i ]] ]
    est <- circda::circ.da(utest, utrain, inatrain, rads = TRUE, type = type)
    per[i, ] <- Rfast::colmeans( est == inatest )
  }

  perf <- Rfast::colmeans(per)
  names(perf) <- colnames(per)
  perf
}
