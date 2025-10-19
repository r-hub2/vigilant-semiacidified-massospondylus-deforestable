

#' @importMethodsFrom terra rast
NULL

#' @importFrom stats cov
#' @importFrom StableEstim KoutParametersEstim ComplexCF
NULL

#' @useDynLib deforestable, .registration = TRUE
NULL

# Cleans up after unloading
.onUnload <- function (libpath) {
  library.dynam.unload("deforestable", libpath)
}

#' @importFrom Rcpp evalCpp
NULL
