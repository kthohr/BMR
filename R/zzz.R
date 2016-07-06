
.onLoad <- function(libname, pkgname) {
    suppressWarnings(Rcpp::loadRcppModules())
}

#loadModule("bvarm_module",TRUE)
