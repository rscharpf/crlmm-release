# Loading required libraries
THISPKG <- "crlmm"

.onAttach <- function(libname, pkgname) {
	packageStartupMessage("Welcome to crlmm version ", packageDescription(THISPKG, fields="Version"))
}

.onUnload <- function( libpath ){
	library.dynam.unload(THISPKG, libpath)
}
.crlmmPkgEnv <- new.env(parent=emptyenv())
