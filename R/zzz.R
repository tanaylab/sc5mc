.onLoad <- function(libname, pkgname) {
	logging::basicConfig()
}

.onUnload <- function (libpath) {
  library.dynam.unload("sc5mc", libpath)
}
