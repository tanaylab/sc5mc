.onLoad <- function(libname, pkgname) {
	logging::basicConfig()
	tgconfig::register_params(system.file('config/metacell_params.yaml', package='sc5mc'), package='sc5mc')
}

.onUnload <- function (libpath) {
  library.dynam.unload("sc5mc", libpath)
}
