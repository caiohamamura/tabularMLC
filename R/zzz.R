.onUnload <- function (libpath) {
  library.dynam.unload("tabularMLC", libpath)
  invisible()
}