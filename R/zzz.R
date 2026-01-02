.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("This is PTMods version",
                              packageVersion("PTMods")))
}
