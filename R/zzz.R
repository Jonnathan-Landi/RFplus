.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, version))
  packageStartupMessage("Type RFplusNews() to see new features/changes/bug fixes.")
}
