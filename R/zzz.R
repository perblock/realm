#' @importFrom utils packageVersion packageDescription
.onAttach <- function(libname, pkgname) {
  if (!interactive()) return()
  base::packageStartupMessage(
    pkgname, ": version ", utils::packageVersion(pkgname)
  )
}
