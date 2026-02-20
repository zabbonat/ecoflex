# =============================================================================
# ecoflex: Package startup
# =============================================================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "ecoflex ", utils::packageVersion("ecoflex"),
    " - Flexible Econometric Models"
  )
}
