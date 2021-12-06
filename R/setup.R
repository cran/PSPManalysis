.onAttach <- function(libname, pkgname) {
  msg <- paste0("\nWelcome to the PSPManalysis package (version ", packageVersion("PSPManalysis"), ") for anlyzing\nphysiologically structured population models\n\n")
  msg <- c(msg, "Explore the demos (shown by demo(package=\"PSPManalysis\") to get an overview\n")
  msg <- c(msg, "Also check out the help pages of the exported functions:\n\n")
  msg <- c(msg, "PSPMdemo, PSPMecodyn, PSPMequi, PSPMevodyn, PSPMind, PSPMclean, PSPMhelp, csbread and showpspm\n\n")
  msg <- c(msg, "For detailed information consult the manual using PSPMhelp()\n")

  packageStartupMessage(msg)
}
