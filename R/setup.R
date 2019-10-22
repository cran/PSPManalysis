.onAttach <- function(libname, pkgname) {
  msg <- "\nWelcome to the PSPManalysis package for anlyzing physiologically structured population models\n"
  msg <- c(msg, "Explore the demos (shown by demo(package=\"PSPManalysis\") to get an overview\n")
  msg <- c(msg, "Also check out the help pages of the exported functions:\n\n")
  msg <- c(msg, "PSPMdemo, PSPMecodyn, PSPMequi, PSPMevodyn, PSPMind, PSPMclean, PSPMhelp, csbread and showpspm\n\n")
  msg <- c(msg, "For detailed information consult the manual using PSPMhelp(\"html\") or PSPMhelp(\"pdf\")\n")

  packageStartupMessage(msg)
}
