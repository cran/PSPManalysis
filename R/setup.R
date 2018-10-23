.onAttach <- function(libname, pkgname) {
  # Checking build tools setup using the newer R package 'pkgbuild'
  result <- try(pkgbuild::rcmd_build_tools("SHLIB"), silent = TRUE)
  if (is.list(result) && ('status' %in% names(result))) {
    if (result$status != 0) {
      if (.Platform$OS.type == "windows") {
        stop("\n\nCommand \"R CMD SHLIB\" failed. Compilation of any source file will hence be unsuccessful.\nInstall Rtools from https://cran.r-project.org/bin/windows/Rtools and allow the installer to adjust the PATH variable.\n\n")
      }
      else {
        stop("\n\nCommand \"R CMD SHLIB\" failed. Compilation of any source file will will hence be unsuccessful.\nPlease install an appropriare C compiler (gcc or clang). On Mac OS execute 'xcode-select --install'\n\n")
      }
    }
  }

  msg <- "\nWelcome to the PSPManalysis package for anlyzing physiologically structured population models\n"
  msg <- c(msg, "Explore the demos (shown by demo(package=\"PSPManalysis\") to get an overview\n")
  msg <- c(msg, "Also check out the help pages of the exported functions:\n\n")
  msg <- c(msg, "PSPMdemo, PSPMecodyn, PSPMequi, PSPMevodyn, PSPMind, PSPMclean, PSPMhelp, csbread and showpspm\n\n")
  msg <- c(msg, "For detailed information consult the manual using PSPMhelp(\"html\") or PSPMhelp(\"pdf\")\n")

  packageStartupMessage(msg)
}
