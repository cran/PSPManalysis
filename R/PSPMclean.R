#' Deletes on request all files produced by the PSPManalysis package.
#'
#' \code{PSPMclean} deletes all PSPManalysis result files (default) and/or
#' all executables (hit 'F') in the current directory.
#'
#' @param str Character (optional). Only valid argument is 'F'. If not or wrongly
#' specified the user will be asked whether to do a full clean up or whether to
#' quit the clean up.
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' PSPMclean()
#'
#' PSPMclean("F")
#' }
#'
#' @export
PSPMclean <- function(str = NULL) {

  if ((!length(str)) || ((str != "f") && (str != "F"))) {
    cat("\nDelete all PSPManalysis result files (default) and/or all executables (hit 'F') in the current directory?\n\n")
    str = readline("Press 'Q' to abort, 'F' for a full clean up, any other key to only clean the result files: ")
  }

  if ((str == "q") || (str == "Q")) {
    cat("\nDelete operation aborted\n\n")
  }
  else
  {
    if ((str == "f") || (str == "F")) {
      cat("\nDeleting all executables in current and temporary directory")
      for (i in dir(".", "*demo\\.so")) file.remove(paste0("./", i))
      for (i in dir(".", "*ind\\.so")) file.remove(paste0("./", i))
      for (i in dir(".", "*equi\\.so")) file.remove(paste0("./", i))
      for (i in dir(".", "*ecodyn\\.so")) file.remove(paste0("./", i))
      for (i in dir(".", "*evodyn\\.so")) file.remove(paste0("./", i))

      for (i in dir(".", "*demo\\.dll")) file.remove(paste0("./", i))
      for (i in dir(".", "*ind\\.dll")) file.remove(paste0("./", i))
      for (i in dir(".", "*equi\\.dll")) file.remove(paste0("./", i))
      for (i in dir(".", "*ecodyn\\.dll")) file.remove(paste0("./", i))
      for (i in dir(".", "*evodyn\\.dll")) file.remove(paste0("./", i))

      for (i in dir(tempdir(), "*demo\\.so")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*ind\\.so")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*equi\\.so")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*ecodyn\\.so")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*evodyn\\.so")) file.remove(paste0(tempdir(), "/", i))

      for (i in dir(tempdir(), "*demo\\.dll")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*ind\\.dll")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*equi\\.dll")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*ecodyn\\.dll")) file.remove(paste0(tempdir(), "/", i))
      for (i in dir(tempdir(), "*evodyn\\.dll")) file.remove(paste0(tempdir(), "/", i))
    }
    cat("\nDeleting all PSPManalysis result files in current directory\n\n")
    for (i in c("PGR", "EQ", "BP", "BPE", "LP", "ESS", "PIP", "ECODYN", "EVODYN", "IND")) {
      allres = dir(".", paste0("*-", i, "-[0-9][0-9][0-9][0-9]\\.bif"))
      for (j in allres) file.remove(j)
      allres = dir(".", paste0("*-", i, "-[0-9][0-9][0-9][0-9]\\.csb"))
      for (j in allres) file.remove(j)
      allres = dir(".", paste0("*-", i, "-[0-9][0-9][0-9][0-9]\\.err"))
      for (j in allres) file.remove(j)
      allres = dir(".", paste0("*-", i, "-[0-9][0-9][0-9][0-9]\\.out"))
      for (j in allres) file.remove(j)
    }
  }
}
