#' Opens the PSPManalysis manual
#'
#' \code{PSPMhelp} opens the manual of the the PSPManalysis package in html format.
#'
#' The manual is created in bookdown format. A PDF version can be downloaded via the 
#' PDF icon in the menu bar.
#' 
#' @return None.
#'
#' @examples
#' \dontrun{
#' PSPMhelp()
#' }
#'
#' @importFrom rstudioapi viewer
#' @importFrom utils unzip
#' @export
PSPMhelp <- function ()
{
  oldwd <- getwd()
  tempDir <- tempdir()
  unlink(paste0(tempDir, "/manual"), recursive = TRUE)
  dir.create(paste0(tempDir, "/manual"))
  setwd(paste0(tempDir, "/manual"))
  unzip(paste0(system.file("manual", package = "PSPManalysis"), "/PSPManalysis-manual.zip"))
  setwd(oldwd)
  htmlFile <- file.path(tempDir, "manual/index.html")
  rstudioapi::viewer(htmlFile, height="maximize")
}
