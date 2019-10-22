#' Opens the PSPManalysis manual
#'
#' \code{PSPMhelp} opens the manual of the the PSPManalysis package in either html or pdf format.
#'
#' @param type String (optional). Either "pdf" or "html". Default if "html"
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' PSPMhelp()
#'
#' PSPMhelp("pdf")
#' }
#'
#' @export
PSPMhelp <- function (type = c("html", "pdf"))
{
  if (!missing(type) && ((type == "pdf") || (type == "PDF"))) vignette("PSPManalysis-pdf")
  else vignette("PSPManalysis")
}

