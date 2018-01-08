#'
#' Shows the model definition file of one of the example models provided with PSPManalysis
#'
#' \code{showpspm} displays the file contents of one of the physiologically structured
#' population models that is provided as an example.
#'
#'   showpspm(modelname = NULL)
#'
#' @param  modelname  (string)
#' \preformatted{}
#'               Name of the example model to be displayed.
#'
#' @examples
#' showpspm("Medfly.R")
#'
#' @export
showpspm <- function(modelname = NULL) {
  oldwd <- getwd()
  modeldir <- system.file("Models", package="PSPManalysis")
  setwd(modeldir)

  if ((!length(modelname)) || (!nchar(modelname)) || (!file.exists(modelname))) {
    cat("\nAvailable example models:\n\n")
    allmodels <- list.files(".", ".[hR]")
    cat(" ")
    cat(paste0(allmodels, sep="\n"))
    cat("\nYou have to specify one of the above model names\n\n")
  }
  else {
    file.show(modelname)
  }

  setwd(oldwd)
}
