#' Computes the individual life history of a physiologically structured population model in a given environment
#'
#' \code{PSPMind} is a utility function to compute the individual life history as
#' defined by a physiologically structured population model, given a specific
#' set of values for the environmental variables
#'
#'   output <- PSPMind(modelname = NULL, environment = NULL, parameters = NULL, options = NULL,
#'                     clean = FALSE, force = FALSE, debug = FALSE, silent = FALSE)
#'
#' @param   modelname  (string, required)
#' \preformatted{}
#'               Basename of the file with the model specification. The file
#'               should have an extension ".h". For example, the model "PNAS2002"
#'               is specified in the file "PNAS2002.h". If the model is specified 
#'               in R include the .R extension explicitly, i.e. specify the model
#'               name as "PNAS2002.R"
#'
#' @param  environment  (row vector, required)
#' \preformatted{}
#'                 Vector of length ENVIRON_DIM (set in the model program
#'                 file; This is the length of the variable 'EnvironmentState' if the
#'                 model is specified in R), specifying the value of the environmental variables
#'                 at which to calculate the individual life history.
#'                 This vector can also be extended with values of the birth rates for all structured
#'                 populations in the model, which would scale the output of the model with these
#'                 birth rates.
#'
#' @param   parameters (row vector, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector of length PARAMETER_NR (set in the model program
#'               file; This is the length of the variable 'DefaultParameters' if the
#'               model is specified in R), specifying the values for the model
#'               parameters to use in the computation. Vectors of other lengths,
#'               including an empty vector will be ignored.
#'
#' @param  options      (row vector of strings, optional, can be left equal to its default NULL)
#' \preformatted{}
#'                 Vector with a pair of strings, consisting of an option name and a value (for
#'                 example c("isort", "1")). The only possible option name and its values is:
#' \preformatted{}
#'
#'                 \verb{"isort", "<index>"}: Index of i-state variable to use as
#'                                      ruling variable for sorting the
#'                                      structured populations
#'
#' @param  clean        (Boolean, optional argument)
#' \preformatted{}
#'                 Specify clean = TRUE as argument to remove all the result files
#'                 of the model before the computation
#'
#' @param  force        (Boolean, optional argument)
#' \preformatted{}
#'                 Specify force = TRUE as argument to force a rebuilding of the model
#'                 before the computation
#'
#' @param  debug        (Boolean, optional argument)
#'                 Specify debug = TRUE as argument to compile the model in verbose
#'                 mode and with debugging flag set
#'
#' @param   silent     (Boolean, optional argument)
#' \preformatted{}
#'               Specify silent = TRUE as argument to suppress reporting of compilation
#'               commands and results on the console
#'
#' @return  The output is a structure with the population state as normally stored in the
#'   .csb output file of \code{\link{PSPMdemo}}, \code{\link{PSPMequi}}, \code{\link{PSPMecodyn}} and
#' \code{\link{PSPMevodyn}}.
#'
#' @examples
#' PSPMind("PNAS2002_5bs", c(1.30341E-05, 3.84655E-05, 4.00802E-06), 
#'         options = c("isort", "1"), clean=TRUE, force=TRUE)
#'
#' @export
PSPMind <- function(modelname = NULL, environment = NULL, parameters = NULL, options = NULL, clean = FALSE, force = FALSE, debug = FALSE, silent = FALSE) {

  Oldwd = model.Name = Rmodel = Varlist = Funlist = libfile.Basename = NULL;
  libfile.Fullname = buildSO("PSPMind", modelname, debug, force, silent)
  setwd(Oldwd)

  if (!file.exists(libfile.Fullname)) stop(paste0("\nExecutable ", libfile.Basename, " not found! Computation aborted.\n"))

  if ((!length(environment)) || (!is.double(environment))) stop('Environmental values should be a vector with double values')
  if ((length(parameters)) && (!is.double(parameters))) stop('If specified parameter values should be a vector with double values')
  if ((length(options)) && (!is.character(options))) stop('If specified options should be an array with strings')

  if (clean) {
    outlist=list.files(pattern=paste0(model.Name, "-.*-.*.", "[cemo][srau][brt]"))
    if (debug)
      cat("\nCleaning :", outlist, "\n")
    for (i in outlist) file.remove(i)
  }

  dyn.load(libfile.Fullname)
  cout <- .Call("PSPMind", model.Name, environment, parameters, options, PACKAGE=paste0(model.Name, "ind"))
  dyn.unload(libfile.Fullname)

  if (Rmodel == 1) {
    rm(list = Filter( exists, Varlist ), envir = .GlobalEnv )
    rm(list = Filter( exists, Funlist ), envir = .GlobalEnv )
  }

  state = NULL
  if (exists("cout")) {
    outfile.name = paste0(cout, ".csb")
    if (file.exists(outfile.name)) {
      state <- csbread(outfile.name, 1)
    }
  }

  setwd(Oldwd)
  return (state)
}
