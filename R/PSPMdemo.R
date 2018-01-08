#' Demographic analysis of a structured population model
#'
#' \code{PSPMdemo} computes the population growth rate of a physiologically structured
#' population model and its sensitivities with respect to all model parameters.
#' \code{PSPMdemo} either carries out these computation for a single parameter set or
#' varies one of the parameters over a range of values specified by the user
#'
#'   output <- PSPMdemo(modelname = NULL, curvepars = NULL, parameters = NULL, options = NULL,
#'                      clean = FALSE, force = FALSE, debug = FALSE)
#'
#'
#' @param   modelname  (string, required)
#' \preformatted{}
#'               Basename of the file with model specification. The file
#'               should have extension ".h". For example, the model "Medfly"
#'               is specified in the file "Medfly.h". If the model is specified in R
#'               include the .R extension explicitly, i.e. specify the model name
#'               as "Medfly.R"
#'
#' @param   curvepars  (row vector, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector of length 5, specifying:
#' \preformatted{}
#'               \verb{curvepars[1]}: the index of the parameter to vary
#'               (in case the model is specified in R, this can be
#'               a string with the name of the parameter as specified
#'               in the variable 'DefaultParameters')
#' \preformatted{}
#'               \verb{curvepars[2]}: the initial value of the parameter
#' \preformatted{}
#'               \verb{curvepars[3]}: the step size in the parameter value
#' \preformatted{}
#'               \verb{curvepars[4]}: lower threshold, below which value of the
#'               parameter the computation stops
#' \preformatted{}
#'               \verb{curvepars[5]}: upper threshold, above which value of the
#'               parameter the computation stops
#'
#' @param   parameters (row vector, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector of length PARAMETER_NR (set in the model program
#'               file; This is the length of the variable 'DefaultParameters' if the
#'               model is specified in R), specifying the values for the model
#'               parameters to use in the computation. Vectors of other lengths,
#'               including an empty vector will be ignored.
#'
#' @param   options    (row vector of strings, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector with pairs of strings, consisting of an option name and a value (for
#'               example c("isort", "1")) or single options (i.e. c("test")).
#'               Possible option names and their values are:
#' \preformatted{}
#'               \verb{"isort", "<index>"}: Index of i-state variable to use as
#'                      ruling variable for sorting the
#'                      structured populations
#' \preformatted{}
#'               \verb{"test"}: Perform only a single integration over
#'                      the life history, reporting dynamics
#'                      of survival, R0, i-state and
#'                      interaction variables
#'
#' @param   clean      (Boolean, optional argument)
#' \preformatted{}
#'               Specify clean = TRUE as argument to remove all the result files
#'               of the model before the computation
#'
#' @param   force      (Boolean, optional argument)
#' \preformatted{}
#'               Specify force = TRUE as argument to force a rebuilding of the model
#'               before the computation
#'
#' @param   debug      (Boolean, optional argument)
#' \preformatted{}
#'               Specify debug = TRUE as argument to compile the model in verbose
#'               mode and with debugging flag set
#'
#' @return  The output is a list containing the following elements:
#' \preformatted{}
#' \verb{curvepoints}: Matrix with output for all computed points along the curve
#' \preformatted{}
#' \verb{curvedesc}: Column vector with strings, summarizing the numerical details
#'              of the computed curve (i.e., initial point, parameter values,
#'              numerical settings used).
#'
#' @examples
#' PSPMdemo("Medfly", c(2, 11, 0.1, 11, 16))
#'
#' @import utils
#' @export
PSPMdemo <- function(modelname = NULL, curvepars = NULL, parameters = NULL, options = NULL, clean = FALSE, force = FALSE, debug = FALSE) {

  Oldwd = model.Name = Rmodel = Varlist = Funlist = libfile.Basename = DefaultParameters = NULL;

  libfile.Fullname = buildSO("PSPMdemo", modelname, debug, force)
  setwd(Oldwd)

  if (!file.exists(libfile.Fullname)) stop(paste0("\nExecutable ", libfile.Basename, " not found! Computation aborted.\n"))

  if (Rmodel == 1) {
    if ((is.character(curvepars[1])) && (length((1:length(DefaultParameters))[curvepars[1]==names(DefaultParameters)]) == 1)) {
      curvepars <- c(as.integer((1:length(DefaultParameters))[curvepars[1]==names(DefaultParameters)])-1,
                 as.double(curvepars[2]), as.double(curvepars[3]), as.double(curvepars[4]), as.double(curvepars[5]))
    }
  }

  if ((length(curvepars))  && ((length(curvepars) != 5) || (!is.double(curvepars)))) stop('If specified the curve parameter argument should be a vector of length 5 with double values')
  if ((length(parameters)) && (!is.double(parameters))) stop('If specified parameter values should be a vector with double values')
  if ((length(options)) && (!is.character(options))) stop('If specified options should be an array with strings')

  if (clean) {
    outlist=list.files(pattern=paste0(model.Name, "-.*-.*.", "[cemo][srau][brt]"))
    if (debug) cat("\nCleaning :", outlist, "\n")
    for (i in outlist) file.remove(i)
  }

  dyn.load(libfile.Fullname)
  cout <- .Call("PSPMdemo", model.Name, curvepars, parameters, options, PACKAGE=paste0(model.Name, "demo"))
  dyn.unload(libfile.Fullname)

  if (Rmodel == 1) {
    rm(list = Filter( exists, Varlist ), envir = .GlobalEnv )
    rm(list = Filter( exists, Funlist ), envir = .GlobalEnv )
  }

  desc = data = NULL
  if (exists("cout")) {
    outfile.name = paste0(cout, ".out")
    if (file.exists(outfile.name)) {
      desc <- readLines(outfile.name)
      data <- as.matrix(read.table(text=desc, blank.lines.skip = TRUE, fill=TRUE))
      desc <- desc[grepl("^#", desc)]
      cat(desc, sep='\n')
    }
  }

  setwd(Oldwd)
  if (length(desc) || length(data)) {
    output = list(curvedesc = desc, curvepoints = data)
    return(output)
  }
}


