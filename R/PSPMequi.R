#' Bifurcation analysis of a structured population model
#'
#' \code{PSPMequi} computes bifurcation curves for a physiologically structured
#' population model as a function of one or two parameters and detects bifurcation
#' points along these curves.
#' When computing equilibrium curves of a physiologically structured population model
#' as a function of one parameter, \code{PSPMequi} can detect transcritical
#' bifurcation points (branching points, BP) of both the structured populations as
#' well as environment variables (BPE), limit points (LP) in the equilibrium curve
#' and evolutionary stationary points (ESS). The location of these bifurcation points
#' can subsequently be computed as a function of second parameter.
#' In addition \code{PSPMequi} can compute the pairwise invasion plot (PIP) as a
#' function of the resident and a mutant value for one evolving parameter.
#'
#'   output <- PSPMequi(modelname = NULL, biftype = NULL, startpoint = NULL, stepsize = NULL, parbnds = NULL,
#'                      parameters = NULL, options = NULL, clean = FALSE, force  = FALSE, debug  = FALSE)
#'
#' @param   modelname  (string, required)
#' \preformatted{}
#'               Basename of the file with model specification. The file
#'               should have extension ".h". For example, the model "PNAS2002"
#'               is specified in the file "PNAS2002.h". If the model is specified in R
#'               include the .R extension explicitly, i.e. specify the model name
#'               as "PNAS2002.R"
#'
#' @param   biftype    (string, required)
#' \preformatted{}
#'               Type of bifurcation to compute: "BP", "BPE", "EQ", "LP", "ESS" or "PIP"
#'
#' @param   startpoint (row vector, required)
#' \preformatted{}
#'               The initial point from which to start the continuation of
#'               the curve
#'
#' @param   stepsize   (double value, required)
#' \preformatted{}
#'               Value of the step size in the first bifurcation parameter
#'
#' @param   parbnds    (row vector of length n*3, required)
#' \preformatted{}
#'               Vector of length n*3 with n=1 for EQ continuation, n=2 for
#'               "BP", "BPE", "LP" and "PIP" continuation and n>=2 for "ESS" continuation.
#'               Each triples specifies:
#'
#' \preformatted{}
#'               \verb{parbnds[1]}: the index of the first bifurcation parameter
#'                           (in case the model is specified in R, this can be
#'                           a string with the name of the parameter as specified
#'                           in the variable 'DefaultParameters')
#' \preformatted{}
#'               \verb{parbnds[2]}: lower threshold, below which value of the
#'                           first bifurcation parameter the computation stops
#' \preformatted{}
#'               \verb{parbnds[3]}: upper threshold, above which value of the
#'                           first bifurcation parameter the computation stops
#' \preformatted{}
#'               In case of two-parameter bifurcations:
#' \preformatted{}
#'               \verb{parbnds[4]}: the index of the second bifurcation parameter
#'                           (in case the model is specified in R, this can be
#'                            a string with the name of the parameter as specified
#'                            in the variable 'DefaultParameters')
#' \preformatted{}
#'               \verb{parbnds[5]}: lower threshold, below which value of the
#'                           second bifurcation parameter the computation stops
#' \preformatted{}
#'               \verb{parbnds[6]}: upper threshold, above which value of the
#'                           second bifurcation parameter the computation stops
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
#'               example c("popBP", "1")) or single options (i.e. c("test")).
#'               Possible option names and their values are:
#' \preformatted{}
#'               \verb{"envBP", "<index>"}: Index of environment variable, of which
#'                                    to continue the transcritical bifurcation
#' \preformatted{}
#'               \verb{"popBP", "<index>"}: Index of structured population, of which
#'                                    to continue the transcritical bifurcation
#' \preformatted{}
#'               \verb{"popEVO", "<index>"}: Index of structured population, for
#'                                    which to compute the selection gradient or
#'                                    perform ESS or PIP continuation
#' \preformatted{}
#'               \verb{"parEVO", "<index>"}: Index of parameter, for which to compute
#'                                    which to compute the selection gradient
#' \preformatted{}
#'               \verb{"envZE", "<index>"}: Index of environment variable in
#'                                    trivial equilibrium (can be used
#'                                    multiple times)
#' \preformatted{}
#'               \verb{"popZE", "<index>"}: Index of structured population in
#'                                    trivial equilibrium (can be used
#'                                    multiple times)
#' \preformatted{}
#'               \verb{"isort", "<index>"}: Index of i-state variable to use as
#'                                    ruling variable for sorting the
#'                                    structured populations
#' \preformatted{}
#'               \verb{"noBP"}: Do not check for branching points while
#'                                    computing equilibrium curves
#' \preformatted{}
#'               \verb{"noLP"}: Do not check for limit points while
#'                                    computing equilibrium curves
#' \preformatted{}
#'               \verb{"single"}: Only compute the first point of the
#'                                    solution curve, do not continue the curve
#' \preformatted{}
#'               \verb{"test"}: Perform only a single integration over
#'                                    the life history, reporting dynamics
#'                                    of survival, R0, i-state and
#'                                    interaction variables
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
#'   \verb{curvepoints}: Matrix with output for all computed points along the curve
#' \preformatted{}
#'   \verb{curvedesc}:   Column vector with strings, summarizing the numerical details
#'                of the computed curve (i.e., initial point, parameter values,
#'                numerical settings used)
#' \preformatted{}
#'   \verb{bifpoints}:   Matrix with the located bifurcation points along the curve
#' \preformatted{}
#'   \verb{biftypes}:    Column vector of strings, containing a description of the
#'                type of bifurcation for each of the located bifurcation points
#'
#' @examples
#' PSPMequi("Indet_growth", "EQ", c(1, 0.22, 0), -0.1, c(6, 0.7, 1.5),
#'          options = c("popEVO", "0"))
#'
#' @import utils
#' @export
PSPMequi <- function(modelname = NULL, biftype = NULL, startpoint = NULL, stepsize = NULL, parbnds = NULL, parameters = NULL,
                     options = NULL, clean = FALSE, force = FALSE, debug = FALSE) {

  Oldwd = model.Name = Rmodel = Varlist = Funlist = libfile.Basename = DefaultParameters = NULL;
  libfile.Fullname = buildSO("PSPMequi", modelname, debug, force)
  setwd(Oldwd)

  if (!file.exists(libfile.Fullname)) stop(paste0("\nExecutable ", libfile.Basename, " not found! Computation aborted.\n"))

  if (Rmodel == 1) {
    parbnds2 <- NULL
    for (i in seq(1, length(parbnds), 3)) {
      if ((is.character(parbnds[i])) && (length((1:length(DefaultParameters))[parbnds[i]==names(DefaultParameters)]) == 1)) {
        parbnds2 <- c(parbnds2, as.integer((1:length(DefaultParameters))[parbnds[i]==names(DefaultParameters)])-1, as.double(parbnds[i+1]), as.double(parbnds[i+2]))
      }
      else
      {
        parbnds2 <- c(parbnds2, as.integer(parbnds[i]), as.double(parbnds[i+1]), as.double(parbnds[i+2]))
      }
    }
    parbnds <- parbnds2
  }

  if ((!length(startpoint)) || (!is.double(startpoint))) stop('Starting values should be a vector with double values')
  if ((length(stepsize) != 1) || (!is.double(stepsize))) stop('Step size argument should be a single double value')
  if (((length(parbnds) %% 3) != 0) || (!is.double(parbnds))) stop('Parameter bounds values should be a vector of length n*3 with double values')
  if ((length(parameters)) && (!is.double(parameters))) stop('If specified parameter values should be a vector with double values')
  if ((length(options)) && (!is.character(options))) stop('If specified options should be an array with strings')

  if (clean) {
    outlist=list.files(pattern=paste0(model.Name, "-.*-.*.", "[bcemo][israu][fbrt]"))
    if (debug) cat("\nCleaning :", outlist, "\n")
    for (i in outlist) file.remove(i)
  }

  dyn.load(libfile.Fullname)
  cout <- .Call("PSPMequi", model.Name, biftype, startpoint, stepsize, parbnds, parameters, options, PACKAGE=paste0(model.Name, "equi"))
  dyn.unload(libfile.Fullname)

  if (Rmodel == 1) {
    rm(list = Filter( exists, Varlist ), envir = .GlobalEnv )
    rm(list = Filter( exists, Funlist ), envir = .GlobalEnv )
  }

  desc = data = bifpoints = biftypes = NULL
  if (exists("cout")) {
    outfile.name = paste0(cout, ".out")
    if (file.exists(outfile.name) && (file.info(outfile.name)$size > 0)) {
      desc <- readLines(outfile.name)
      data <- as.matrix(read.table(text=desc, blank.lines.skip = TRUE, fill=TRUE))
      desc <- desc[grepl("^#", desc)]
      cat(desc, sep='\n')
    }

    biffile.name = paste0(cout, ".bif")
    if (file.exists(biffile.name) && (file.info(biffile.name)$size > 0)) {
      bifinput <- readLines(biffile.name)
      bifpoints <- as.matrix(read.table(text=bifinput, blank.lines.skip = TRUE, comment.char='*', fill=TRUE))
      biftypes = gsub("^.*\\*\\*\\*\\*\\s+|\\s+\\*\\*\\*\\*.*$", "", bifinput)
    }
  }

  setwd(Oldwd)
  if (length(desc) || length(data) || length(bifpoints) || length(biftypes)) {
    if (length(bifpoints) || length(biftypes)) {
      output = list(curvedesc = desc, curvepoints = data, bifpoints = bifpoints, biftypes = biftypes)
    }
    else output = list(curvedesc = desc, curvepoints = data)
    return(output)
  }
}
