#'  Evolutionary dynamics for a structured population model computed following the canonical equation
#'
#' \code{PSPMevodyn} computes the dynamics of a physiologically structured population model
#' over evolutionary time for an arbitrary number of evolving parameters. The
#' evolutionary trajectory of these evolving parameters is determined by the
#' canonical equation of Adaptive Dynamics, which is solved using a simple
#' Euler integration scheme.
#'
#'   output <- PSPMevodyn(modelname = NULL, startpoint = NULL, curvepars = NULL, evopars = NULL, covars = NULL,
#'                        parameters = NULL, options = NULL, clean = FALSE, force = FALSE, debug = FALSE)
#'
#' @param  modelname  (string, required)
#' \preformatted{}
#'               Basename of the file with model specification. The file
#'               should have extension ".h". For example, the model "PNAS2002"
#'               is specified in the file "PNAS2002.h". If the model is specified in R
#'               include the .R extension explicitly, i.e. specify the model name
#'               as "PNAS2002.R"
#'
#' @param  startpoint (row vector, required)
#' \preformatted{}
#'               The initial point from which to start the simulation of the dynamics over
#'               evolutionary time, including the initial values of the evolving parameters
#'
#' @param  curvepars  (row vector of length 2, required)
#' \preformatted{}
#'               Vector of length 2 specifying:
#'
#' \preformatted{}
#'               \verb{curvepars[1]}: the maximum step size in evolutionary time during
#'                             the integration of the canonical equation
#' \preformatted{}
#'               \verb{curvepars[2]}: the maximum evolutionary time at which to stop
#'                             the integration of the canonical equation
#'
#' @param  evopars    (row vector of length n*3, required)
#' \preformatted{}
#'               Vector of length n*3 specifying:
#'
#' \preformatted{}
#'               \verb{evopars[1]}: the index of the first evolution parameter
#'                           (in case the model is specified in R, this can be
#'                            a string with the name of the parameter as specified
#'                            in the variable 'DefaultParameters')
#' \preformatted{}
#'               \verb{evopars[2]}: lower threshold, below which value of the
#'                           first evolution parameter the computation stops
#' \preformatted{}
#'               \verb{evopars[3]}: upper threshold, above which value of the
#'                           first evolution parameter the computation stops
#' \preformatted{}
#'               ......
#' \preformatted{}
#'
#'               \verb{evopars[n*3-2]}: the index of the last evolution parameter
#'                               (in case the model is specified in R, this can be
#'                                a string with the name of the parameter as specified
#'                                in the variable 'DefaultParameters')
#' \preformatted{}
#'               \verb{evopars[n*3-1]}: lower threshold, below which value of the
#'                               last evolution parameter the computation stops
#' \preformatted{}
#'               \verb{evopars[n*3]}:   upper threshold, above which value of the
#'                               last evolution parameter the computation stops
#'
#' @param  covars     (row vector or matrix, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector of length N*N or NxN matrix, where N is the number of evolving
#'               parameters. The vector or matrix elements specify the values of the
#'               covariance matrix in the selection gradients. Vectors of other lengths,
#'               including an empty vector will be ignored.
#'
#' @param  parameters (row vector, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector of length PARAMETER_NR (set in the model program
#'               file), specifying the values for the model parameters to
#'               use in the computation. Vectors of other lengths, including
#'               an empty vector will be ignored.
#'
#' @param  options    (row vector of strings, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector with pairs of strings, consisting of an option name and a value (for
#'               example c("popEVO", "1")) or single options (i.e. c("test")).
#'               Possible option names and their values are:
#'
#' \preformatted{}
#'               \verb{"popEVO", "<index>"}: Index of structured population, for
#'                                    which to compute the selection gradient or
#'                                    perform ESS or PIP continuation
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
#'               \verb{"test"}: Perform only a single integration over
#'                                    the life history, reporting dynamics
#'                                    of survival, R0, i-state and
#'                                    interaction variables
#'
#' @param  clean      (Boolean, optional argument)
#' \preformatted{}
#'               Specify clean = TRUE as argument to remove all the result files
#'               of the model before the computation
#'
#' @param  force      (Boolean, optional argument)
#' \preformatted{}
#'               Specify force = TRUE as argument to force a rebuilding of the model
#'               before the computation
#'
#' @param  debug      (Boolean, optional argument)
#' \preformatted{}
#'               Specify debug = TRUE as argument to compile the model in verbose
#'               mode and with debugging flag set
#'
#' @return  The output is a list containing the following elements:
#'
#' \preformatted{}
#'   \verb{curvepoints}: Matrix with output for all computed points along the curve
#'
#' \preformatted{}
#'   \verb{curvedesc}:   Column vector with strings, summarizing the numerical details
#'                of the computed curve (i.e., initial point, parameter values,
#'                numerical settings used)
#'
#' @examples
#' PSPMevodyn("Indet_growth", c(0.22, 0.03554, 1.0), c(0.05, 2), c(6, 0.5, 1.5), 
#'            options = c("popEVO", "0"))
#'
#' @import utils
#' @export
PSPMevodyn <- function(modelname = NULL, startpoint = NULL, curvepars = NULL, evopars = NULL, covars = NULL, parameters = NULL, options = NULL, clean = FALSE, force = FALSE, debug = FALSE) {

  Oldwd = model.Name = Rmodel = Varlist = Funlist = libfile.Basename = DefaultParameters = NULL;
  libfile.Fullname = buildSO("PSPMevodyn", modelname, debug, force)
  setwd(Oldwd)

  if (!file.exists(libfile.Fullname)) stop(paste0("\nExecutable ", libfile.Basename, " not found! Computation aborted.\n"))

  if (Rmodel == 1) {
    evopars2 <- NULL
    for (i in seq(1, length(evopars), 3)) {
      if ((is.character(evopars[i])) && (length((1:length(DefaultParameters))[evopars[i]==names(DefaultParameters)]) == 1)) {
        evopars2 <- c(evopars2, as.integer((1:length(DefaultParameters))[evopars[i]==names(DefaultParameters)])-1, as.double(evopars[i+1]), as.double(evopars[i+2]))
      }
      else
      {
        evopars2 <- c(evopars2, as.integer(evopars[i]), as.double(evopars[i+1]), as.double(evopars[i+2]))
      }
    }
    evopars <- evopars2
  }

  if ((!length(startpoint)) || (!is.double(startpoint))) stop('Starting values should be a vector with double values')
  if ((length(curvepars) != 2) || (!is.double(curvepars))) stop('Curve parameters should be a vector of length 2 with double values')
  if (((length(evopars) %% 3) != 0) || (!is.double(evopars))) stop('Evolutionary parameter values should be a vector of length n*3 with double values')
  if ((length(covars)) && (!is.double(covars))) stop('If specfied covariance values should be a vector or matrix with double values')
  if ((length(parameters)) && (!is.double(parameters))) stop('If specified parameter values should be a vector with double values')
  if ((length(options)) && (!is.character(options))) stop('If specified options should be an array with strings')

  if (clean) {
    outlist=list.files(pattern=paste0(model.Name, "-.*-.*.", "[bcemo][israu][fbrt]"))
    if (debug)
      cat("\nCleaning :", outlist, "\n")
    for (i in outlist) file.remove(i)
  }

  dyn.load(libfile.Fullname)
  # Turn the covariance matrix into a vector if necessary
  if (length(covars) > 0) covars = c(t(covars))
  cout <- .Call("PSPMevodyn", model.Name, startpoint, curvepars, evopars, covars, parameters, options, PACKAGE=paste0(model.Name, "evodyn"))
  dyn.unload(libfile.Fullname)

  if (Rmodel == 1) {
    rm(list = Filter( exists, Varlist ), envir = .GlobalEnv )
    rm(list = Filter( exists, Funlist ), envir = .GlobalEnv )
  }

  desc = data = NULL
  if (exists("cout")) {
    outfile.name = paste0(cout, ".out")
    if (file.exists(outfile.name) && (file.info(outfile.name)$size > 0)) {
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
  else return
}

