#'
#' Ecological dynamics of a structured population model computed using the Escalator Boxcar Train
#'
#' \code{PSPMecodyn} computes the dynamics of a physiologically structured population model
#' starting from an environmental and population state that is computed with \code{\link{PSPMequi}}.
#' If starting from an arbitrary state is required, the list specifying the initial state should
#' have the same layout as produced by \code{\link{PSPMequi}}.
#'
#'   output <- PSPMecodyn(modelname = NULL, startstate = NULL, timepars = NULL, bifpars = NULL,
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
#' @param  startstate (list, required)
#' \preformatted{}
#'               The initial environmental and population state from which to start the
#'               simulation of the dynamics. This list should have the identical layout
#'               as a list returned by the function csbread().
#'               As a minimum, the list should contain a vecgtor 'Environment' specifying
#'               the initial values of the environmental variables, and a matrix 'Pop00'
#'               (assuming there is only a single population in the model), which specifies
#'               on each row the number and individual state variables of a cohort of
#'               while the different rows specify all the cohorts in the population.
#'
#' @param  timepars   (row vector of length 4, required)
#' \preformatted{}
#'               Vector of length 4 specifying the settings for the time integration:
#' \preformatted{}
#'              \verb{timepars[1]}: Cohort cycle time, i.e. time interval between starts of new
#'              boundary cohorts
#' \preformatted{}
#'              \verb{timepars[2]}: Output time interval, i.e. time interval between data output
#'              to .out file
#' \preformatted{}
#'              \verb{timepars[3]}: State output interval, i.e. time interval between complete
#'              state output to .csb file
#' \preformatted{}
#'              \verb{timepars[4]}: Maximum integration time, i.e. maximum time value until which
#'              to continue the integration
#'
#' @param  bifpars    (row vector of length 6, optional)
#' \preformatted{}
#'               Vector of length 6 specifying the settings for the bifurcation settings. If
#'               not specified a normal time integration is carried out.
#' \preformatted{}
#'               \verb{bifpars[1]}: Index of the bifurcation parameter
#' \preformatted{}
#'               \verb{bifpars[2]}: Starting value of the bifurcation parameter
#' \preformatted{}
#'               \verb{bifpars[3]}: Step size in the bifurcation parameter
#' \preformatted{}
#'               \verb{bifpars[4]}: Final value of the bifurcation parameter
#' \preformatted{}
#'               \verb{bifpars[5]}: Period of producing data output during each bifurcation interval
#' \preformatted{}
#'               \verb{bifpars[6]}: Period of producing state output during each bifurcation interval
#'
#' @param   parameters (row vector, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector of length PARAMETER_NR (set in the model program
#'               file; This is the length of the variable 'DefaultParameters' if the
#'               model is specified in R), specifying the values for the model
#'               parameters to use in the computation. Vectors of other lengths,
#'               including an empty vector will be ignored.
#'
#' @param  options    (row vector of strings, optional, can be left equal to its default NULL)
#' \preformatted{}
#'               Vector with pairs of strings, consisting of an option name and a value (for
#'               example c("info", "1")).
#'               Possible option names and their values are:
#' \preformatted{}
#'               \verb{"info", "<index>"}:  Level of performance information on the DOPRI5
#'                       integrator written to .err file (1, 2, 3 or 4)
#' \preformatted{}
#'               \verb{"report", "<index>"}:  Interval between reporting of data output to
#'                       console ( > 0)
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
#' \preformatted{}
#'              \verb{curvepoints}: Matrix with output for all computed points along the curve
#' \preformatted{}
#'              \verb{curvedesc}:   Column vector with strings, summarizing the numerical details
#'              of the computed curve (i.e., initial point, parameter values,
#'              numerical settings used)
#' @examples
#' initstate <- list(Environment = c(1.561276e-04, 1.270327e-04, 4.008016e-06), 
#'                   Pop00 = matrix(c(0.001, 0, 7.0, 1.0E-5, 300, 111), ncol = 3, byrow = TRUE))
#' PSPMecodyn("PNAS2002", initstate, c(1, 1, 10, 100))
#'
#' @import utils
#' @export
PSPMecodyn <- function(modelname = NULL, startstate = NULL, timepars = NULL, bifpars = NULL, parameters = NULL, options = NULL, clean = FALSE, force = FALSE, debug = FALSE) {

  if ((!length(modelname)) || (!nchar(modelname))) stop("You have to specify a model name")

  pkgmodeldir <- system.file("Models", package="PSPManalysis")
  modelfound <- 0
  for (modeldir in c(".", pkgmodeldir))
    {
      if (modeldir == pkgmodeldir)
        {
          str <- readline(paste0("\nModel file not found in current directory\nContinue searching for model file in package directory? [y(es)/q(uit)] : "))
          if ((str != '') && (str != 'y') && (str != 'Y')) break
        }
      if ((regexpr("\\.h", modelname) == (nchar(modelname)-1)) && file.exists(paste0(modeldir, "/", modelname)))
        {
          hfile.fullname <- paste0(modeldir, "/", modelname)
          rmodel <- 0
          modelfound <- 1
          break
        }
      if ((regexpr("\\.R", modelname) == (nchar(modelname)-1)) && file.exists(paste0(modeldir, "/", modelname)))
        {
          hfile.fullname <- paste0(modeldir, "/", modelname)
          rmodel <- 1
          modelfound <- 1
          break
        }
      if (file.exists(paste0(modeldir, "/", modelname, ".h")))
        {
          hfile.fullname <- paste0(modeldir, "/", modelname, ".h")
          rmodel <- 0
          modelfound <- 1
          break
        }
      if (file.exists(paste0(modeldir, "/", modelname, ".R")))
        {
          hfile.fullname <- paste0(modeldir, "/", modelname, ".R")
          rmodel <- 1
          modelfound <- 1
          break
        }
    }

  if (modelfound != 1) {
    if (regexpr("\\.[hR]", modelname) == (nchar(modelname)-1))
      stop(paste('No model source file named', modelname, 'can be found', sep=' '))
    else
      stop(paste('No model source file named', paste(modelname, ".h", sep=''), 'or', paste(modelname, ".R", sep=''), 'can be found', sep=' '))
  }

  if (rmodel == 1) {
    varlist <- c("PSPMdimensions", "NumericalOptions", "EnvironmentState", "DefaultParameters")
    funlist <- c("StateAtBirth", "LifeHistoryRates", "LifeStageEndings", "DiscreteChanges", "EnvEqui")
    rm(list = Filter(exists, varlist), envir = .GlobalEnv)
    rm(list = Filter(exists, funlist), envir = .GlobalEnv)
    source(normalizePath(hfile.fullname))
    pspmdims <- get("PSPMdimensions", envir = .GlobalEnv)
    PopulationNr <- pspmdims["PopulationNr"];
    IStateDimension <- pspmdims["IStateDimension"];
    LifeHistoryStages <- pspmdims["LifeHistoryStages"];
    ImpactDimension <- pspmdims["ImpactDimension"];

    EnvironmentDim <- length(get("EnvironmentState", envir = .GlobalEnv))
    ParameterNr <- length(get("DefaultParameters", envir = .GlobalEnv))
  }

  if ((!length(startstate)) || (!is.list(startstate))) stop('Starting state should be a list with the same layout as returned by csbread()')
  if ((length(timepars) != 4) || (!is.double(timepars))) stop('Time settings should be a vector of length 4 with double values')
  if ((length(bifpars)) && ((!is.double(bifpars)) || (length(bifpars) != 6))) stop('If specfied bifurcation settings should be a vector of length 6 with double values')
  if ((length(parameters)) && (!is.double(parameters))) stop('If specified parameter values should be a vector with double values')
  if ((length(options)) && (!is.character(options))) stop('If specified options should be an array with strings')

  hfile.dirname <- dirname(normalizePath(hfile.fullname))
  hfile.basename <- basename(normalizePath(hfile.fullname))
  hfile.abspath <- normalizePath(paste0(hfile.dirname, '/', hfile.basename))
  if (rmodel == 1) model.Name <- sub("\\.R", "", hfile.basename) else model.Name <- sub("\\.h", "", hfile.basename)

  oldwd <- getwd()
  PSPMsrcdir.name <- normalizePath(system.file("C", package="PSPManalysis"))
  tmpdir = normalizePath(tempdir())

  if (clean) {
    outlist <- list.files(pattern=paste0(model.Name, "-.*-.*.", "[bcemo][israu][fbrt]"))
    if (debug) cat("\nCleaning :", outlist, "\n")
    for (i in outlist) file.remove(i)
  }

  libfile.basename <- paste0(model.Name, "ecodyn", .Platform$dynlib.ext)
  libfile.fullname <- paste0(tmpdir, '/', libfile.basename)

  if (is.loaded(libfile.fullname)) dyn.unload(libfile.fullname)

  if (file.exists(libfile.fullname)) {
    # The comma in front of mtime is needed to get output as class POSIXct
    build.flag <- (difftime(file.info(hfile.fullname)[,'mtime'], file.info(libfile.fullname)[,'mtime'], units = "secs") >= 0)
  }
  else build.flag <- TRUE

  if (build.flag || force) {
    setwd(tmpdir)

    srclist <- c("PSPMecodyn.c", "escbox/ebtcohrt.c", "escbox/ebtdopri5.c", "escbox/ebtutils.c")
    inclist <- c("escbox/ebtbifstats.c", "escbox/ebtcohrt.h", "escbox/ebtdopri5.h", "escbox/ebtutils.h", "escbox/ebtmain.h", "escbox/escbox.h")
    objlist <- c("PSPMecodyn.o", "ebtcohrt.o", "ebtdopri5.o", "ebtutils.o")
    if (file.exists(libfile.fullname)) file.remove(libfile.fullname)
    file.remove(Filter(file.exists, objlist))

    cat("\nBuilding executable", libfile.basename, "using sources from", PSPMsrcdir.name, "...\n\n")

    # Check whether the source files are to be found and copy them to the temporary directory
    for (i in srclist) {
      fname <- paste0(PSPMsrcdir.name, '/', i)
      if (!file.exists(fname)) {
        setwd(oldwd)
        stop(paste("\nFile:", fname, "not found!\n Reinstall the PSPManalysis package\n\n"))
      }
      else file.copy(fname, tmpdir, overwrite = TRUE)
    }

    for (i in inclist) {
      fname <- paste0(PSPMsrcdir.name, '/', i)
      if (!file.exists(fname)) {
        setwd(oldwd)
        stop(paste("\nFile:", fname, "not found!\nReload the PSPManalysis.R script with source() and make sure it is in the same directory as all other PSPManalysis C source files\n\n"))
      }
    }

    # Copy the header file to the temporary directory
    file.copy(hfile.abspath, tmpdir, overwrite = TRUE)

    # Get the model dimension settings from the source file
    if (rmodel == 0) {
      con <- file(hfile.abspath, "r")
      PopulationNr <- -1;
      EnvironmentDim <- -1;
      IStateDimension <- -1;
      LifeHistoryStages <- -1;
      ImpactDimension <- -1;
      ParameterNr <- -1;
      while ( TRUE ) {
        line <- readLines(con, n = 1)
        if ( length(line) == 0 ) {
          break
        }
        txt <- "^#define\\s+POPULATION_NR\\s+"
        startpos <- regexpr(paste0(txt, "[1-9][0-9]*"), line)
        if (startpos != -1) {
          numstr <- sub(txt, "", substr(line,startpos,attr(startpos, "match.length")))
          PopulationNr <- strtoi(numstr);
        }
        txt <- "^#define\\s+ENVIRON_DIM\\s+"
        startpos <- regexpr(paste0(txt, "[1-9][0-9]*"), line)
        if (startpos != -1) {
          numstr <- sub(txt, "", substr(line,startpos,attr(startpos, "match.length")))
          EnvironmentDim <- strtoi(numstr);
        }
        txt <- "^#define\\s+I_STATE_DIM\\s+"
        startpos <- regexpr(paste0(txt, "[1-9][0-9]*"), line)
        if (startpos != -1) {
          numstr <- sub(txt, "", substr(line,startpos,attr(startpos, "match.length")))
          IStateDimension <- strtoi(numstr);
        }
        txt <- "^#define\\s+STAGES\\s+"
        startpos <- regexpr(paste0(txt, "[1-9][0-9]*"), line)
        if (startpos != -1) {
          numstr <- sub(txt, "", substr(line,startpos,attr(startpos, "match.length")))
          LifeHistoryStages <- strtoi(numstr);
        }
        txt <- "^#define\\s+PARAMETER_NR\\s+"
        startpos <- regexpr(paste0(txt, "[1-9][0-9]*"), line)
        if (startpos != -1) {
          numstr <- sub(txt, "", substr(line,startpos,attr(startpos, "match.length")))
          ParameterNr <- strtoi(numstr);
        }
        txt <- "^#define\\s+INTERACT_DIM\\s+"
        startpos <- regexpr(paste0(txt, "[1-9][0-9]*"), line)
        if (startpos != -1) {
          numstr <- sub(txt, "", substr(line,startpos,attr(startpos, "match.length")))
          ImpactDimension <- strtoi(numstr);
        }
      }
      close(con)
    }

    if (PopulationNr < 1) stop(paste0("POPULATION_NR equals ", PopulationNr, ", but should be larger than 0!\n\n"))
    if (EnvironmentDim < 1) stop(paste0("ENVIRON_DIM equals ", EnvironmentDim, ", but should be larger than 0!\n\n"))
    if (IStateDimension < 1) stop(paste0("I_STATE_DIM equals ", IStateDimension, ", but should be larger than 0!\n\n"))
    if (LifeHistoryStages < 1) stop(paste0("STAGES equals ", LifeHistoryStages, ", but should be larger than 0!\n\n"))
    if (ParameterNr < 1) stop(paste0("PARAMETER_NR equals ", ParameterNr, ", but should be larger than 0!\n\n"))
    if (ImpactDimension < 1) stop(paste0("INTERACT_DIM equals ", ImpactDimension, ", but should be larger than 0!\n\n"))

    # Construct the command line
    buildargs <- paste0("--output=\"", libfile.basename, "\"")
    buildargs <- paste(buildargs, "PSPMecodyn.c ebtcohrt.c ebtdopri5.c ebtutils.c", sep=" ")

    # Define the basic compilation flags
    cppflags <- "-DR_PACKAGE"
    if (exists("CFLAGS")) cppflags <- paste0(cppflags, " ", get("CFLAGS"))
    if (debug) cppflags <- paste0(cppflags, " -DDEBUG=1 -g -Wall")
    cppflags <- paste0(cppflags, " -I.", " -I\"", PSPMsrcdir.name, "\" -I\"", PSPMsrcdir.name, "/escbox\" ")

    hasfftw <- FALSE
    if (hasfftw) cppflags <- paste0(cppflags, " -DHAS_FFTW3=1")

    # Define the model-specific flags
    modelflags <- paste0(" -DRFUNCTIONS=", rmodel)
    if (rmodel == 0) modelflags <- paste0(modelflags, " -DPROBLEMHEADER=\"", hfile.basename, "\"")

    modelflags <- paste0(modelflags, " -DENVIRON_DIM=", as.integer(EnvironmentDim))
    modelflags <- paste0(modelflags, " -DPOPULATION_NR=", as.integer(PopulationNr))
    modelflags <- paste0(modelflags, " -DSTAGES=", as.integer(LifeHistoryStages))
    modelflags <- paste0(modelflags, " -DI_STATE_DIM=", as.integer(IStateDimension))
    modelflags <- paste0(modelflags, " -DINTERACT_DIM=", as.integer(ImpactDimension))
    modelflags <- paste0(modelflags, " -DPARAMETER_NR=", as.integer(ParameterNr))

    # Define statements in case of a model defined in R
    if (rmodel == 1) {
      if (exists("DiscreteChanges", envir = .GlobalEnv)) modelflags <- paste0(modelflags, " -DDISCRETECHANGES=1")
      else modelflags <- paste0(modelflags, " -DDISCRETECHANGES=0")
      if (exists("NumericalOptions", envir = .GlobalEnv)) {
        numopts <- get("NumericalOptions", envir = .GlobalEnv)
        for (i in 1:length(numopts)) modelflags <- paste0(modelflags, " -D", names(numopts)[i], "=", as.character(numopts[i]))
      }
    }
    
    if (hasfftw) buildenv <- c(PKG_CPPFLAGS = paste(cppflags, modelflags), PKG_LIBS = "-lm $(FFTW3_LIBS)")
    else buildenv <- c(PKG_CPPFLAGS = paste(cppflags, modelflags), PKG_LIBS = "-lm")

    if (debug) {
      cat("Command                :  R CMD SHLIB\n")
      cat(paste0("Arguments              :  ", paste(buildargs, sep=" ", collapse=" "), "\n"))
      cat("Environment variables  :  \n")
      for (i in 1:length(buildenv)) {
          cat(sprintf("\t%14s = %s\n", names(buildenv)[i], buildenv[i]))
        }
      cat("\n")
    }

    # Compilation steps using the newer R package 'pkgbuild'
    # result <- try(pkgbuild::rcmd_build_tools("SHLIB", cmdargs=buildargs, env=buildenv, wd=tmpdir), silent = TRUE)
    # if (is.list(result) && ('status' %in% names(result))) {
    #   if (result$status != 0) {
    #     setwd(oldwd)
    #     cat(result$stdout)
    #     cat(result$stderr)
    #     stop(paste0("\nCompilation of ", libfile.basename, " failed!\n"))
    #   }
    #   else cat(result$stderr)
    # }
    # else {
    # Compilation steps using the older R package 'devtools'
    result <- devtools::RCMD("SHLIB", options=buildargs, path = tmpdir, env_vars=buildenv)
    if ((is.integer(result) && (result != 0)) || (is.logical(result) && (!result))) {
      setwd(oldwd)
      stop(paste0("\nCompilation of ", libfile.basename, " failed!\n"))
    }
    # }

    if (!file.exists(libfile.fullname)) {
        setwd(oldwd)
        stop(paste0("\nCompilation succeeded, but file ", libfile.basename, " not found!\n"))
      }
    else cat("\nCompilation of ", libfile.basename, " succeeded!\n\n")
  }
  else cat("Dynamic library file", libfile.basename, "is up-to-date\n\n");

  setwd(oldwd)
  dyn.load(libfile.fullname)
  cout <- .Call("PSPMecodyn", model.Name, startstate, timepars, bifpars, parameters, options, PACKAGE=paste0(model.Name, "ecodyn"))
  dyn.unload(libfile.fullname)

  if (rmodel == 1) {
    rm(list = Filter( exists, varlist ), envir = .GlobalEnv )
    rm(list = Filter( exists, funlist ), envir = .GlobalEnv )
  }

  desc <- data <- NULL
  if (exists("cout")) {
    outfile.name <- paste0(cout, ".out")
    if (file.exists(outfile.name) && (file.info(outfile.name)$size > 0)) {
      desc <- readLines(outfile.name)
      data <- as.matrix(read.table(text=desc, blank.lines.skip = TRUE, fill=TRUE))
      desc <- desc[grepl("^#", desc)]
      cat(desc, sep='\n')
    }
  }

  setwd(oldwd)
  if (length(desc) || length(data)) {
    output <- list(curvedesc = desc, curvepoints = data)
    return(output)
  }
  else return
}


