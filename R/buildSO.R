buildSO <- function(PSPMmodule = NULL, modelname = NULL, debug = FALSE, force = FALSE, silent = FALSE) {

  oldwd <- getwd()
  eval.parent(substitute(Oldwd<-oldwd))

  if ((!length(modelname)) || (!nchar(modelname))) stop("You have to specify a model name")

  pkgmodeldir <- system.file("Models", package="PSPManalysis")
  modelfound <- 0
  modeldname <- dirname(modelname)
  modelbname <- basename(modelname)
  for (modeldir in c(modeldname, pkgmodeldir))
    {
      if (modeldir == pkgmodeldir)
        {
          str = readline(paste0("\nModel file not found in current directory\nContinue searching for model file in package directory? [y(es)/q(uit)] : "))
          if ((str != '') && (str != 'y') && (str != 'Y')) break
        }
      if ((regexpr("\\.h", modelbname) == (nchar(modelbname)-1)) && file.exists(paste0(modeldir, "/", modelbname)))
        {
          hfile.fullname <- paste0(modeldir, "/", modelbname)
          rmodel <- 0
          modelfound <- 1
          break
        }
      if ((regexpr("\\.R", modelbname) == (nchar(modelbname)-1)) && file.exists(paste0(modeldir, "/", modelbname)))
        {
          hfile.fullname <- paste0(modeldir, "/", modelbname)
          rmodel <- 1
          modelfound <- 1
          break
        }
      if (file.exists(paste0(modeldir, "/", modelbname, ".h")))
        {
          hfile.fullname <- paste0(modeldir, "/", modelbname, ".h")
          rmodel <- 0
          modelfound <- 1
          break
        }
      if (file.exists(paste0(modeldir, "/", modelbname, ".R")))
        {
          hfile.fullname <- paste0(modeldir, "/", modelbname, ".R")
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
    rm(list = Filter( exists, varlist ), envir = .GlobalEnv )
    rm(list = Filter( exists, funlist ), envir = .GlobalEnv )

    source(hfile.fullname)
  }

  PSPMsrcdir.name <- normalizePath(system.file("C", package="PSPManalysis"))
  tmpdir <- normalizePath(tempdir())

  hfile.dirname <- dirname(normalizePath(hfile.fullname))
  hfile.basename <- basename(normalizePath(hfile.fullname))
  hfile.abspath <- normalizePath(paste0(hfile.dirname, '/', hfile.basename))
  if (rmodel == 1) model.name <- sub("\\.R", "", hfile.basename) else model.name <- sub("\\.h", "", hfile.basename)

  libfile.basename <- paste0(model.name, sub("PSPM", "", PSPMmodule), .Platform$dynlib.ext)
  libfile.fullname <- paste0(tmpdir, '/', libfile.basename)

  if (is.loaded(libfile.fullname)) dyn.unload(libfile.fullname)

  # Check whether the executable is up-to-date
  build.flag <- !file.exists(libfile.fullname);

  if (!build.flag) {
    # The comma in front of mtime is needed to get output as class POSIXct
    build.flag <- (difftime(file.info(hfile.fullname)[,'mtime'], file.info(libfile.fullname)[,'mtime'], units = "secs") >= 0)
  }

  if (build.flag || force) {
    # Switch to the temporary directory for building
    setwd(tmpdir)

    srclist <- c(paste0(PSPMmodule, ".c"), "biftest.c", "curve.c", "io.c")
    objlist <- c(paste0(PSPMmodule, ".o"), "biftest.o", "curve.o", "io.o")

    if (file.exists(libfile.fullname)) file.remove(libfile.fullname)
    if (file.exists(libfile.basename)) file.remove(libfile.basename)
    if (file.exists(hfile.basename)) file.remove(hfile.basename)
    if (file.exists(paste0(PSPMmodule, ".o"))) file.remove(paste0(PSPMmodule, ".o"));

    if (force) file.remove(Filter(file.exists, srclist))
    if (force) file.remove(Filter(file.exists, objlist))

    if (!silent) cat("\nBuilding executable", libfile.basename, "using sources from", PSPMsrcdir.name, "...\n\n")

    # Check whether the source files are to be found and copy them to the temporary directory
    for (i in srclist) {
      fname <- paste0(PSPMsrcdir.name, "/", i)
      if (!file.exists(fname)) {
        setwd(oldwd)
        stop(paste("\nFile:", fname, "not found! Reinstall the PSPManalysis package\n\n"))
      }
      else file.copy(fname, tmpdir, overwrite = TRUE)
    }

    # Copy the header file to the temporary directory
    file.copy(hfile.abspath, tmpdir, overwrite = TRUE)

    # Construct the command line
    buildargs <- c(paste0("--output=\'", libfile.basename, "\'"), srclist)

    # Define the basic compilation flags
    cppflags <- "-DR_PACKAGE"
    if (exists("CFLAGS")) cppflags <- paste0(cppflags, " ", get("CFLAGS"))
    if (debug) cppflags <- paste0(cppflags, " -DDEBUG=1 -g -Wall")
    else cppflags <- paste0(cppflags, " -Wno-format-overflow")
    cppflags <- paste0(cppflags, " -I.", " -I\"", PSPMsrcdir.name, "\"")

    # Define the model-specific flags
    modelflags <- paste0(" -DRFUNCTIONS=", rmodel)
    if (rmodel == 0) modelflags <- paste0(modelflags, " -DPROBLEMHEADER=\"", hfile.basename, "\"")

    # Define statements in case of a model defined in R
    if (rmodel == 1) {
      pspmdims <- get("PSPMdimensions", envir = .GlobalEnv)
      modelflags <- paste0(modelflags, " -DPOPULATION_NR=", as.integer(pspmdims["PopulationNr"]))
      modelflags <- paste0(modelflags, " -DSTAGES=", as.integer(pspmdims["LifeHistoryStages"]))
      modelflags <- paste0(modelflags, " -DI_STATE_DIM=", as.integer(pspmdims["IStateDimension"]))
      if (is.element('ImpactDimension', names(pspmdims))) modelflags <- paste0(modelflags, " -DINTERACT_DIM=", as.integer(pspmdims["ImpactDimension"]))
      if (exists("EnvironmentState")) modelflags <- paste0(modelflags, " -DENVIRON_DIM=", length(get("EnvironmentState", envir = .GlobalEnv)))
      modelflags <- paste0(modelflags, " -DPARAMETER_NR=", length(get("DefaultParameters", envir = .GlobalEnv)))
      if (exists("DiscreteChanges", envir = .GlobalEnv)) modelflags <- paste0(modelflags, " -DDISCRETECHANGES=1")
      else modelflags <- paste0(modelflags, " -DDISCRETECHANGES=0")
      if (exists("NumericalOptions", envir = .GlobalEnv)) {
        numopts <- get("NumericalOptions", envir = .GlobalEnv)
        for (i in 1:length(numopts)) modelflags <- paste0(modelflags, " -D", names(numopts)[i], "=", as.character(numopts[i]))
      }
    }

    # Define the environment variables to be set
    buildenv <- c(PKG_CFLAGS = paste(cppflags, modelflags), PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS)")

    if (debug) {
      cat("Current working directory: ", getwd(), "\n")
      cat("Command                  :  R CMD SHLIB\n")
      cat(paste0("Arguments                :  ", paste(buildargs, sep=" ", collapse=" "), "\n"))
      cat("Environment variables    :  \n")
      for (i in 1:length(buildenv)) {
          cat(sprintf("\t%14s = %s\n", names(buildenv)[i], buildenv[i]))
        }
      cat("\n")
    }

    # Compilation steps using the newer R package 'pkgbuild'
    result <- try(pkgbuild::rcmd_build_tools("SHLIB", cmdargs=buildargs, env=buildenv, wd=tmpdir), silent = TRUE)
    if (is.list(result) && ('status' %in% names(result))) {
      if (result$status != 0) {
        setwd(oldwd)
        cat(result$stdout)
        cat(result$stderr)
        stop(paste0("\nCompilation of ", libfile.basename, " failed!\n"))
      }
      else cat(result$stderr)
    }

    if (!silent) cat("\nCompilation of ", libfile.basename, " succeeded!\n\n")
  }
  else if (!silent) cat("\nDynamic library file", libfile.basename, "is up-to-date\n");

  setwd(oldwd)

  eval.parent(substitute(model.Name<-model.name))
  eval.parent(substitute(Rmodel<-rmodel))
  eval.parent(substitute(libfile.Basename<-libfile.basename))
  if (rmodel == 1) {
    eval.parent(substitute(Varlist<-varlist))
    eval.parent(substitute(Funlist<-funlist))
  }

  return(libfile.fullname)
}

