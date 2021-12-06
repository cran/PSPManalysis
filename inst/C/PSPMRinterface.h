/*
  PSPMRinterface.h -  Header file containing all the function that take care of interacting with R,
                      in case the life history model is specified by R functions.

  Last modification: AMdR - Oct 14, 2020
*/

#include <R.h>
#include <Rdefines.h>

/*

 *==============================  Timing of PSPMequi  ================================================================================

  source("../PSPManalysis.R")

  system.time(PSPMequi("PNAS2002.R", "EQ", c(0.000253602, 8.8569E-06, 0, 4.00802E-06, 2.33334E-06), -0.1, c(1, 0, 0.0004), NULL, NULL, clean = TRUE, force = TRUE))

  Without compilation using cmpfun():
 
      user  system elapsed 
   280.812   1.923 283.578 

  With compilation using cmpfun():
 
     user  system elapsed 
  270.097   1.701 272.532 
  281.084   1.169 283.046 

  With compilation using cmpfun() and enableJIT(3):
 
     user  system elapsed 
  275.394   1.611 277.693 

  C-specified model:

  system.time(PSPMequi("PNAS2002", "EQ", c(0.000253602, 8.8569E-06, 0, 4.00802E-06, 2.33334E-06), -0.1, c(1, 0, 0.0004), NULL, NULL, clean = TRUE, force = TRUE))

   user  system elapsed 
  1.464   0.159   1.665 

 *==============================  Timing of PSPMecodyn  ==============================================================================

  source("../PSPManalysis.R")
  init <- csbread("Checks/Initstate", 1)

  system.time(PSPMecodyn("PNAS2002.R", init, c(1, 1, 0, 500), clean = TRUE, force = TRUE))

     user  system elapsed 
  675.784   5.742 687.237
  677.218   3.112 682.694 
 
  system.time(PSPMecodyn("PNAS2002", init, c(1, 1, 0, 500), clean = TRUE, force = TRUE))

    user  system elapsed 
  1.910   0.247   2.478 

 */

/*
 *====================================================================================================================================
 * Static variables specifically needed when the model is specified in R
 *====================================================================================================================================
 */

static long int                   N_Protected = 0L;                                 // initialize this with zero at the first time

static SEXP                       R_BirthStates;
static int                        *R_BirthStatesPnt;
static SEXP                       R_lifestage;
static int                        *R_lifestagePnt;

#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
static SEXP                       R_E, R_EnvType, R_EnvNames;
static double                     *R_EPnt;
#endif
static SEXP                       R_Parms, R_DefParms, R_ParNames;
static double                     *R_ParmsPnt;
static SEXP                       R_istate;
static double                     *R_istatePnt;
static SEXP                       R_bstate;
static double                     *R_bstatePnt;
static SEXP                       R_BirthStateNr;
static int                        *R_BirthStateNrPnt;
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
static SEXP                       R_I;
static double                     *R_IPnt = NULL;
#endif

static SEXP                       R_StateAtBirth;
static SEXP                       R_LifeHistoryRates;
static SEXP                       R_LifeStageEndings;
#if (DISCRETECHANGES == 1)
static SEXP                       R_DiscreteChanges;
#endif
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
static SEXP                       R_EnvEqui;
#endif

void my_unprotect(int n)
{
  UNPROTECT(n);
  N_Protected -= n;
  return;
}


void Rinterface_Init()
{
  int i;

  N_Protected = 0L;

  PROTECT(R_DefParms = findVar(install("DefaultParameters"), R_GlobalEnv));
  N_Protected++;
  PROTECT(R_ParNames = getAttrib(R_DefParms, R_NamesSymbol));
  N_Protected++;
  for (i = 0; i < ParameterNr; i++) parameternames[i] = CHAR(STRING_ELT(R_ParNames, i));
  memcpy(parameter, REAL(R_DefParms), ParameterNr*sizeof(double));

  PROTECT(R_Parms = allocVector(REALSXP, ParameterNr));
  N_Protected++;
  setAttrib(R_Parms, R_NamesSymbol, R_ParNames);
  R_ParmsPnt = REAL(R_Parms);

#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  PROTECT(R_EnvType = findVar(install("EnvironmentState"), R_GlobalEnv));
  N_Protected++;
  PROTECT(R_EnvNames = getAttrib(R_EnvType, R_NamesSymbol));
  N_Protected++;
  for (i = 0; i < EnvironDim; i++)
    {
      if (!strcmp(CHAR(STRING_ELT(R_EnvType, i)), "GENERALODE"))
        EnvironmentType[i] = GENERALODE;
      else if (!strcmp(CHAR(STRING_ELT(R_EnvType, i)), "PERCAPITARATE"))
        EnvironmentType[i] = PERCAPITARATE;
      else if (!strcmp(CHAR(STRING_ELT(R_EnvType, i)), "POPULATIONINTEGRAL"))
        EnvironmentType[i] = POPULATIONINTEGRAL;
      else
        error("\nEnvironmental variable %d not defined as \"GENERALODE\", \"PERCAPITARATE\", or \"POPULATIONINTEGRAL\"!\n\n", i);
    }

  PROTECT(R_E = allocVector(REALSXP, EnvironDim));
  N_Protected++;
  setAttrib(R_E, R_NamesSymbol, R_EnvNames);
  R_EPnt = REAL(R_E);
#endif

  PROTECT(R_BirthStates = allocVector(INTSXP, PopulationNr));
  N_Protected++;
  R_BirthStatesPnt = INTEGER(R_BirthStates);

  PROTECT(R_BirthStateNr = allocVector(INTSXP, 1));
  N_Protected++;
  R_BirthStateNrPnt = INTEGER(R_BirthStateNr);

#if (POPULATION_NR > 1)
  PROTECT(R_lifestage = allocVector(INTSXP, PopulationNr));
  N_Protected++;
  PROTECT(R_istate = allocMatrix(REALSXP, PopulationNr, IStateDim));
  N_Protected++;
  PROTECT(R_bstate = allocMatrix(REALSXP, PopulationNr, IStateDim));
  N_Protected++;
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  PROTECT(R_I = allocMatrix(REALSXP, PopulationNr, InteractDim));
  N_Protected++;
#endif
#else
  PROTECT(R_lifestage = allocVector(INTSXP, PopulationNr));
  N_Protected++;
  PROTECT(R_istate = allocVector(REALSXP, IStateDim));
  N_Protected++;
  PROTECT(R_bstate = allocVector(REALSXP, IStateDim));
  N_Protected++;
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  PROTECT(R_I = allocVector(REALSXP, InteractDim));
  N_Protected++;
#endif
#endif
  R_lifestagePnt  = INTEGER(R_lifestage);
  R_istatePnt     = REAL(R_istate);
  R_bstatePnt     = REAL(R_bstate);
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  R_IPnt = REAL(R_I);
#endif

  R_StateAtBirth = PROTECT(allocList(3));
  N_Protected++;
  SET_TYPEOF(R_StateAtBirth, LANGSXP);
  SETCAR(R_StateAtBirth, install("StateAtBirth"));

  R_LifeHistoryRates = PROTECT(allocList(7));
  N_Protected++;
  SET_TYPEOF(R_LifeHistoryRates, LANGSXP);
  SETCAR(R_LifeHistoryRates, install("LifeHistoryRates"));

  R_LifeStageEndings = PROTECT(allocList(7));
  N_Protected++;
  SET_TYPEOF(R_LifeStageEndings, LANGSXP);
  SETCAR(R_LifeStageEndings, install("LifeStageEndings"));

#if (DISCRETECHANGES == 1)
  R_DiscreteChanges = PROTECT(allocList(7));
  N_Protected++;
  SET_TYPEOF(R_DiscreteChanges, LANGSXP);
  SETCAR(R_DiscreteChanges, install("DiscreteChanges"));
#endif

#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  R_EnvEqui = PROTECT(allocList(4));
  N_Protected++;
  SET_TYPEOF(R_EnvEqui, LANGSXP);
  SETCAR(R_EnvEqui, install("EnvEqui"));
#endif

  return;
}


void Rinterface_End(void)
{
  UNPROTECT((int)N_Protected);
  return;
}


/*
 *====================================================================================================================================
 *  SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
 *====================================================================================================================================
 */

/*
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 */

void SetBirthStates(int BirthStates[PopulationNr], double E[])
{
  int  nlen;
  SEXP Result;                                                                      // Protected
  SEXP R_pnt;                                                                       // Unprotected

  R_pnt = R_StateAtBirth;
  R_pnt = CDR(R_pnt);
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  memcpy(R_EPnt, E, EnvironDim*sizeof(double));
  SETCAR(R_pnt, R_E);
#else
  SETCAR(R_pnt, R_NilValue);
#endif
  R_pnt = CDR(R_pnt);
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(Result = eval(R_StateAtBirth, R_GlobalEnv));
  N_Protected++;

  nlen = length(Result);
#if (POPULATION_NR == 1)
  if (nlen == IStateDim)
    {
      BirthStates[0] = 1;
      setAttrib(R_istate, R_NamesSymbol, getAttrib(Result, R_NamesSymbol));
      setAttrib(R_bstate, R_NamesSymbol, getAttrib(Result, R_NamesSymbol));
    }
  else if ((nlen % IStateDim == 0) && isMatrix(Result) && (ncols(Result) == IStateDim))
    {
      BirthStates[0] = nrows(Result);
      setAttrib(R_istate, R_NamesSymbol, GetColNames(Result));
      setAttrib(R_bstate, R_NamesSymbol, GetColNames(Result));
    }
  else
    error("\nStateAtBirth() should return a vector with %d elements or a matrix with %d columns\n and a number of rows equal to the number of "
          "possible states at birth!\n\n",
          IStateDim, IStateDim);
#else
  int  pp;
  SEXP t, names;

  if (isMatrix(Result) && (ncols(Result) == IStateDim) && (nrows(Result) == PopulationNr))
    {
      for (pp = 0; pp < PopulationNr; pp++) BirthStates[pp] = 1;
      setAttrib(R_istate, R_DimNamesSymbol, getAttrib(Result, R_DimNamesSymbol));
      setAttrib(R_bstate, R_DimNamesSymbol, getAttrib(Result, R_DimNamesSymbol));
    }
  else if (isArray(Result) && (nlen % (PopulationNr*IStateDim) == 0))
    {
      t = getAttrib(Result, R_DimSymbol);
      if (TYPEOF(t) == INTSXP && (LENGTH(t) == 3) && (INTEGER(t)[0] == PopulationNr) && (INTEGER(t)[2] == IStateDim))
        {
          for (pp = 0; pp < PopulationNr; pp++) BirthStates[pp] = INTEGER(t)[1];
          PROTECT(names = allocVector(VECSXP, 2));
          N_Protected++;

          SET_VECTOR_ELT(names, 0, VECTOR_ELT(getAttrib(Result, R_DimNamesSymbol), 0));
          SET_VECTOR_ELT(names, 1, VECTOR_ELT(getAttrib(Result, R_DimNamesSymbol), 2));

          setAttrib(R_istate, R_DimNamesSymbol, names);
          setAttrib(R_bstate, R_DimNamesSymbol, names);
          my_unprotect(1);
        }
      else
        error("\nStateAtBirth() should return a 2- or 3-dimensional array with the first dimension of length %d and the last dimension of length %d!\n\n",
              PopulationNr, IStateDim);
    }
  else
    error("\nStateAtBirth() should return a 2- or 3-dimensional array with the first dimension of length %d and the second dimension of length %d!\n\n",
          PopulationNr, IStateDim);
#endif

  my_unprotect(1);

  return;
}


/*
 * Specify all the possible states at birth for all individuals in all
 * structured populations in the problem. BirthStateNr represents the index of
 * the state of birth to be specified. Each state at birth should be a single,
 * constant value for each i-state variable.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 */

void StateAtBirth(double *istate[PopulationNr], int BirthStateNr, double E[])
{
  int     ii, pp;
  SEXP    Result;                                                                   // Protected
  SEXP    R_pnt;                                                                    // Unprotected
  double  *respnt;

  R_pnt = R_StateAtBirth;
  R_pnt = CDR(R_pnt);
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  memcpy(R_EPnt, E, EnvironDim*sizeof(double));
  SETCAR(R_pnt, R_E);
#else
  SETCAR(R_pnt, R_NilValue);
#endif
  R_pnt = CDR(R_pnt);
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(Result = eval(R_StateAtBirth, R_GlobalEnv));
  N_Protected++;
  respnt = REAL(coerceVector(Result, REALSXP));

  for (pp = 0; pp < PopulationNr; pp++)
    {
      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < IStateDim; ii++) 
        *(istate[pp] + ii) = *(respnt + ii*PopulationNr*MaxStatesAtBirth  + PopulationNr*BirthStateNr + pp);
    }

  my_unprotect(1);

  return;
}


/*
 * Specify the threshold determining the end point of each discrete life
 * stage in individual life history as function of the i-state variables and
 * the individual's state at birth for all populations in every life stage.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 */

// Some array dimensions have to be explicitly defined as POPULATION_NR to avoid errors with access attributes when using gcc11

void IntervalLimit(int lifestage[PopulationNr], double *istate[PopulationNr], double *bstate[POPULATION_NR], int BirthStateNr, double E[],
                   double limit[POPULATION_NR])
{
  int  ii, pp;
  SEXP R_limit;                                                                     // Protected
  SEXP R_pnt;                                                                       // Unprotected

  for (pp = 0; pp < PopulationNr; pp++)
    {
      R_lifestagePnt[pp] = lifestage[pp] + 1;

      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < IStateDim; ii++) 
        {
          *(R_istatePnt + ii*PopulationNr + pp) = *(istate[pp] + ii);
          *(R_bstatePnt + ii*PopulationNr + pp) = *(bstate[pp] + ii);
        }
    }
  *R_BirthStateNrPnt = BirthStateNr + 1;
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  memcpy(R_EPnt, E, EnvironDim*sizeof(double));
#endif
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));

  R_pnt = R_LifeStageEndings;
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_lifestage);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_istate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_bstate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_BirthStateNr);
  R_pnt = CDR(R_pnt);
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  SETCAR(R_pnt, R_E);
#else
  SETCAR(R_pnt, R_NilValue);
#endif
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(R_limit = eval(R_LifeStageEndings, R_GlobalEnv));
  N_Protected++;
  
#if (POPULATION_NR == 1)
  if (length(R_limit) != 1) error("\nMaturation through the stage not specified as a single numeric value!\n\n");
#else
  if (length(R_limit) != PopulationNr)
    error("\nMaturation through the stage not specified as a vector with %d elements!\n\n", PopulationNr);
#endif
  memcpy(limit, REAL(R_limit), PopulationNr*sizeof(double));
  my_unprotect(1);

  return;
}


/*
 * Specify the development, fecundity and mortality of individuals and its impact on
 * its environment as a function of the i-state variables and the individual's state
 * at birth for all populations in every life stage.
 *
 * This routine is only used in EBT simulations
 */

#if (defined(PSPMECODYN))

// Some array dimensions have to be explicitly defined as POPULATION_NR to avoid errors with access attributes when using gcc11

void SetLifeHistoryRates(int lifestage[PopulationNr], double *istate[PopulationNr],
                         double *bstate[POPULATION_NR], int BirthStateNr, double E[],
                         double development[POPULATION_NR][I_STATE_DIM],
                         double *fecundity[POPULATION_NR],
                         double mortality[POPULATION_NR],
                         double impact[POPULATION_NR][INTERACT_DIM])
{
  int     bb, ii, jj, pp, nprotect = 0;
  SEXP    Result, names, R_pnt, R_dev, R_fec, R_mort, R_imp;                                             // Protected
  int     devdone = 0, fecdone = 0, mordone = 0, impdone = 0;

  for (pp = 0; pp < PopulationNr; pp++)
    {
      R_lifestagePnt[pp] = lifestage[pp] + 1;

      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < IStateDim; ii++) 
      {
        *(R_istatePnt + ii*PopulationNr + pp) = *(istate[pp] + ii);
        *(R_bstatePnt + ii*PopulationNr + pp) = *(bstate[pp] + ii);
      }
  }
  *R_BirthStateNrPnt = BirthStateNr + 1;
  memcpy(R_EPnt, E, EnvironDim*sizeof(double));
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));

  R_pnt = R_LifeHistoryRates;
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_lifestage);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_istate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_bstate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_BirthStateNr);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_E);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(Result = eval(R_LifeHistoryRates, R_GlobalEnv));
  N_Protected++;
  nprotect++;
  PROTECT(names = getAttrib(Result, R_NamesSymbol));
  N_Protected++;
  nprotect++;

  // Get the pointers to the list components
  for (jj = 0; jj < length(Result); jj++)
    {
      if (development && (!strcmp(CHAR(STRING_ELT(names, jj)), "development")))
        {
          PROTECT(R_dev = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
#if (POPULATION_NR == 1)
          if (length(R_dev) != IStateDim) error("\nDevelopment not specified as a vector %d elements!\n\n", IStateDim);
#else
          if (length(R_dev) != (PopulationNr*IStateDim))
            error("\nDevelopment not specified as a matrix with %d rows and %d columns!\n\n", PopulationNr, IStateDim);
#endif
          for (ii = 0; ii < IStateDim; ii++)
            for (pp = 0; pp < PopulationNr; pp++) 
              development[pp][ii] = *(REAL(R_dev) + ii*PopulationNr + pp);
          devdone = 1;
        }
      else if (fecundity && (!strcmp(CHAR(STRING_ELT(names, jj)), "fecundity")))
        {
          PROTECT(R_fec = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
#if (POPULATION_NR == 1)
          if (length(R_fec) != MaxStatesAtBirth) error("\nFecundity not specified as a vector with %d elements!\n\n", MaxStatesAtBirth);
#else
          if (length(R_fec) != (PopulationNr*MaxStatesAtBirth))
            error("\nFecundity not specified as a matrix with %d rows and %d columns!\n\n", PopulationNr, MaxStatesAtBirth);
#endif
          for (bb = 0; bb < MaxStatesAtBirth; bb++)
            for (pp = 0; pp < PopulationNr; pp++) 
              fecundity[pp][bb] = *(REAL(R_fec) + bb*PopulationNr + pp);
          fecdone = 1;
        }
      else if (mortality && (!strcmp(CHAR(STRING_ELT(names, jj)), "mortality")))
        {
          PROTECT(R_mort = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
          if (length(R_mort) != PopulationNr) error("\nMortality not specified as a vector with %d elements!\n\n", PopulationNr);
          memcpy(mortality, REAL(R_mort), PopulationNr*sizeof(double));
          mordone = 1;
        }
      else if (impact && (!strcmp(CHAR(STRING_ELT(names, jj)), "impact")))
        {
          PROTECT(R_imp = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
#if (POPULATION_NR == 1)
          if (length(R_imp) != InteractDim) error("\nImpact not specified as a vector %d elements!\n\n", InteractDim);
#else
          if (length(R_imp) != (PopulationNr*InteractDim))
            error("\nImpact not specified as a matrix with %d rows and %d columns!\n\n", PopulationNr, InteractDim);
#endif
          for (ii = 0; ii < InteractDim; ii++)
            for (pp = 0; pp < PopulationNr; pp++) 
              impact[pp][ii] = *(REAL(R_imp) + ii*PopulationNr + pp);
          impdone = 1;
        }
    }

  my_unprotect(nprotect);

  if (development && !devdone)
    error("\nModel function LifeHistoryRates() does not return a list element 'development'!\n\n");
  if (fecundity && !fecdone)
    error("\nModel function LifeHistoryRates() does not return a list element 'fecundity'!\n\n");
  if (mortality && !mordone)
    error("\nModel function LifeHistoryRates() does not return a list element 'mortality'!\n\n");
  if (impact && !impdone)
    error("\nModel function LifeHistoryRates() does not return a list element 'impact'!\n\n");

  return;
}

#else

/*
 * Specify the development, fecundity and mortality of individuals and its impact on
 * its environment as a function of the i-state variables and the individual's state
 * at birth for all populations in every life stage.
 *
 * Also assign these rates to the appropriate ODEs (not used in EBT simulations)
 */

// Some array dimensions have to be explicitly defined as POPULATION_NR to avoid errors with access attributes when using gcc11

void SetDerivatives(const int totalOdeDim, const int *birthStateNr, const int BirthStateNr, int populationStart[PopulationNr],
                    int lifeStage[PopulationNr], double *bstate[POPULATION_NR], double age, double *y, double *dyda, double *discard)
{
  int     ii, pp;
  double  *iStatePnt[PopulationNr];
  
  double  *development = NULL, *fecundity = NULL, *mortality = NULL;
#if (defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  double  *impact = NULL;
  SEXP    R_imp;                                                                    // Protected
#endif
  double  *derpnt;
  int     jj, nprotect = 0;
  double  survival;
#if (defined(PSPMDEMO) && (PULSED == 0))
  double  cumrepro;
  SEXP    R_fec;                                                                    // Protected
#elif (defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  SEXP    R_fec;                                                                    // Protected
#endif
  SEXP    Result, names, R_dev, R_mort;                                             // Protected
  SEXP    R_pnt;                                                                    // Unprotected

  memset(dyda, 0, totalOdeDim*sizeof(double));
  for (pp = 0; pp < PopulationNr; pp++)
    {
      R_lifestagePnt[pp] = lifeStage[pp] + 1;
      iStatePnt[pp]      = (populationStart[pp] < 0) ? bstate[pp] : (y + populationStart[pp]);

      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < IStateDim; ii++) 
      {
        *(R_istatePnt + ii*PopulationNr + pp) = *(iStatePnt[pp] + ii);
        *(R_bstatePnt + ii*PopulationNr + pp) = *(bstate[pp] + ii);
      }
  }
  *R_BirthStateNrPnt = BirthStateNr + 1;
#if (defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  memcpy(R_EPnt, Evar, EnvironDim*sizeof(double));
#endif
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));

  R_pnt = R_LifeHistoryRates;
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_lifestage);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_istate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_bstate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_BirthStateNr);
  R_pnt = CDR(R_pnt);
#if (defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  SETCAR(R_pnt, R_E);
#else
  SETCAR(R_pnt, R_NilValue);
#endif
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(Result = eval(R_LifeHistoryRates, R_GlobalEnv));
  N_Protected++;
  nprotect++;
  PROTECT(names = getAttrib(Result, R_NamesSymbol));
  N_Protected++;
  nprotect++;

  // Get the pointers to the list components
  for (jj = 0; jj < length(Result); jj++)
    {
      if (!strcmp(CHAR(STRING_ELT(names, jj)), "development"))
        {
          PROTECT(R_dev = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
#if (POPULATION_NR == 1)
          if (length(R_dev) != IStateDim) error("\nDevelopment not specified as a vector %d elements!\n\n", IStateDim);
#else
          if (length(R_dev) != (PopulationNr*IStateDim))
            error("\nDevelopment not specified as a matrix with %d rows and %d columns!\n\n", PopulationNr, IStateDim);
#endif
          development = REAL(R_dev);
        }
#if ((defined(PSPMDEMO) && (PULSED == 0)) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
      else if (!strcmp(CHAR(STRING_ELT(names, jj)), "fecundity"))
        {
          PROTECT(R_fec = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
#if (POPULATION_NR == 1)
          if (length(R_fec) != MaxStatesAtBirth) error("\nFecundity not specified as a vector with %d elements!\n\n", MaxStatesAtBirth);
#else
          if (length(R_fec) != (PopulationNr*MaxStatesAtBirth))
            error("\nFecundity not specified as a matrix with %d rows and %d columns!\n\n", PopulationNr, MaxStatesAtBirth);
#endif
          fecundity = REAL(R_fec);
        }
#endif
      else if (!strcmp(CHAR(STRING_ELT(names, jj)), "mortality"))
        {
          PROTECT(R_mort = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
          if (length(R_mort) != PopulationNr) error("\nMortality not specified as a vector with %d elements!\n\n", PopulationNr);
          mortality = REAL(R_mort);
        }
#if (defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
      else if (!strcmp(CHAR(STRING_ELT(names, jj)), "impact"))
        {
          PROTECT(R_imp = coerceVector(VECTOR_ELT(Result, jj), REALSXP));
          N_Protected++;
          nprotect++;
#if (POPULATION_NR == 1)
          if (length(R_imp) != InteractDim) error("\nImpact not specified as a vector %d elements!\n\n", InteractDim);
#else
          if (length(R_imp) != (PopulationNr*InteractDim))
            error("\nImpact not specified as a matrix with %d rows and %d columns!\n\n", PopulationNr, InteractDim);
#endif
          impact = REAL(R_imp);
        }
#endif
    }

  if (!development)
    error("\nModel function LifeHistoryRates() does not return a list element 'development'!\n\n");
#if ((defined(PSPMDEMO) && (PULSED == 0)) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  if (!fecundity)
    error("\nModel function LifeHistoryRates() does not return a list element 'fecundity'!\n\n");
#endif
  if (!mortality)
    error("\nModel function LifeHistoryRates() does not return a list element 'mortality'!\n\n");

#if (defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  if (!impact)
    error("\nModel function LifeHistoryRates() does not return a list element 'impact'!\n\n");
#endif

  // Now assign the derivatives
  for (pp = 0; pp < PopulationNr; pp++)
    {
      if (populationStart[pp] < 0) continue;

      derpnt = (dyda + populationStart[pp]);

      // R has column-wise storage of matrix elements !!!
      for (jj = 0; jj < IStateDim; jj++, derpnt++)
        *derpnt = development[jj*PopulationNr + pp];
#if defined(PSPMDEMO)
      *derpnt = -(mortality[pp] + PGRvar[pp]);
      derpnt++;
#else
      *derpnt = -(mortality[pp]);
      derpnt++;
#endif
      survival = exp(iStatePnt[pp][IStateDim]);
// Here we assume that the fecundity matrix returned is a rectangular matrix with equal number of
// columns (MaxStatesAtBirth) for all populations. Only the first birthStateNr[pp] entries for each
// population pp are used. Hence it is assumed that the user fills up the possibly remaining entries
// (MaxStatesAtBirth - birthStateNr[pp]) with nonsense values
#if (defined(PSPMDEMO) && (PULSED == 0))
      cumrepro = 0.0;
#endif
      for (jj = 0; jj < birthStateNr[pp]; jj++, derpnt++)
        {
          // R has column-wise storage of matrix elements !!!
          *derpnt = fecundity[jj*PopulationNr + pp]*survival;
#if (defined(PSPMDEMO) && (PULSED == 0))
          cumrepro += age*(*derpnt);
#endif
        }

#if (defined(PSPMDEMO) && (PULSED == 0))
      *derpnt = cumrepro;                                                           // Generation time, last element
#else
      // R has column-wise storage of matrix elements !!!
      for (jj = 0; jj < InteractDim; jj++, derpnt++)
        *derpnt = impact[jj*PopulationNr + pp]*survival;

#if (defined(PSPMEQUI) || defined(PSPMEVODYN)) && (FULLSTATEOUTPUT > 0)
      if (DoStateOutput)
        {
          for (jj = 0; jj < IStateDim; jj++, derpnt++)
            *derpnt = iStatePnt[pp][jj]*survival;                                    // Average i-state in cohort
          *derpnt   = survival;                                                     // Number of individuals in cohort is last
        }
#endif
#endif
    }

  my_unprotect(nprotect);

  return;
}

#endif // defined(PSPMECODYN)

/*
 * Specify the possible discrete changes (jumps) in the individual state
 * variables when ENTERING the stage specified by 'lifestage[]'.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 */

// Some array dimensions have to be explicitly defined as POPULATION_NR to avoid errors with access attributes when using gcc11

void DiscreteChanges(int lifestage[PopulationNr], double *istate[PopulationNr], double *bstate[POPULATION_NR], int BirthStateNr, double E[])
{
#if (DISCRETECHANGES == 1)
  int   ii, pp;

  SEXP  Result, Result2;                                                            // Protected
  SEXP  R_pnt;                                                                      // Unprotected

  for (pp = 0; pp < PopulationNr; pp++)
    {
      R_lifestagePnt[pp] = lifestage[pp] + 1;

      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < IStateDim; ii++) 
      {
        *(R_istatePnt + ii*PopulationNr + pp) = *(istate[pp] + ii);
        *(R_bstatePnt + ii*PopulationNr + pp) = *(bstate[pp] + ii);
      }
    }
  *R_BirthStateNrPnt = BirthStateNr + 1;
#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))
  memcpy(R_EPnt, E, EnvironDim*sizeof(double));
#endif
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));

  R_pnt = R_DiscreteChanges;
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_lifestage);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_istate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_bstate);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_BirthStateNr);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_E);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(Result = eval(R_DiscreteChanges, R_GlobalEnv));
  N_Protected++;
  PROTECT(Result2 = coerceVector(Result, REALSXP));
  N_Protected++;

  // Only copy values for those populations that really change stage 
  for (pp = 0; pp < PopulationNr; pp++)
    {
      if (lifestage[pp] < 1) continue;

      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < IStateDim; ii++)
        *(istate[pp] + ii) = *(REAL(Result2) + ii*PopulationNr + pp);
    }
  my_unprotect(2);

#endif

  return;
}

#if (defined(PSPMECODYN) || defined(PSPMEQUI) || defined(PSPMEVODYN) || defined(PSPMIND))

void EnvEqui(double E[], double I[PopulationNr][InteractDim], double condition[EnvironDim])
{
  int   ii, pp;
  SEXP  Result, Result2;                                                             // Protected
  SEXP  R_pnt;                                                                      // Unprotected

  for (pp = 0; pp < PopulationNr; pp++)
    {
      // R has column-wise storage of matrix elements !!!
      for (ii = 0; ii < InteractDim; ii++)
        *(R_IPnt + ii*PopulationNr + pp) = I[pp][ii];
    }
  memcpy(R_EPnt, E, EnvironDim*sizeof(double));
  memcpy(R_ParmsPnt, parameter, ParameterNr*sizeof(double));

  R_pnt = R_EnvEqui;
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_I);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_E);
  R_pnt = CDR(R_pnt);
  SETCAR(R_pnt, R_Parms);
  R_pnt = CDR(R_pnt);

  PROTECT(Result = eval(R_EnvEqui, R_GlobalEnv));
  N_Protected++;
  PROTECT(Result2 = coerceVector(Result, REALSXP));
  N_Protected++;

  memcpy(condition, REAL(Result2), EnvironDim*sizeof(double));

  my_unprotect(2);

  return;
}

#endif

/*==================================================================================================================================*/
