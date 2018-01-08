/***
  NAME
    PSPMind

  PURPOSE
    Generic, problem-independent specification for problems with
    POPULATION_NR structured population, whose individuals are characterized
    by I_STATE_DIM state variables, that interact with ENVIRON_DIM
    environment variables. All problem-specific life-history functions are
    specified in an include file

    Copyright (C) 2015, Andre M. de Roos, University of Amsterdam

    This file is part of the PSPManalysis software package.

    PSPManalysis is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    PSPManalysis is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PSPManalysis. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Dec 05, 2017
***/

#define PSPMIND                   1

#if (!defined(RFUNCTIONS) || (RFUNCTIONS != 1))
#define RFUNCTIONS                0
#endif
#if (!defined(MFUNCTIONS) || (MFUNCTIONS != 1))
#define MFUNCTIONS                0
#endif

#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
double                            parameter[PARAMETER_NR];
const char                        *parameternames[PARAMETER_NR];
int                               EnvironmentType[ENVIRON_DIM];
#endif

#include "globals.h"

/*
 *====================================================================================================================================
 *	Import the population and environment dimension settings
 *====================================================================================================================================
 */

#define Survival(p)               (exp(istate[p][IStateDim]))
#define SetSurvival(p, s)         istate[p][IStateDim] = (((s) >= 0) && ((s) < exp(istate[p][IStateDim]))) ? (log(max((s), DBL_EPSILON))) : (istate[p][IStateDim])


#if ((RFUNCTIONS != 1) && (MFUNCTIONS != 1))
#if defined(PROBLEMHEADER)                                                          // Include header file
#define HEADERNAME <PROBLEMHEADER>
#include HEADERNAME
#else
#error No header file defined!
#endif
#endif

#include "defaults.h"

#if !defined(ENVIRON_DIM) || (ENVIRON_DIM < 1)
#error Simulating evolutionary dynamics requires ENVIRON_DIM to be larger than 0
#endif

#if !defined(INTERACT_DIM) || (INTERACT_DIM < 1)
#error INTERACT_DIM should be defined larger than 0
#endif

#if !defined(PARAMETER_NR) || (PARAMETER_NR < 3)
#error PARAMETER_NR should be defined larger than 2
#endif


/*
 *====================================================================================================================================
 *		Definition of global variables and parameters
 *====================================================================================================================================
 */

// Global dimension variables
static int                        BirthStateNr[POPULATION_NR];
static double                     Beq[POPULATION_NR];

// These are the variables to solve for
static double                     Evar[ENVIRON_DIM];

// Global variables to hold variables shared among routines
static int                        LastMemAllocated = 0;

// Global pointers into the heap
static double                     *RightEigenvecMem = NULL;

static double                     *BirthStateMem = NULL;
static double                     *PopDensMem     = NULL;
static int                        *CohortsMem     = NULL;
static double                     *CohortLimitMem = NULL;

// Global flags to tailor execution
const int                         TestRun       = 0;
static int                        DoStateOutput = 0;
static int                        SortIndex     = 0;

// Global variables for other purposes

static char                       progname[MAXPATHLEN];
static double                     LogMinSurvival;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
static char                       parstring[MAX_STR_LEN];
static char                       optstring[MAX_STR_LEN];
static char                       evarstring[MAX_STR_LEN];
#endif


/*
 *====================================================================================================================================
 *              Include header files with generic routines
 *====================================================================================================================================
 */

#include "memory.h"
#include "lifehistory.h"
#include "dopri5.h"

/*==================================================================================================================================*/

int Equation()

{
  int     i, j, b = 0, p, retval = SUCCES;
  int     MaxCohortDim  = CohortDim + 1;
  double  *NextGenMatrix = NULL, *FinalIstateMem = NULL;
  double  R0[PopulationNr], norm;

  /*
   *===========================================================================
   * Set the dimensions
   *===========================================================================
   */
  SetBirthStates(BirthStateNr, Evar);

  MaxStatesAtBirth = 0;
  for (p = 0; p < PopulationNr; p++)
    {
      BirthStateNr[p]  = max(BirthStateNr[p], 1);
      MaxStatesAtBirth = max(MaxStatesAtBirth, BirthStateNr[p]);
      MaxCohortDim     = CohortDim + MaxStatesAtBirth;
      PopDensCohortDim = MaxCohortDim;
    }

  // Allocate the local memory and initialize it to 0
  NextGenMatrix  = calloc(MaxStatesAtBirth*MaxStatesAtBirth, sizeof(double));
  FinalIstateMem = calloc(MaxStatesAtBirth*PopulationNr*MaxCohortDim, sizeof(double));
  retval         = AllocateHeapMemory();

  if ((retval != SUCCES) || (!NextGenMatrix) || (!FinalIstateMem))
    {
      if (FinalIstateMem) free(FinalIstateMem);
      if (NextGenMatrix) free(NextGenMatrix);
      FreeHeapMemory();
      return ReportMemError("Equation");
    }


  /*
   *===========================================================================
   * The life history integration loop. Integration is carried out for each
   * state at birth separately. The states at birth are processed from the
   * highest to the lowest, because not all populations may have an equal
   * number of state at birth.
   *===========================================================================
   */

  retval        = SUCCES;
  DoStateOutput = 0;

#if (defined(OPENMP) && (RFUNCTIONS != 1) && (MFUNCTIONS != 1))                     // If using R-defined or Matlab-defined model this parallelization is impossible
#pragma omp parallel for if (MaxStatesAtBirth > 1)                                  // Only use threading with multiple states at birth
  for (b = 0; b < MaxStatesAtBirth; b++)
    if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;
#else
  for (b = 0; b < MaxStatesAtBirth; b++)
    if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;
#endif

  if (retval == FAILURE)
    {
      if (FinalIstateMem) free(FinalIstateMem);
      if (NextGenMatrix) free(NextGenMatrix);
      return FAILURE;
    }

  retval        = SUCCES;
  DoStateOutput = 1;

  memset(BirthStateMem, 0, PopulationNr*MaxStatesAtBirth*IStateDim*sizeof(double));
  memset(PopDensMem,    0, PopulationNr*MaxStatesAtBirth*PopDensCohortDim*CohortNr*sizeof(double));
  memset(CohortsMem,    0, PopulationNr*MaxStatesAtBirth*sizeof(int));

#ifdef OPENMP
#pragma omp parallel for if (MaxStatesAtBirth > 1)                                  // Only use threading with multiple states at birth
  for (b = 0; b < MaxStatesAtBirth; b++)
    if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;
#else                                                                               // Repeat to stop complaints about ignored #pragma
  for (b = 0; b < MaxStatesAtBirth; b++)
    if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;
#endif

  if (retval == FAILURE)
    {
      if (FinalIstateMem) free(FinalIstateMem);
      if (NextGenMatrix) free(NextGenMatrix);
      return FAILURE;
    }

  /*
   *===========================================================================
   * Compute the final values of the fixed point equation F(y)=0,
   *===========================================================================
   */
  for (p = 0; p < CurPopulationNr; p++)
    {
      if (BirthStateNr[p] == 1)
        {
          R0[p] = FinalIstate(0, p, IStateDim + 1);
          RightEigenvec(p, 0) = Beq[p];
          for (i = 0; i < CohortNr; i++)
            PopDens(p, 0, IStateDim, i) *= Beq[p];
 
        }
      else
        {
          // Notice that in the next generation matrix the first index is the index of the offspring, the second index
          // is the index of the parent. This means that row b of the next generation represent all offspring with state
          // at birht #b produced by all the different parent types with state at birth #j.
          // Hence, the data in Istate[][][] have to be transposed, because here the first index is the parent index
          // whereas the second index is the type of offspring they produce
          for (b = 0; b < BirthStateNr[p]; b++)
            for (j = 0; j < BirthStateNr[p]; j++) NextGenMatrix[b*BirthStateNr[p] + j] = FinalIstate(j, p, IStateDim + 1 + b);

          if (Eigenval(BirthStateNr[p], NextGenMatrix, 0, R0 + p, DOMINANT, &(RightEigenvec(p, 0)), NULL, RHSTOL) != SUCCES)
            {
              ErrorMsg(__FILE__, __LINE__, "Computation of dominant eigenvalue failed for population %d!", p);
              if (FinalIstateMem) free(FinalIstateMem);
              if (NextGenMatrix) free(NextGenMatrix);
              return FAILURE;
            }

          // Scale the right eigenvector such that the elements sum to 1
          for (b = 0, norm = 0; b < BirthStateNr[p]; b++)
            {
              if (fabs(RightEigenvec(p, b)) < Odesolve_Func_Tol)
                RightEigenvec(p, b) = 0.0;
              else
                norm += RightEigenvec(p, b);
            }
          SCAL(BirthStateNr[p], Beq[p]/norm, &(RightEigenvec(p, 0)), 1);

          // Check the validity of the right eigenvector
          for (b = 0; b < BirthStateNr[p]; b++)
            {
              if (RightEigenvec(p, b) < 0)
                {
                  ErrorMsg(__FILE__, __LINE__, "Negative and positive elements in right eigenvector of population %d!", p);
                  if (FinalIstateMem) free(FinalIstateMem);
                  if (NextGenMatrix) free(NextGenMatrix);
                  return FAILURE;
                }
            }
        }
    }


  STDOUT("\n\n\n%31s", "");
  for (i = 0; i < IStateDim; i++) STDOUT("%12s%2d]", "Istate[", i);
  STDOUT("       Survival%15s", "R0");
  for (i = 0; i < InteractDim; i++) STDOUT("        Impact[%2d]", i);

  for (p = 0; p < CurPopulationNr; p++)
    {
      for (b = 0; b < BirthStateNr[p]; b++)
        {
          STDOUT("\nPop. #%2d - Bstate %2d - (Final):", p, b);
          for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", FinalIstate(b, p, i));
          STDOUT("%15.6G", exp(FinalIstate(b, p, IStateDim)));
          STDOUT("%15.6G", ASUM(BirthStateNr[p], FinalIstatePnt(b, p, (IStateDim + 1)), 1));
          for (i = 0; i < InteractDim; i++) STDOUT("%18.6G", FinalIstate(b, p, (IStateDim + 1) + BirthStateNr[p] + i));
        }
      if (BirthStateNr[p] > 1)
        {
          STDOUT("\n\n\n");
          STDOUT("|");
          for (j = 0; j < (15*BirthStateNr[p] - 26)/2; j++) STDOUT("-");
          STDOUT("  Next generation matrix  ");
          for (j = 0; j < (15*BirthStateNr[p] - 25)/2; j++) STDOUT("-");
          STDOUT("|\n");
          for (b = 0; b < BirthStateNr[p]; b++)
            {
              for (j = 0; j < BirthStateNr[p]; j++)
                {
                  STDOUT("%15.6G", NextGenMatrix[b*BirthStateNr[p] + j]);
                }
              STDOUT(" |");
              STDOUT("%15.6G", ASUM(BirthStateNr[p], NextGenMatrix + b*BirthStateNr[p], 1));
              STDOUT("\n");
            }
          for (j = 0; j < BirthStateNr[p]; j++) STDOUT("  -------------");
          STDOUT("\n");
          for (j = 0; j < BirthStateNr[p]; j++) STDOUT("%15.6G", ASUM(BirthStateNr[p], NextGenMatrix + j, BirthStateNr[p]));

          STDOUT("\n\n\nEig(M) : %12.6E", R0[p]);

          STDOUT("\n\n\nStable birth distribution\n");
          for (b = 0; b < BirthStateNr[p]; b++) STDOUT("%14.6G  ", RightEigenvec(p, b));
        }
    }

  STDOUT("\n\n");
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  if (FinalIstateMem) free(FinalIstateMem);
  if (NextGenMatrix) free(NextGenMatrix);

  return retval;
}


/*==================================================================================================================================*/

int DefineOutput(double *x, double *output)

{
  // Need this stub here to prevent the linking stage to fail
  return SUCCES;
}


/*
 *====================================================================================================================================
 *  Implementation of the routine that computes the individual life history
 *====================================================================================================================================
 */

void ComputeLifeHistory(const int argc, char **argv)
{
  register int  i;
  char          csbname[MAX_STR_LEN];
  struct stat   buffer;

#if defined(R_PACKAGE)
  STDOUT("\n");
#else
  fprintf(stderr, "\n");
#endif

  i = 0;
  while (1)
    {
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      sprintf(csbname, "%s-%s-%04d.mat", progname, "IND", i);
#else
      sprintf(csbname, "%s-%s-%04d.csb", progname, "IND", i);
#endif
      if (stat(csbname, &buffer)) break;
      i++;
    }
  sprintf(runname, "%s-%s-%04d", progname, "IND", i);

  Equation();
  WriteStateToFile(2);

  STDOUT("\n");
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif
  FreeHeapMemory();

  return;
}


/*==================================================================================================================================*/

void InitialiseVars(void)

{
  int ii;
#if (defined(_MSC_VER) && (_MSC_VER < 1500)) || (defined(R_PACKAGE) && defined(_WIN32))
  (void)_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  // Initialize some variables
  strcpy(runname, "");

  // Get the machine precisions
  epsMach = dlamch("Epsilon");

  CurPopulationNr   = PopulationNr;
  CohortDim         = IStateDim + 1 + InteractDim;
  Bifparone         = -1;
  
  CurveType         = IND;

  eVarPntr      = Evar;
  parPntr       = parameter;

  for (ii = 0; ii < PopulationNr; ii++) Beq[ii] = 1.0;

  // The following variables can be modified by the user with optional #define statements
  LogMinSurvival      = log(MIN_SURVIVAL);
  CohortNr            = COHORT_NR + Stages;
  
  Odesolve_Init_Step  = ODESOLVE_INIT_STEP;
  Odesolve_Fixed_Step = ODESOLVE_FIXED_STEP;
  Odesolve_Min_Step   = ODESOLVE_MIN_STEP;
  Odesolve_Max_Step   = ODESOLVE_MAX_STEP;
  Odesolve_Abs_Err    = ODESOLVE_ABS_ERR;
  Odesolve_Rel_Err    = ODESOLVE_REL_ERR;
  Odesolve_Func_Tol   = ODESOLVE_FUNC_TOL;

  return;
}


/*
 *====================================================================================================================================
 *  UNIX shell interface function main() and supporting functions.
 *====================================================================================================================================
 */

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)

static void Usage(char *progname)

{
  int   i;
  char  tmpstr[MAX_STR_LEN];
  char  desc[MAX_STR_LEN];
  char  varstr[MAX_STR_LEN];

  strcpy(desc, "Aim:\tSimulating the life history of a single individual of a structured population given a specific environmental state");
  strcpy(varstr, "");
  for (i = 0; i < EnvironDim; i++)
    {
      sprintf(tmpstr, " E[%d]", i);
      strcat(varstr, tmpstr);
    }
  strcat(varstr, " [");
  for (i = 0; i < PopulationNr; i++)
    {
      sprintf(tmpstr, " b[%d]", i);
      strcat(varstr, tmpstr);
    }  
  strcat(varstr, "]");

  fprintf(stderr, "Usage:\t%s [<options>]%s", progname, varstr);
  fprintf(stderr, "\n\n%s\n\n", desc);
  fprintf(stderr, "Possible options are:\n\n");
  fprintf(stderr, "\t-isort   <index>  : Index of i-state variable to use as ruling variable for sorting the structured populations\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "\n%s, Copyright (C) 2015, Andre M. de Roos, University of Amsterdam\n\n", progname);
  fprintf(stderr, "This program comes with ABSOLUTELY NO WARRANTY; without even the implied warranty of\n");
  fprintf(stderr, "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public\n");
  fprintf(stderr, "License (<http://www.gnu.org/philosophy/why-not-lgpl.html>) for more details\n\n");

  exit(1);                                                                          // Only executed when in command-line mode

  return;
}


/*==================================================================================================================================*/

int main(int argc, char **argv)

{
  int   evarnr  = 0, tmpint;
  char  **argpnt1 = NULL;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = ENVIRON_DIM;
  InteractDim  = INTERACT_DIM;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

  if (argc < 2) Usage(argv[0]);

  argpnt1 = argv + 1;

  while (*argpnt1)
    {
      if (!strcmp(*argpnt1, "-?") || !strcmp(*argpnt1, "--help"))
        {
          Usage(argv[0]);
        }
      else if (!strcmp(*argpnt1, "-isort"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo index of i-state variable specified for argument -isort!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= IStateDim))
            {
              fprintf(stderr, "\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, IStateDim);
              Usage(argv[0]);
            }
          SortIndex = tmpint;
        }
      else if ((!strncmp(*argpnt1, "--", 2)))
        {
          fprintf(stderr, "\nUnknown command line option: %s\n", *argpnt1);
          Usage(argv[0]);
        }
      else if ((!strncmp(*argpnt1, "-", 1)) && isalpha(*(*argpnt1 + 1)))
        {
          fprintf(stderr, "\nUnknown command line option: %s\n", *argpnt1);
          Usage(argv[0]);
        }
      else
        {
          if (evarnr < EnvironDim)
            Evar[evarnr] = atof(*argpnt1);
          else if (evarnr < (EnvironDim + PopulationNr))
            Beq[evarnr - EnvironDim] = atof(*argpnt1);
          else
            {
              fprintf(stderr, "\nToo many command-line arguments.\n");
              fprintf(stderr, "The vector with environmental variables and (optional) birth rates should have length %d or %d.\n\n", 
                      EnvironDim, EnvironDim + PopulationNr);
              Usage(argv[0]);
            }
          evarnr++;
        }
      argpnt1++;
    }

  strcpy(progname, argv[0]);
  progname[strlen(progname) - 3] = '\0';                                            // Cut off the 'ind' appendix

  ComputeLifeHistory(argc, argv);

  return 0;
}


/*
 *====================================================================================================================================
 *	Matlab interface and exit/cleanup function
 *====================================================================================================================================
 */

#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)

#if defined(MATLAB_MEX_FILE)
#include "mat.h"

extern MATFile                    *pmat;
#endif

static void CloseStreams(void)
{
  if (pmat) matClose(pmat);

  FreeHeapMemory();

  return;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t                          nrows, ncols;
  double                          *inputVals, tmpdouble;
  int                             tmpint, rind, irhs;
  const mxArray                   *cell_element_ptr;
  char                            optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN], varstr[MAX_STR_LEN];

  mwIndex                         i, j;
  size_t                          total_num_of_cells, buflen;
  int                             status;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = ENVIRON_DIM;
  InteractDim  = INTERACT_DIM;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

  // check for proper number of arguments
  if (nrhs != 3)
    mexErrMsgIdAndTxt("MATLAB:PSPMind:nrhs", "\nIncorrect number of command-line arguments.\n\nUse: %s(%s)\n\n%22s: %s\n%22s: %s\n%22s: %s",
                      mexFunctionName(), "<environment>, <parameters>, <options>", "<environment>",
                      "Values of the environmental variables", "<parameters>",
                      "Array of parameter values to use (empty array or of same length as parameter array)", "<options>", "Possible option: isort");

  // check for proper number of output variables
  if (nlhs != 1) mexErrMsgIdAndTxt("MATLAB:PSPMind:nlhs", "\nA single output argument is required!");

#if (MFUNCTIONS == 1)
  Minterface_Init();
#endif

  //============================== Process the options argument ======================================================================
  // Extract the contents of MATLAB cell into the C array
  irhs = 2;
  if (!mxIsCell(prhs[irhs])) mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options", "\nOptions should be specified as a cell array!");

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(optstring, "{");
  for (i = 0; i < total_num_of_cells; i++)
    {
      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optname, buflen);
      if (strcmp(optname, "isort")) mexErrMsgIdAndTxt("MATLAB:PSPMind:options", "\nIllegal option %s!", optname);

      if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMind:options", "\nNo value specified for option %s!", optname);

      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optval, buflen);
      if (status) mexErrMsgIdAndTxt("MATLAB:PSPMind:options", "\nError in retrieving value for option %s!", optname);

      // optstring still equal to "{"
      if (strlen(optstring) == 1)
        strcat(optstring, "'");
      else
        strcat(optstring, ", '");
      strcat(optstring, optname);
      strcat(optstring, "', '");
      strcat(optstring, optval);
      strcat(optstring, "'");

      tmpint = atoi(optval);
      if ((tmpint < 0) || (tmpint >= IStateDim))
        mexErrMsgIdAndTxt("MATLAB:PSPMind:options",
                          "\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!", tmpint,
                          IStateDim);
      SortIndex = tmpint;
    }
  strcat(optstring, "}");

  //============================== Process the parameters argument ===================================================================

  irhs  = 1;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols == ParameterNr) && (nrows == 1))
    memcpy(parameter, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:PSPMind:parameters", "\nParameter argument ignored as it is not a row vector of length %d!", ParameterNr);

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(parstring, "[");
  for (i = 0; i < total_num_of_cells; i++)
    {
      if (i) strcat(parstring, " ");
      memcpy(&tmpdouble, mxGetPr(prhs[irhs]) + i, mxGetElementSize(prhs[irhs]));
      sprintf(tmpstr, "%.6G", tmpdouble);
      strcat(parstring, tmpstr);
    }
  strcat(parstring, "]");

  //============================== Process the initial point argument ================================================================

  irhs = 0;
  memset((void *)Evar, 0, EnvironDim*sizeof(double));
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols == EnvironDim) && (nrows == 1))
    memcpy(Evar, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  else if ((ncols == (EnvironDim + PopulationNr)) && (nrows == 1))
    {
      memcpy(Evar, mxGetPr(prhs[irhs]), EnvironDim*mxGetElementSize(prhs[irhs]));
      memcpy(Beq,  mxGetPr(prhs[irhs]) + EnvironDim, PopulationNr*mxGetElementSize(prhs[irhs]));
    }
  else if (ncols)
    mexErrMsgIdAndTxt("MATLAB:PSPMind:environment", "\nEnvironmental values argument has to be a row vector of length %d or %d!", 
                      EnvironDim, EnvironDim + PopulationNr);

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(evarstring, "[");
  for (i = 0; i < total_num_of_cells; i++)
    {
      if (i) strcat(evarstring, " ");
      memcpy(&tmpdouble, mxGetPr(prhs[irhs]) + i, mxGetElementSize(prhs[irhs]));
      sprintf(tmpstr, "%.6G", tmpdouble);
      strcat(evarstring, tmpstr);
    }
  strcat(evarstring, "]");

  //============================= Get the program name ===============================================================================
  // Get the name of the mex file
  strcpy(progname, mexFunctionName());
  progname[strlen(progname) - 3] = '\0';                                            // Cut off the 'ind' appendix

  mexAtExit(CloseStreams);

  // call the computational routine
  ComputeLifeHistory(0, NULL);
  FreeHeapMemory();

  plhs[0] = mxCreateString(runname);

#if (MFUNCTIONS == 1)
  Minterface_End();
#endif

  return;
}


/*
 *====================================================================================================================================
 *	R interface function
 *====================================================================================================================================
 */

#elif defined(R_PACKAGE)

SEXP PSPMind(SEXP moduleName, SEXP evarVals, SEXP parVals, SEXP optVals)

{
  int  i, ncols, tmpint;
  char optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN];
  SEXP resfil;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = ENVIRON_DIM;
  InteractDim  = INTERACT_DIM;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

#if (RFUNCTIONS == 1)
  Rinterface_Init();
#endif

  //============================== Process the options argument ======================================================================

  if (length(optVals) && (!isString(optVals)))
    error("\nOptions should either be specified as NULL or as a vector of strings\n\n");
  strcpy(optstring, "");

  ncols =length(optVals);
  for (i = 0; i < ncols; i++)
    {
      strcpy(optname, CHAR(STRING_ELT(optVals, i)));
      if (strcmp(optname, "isort")) error("\nIllegal option %s!\n\n", optname);

      if (!(++i < ncols)) error("\nNo value specified for option %s!\n\n", optname);

      strcpy(optval, CHAR(STRING_ELT(optVals, i)));
      
      // optstring still equal to ""
      if (!strlen(optstring))
        strcat(optstring, "c(\"");
      else
        strcat(optstring, ", \"");
      strcat(optstring, optname);
      strcat(optstring, "\", \"");
      strcat(optstring, optval);
      strcat(optstring, "\"");

      tmpint = atoi(optval);
      if ((tmpint < 0) || (tmpint >= IStateDim))
        error("\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
              IStateDim);
      SortIndex = tmpint;
    }

  if (strlen(optstring))
    strcat(optstring, ")");
  else
    strcpy(optstring, "NULL");

  //============================== Process the parameters argument ===================================================================

  ncols = length(parVals);
  if (isReal(parVals) && (ncols == ParameterNr))
    memcpy(parameter, REAL(parVals), ncols*sizeof(double));
  else if (ncols)
    warning("\nParameter argument ignored as it is not a row vector of appropriate length\n\n");

  if (ncols)
    {
      strcpy(parstring, "c(");
      for (i = 0; i < ncols; i++)
        {
          if (i) strcat(parstring, ", ");
          sprintf(tmpstr, "%.6G", REAL(parVals)[i]);
          strcat(parstring, tmpstr);
        }
      strcat(parstring, ")");
    }
  else
    strcpy(parstring, "NULL");

  //============================== Process the environmental values argument =========================================================

  ncols = length(evarVals);
  if (isReal(evarVals) && (ncols == EnvironDim))
    memcpy(Evar, REAL(evarVals), EnvironDim*sizeof(double));
  else if (isReal(evarVals) & (ncols == (EnvironDim + PopulationNr)))
    {
      memcpy(Evar, REAL(evarVals), EnvironDim*sizeof(double));
      memcpy(Beq,  REAL(evarVals) + EnvironDim, PopulationNr*sizeof(double));
    }
  else
    error("\nVector argument with environmental variables and birth rates has to be a row vector of appropriate length\n\n");

  if (ncols)
    strcpy(evarstring, "c(");
  else
    strcpy(evarstring, "NULL");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(evarstring, ", ");
      sprintf(tmpstr, "%.6G", REAL(evarVals)[i]);
      strcat(evarstring, tmpstr);
    }
  if (ncols) strcat(evarstring, ")");

  //============================= Get the program name ===============================================================================
  // Get the name of the module file
  if (isString(moduleName) && (length(moduleName) == 1))
    strcpy(progname, CHAR(STRING_ELT(moduleName, 0)));
  else
    error("\nModel name argument must be a single string\n\n");

  // call the computational routine
  ComputeLifeHistory(0, NULL);

  FreeHeapMemory();

#if (RFUNCTIONS == 1)
  Rinterface_End();
#endif

  PROTECT(resfil = mkString(runname));

  UNPROTECT(1);

  return resfil;
}

/*==================================================================================================================================*/
#endif
