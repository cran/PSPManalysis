/***
  NAME
    PSPMevodyn

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

    Last modification: AMdR - Jan 19, 2023
***/

#define PSPMEVODYN                1

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
#define Birthrate(p)              (birthRatePntr[p])


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
 *		Definition of problem dimensions
 *====================================================================================================================================
 */

#undef PULSED

/*
 *====================================================================================================================================
 *		Definition of global variables and parameters
 *====================================================================================================================================
 */

// Global dimension variables
static int                        BirthStateNr[POPULATION_NR];
static int                        MaxPntDim;

// These are the variables to solve for
static double                     Evar[ENVIRON_DIM];
static double                     EnvEquiCondition[ENVIRON_DIM];

static int                        EnvTrivEqui[ENVIRON_DIM];

static double                     R0[POPULATION_NR];
static double                     Beq[POPULATION_NR];
static double                     InteractVars[POPULATION_NR][INTERACT_DIM];

static int                        PopTrivEqui[POPULATION_NR];
static int                        R0ResIndex[POPULATION_NR];

// Global variables to hold variables shared among routines
static int                        LastMemAllocated = 0;
static int                        pntdim;

// Global pointers into the heap
static double                     *RightEigenvecMem = NULL;

#if (FULLSTATEOUTPUT > 0)
static double                     *BirthStateMem = NULL;
static double                     *PopDensMem    = NULL;
static int                        *CohortsMem    = NULL;

static double                     *CohortLimitMem = NULL;
#if (FULLSTATEOUTPUT == 1)
static double                     CohortMin[POPULATION_NR], CohortMax[POPULATION_NR];
#endif
#endif

// Global flags to tailor execution
static int                        TestRun       = 0;
static int                        DoStateOutput = 0;
static int                        SortIndex     = 0;
static int                        ReportLevel   = 1;
static int                        evoParsIndex[PARAMETER_NR];
static int                        evoPopIndex[PARAMETER_NR];

// Global variables for other purposes

static char                       progname[MAXPATHLEN];
static double                     LogMinSurvival;
static double                     point[ENVIRON_DIM + POPULATION_NR];
static double                     pntmin[ENVIRON_DIM + POPULATION_NR];
static double                     evoParMin[PARAMETER_NR];
static double                     evoParMax[PARAMETER_NR];
static double                     Covariances[PARAMETER_NR*PARAMETER_NR];
static double                     evoTime, evoMaxTime;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
static char                       parstring[MAX_STR_LEN];
static char                       optstring[MAX_STR_LEN];
static char                       pntstring[MAX_STR_LEN];
static char                       evotimestring[MAX_STR_LEN];
static char                       evoparstring[MAX_STR_LEN];
static char                       covmatstring[MAX_STR_LEN];
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

int Equation(double *argument, double *result)

{
  int     indx = 0, e, i, j, b = 0, p, retval = SUCCES;
  int     MaxCohortDim = CohortDim + 1 + InteractDim;
  double  norm;
  double  *NextGenMatrix = NULL, *FinalIstateMem = NULL;
  double  *respntr;

  /*
   *===========================================================================
   * Map current estimate of solution to global variables
   *===========================================================================
   */
  for (e = 0; e < EnvironDim; e++)
    {
      if (EnvTrivEqui[e])
        Evar[e] = 0.0;
      else
        {
          Evar[e] = argument[indx]*pnt_scale[indx];
          indx++;
        }
    }

  for (p = 0; p < CurPopulationNr; p++)
    {
      if (PopTrivEqui[p])
        Beq[p] = 0.0;
      else
        {
          Beq[p] = argument[indx]*pnt_scale[indx];
          indx++;
        }
    }

  /*
   *===========================================================================
   * Set the dimensions
   *===========================================================================
   */
  CurPopulationNr = PopulationNr;
  SetBirthStates(BirthStateNr, Evar);

  MaxStatesAtBirth = 0;
  for (p = 0; p < CurPopulationNr; p++)
    {
      BirthStateNr[p]  = max(BirthStateNr[p], 1);
      MaxStatesAtBirth = max(MaxStatesAtBirth, BirthStateNr[p]);
      MaxCohortDim     = CohortDim + MaxStatesAtBirth + InteractDim;
    }

  // Allocate the local memory and initialize it to 0
  NextGenMatrix  = calloc(MaxStatesAtBirth*MaxStatesAtBirth, sizeof(double));
  FinalIstateMem = calloc(MaxStatesAtBirth*CurPopulationNr*MaxCohortDim, sizeof(double));
  retval         = AllocateHeapMemory();

  if ((retval != SUCCES) || (!NextGenMatrix) || (!FinalIstateMem))
    {
      if (FinalIstateMem) free(FinalIstateMem);
      if (NextGenMatrix) free(NextGenMatrix);
      FreeHeapMemory();
      return ReportMemError("Equation");
    }

#if (FULLSTATEOUTPUT > 0)
  memset(BirthStateMem, 0, CurPopulationNr*MaxStatesAtBirth*IStateDim*sizeof(double));
  memset(PopDensMem, 0, CurPopulationNr*MaxStatesAtBirth*PopDensCohortDim*CohortNr*sizeof(double));
  memset(CohortsMem, 0, CurPopulationNr*MaxStatesAtBirth*sizeof(int));
#endif

  /*
   *===========================================================================
   * Testing output
   *===========================================================================
   */

  if (TestRun)
    {
      STDOUT("\n");
      for (e = 0; e < EnvironDim; e++) STDOUT("\nEnvironment variable #%d:     %15.6G", e, Evar[e]);
      for (p = 0; p < CurPopulationNr; p++) STDOUT("\nBirth rate of population #%d: %15.6G", p, Beq[p]);
      for (i = 0; i < evoParsDim; i++) STDOUT("\nParameter #%d:                %15.6G", evoParsIndex[i], parameter[evoParsIndex[i]]);

      STDOUT("\n\n\n%40s", "");
      for (i = 0; i < IStateDim; i++) STDOUT("%12s%2d]", "Istate[", i);

      STDOUT("       Survival%15s", "R0");
      for (i = 0; i < InteractDim; i++) STDOUT("     Impact[%2d]", i);
      STDOUT("\n");
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  /*
   *===========================================================================
   * The life history integration loop. Integration is carried out for each
   * state at birth separately. The states at birth are processed from the
   * highest to the lowest, because not all populations may have an equal
   * number of state at birth.
   *===========================================================================
   */

  retval = SUCCES;

  int openmpDone = 0;
#if (defined(OPENMP) && (RFUNCTIONS != 1) && (MFUNCTIONS != 1))                     // If using R-defined or Matlab-defined model this parallelization is impossible
  if (!TestRun)
    {
#pragma omp parallel for if (MaxStatesAtBirth > 1)                                  // Only use threading with multiple states at birth
      for (b = 0; b < MaxStatesAtBirth; b++)
        if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;
      openmpDone = 1;
    }
#endif
  if (!openmpDone)
    {
      for (b = 0; b < MaxStatesAtBirth; b++)
        if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;
    }

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
          RightEigenvec(p, 0) = 1.0;
          for (i               = 0; i < InteractDim; i++)
            InteractVars[p][i] = Beq[p]*FinalIstate(0, p, CohortDim + 1 + i);       // Multiply all measures with the birth rate
#if (FULLSTATEOUTPUT > 0)
          CohortLimit(0, p) = (FinalIstate(0, p, SortIndex) - BirthState(p, 0, SortIndex))/COHORT_NR;
#endif
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
          SCAL(BirthStateNr[p], 1.0/norm, &(RightEigenvec(p, 0)), 1);

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

          for (i = 0; i < InteractDim; i++)
            {
              InteractVars[p][i] = 0.0;
              for (b = 0; b < BirthStateNr[p]; b++) InteractVars[p][i] += RightEigenvec(p, b)*FinalIstate(b, p, CohortDim + BirthStateNr[p] + i);
              InteractVars[p][i] *= Beq[p];                                         // Multiply all measures with the birth rate
            }
#if (FULLSTATEOUTPUT == 1)
          CohortMin[p] = SAFETY*DBL_MAX;
          CohortMax[p] = -SAFETY*DBL_MAX;
          for (b = 0; b < BirthStateNr[p]; b++)
            {
              if (RightEigenvec(p, b))
                {
                  CohortMin[p] = min(CohortMin[p], BirthState(p, b, SortIndex));
                  CohortMax[p] = max(CohortMax[p], FinalIstate(b, p, SortIndex));
                }
            }
          for (b = 0; b < BirthStateNr[p]; b++) CohortLimit(b, p) = (CohortMax[p] - CohortMin[p])/COHORT_NR;
#elif (FULLSTATEOUTPUT == 2)
          for (b = 0; b < BirthStateNr[p]; b++) CohortLimit(b, p) = (FinalIstate(b, p, SortIndex) - BirthState(p, b, SortIndex))/COHORT_NR;
#endif
        }
    }

  EnvEqui(Evar, InteractVars, EnvEquiCondition);

  if (TestRun)
    {
      STDOUT("\n\n\n%31s", "");
      for (i = 0; i < IStateDim; i++) STDOUT("%12s%2d]", "Istate[", i);
      STDOUT("       Survival%15s", "R0");
      for (i = 0; i < InteractDim; i++) STDOUT("   InteractVar[%2d]", i);

      for (p = 0; p < CurPopulationNr; p++)
        {
          for (b = 0; b < BirthStateNr[p]; b++)
            {
              STDOUT("\nPop. #%2d - Bstate %2d - (Final):", p, b);
              for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", FinalIstate(b, p, i));
              STDOUT("%15.6G", exp(FinalIstate(b, p, IStateDim)));
              STDOUT("%15.6G", SUM(BirthStateNr[p], FinalIstatePnt(b, p, CohortDim), 1));
              for (i = 0; i < InteractDim; i++) STDOUT("%18.6G", Beq[p]*FinalIstate(b, p, CohortDim + BirthStateNr[p] + i));
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
                  STDOUT("%15.6G", SUM(BirthStateNr[p], NextGenMatrix + b*BirthStateNr[p], 1));
                  STDOUT("\n");
                }
              for (j = 0; j < BirthStateNr[p]; j++) STDOUT("  -------------");
              STDOUT("\n");
              for (j = 0; j < BirthStateNr[p]; j++) STDOUT("%15.6G", SUM(BirthStateNr[p], NextGenMatrix + j, BirthStateNr[p]));

              STDOUT("\n\n\nEig(M) : %12.6E", R0[p]);

              STDOUT("\n\n\nStable birth distribution\n");
              for (b = 0; b < BirthStateNr[p]; b++) STDOUT("%14.6G  ", RightEigenvec(p, b));
            }
        }

      STDOUT("\n\n");
      for (e = 0; e < EnvironDim; e++) STDOUT("\nEquilibrium condition environment variable %2d:\t\t%18.6G", e, EnvEquiCondition[e]);
      STDOUT("\n");
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  /*
   * Now assign all the result values that determine the fixed point F(y) = 0.
   *
   * Remember the first environment variable in evolutionary dynamic simulation has index 0 in pnt_scale
   */
  indx    = 0;
  respntr = result;
  for (e = 0; e < EnvironDim; e++)
    {
      if (EnvTrivEqui[e]) continue;
      if (EnvironmentType[e] == GENERALODE)
        *respntr = EnvEquiCondition[e]/pnt_scale[indx];
      else if (EnvironmentType[e] == POPULATIONINTEGRAL)
        *respntr = (Evar[e] - EnvEquiCondition[e])/pnt_scale[indx];
      else
        *respntr = EnvEquiCondition[e];
      respntr++;
      indx++;
    }
  for (p = 0; p < CurPopulationNr; p++)
    {
      if (PopTrivEqui[p]) continue;
      *respntr = R0[p] - 1.0;
      respntr++;
    }

  if (FinalIstateMem) free(FinalIstateMem);
  if (NextGenMatrix) free(NextGenMatrix);
  return retval;
}


/*==================================================================================================================================*/

int DefineOutput(double *x, double *output)

{
  int     outnr = 0, i, j;
  double  result[MaxPntDim];

#if (FULLSTATEOUTPUT > 0)
  DoStateOutput = 1;
  Equation(x, result);
  DoStateOutput = 0;
  WriteStateToFile(FULLSTATEOUTPUT);
#else
  Equation(x, result);
#endif

  output[outnr++] = evoTime;
  for (i = 0; i < EnvironDim; i++) output[outnr++] = Evar[i];
  for (i = 0; i < PopulationNr; i++) output[outnr++] = Beq[i];
  for (i = 0; i < evoParsDim; i++) output[outnr++] = parameter[evoParsIndex[i]];

  // Also add all interaction variables
  for (i = 0; i < PopulationNr; i++)
    for (j = 0; j < InteractDim; j++) output[outnr++] = InteractVars[i][j];

  for (i = 0; i < EnvironDim; i++)
    if (EnvironmentType[i] == PERCAPITARATE) output[outnr++] = EnvEquiCondition[i];

  for (i = 0; i < CurPopulationNr; i++) output[outnr++] = R0[i];

  // Should we indeed add how accurate the solution was?
  output[outnr++] = NRM2(pntdim, result, 1);

  return outnr;
}


/*
 *====================================================================================================================================
 *  Implementation of the routine that computes the entire curve
 *====================================================================================================================================
 */

void ComputeCurve(const int argc, char **argv)
{
  register int  i, j, colnr;
  int           pntnr = 0, outmax, CurveEnd = 0, Converged, AtBoundary[ParameterNr], TotalAtBoundary;
  int           cycles, last = 1, retval = 0;
  double        oldpars[ParameterNr], evoparvec[ParameterNr], selectdiff[ParameterNr];
  double        oldpoint[MaxPntDim], evopntvec[MaxPntDim];
  double        Jac[MaxPntDim*MaxPntDim], JacCopy[MaxPntDim*MaxPntDim], dFdp[ParameterNr*MaxPntDim];
  double        y[MaxPntDim], rhs[MaxPntDim], pardif, oldcurvestep;
  char          csbname[3*MAX_STR_LEN], errname[3*MAX_STR_LEN], outname[3*MAX_STR_LEN];
  char          tmpstr[MAX_STR_LEN];
  struct stat   buffer;

#if defined(R_PACKAGE)
  STDOUT("\n");
#else
  fprintf(stderr, "\n");
#endif
  (void)SetScales(point, pntdim);

  if (TestRun)
    {
      double rhs[MaxPntDim], rhsnorm;

#if (defined(R_PACKAGE))
      STDOUT("\n\nExecuting : ");
      STDOUT("PSPMevodyn(\"%s\", %s, %s, %s, %s, %s, %s)", progname, pntstring, evotimestring, evoparstring, covmatstring, parstring, optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      STDOUT("\n\nExecuting : ");
      STDOUT("PSPMevodyn('%s', %s, %s, %s, %s, %s, %s)", progname, pntstring, evotimestring, evoparstring, covmatstring, parstring, optstring);
#else
      STDOUT("Executing : ");
      for (i = 0; i < argc; i++) STDOUT("%s ", argv[i]);
#endif
      STDOUT("\n\n");

      STDOUT("Parameter values  : \n");
      for (i = 0; i < ParameterNr; i++)
        {
          if (!(i % 3)) STDOUT("\n");
          STDOUT("\t%-10s:", parameternames[i]);
          STDOUT("  %-13G", parameter[i]);
        }
      STDOUT("\n");
      fflush(NULL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif

      errfile = stderr;
      Equation(point, rhs);

      // Compute rhsnorm
      rhsnorm = NRM2(pntdim, rhs, 1);
      rhsnorm = rhsnorm/(1.0 + rhsnorm);

      STDOUT("\n\nNorm of fixed point conditions: %18.6G", rhsnorm);

      STDOUT("\n\n");
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif

      FreeHeapMemory();
      return;
    }

  if (strlen(runname))
    {
      snprintf(errname, sizeof(errname), "%s.err", runname);
      snprintf(outname, sizeof(outname), "%s.out", runname);
    }
  else
    {
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)
      strcpy(progname, argv[0]);
      progname[strlen(progname) - 6] = '\0';                                        // Cut off the 'evodyn' appendix
#endif
      i = 0;
      while (1)
        {
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
          snprintf(csbname, sizeof(csbname), "%s-%s-%04d.mat", progname, "EVODYN", i);
#else
          snprintf(csbname, sizeof(csbname), "%s-%s-%04d.csb", progname, "EVODYN", i);
#endif
          snprintf(errname, sizeof(errname), "%s-%s-%04d.err", progname, "EVODYN", i);
          snprintf(outname, sizeof(outname), "%s-%s-%04d.out", progname, "EVODYN", i);
          if (stat(csbname, &buffer) && stat(errname, &buffer) && stat(outname, &buffer)) break;
          i++;
        }
      snprintf(runname, sizeof(runname), "%s-%s-%04d", progname, "EVODYN", i);
    }

  errfile = fopen(errname, "w");
  outfile = fopen(outname, "w");

  if (outfile)
    {
      fprintf(outfile, "#\n# Executing : ");

#if defined(R_PACKAGE)
      fprintf(outfile, "PSPMevodyn(\"%s\", %s, %s, %s, %s, %s, %s)", progname, pntstring, evotimestring, evoparstring, covmatstring, parstring,
              optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      fprintf(outfile, "PSPMevodyn('%s', %s, %s, %s, %s, %s, %s)", progname, pntstring, evotimestring, evoparstring, covmatstring, parstring,
              optstring);
#else
      for (i = 0; i < argc; i++) fprintf(outfile, "%s ", argv[i]);
#endif
      fprintf(outfile, "\n#\n");

      fprintf(outfile, "# Parameter values  : \n#");
      for (i = 0; i < ParameterNr; i++)
        {
          if (!(i % 3)) fprintf(outfile, "\n# ");
          fprintf(outfile, "\t%-10s:", parameternames[i]);
          fprintf(outfile, "  %-13G", parameter[i]);
        }
      fprintf(outfile, "\n#\n");
      
      for (i = 0; i < evoParsDim; i++)
        {
          fprintf(outfile, "# Index and name of evolution parameter #%d                       : %d (%s)\n", i + 1, evoParsIndex[i],
                  parameternames[evoParsIndex[i]]);
          fprintf(outfile, "# Index of structured population for which parameter #%d evolves  : %d\n", i + 1, evoPopIndex[i]);
        }
      fprintf(outfile, "#\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
      colnr = 0;
#else
      colnr = 1;
#endif
      snprintf(tmpstr, sizeof(tmpstr), "%d:Evol.time", colnr++);
      fprintf(outfile, "#%13s", tmpstr);
      for (i = 0; i < EnvironDim; i++)
        {
          snprintf(tmpstr, sizeof(tmpstr), "%d:E[%d]", colnr++, i);
          fprintf(outfile, "%16s", tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          snprintf(tmpstr, sizeof(tmpstr), "%d:b[%d]", colnr++, i);
          fprintf(outfile, "%16s", tmpstr);
        }
      for (i = 0; i < evoParsDim; i++)
        {
          snprintf(tmpstr, sizeof(tmpstr), "%d:%s", colnr++, parameternames[evoParsIndex[i]]);
          fprintf(outfile, "%16s", tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        for (j = 0; j < InteractDim; j++)
          {
            snprintf(tmpstr, sizeof(tmpstr), "%d:I[%d][%d]", colnr++, i, j);
            fprintf(outfile, "%16s", tmpstr);
          }
      for (i = 0; i < EnvironDim; i++)
        if (EnvironmentType[i] == PERCAPITARATE)
          {
            snprintf(tmpstr, sizeof(tmpstr), "%d:pcgE[%d]", colnr++, i);
            fprintf(outfile, "%16s", tmpstr);
          }
      for (i = 0; i < CurPopulationNr; i++)
        {
          snprintf(tmpstr, sizeof(tmpstr), "%d:R0[%d]", colnr++, i);
          fprintf(outfile, "%16s", tmpstr);
        }
      snprintf(tmpstr, sizeof(tmpstr), "%d:RHS norm\n", colnr++);
      fprintf(outfile, "%17s", tmpstr);
      fflush(outfile);
    }

  // Continue the curve
  while (evoTime < evoMaxTime)
    {
      // Compute fixed point with new varied parameter
      cycles       = 0;
      oldcurvestep = curvestep;
      while (Stepreduce <= MAX_STEPREDUCE)
        {
          retval = FindPoint(pntdim, point, NULL, NULL, DYTOL, RHSTOL, MAXITER, Equation);

          if (retval == SUCCES) break;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif

          NumProcError(__FILE__, __LINE__, retval);

          if (!pntnr)
            {
              ErrorMsg(__FILE__, __LINE__, "No convergence while locating first point of curve");
              STDOUT("\n");
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
              FreeHeapMemory();
              return;
            }
          else if (!Stepchange)                                                     // If unsuccesfull while repeating point return
            {
              ErrorMsg(__FILE__, __LINE__, "Failed to locate a solution point after scaling");
              STDOUT("\n");
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
              FreeHeapMemory();
              return;
            }

          // Generate prediction of solution point with smaller step size
          cycles++;
          Stepreduce *= 2;

          COPY(ParameterNr, oldpars, 1, parameter, 1);
          COPY(pntdim, oldpoint, 1, point, 1);

          AXPY(ParameterNr, (curvestep/Stepreduce), evoparvec, 1, parameter, 1);
          AXPY(pntdim, (curvestep/Stepreduce), evopntvec, 1, point, 1);

          ReportMsg("\n\nPrediction :\t");
          for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
          for (i = 0; i < evoParsDim; i++) ReportMsg("%16.8E  ", parameter[evoParsIndex[i]]);
          ReportMsg("\n");
        }

      // If unsuccesfull return
      if (retval != SUCCES)
        {
          ErrorMsg(__FILE__, __LINE__, "Failed to locate a solution point");
          break;
        }

      // Increase step when both this and previous point were located with
      // the current step size
      if ((last && !cycles) && (Stepreduce > 1)) Stepreduce /= 2;
      last = (!cycles);

      // Scale the point vector anew if necessary and redo the current point
      if ((retval = SetScales(point, pntdim)))
        {
          ReportMsg("\n\nVariable %d rescaled!\n", retval);
          Stepchange = 0;
        }
      // Otherwise generate output and predict new point on the curve
      else
        {
          evoTime += oldcurvestep/Stepreduce;
          curvestep = min(curvestep, Maxcurvestep);

          // Report on located point
          ReportMsg("\nNew point :\t");
          ReportMsg("%16f  ", evoTime);
          for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
          for (i = 0; i < evoParsDim; i++) ReportMsg("%16.8E  ", parameter[evoParsIndex[i]]);
          ReportMsg("\n");

          if (((pntnr + 1) % ReportLevel) == 0)
            {
#if (defined(R_PACKAGE))
              STDOUT("%9f", evoTime);
              for (i = 0; i < pntdim; i++) STDOUT(",%15.8E", point[i]*pnt_scale[i]);
              for (i = 0; i < evoParsDim; i++) STDOUT(",%15.8E", parameter[evoParsIndex[i]]);
#else
              STDOUT("%9f  ", evoTime);
              for (i = 0; i < pntdim; i++) STDOUT("%16.8E", point[i]*pnt_scale[i]);
              for (i = 0; i < evoParsDim; i++) STDOUT("%16.8E", parameter[evoParsIndex[i]]);
#endif
              STDOUT("\n");
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
            }

          // Compute the new direction of the evolutionary step
          COPY(pntdim, point, 1, y, 1);
          memset((void *)AtBoundary, 0, ParameterNr*sizeof(int));
          memset((void *)evoparvec, 0, ParameterNr*sizeof(double));
          memset((void *)selectdiff, 0, ParameterNr*sizeof(double));
          memset((void *)evopntvec, 0, MaxPntDim*sizeof(double));
          Jacobian(pntdim, y, pntdim, Jac, Equation, CENTRAL);
          for (i = 0, TotalAtBoundary = 0; i < evoParsDim; i++)
            {
              // I can not use SelectionGradient here, as I need the derivatives of all components of
              // F(y) w.r.t. the evolutionary parameter
              pardif = max(fabs(Jacobian_Step*(parameter[evoParsIndex[i]])), Jacobian_Min_Step);
              if (FastNumerics == 1)
                {
                  if (CentralDerivative(pntdim, Equation, y, rhs, parameter + evoParsIndex[i], pardif, rhs, dFdp + i*pntdim, 1) == FAILURE)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                      FreeHeapMemory();
                      return;
                    }
                }
              else
                {
                  for (j = 0; j < pntdim; j++)
                    {
                      if (CentralDerivative(pntdim, Equation, y, rhs, parameter + evoParsIndex[i], pardif, rhs + j, dFdp + i*pntdim + j, 0) == FAILURE)
                        {
                          ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                          FreeHeapMemory();
                          return;
                        }
                    }
                }

              selectdiff[i] =*(dFdp + i*pntdim + R0ResIndex[evoPopIndex[i]]);

              // Solve the linear system dFdp + dFdE*dEdp = 0
              COPY(pntdim*pntdim, Jac, 1, JacCopy, 1);
              SCAL(pntdim, -1.0, dFdp + i*pntdim, 1);
              if (SolveLinearSystem(pntdim, JacCopy, dFdp + i*pntdim, DYTOL) != SUCCES) continue;

              // Stop the evolution in a parameter when it has reached the limit of its domain and
              // the selection differential is pointing outward. The selection differential is set to 0
              // here to avoid effects through the variance/covariance matrix (see next loop)
              if (((parameter[evoParsIndex[i]] >= (1.0 - DBL_EPSILON)*evoParMax[i]) && (selectdiff[i] > 0)) ||
                  ((parameter[evoParsIndex[i]] <= (1.0 + DBL_EPSILON)*evoParMin[i]) && (selectdiff[i] < 0)))
                {
                  AtBoundary[i] = 1;
                  TotalAtBoundary++;
                  selectdiff[i] = 0.0;
                }
            }
          for (i = 0; i < evoParsDim; i++)
            {
              if (AtBoundary[i])
                evoparvec[evoParsIndex[i]] = 0.0;
              else
                {
                  // The change in an evolutionary parameter is a linear combination of the selection differentials
                  // in all evolutionary parameters, where the weighing factors are given the corresponding row in the
                  // variance/covariance matrix
                  for (j = 0; j < evoParsDim; j++) evoparvec[evoParsIndex[i]] += Covariances[i*evoParsDim + j]*selectdiff[j];
                  evoparvec[evoParsIndex[i]] *= Beq[evoPopIndex[i]];

                  // The final estimate of the next solution point is a linear combination of all tangent vectors
                  // for the different evolutionary parameters
                  AXPY(pntdim, evoparvec[evoParsIndex[i]], dFdp + i*pntdim, 1, evopntvec, 1);
                }
            }

          // Additional call to reset global variables
          Equation(y, rhs);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif

          // Signal curve stop if one of the components has become negative or all of the evolutionary parameters are out of bounds
          CurveEnd = CurveEnd || (TotalAtBoundary == evoParsDim);
          for (i                     = 0; i < pntdim; i++)
            if (pntnr > 20) CurveEnd = CurveEnd || (point[i]*pnt_scale[i] < pntmin[i] - epsMach);

          // Generate output: invoked after setting CurveEnd to allow for output of last solution point on branch
          outmax = DefineOutput(point, Output);
          if (outmax)
            {
              fprintf(outfile, "%14f", evoTime);
              PrettyPrintArray(outfile, outmax - 1, Output + 1);
            }

          for (i = 0, Converged = 1; i < evoParsDim; i++)
            Converged = Converged && (fabs(oldpars[evoParsIndex[i]] - parameter[evoParsIndex[i]]) <
                                      0.5*ESSTOL*(oldcurvestep/Stepreduce)*fabs(oldpars[evoParsIndex[i]] + parameter[evoParsIndex[i]]));
          // End the continuation at end of curve after generation of last output point
          if (CurveEnd || Converged) break;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif

          // Store info on current solution point
          COPY(ParameterNr, parameter, 1, oldpars, 1);
          COPY(pntdim, point, 1, oldpoint, 1);
          Stepchange = 1;

          // Report on selection gradient and stepsize
          for (i = 0; i < evoParsDim; i++)
            {
              ReportMsg("Selection gradient in parameter %2d: %.8G\n", evoParsIndex[i], evoparvec[evoParsIndex[i]]);
              ReportMsg("Realized step in parameter      %2d: %.8G\n", evoParsIndex[i], curvestep*evoparvec[evoParsIndex[i]]/Stepreduce);
            }

          AXPY(ParameterNr, (curvestep/Stepreduce), evoparvec, 1, parameter, 1);
          AXPY(pntdim, (curvestep/Stepreduce), evopntvec, 1, point, 1);

          // Prediction
          ReportMsg("\n\nPrediction :\t");
          for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
          for (i = 0; i < evoParsDim; i++) ReportMsg("%16.8E  ", parameter[evoParsIndex[i]]);
          ReportMsg("\n");

          pntnr++;
        }
      fflush(NULL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

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
  int         i;

#if (defined(_MSC_VER) && (_MSC_VER < 1500)) || (defined(R_PACKAGE) && defined(_WIN32) && !(defined(_UCRT) || ((__MSVCRT_VERSION__ >= 0x1400) || (__MSVCRT_VERSION__ >= 0xE00 && __MSVCRT_VERSION__ < 0x1000))))
  (void)_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  // Initialize some variables
  errfile    = NULL;
  outfile    = NULL;
  Stepchange = 0;
  Stepreduce = 1;
  strcpy(runname, "");

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  CtrlCPressed = 0;
#endif

  // Get the machine precisions
  epsMach = dlamch("Epsilon" FCONE);

  MaxPntDim         = EnvironDim + PopulationNr;
  CurPopulationNr   = PopulationNr;
  CohortDim         = IStateDim + 1;
  PopDensCohortDim  = CohortDim;
  evoParsDim        = 0;
  pntdim            = MaxPntDim;

  PopBPIndex        = 0;
  EnvBPIndex        = 0;
  Bifparone         = -1;
  Bifpartwo         = -1;

  TestRun       = 0;
  DoStateOutput = 0;
  SortIndex     = 0;
  ReportLevel   = 1;
  
  CurveType         = EVODYN;
  
  eVarPntr          = Evar;
  birthRatePntr     = Beq;
  parPntr           = parameter;
  evoParsIndexPntr  = evoParsIndex;
  timePntr          = &evoTime;
  
  memset(PopTrivEqui, 0, PopulationNr*sizeof(int));
  memset(EnvTrivEqui, 0, EnvironDim*sizeof(int));
#if (ALLOWNEGATIVE)
  for (i = 0; i < MaxPntDim; i++) pntmin[i] = -SAFETY*DBL_MAX;
#else
  for (i = 0; i < MaxPntDim; i++) pntmin[i] = 0.0;
#endif
  for (i = 0; i < ParameterNr; i++) evoParMin[i]    = -SAFETY*DBL_MAX;
  for (i = 0; i < ParameterNr; i++) evoParMax[i]    = SAFETY*DBL_MAX;
  for (i = 0; i < ParameterNr; i++) evoPopIndex[i]  = -1;

  // The following variables can be modified by the user with optional #define statements 
  LogMinSurvival      = log(MIN_SURVIVAL);
  CohortNr            = COHORT_NR + 1;
  
  Jacobian_Min_Step   = JACOBIAN_MIN_STEP;
  Jacobian_Step       = JACOBIAN_STEP;
  Jacobian_Updates    = JACOBIAN_UPDATES;
  FastNumerics        = FASTNUMERICS;

  Odesolve_Init_Step  = ODESOLVE_INIT_STEP;
  Odesolve_Fixed_Step = ODESOLVE_FIXED_STEP;
  Odesolve_Min_Step   = ODESOLVE_MIN_STEP;
  Odesolve_Max_Step   = ODESOLVE_MAX_STEP;
  Odesolve_Abs_Err    = ODESOLVE_ABS_ERR;
  Odesolve_Rel_Err    = ODESOLVE_REL_ERR;
  Odesolve_Func_Tol   = ODESOLVE_FUNC_TOL;

  Time                = 0;

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
  int  i;
  char tmpstr[MAX_STR_LEN];
  char desc[MAX_STR_LEN];
  char varstr[MAX_STR_LEN];

  strcpy(desc, "Aim:\tSimulating evolutionary dynamics of a structured population using the canonical equation");
  strcpy(varstr, "");
  for (i = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i]) continue;
      snprintf(tmpstr, sizeof(tmpstr), " E[%d]", i);
      strcat(varstr, tmpstr);
    }
  for (i = 0; i < PopulationNr; i++)
    {
      if (PopTrivEqui[i]) continue;
      snprintf(tmpstr, sizeof(tmpstr), " b[%d]", i);
      strcat(varstr, tmpstr);
    }
  strcat(varstr, " Par.1 ... Par.n");

  fprintf(stderr, "Usage:\t%s [<options>]%s %s %s", progname, varstr, "<max. evol. time step> <max. evol. time>",
          "<pop. nr.1> <index par.1> <min. par.1> <max. par.1> .... <pop. nr.n> <index par.n> <min. par.n> <max. par.n>");
  fprintf(stderr, "\n\n%s\n\n", desc);
  fprintf(stderr, "Evolutionary time step and maximum evolutionary time have to be positive, all indices of evolutionary ");
  fprintf(stderr, "parameters have to be in the appropriate range\n");
  fprintf(stderr, "Add further quadruple values for each parameter that should evolve\n\n");
  fprintf(stderr, "Possible options are:\n\n");
  fprintf(stderr, "\t-envZE   <index>  : Index of environment variable in trivial equilibrium (can be used multiple times)\n");
  fprintf(stderr, "\t-popZE   <index>  : Index of structured population in trivial equilibrium (can be used multiple times)\n");
  fprintf(stderr, "\t-evoPars <number> : Number of life history parameters of structured population that evolve\n");
  fprintf(stderr, "\t-isort   <index>  : Index of i-state variable to use as ruling variable for sorting the structured populations\n");
  fprintf(stderr, "\t-report  <value>  : Interval of reporting computed output to console. Minimum value of 1 implies output of every point.\n");
  fprintf(stderr, "\t-test             : Perform only a single integration over the life history, reporting dynamics of survival, R0,\n");
  fprintf(stderr, "\t                    i-state and interaction variables\n");
  fprintf(stderr, "\nThe value for -evoPars has to be set on the command-line\n");
  fprintf(stderr, "The values for -envZE and -popZE are undefined unless set on the command-line\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "\n%s, Copyright (C) 2015, Andre M. de Roos, University of Amsterdam\n\n", progname);
  fprintf(stderr, "This program comes with ABSOLUTELY NO WARRANTY; without even the implied warranty of\n");
  fprintf(stderr, "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public\n");
  fprintf(stderr, "License (<http://www.gnu.org/philosophy/why-not-lgpl.html>) for more details\n\n");

  exit(1);                                                                           // Only executed when in command-line mode

  return;
}


/*==================================================================================================================================*/

int main(int argc, char **argv)

{
  register int  i, j, rind;
  int           my_argc, tmpint, cvfdone = 0, epd = 0;
  char          **      argpnt1 = NULL, **argpnt2 = NULL, *my_argv[argc];
  char          *ch;
  double        initpars[PARAMETER_NR];

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = ENVIRON_DIM;
  InteractDim  = INTERACT_DIM;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

  //=========================== Process the command-line arguments ===================================================================
  if (argc < 2) Usage(argv[0]);

  argpnt1 = argv;
  argpnt2 = my_argv;
  my_argc = 0;

  // Store the program name
  *argpnt2 =*argpnt1;
  my_argc++;
  argpnt1++;
  argpnt2++;

  while (*argpnt1)
    {
      if (!strcmp(*argpnt1, "-?") || !strcmp(*argpnt1, "--help"))
        {
          Usage(argv[0]);
        }
      else if (!strcmp(*argpnt1, "-evoPars"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo number specified for the number of evolving life history parameters!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if (tmpint <= 0)
            {
              fprintf(stderr, "\nNumber of evolving life history parameters should be larger than 0!\n");
              Usage(argv[0]);
            }
          evoParsDim = tmpint;
        }
      else if (!strcmp(*argpnt1, "-envZE"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo index of environment variable in boundary (trivial) equilibrium specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            {
              fprintf(stderr, "\nIndex of environment variable in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, EnvironDim);
              Usage(argv[0]);
            }
          if (EnvironmentType[tmpint] == GENERALODE)
            {
              fprintf(stderr, "\nDynamics of environment variable %d not specified as PERCAPITARATE or POPULATIONINTEGRAL!\n", tmpint);
              fprintf(stderr, "Enforcing a boundary (zero-valued) equilibrium for it not possible\n");
              Usage(argv[0]);
            }
          EnvTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(*argpnt1, "-popZE"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo index of population in trivial equilibrium specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            {
              fprintf(stderr, "\nIndex of structured population in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, PopulationNr);
              Usage(argv[0]);
            }
          PopTrivEqui[tmpint] = 1;
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
      else if (!strcmp(*argpnt1, "-report"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo index of i-state variable specified for argument -report!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          ReportLevel = max(tmpint, 1);
        }
      else if (!strcmp(*argpnt1, "-cvf"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nName of CVF file not specified specified!\n");
              Usage(argv[0]);
            }
          (void)strcpy(runname, *argpnt1);
          if ((ch = strstr(runname, ".cvf"))) *ch = '\0';
          cvfdone = ReadCvfFile(runname);
          if (cvfdone != SUCCES)
            {
              fprintf(stderr, "\nCVF file %s not read\n", runname);
              Usage(argv[0]);
            }
        }
      else if ((!strcmp(*argpnt1, "-test")))
        {
          TestRun = 1;
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
          *argpnt2 =*argpnt1;
          my_argc++;
          argpnt2++;
        }
      argpnt1++;
    }

  if (evoParsDim <= 0)
    {
      fprintf(stderr, "\nNumber of evolving life history parameters should be larger than 0! Define with option '-evoPars'.\n");
      Usage(argv[0]);
    }

  // Determine the index of the R0 values of the structured populations in the result vector and the dimension of the point vector
  for (i = 0, rind = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i])
        pntdim--;
      else
        {
          // Switch off negativity check for population integrals
          if (EnvironmentType[i] == POPULATIONINTEGRAL) pntmin[rind + 1] = -SAFETY*DBL_MAX;
          rind++;
        }
    }
  for (i = 0; i < CurPopulationNr; i++)
    {
      if (PopTrivEqui[i])
        {
          pntdim--;
          R0ResIndex[i] = -1;
        }
      else
        R0ResIndex[i] = rind++;
    }

  // Program name, pntdim variables, evoParsDim parameter values, 2 evolutionary time settings, n*4 evolutionary parameter settings
  if ((my_argc - pntdim - evoParsDim - 3) % 4)
    {
      fprintf(stderr, "\nEvolutionary parameter settings are not quadruplets!\n\n");
      Usage(argv[0]);
    }

  // Map all command-line variables into the argument vector
  memset((void *)point, 0, pntdim*sizeof(double));
  for (i = 0, j = 1; i < pntdim; i++, j++) point[i] = atof(my_argv[j]);
  for (i = 0; i < evoParsDim; i++, j++) initpars[i] = atof(my_argv[j]);

  if ((Maxcurvestep = atof(my_argv[j++])) <= 0)
    {
      fprintf(stderr, "\nMaximum evolutionary time step parameter is negative!\n\n");
      Usage(argv[0]);
    }
  curvestep = 0.1*Maxcurvestep;

  if ((evoMaxTime = atof(my_argv[j++])) <= 0)
    {
      fprintf(stderr, "\nMaximum evolutionary simulation time parameter is negative!\n\n");
      Usage(argv[0]);
    }

  epd = 0;
  while (j < my_argc)
    {
      tmpint = atoi(my_argv[j++]);
      if ((tmpint < 0) || (tmpint >= PopulationNr))
        {
          fprintf(stderr,
                  "\nIndex of structured population for evolutionary dynamics simulation (%d) not in the appropriate range (0 <= i < %d)!\n",
                  tmpint, PopulationNr);
          Usage(argv[0]);
        }
      evoPopIndex[epd] = tmpint;

      tmpint = atoi(my_argv[j++]);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        {
          fprintf(stderr, "\nIndex of evolutionary parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint, ParameterNr);
          Usage(argv[0]);
        }
      for (i = 0; i < epd; i++)
        if (tmpint == evoParsIndex[i])
          {
            fprintf(stderr, "\nMultiple specifications for the same evolutionary parameter (%d) not allowed!\n\n", tmpint);
            Usage(argv[0]);
          }
      evoParsIndex[epd] = tmpint;

      parameter[tmpint] = initpars[epd];

      evoParMin[epd] = atof(my_argv[j++]);
      if (parameter[tmpint] < evoParMin[epd])
        {
          fprintf(stderr, "\nValue of evolutionary parameter too small (%G < %G)!\n\n", parameter[tmpint], evoParMin[epd]);
          Usage(argv[0]);
        }
      evoParMax[epd] = atof(my_argv[j++]);
      if (parameter[tmpint] > evoParMax[epd])
        {
          fprintf(stderr, "\nValue of evolutionary parameter too large (%G > %G)!\n\n", parameter[tmpint], evoParMax[epd]);
          Usage(argv[0]);
        }
      epd++;
    }
  if (epd != evoParsDim)
    {
      fprintf(
          stderr,
          "\nNumber of evolving parameters specified with '-evoPars %d' does not correspond to the number of specified parameter triplets (%d)!\n\n",
          evoParsDim, epd);
      Usage(argv[0]);
    }

  memset((void *)Covariances, 0, ParameterNr*ParameterNr*sizeof(double));
  for (i = 0; i < evoParsDim; i++) Covariances[i*evoParsDim + i] = 1.0;

  ComputeCurve(argc, argv);
  
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
  if (errfile) fclose(errfile);
  if (outfile) fclose(outfile);
#if defined(MATLAB_MEX_FILE) && (FULLSTATEOUTPUT > 0)
  if (pmat) matClose(pmat);
#endif

  ResetCurve();
  FreeHeapMemory();

  return;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t          nrows, ncols;
  double          *inputVals, tmpdouble;
  int             tmpint, rind, irhs;
  const mxArray   *cell_element_ptr;
  char            optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN], varstr[MAX_STR_LEN];

  mwIndex         i, j;
  size_t          total_num_of_cells, buflen;
  int             status;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = ENVIRON_DIM;
  InteractDim  = INTERACT_DIM;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

  // check for proper number of arguments
  if (nrhs != 6)
    mexErrMsgIdAndTxt(
        "MATLAB:PSPMevodyn:nrhs",
        "\nIncorrect number of command-line arguments.\n\nUse: %s(%s)\n\n%22s: %s\n%22s: %s\n%22s: %s\n%22s: %s\n%22s: %s\n%22s: %s",
        mexFunctionName(), "<point>, <curve settings>, <evolution settings>, <covariance settings>, <parameters>, <options>", "<point>",
        "Initial point of the computation", "<curve settings>", "Maximum step size and maximum integration time of evolutionary dynamics (2 values)",
        "<evolution settings>", "Index, minimum and maximum of parameters to evolve (N*3 values with N the number of evolving parameters)",
        "<covariance settings>", "Variance-covariance matrix of mutations in vector format. (NxN values with N the number of evolving parameters)",
        "<parameters>", "Array of parameter values to use (empty array or of same length as parameter array)", "<options>",
        "Possible options: envZE, popZE, test or isort");

  // check for proper number of output variables
  if (nlhs != 1) mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:nlhs", "\nA single output argument is required!");

#if (MFUNCTIONS == 1)
  Minterface_Init();
#endif

  //============================== Process the options argument ======================================================================
  // Extract the contents of MATLAB cell into the C array
  irhs = 5;
  if (!mxIsCell(prhs[irhs])) mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options", "\nOptions should be specified as a cell array!");

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(optstring, "{");
  for (i = 0; i < total_num_of_cells; i++)
    {
      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optname, buflen);
      if (!((!strcmp(optname, "envZE")) || (!strcmp(optname, "popZE")) || (!strcmp(optname, "test")) ||
            (!strcmp(optname, "isort")) || (!strcmp(optname, "report"))))
        mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options", "\nIllegal option %s!", optname);

      if (!strcmp(optname, "test"))
        {
          TestRun = 1;

          // optstring still equal to "{"
          if (strlen(optstring) == 1)
            strcat(optstring, "'");
          else
            strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "'");
          continue;
        }

      if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options", "\nNo value specified for option %s!", optname);

      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optval, buflen);
      if (status) mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options", "\nError in retrieving value for option %s!", optname);

      // optstring still equal to "{"
      if (strlen(optstring) == 1)
        strcat(optstring, "'");
      else
        strcat(optstring, ", '");
      strcat(optstring, optname);
      strcat(optstring, "', '");
      strcat(optstring, optval);
      strcat(optstring, "'");

      if (!strcmp(optname, "envZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options",
                              "\nIndex of environment variable in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!",
                              tmpint, EnvironDim);
          if (EnvironmentType[tmpint] == GENERALODE)
            mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options",
                              "\nDynamics of environment variable %d not specified as PERCAPITARATE or POPULATIONINTEGRAL.\n%s!", tmpint,
                              "Enforcing a boundary (zero-valued) equilibrium for it not possible");
          EnvTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(optname, "popZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options",
                              "\nIndex of structured population in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!",
                              tmpint, PopulationNr);
          PopTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(optname, "isort"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= IStateDim))
            mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:options",
                              "\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!",
                              tmpint, IStateDim);
          SortIndex = tmpint;
        }
      else if (!strcmp(optname, "report"))
        {
          tmpint = atoi(optval);
          ReportLevel = max(tmpint, 1);
        }
    }
  strcat(optstring, "}");

  // Determine the index of the R0 values of the structured populations in the result vector and the dimension of the point vector
  for (i = 0, rind = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i])
        pntdim--;
      else
        {
          // Switch off negativity check for population integrals
          if (EnvironmentType[i] == POPULATIONINTEGRAL) pntmin[rind + 1] = -SAFETY*DBL_MAX;
          rind++;
        }
    }
  for (i = 0; i < CurPopulationNr; i++)
    {
      if (PopTrivEqui[i])
        {
          pntdim--;
          R0ResIndex[i] = -1;
        }
      else
        R0ResIndex[i] = rind++;
    }

  //============================== Process the parameters argument ===================================================================

  irhs  = 4;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols == ParameterNr) && (nrows == 1))
    memcpy(parameter, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:PSPMevodyn:parameters", "\nParameter argument ignored as it is not a row vector of length %d!", ParameterNr);

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(parstring, "[");
  for (i = 0; i < total_num_of_cells; i++)
    {
      if (i) strcat(parstring, " ");
      memcpy(&tmpdouble, mxGetPr(prhs[irhs]) + i, mxGetElementSize(prhs[irhs]));
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", tmpdouble);
      strcat(parstring, tmpstr);
    }
  strcat(parstring, "]");

  //============================== Process the evolutionary parameters argument ======================================================

  irhs  = 2;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);

  if (!ncols || (ncols % 4) || (nrows != 1))
    mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evopars",
                      "\nEvolutionary parameter settings needs to be multiple of 4: population number, index, minimum and maximum of the parameter!");

  inputVals = mxGetPr(prhs[irhs]);

  evoParsDim = 0;
  for (i = 0; i < ncols; i += 4)
    {
      tmpint = (int)floor(inputVals[i] + MICRO);
      if ((tmpint < 0) || (tmpint >= PopulationNr))
        mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evopars", "\nIndex of structured population for evolutionary simulations (%d) not in the appropriate range (0 <= i < %d)!\n\n",
                          tmpint, PopulationNr);
      evoPopIndex[evoParsDim] = tmpint;
      
      tmpint = (int)floor(inputVals[i + 1] + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evopars", "\nIndex of evolutionary parameter #%d (%d) not in the appropriate range (0 <= i < %d)!", i,
                          tmpint, ParameterNr);

      for (j = 0; j < evoParsDim; j++)
        if (tmpint == evoParsIndex[j])
          mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evopars", "\nMultiple specifications for the same evolutionary parameter (%d) not allowed!", tmpint);

      evoParsIndex[evoParsDim] = tmpint;
      evoParMin[evoParsDim]    = inputVals[i + 2];
      evoParMax[evoParsDim]    = inputVals[i + 3];
      if (evoParMin[evoParsDim] >= evoParMax[evoParsDim])
        mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evopars", "\nMinimum parameter bound (%G not smaller than maximum (%G)!\n\n", evoParMin[evoParsDim],
                          evoParMax[evoParsDim]);
      evoParsDim++;
    }

  strcpy(evoparstring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(evoparstring, " ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", inputVals[i]);
      strcat(evoparstring, tmpstr);
    }
  strcat(evoparstring, "]");

  //============================== Process the covariance matrix argument ============================================================

  irhs      = 3;
  nrows     = mxGetM(prhs[irhs]);
  ncols     = mxGetN(prhs[irhs]);
  inputVals = mxGetPr(prhs[irhs]);

  memset((void *)Covariances, 0, ParameterNr*ParameterNr*sizeof(double));
  if ((ncols != (evoParsDim*evoParsDim)) || (nrows != 1))
    {
      if (ncols)
        mexWarnMsgIdAndTxt("MATLAB:PSPMevodyn:covariances", "\nCovariance matrix argument ignored as it is not a vector of length %d*%d!", evoParsDim,
                           evoParsDim);
      for (i = 0; i < evoParsDim; i++) Covariances[i*evoParsDim + i] = 1.0;
    }
  else
    {
      for (i = 0; i < evoParsDim; i++)
        for (j = 0; j < evoParsDim; j++) Covariances[i*evoParsDim + j] = inputVals[i*evoParsDim + j];
    }

  strcpy(covmatstring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(covmatstring, " ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", inputVals[i]);
      strcat(covmatstring, tmpstr);
    }
  strcat(covmatstring, "]");

  //============================== Process the evolutionary time stepping argument ===================================================

  irhs  = 1;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);

  if ((ncols != 2) || (nrows != 1))
    mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:bounds",
                      "\nEvolutionary time argument must be a row vector of length 2: maximum evolutionary step size and simulation time!");

  inputVals = mxGetPr(prhs[irhs]);

  if ((Maxcurvestep = inputVals[0]) <= 0) mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evotime", "\nMaximum evolutionary time step parameter is negative!");
  curvestep         = 0.1*Maxcurvestep;
  if ((evoMaxTime = inputVals[1]) <= 0)
    mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:evotime", "\nMaximum evolutionary simulation time parameter is negative!");

  strcpy(evotimestring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(evotimestring, " ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", inputVals[i]);
      strcat(evotimestring, tmpstr);
    }
  strcat(evotimestring, "]");

  //============================== Process the initial point argument ================================================================

  irhs = 0;
  memset((void *)point, 0, pntdim*sizeof(double));
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((nrows != 1) || (ncols != (pntdim + evoParsDim)))
    {
      strcpy(varstr, "[");
      for (i = 0; i < EnvironDim; i++)
        {
          if (EnvTrivEqui[i]) continue;
          snprintf(tmpstr, sizeof(tmpstr), " E[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          if (PopTrivEqui[i]) continue;
          snprintf(tmpstr, sizeof(tmpstr), " b[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < evoParsDim; i++)
        {
          snprintf(tmpstr, sizeof(tmpstr), " Par[%d]", (int)evoParsIndex[i]);
          strcat(varstr, tmpstr);
        }
      strcat(varstr, " ]");
      mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:point", "\nInitial point for evolutionary dynamics must be a row vector of length %d:  %s!",
                        (pntdim + evoParsDim), varstr);
    }
  memcpy(point, mxGetPr(prhs[irhs]), pntdim*mxGetElementSize(prhs[irhs]));

  inputVals = mxGetPr(prhs[irhs]);
  for (i = 0; i < evoParsDim; i++)
    {
      parameter[evoParsIndex[i]] = inputVals[pntdim + i];
      if ((parameter[evoParsIndex[i]] < evoParMin[i]) || (parameter[evoParsIndex[i]] > evoParMax[i]))
        mexErrMsgIdAndTxt("MATLAB:PSPMevodyn:point", "\nValue of evolutionary parameter #%d (%G) not in the appropriate range (%G - %G)!",
                          evoParsIndex[i], parameter[evoParsIndex[i]], evoParMin[i], evoParMax[i]);
    }

  strcpy(pntstring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(pntstring, " ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", inputVals[i]);
      strcat(pntstring, tmpstr);
    }
  strcat(pntstring, "]");

  //============================= Get the program name ===============================================================================
  // Get the name of the mex file
  strcpy(progname, mexFunctionName());
  progname[strlen(progname) - 6] = '\0';                                            // Cut off the 'evodyn' appendix

  mexAtExit(CloseStreams);

  // call the computational routine
  ComputeCurve(0, NULL);
  FreeHeapMemory();

  if (TestRun)
    plhs[0] = mxCreateString("");
  else
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

SEXP PSPMevodyn(SEXP moduleName, SEXP initVals, SEXP evotimeVals, SEXP evoparsVals, SEXP covmatVals, SEXP parVals, SEXP optVals)

{
  int   i, j, ncols, tmpint, rind;
  char  optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN], varstr[MAX_STR_LEN];
  SEXP  resfil;

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
      if (!((!strcmp(optname, "envZE")) || (!strcmp(optname, "popZE")) || (!strcmp(optname, "test")) ||
            (!strcmp(optname, "isort")) || (!strcmp(optname, "report"))))
        error("\nIllegal option %s!\n\n", optname);

      if (!strcmp(optname, "test"))
        {
          TestRun = 1;
          if (!strlen(optstring))
            strcat(optstring, "c(\"");
          else
            strcat(optstring, ", \"");
          strcat(optstring, optname);
          strcat(optstring, "\"");
          continue;
        }

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

      if (!strcmp(optname, "envZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            error("\nIndex of environment variable in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  EnvironDim);
          if (EnvironmentType[tmpint] == GENERALODE)
            error("\nDynamics of environment variable %d not specified as PERCAPITARATE or POPULATIONINTEGRAL!\n%s.\n\n", tmpint,
                  "Enforcing a boundary (zero-valued) equilibrium for it not possible");
          EnvTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(optname, "popZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            error("\nIndex of structured population in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  PopulationNr);
          PopTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(optname, "isort"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= IStateDim))
            error("\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  IStateDim);
          SortIndex = tmpint;
        }
      else if (!strcmp(optname, "report"))
        {
          tmpint = atoi(optval);
          ReportLevel = max(tmpint, 1);
        }
    }

  if (strlen(optstring))
    strcat(optstring, ")");
  else
    strcpy(optstring, "NULL");

  // Determine the index of the R0 values of the structured populations in the result vector and the dimension of the point vector
  for (i = 0, rind = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i])
        pntdim--;
      else
        {
          // Switch off negativity check for population integrals
          if (EnvironmentType[i] == POPULATIONINTEGRAL) pntmin[rind + 1] = -SAFETY*DBL_MAX;
          rind++;
        }
    }
  for (i = 0; i < CurPopulationNr; i++)
    {
      if (PopTrivEqui[i])
        {
          pntdim--;
          R0ResIndex[i] = -1;
        }
      else
        R0ResIndex[i] = rind++;
    }

  //============================== Process the parameters argument ===================================================================

  ncols = length(parVals);
  if (isReal(parVals) && (ncols == ParameterNr))
    memcpy(parameter, REAL(parVals), ncols*sizeof(double));
  else if (ncols)
    warning("\nParameter argument ignored as it is not a row vector of same length as parameter array\n\n");

  if (ncols)
    {
      strcpy(parstring, "c(");
      for (i = 0; i < ncols; i++)
        {
          if (i) strcat(parstring, ", ");
          snprintf(tmpstr, sizeof(tmpstr), "%.6G", REAL(parVals)[i]);
          strcat(parstring, tmpstr);
        }
      strcat(parstring, ")");
    }
  else
    strcpy(parstring, "NULL");

  //============================== Process the evolutionary parameters argument ======================================================

  ncols = length(evoparsVals);
  if (!isReal(evoparsVals) || !ncols || (ncols % 4)) 
    error("\nEvolutionary parameter settings needs to be multiple of 4: population number, index, minimum and maximum of the parameter!\n\n");

  evoParsDim = 0;
  for (i = 0; i < ncols; i += 4)
    {
      tmpint = (int)floor(REAL(evoparsVals)[i] + MICRO);
      if ((tmpint < 0) || (tmpint >= PopulationNr))
        error("\nIndex of structured population for evolutionary simulations (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
              PopulationNr);
      evoPopIndex[evoParsDim] = tmpint;
      
      tmpint = (int)floor(REAL(evoparsVals)[i + 1] + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        error("\nIndex of evolutionary parameter #%d (%d) not in the appropriate range (0 <= i < %d)!\n\n", i, tmpint, ParameterNr);

      for (j = 0; j < evoParsDim; j++)
        if (tmpint == evoParsIndex[j]) error("\nMultiple specifications for the same evolutionary parameter (%d) not allowed!\n\n", tmpint);

      evoParsIndex[evoParsDim] = tmpint;
      evoParMin[evoParsDim]    = REAL(evoparsVals)[i + 2];
      evoParMax[evoParsDim]    = REAL(evoparsVals)[i + 3];
      if (evoParMin[evoParsDim] >= evoParMax[evoParsDim])
        error("\nMinimum parameter bound (%G not smaller than maximum (%G)!\n\n", evoParMin[evoParsDim], evoParMax[evoParsDim]);
      evoParsDim++;
    }

  strcpy(evoparstring, "c(");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(evoparstring, ", ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", REAL(evoparsVals)[i]);
      strcat(evoparstring, tmpstr);
    }
  strcat(evoparstring, ")");

  //============================== Process the covariance matrix argument ============================================================

  memset((void *)Covariances, 0, ParameterNr*ParameterNr*sizeof(double));

  ncols = length(covmatVals);
  if (!isReal(covmatVals) || (ncols != (evoParsDim*evoParsDim)))
    {
      if (ncols) warning("\nCovariance argument ignored as it is not a row vector of length %d*%d\n\n", evoParsDim, evoParsDim);

      for (i = 0; i < evoParsDim; i++) Covariances[i*evoParsDim + i] = 1.0;
    }
  else
    {
      for (i = 0; i < evoParsDim; i++)
        for (j = 0; j < evoParsDim; j++) Covariances[i*evoParsDim + j] = REAL(covmatVals)[i*evoParsDim + j];
    }

  if (ncols)
    strcpy(covmatstring, "c(");
  else
    strcpy(covmatstring, "NULL");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(covmatstring, ", ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", REAL(covmatVals)[i]);
      strcat(covmatstring, tmpstr);
    }
  if (ncols) strcat(covmatstring, ")");

  //============================== Process the evolutionary time stepping argument ===================================================

  ncols = length(evotimeVals);
  if (!isReal(evotimeVals) || (ncols != 2))
    error("\nThe evolutionary time argument must be a row vector of length 2: maximum evolutionary step size and simulation time!\n\n");

  if ((Maxcurvestep = REAL(evotimeVals)[0]) <= 0) error("\nMaximum evolutionary time step parameter is negative!\n\n");
  curvestep         = 0.1*Maxcurvestep;
  if ((evoMaxTime = REAL(evotimeVals)[1]) <= 0) error("\nMaximum evolutionary simulation time parameter is negative!\n\n");

  strcpy(evotimestring, "c(");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(evotimestring, ", ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", REAL(evotimeVals)[i]);
      strcat(evotimestring, tmpstr);
    }
  strcat(evotimestring, ")");

  //============================== Process the initial point argument ================================================================

  memset((void *)point, 0, pntdim*sizeof(double));

  ncols = length(initVals);
  if (!isReal(initVals) || (ncols != (pntdim + evoParsDim)))
    {
      strcpy(varstr, "c(");
      for (i = 0; i < EnvironDim; i++)
        {
          if (EnvTrivEqui[i]) continue;
          snprintf(tmpstr, sizeof(tmpstr), " E[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          if (PopTrivEqui[i]) continue;
          snprintf(tmpstr, sizeof(tmpstr), " b[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < evoParsDim; i++)
        {
          snprintf(tmpstr, sizeof(tmpstr), " Par[%d]", (int)evoParsIndex[i]);
          strcat(varstr, tmpstr);
        }
      strcat(varstr, " )");
      error("\nInitial point for evolutionary dynamics must be a row vector of length %d:  %s!\n\n", pntdim + evoParsDim, varstr);
    }
  memcpy(point, REAL(initVals), pntdim*sizeof(double));

  for (i = 0; i < evoParsDim; i++)
    {
      parameter[evoParsIndex[i]] = REAL(initVals)[pntdim + i];
      if ((parameter[evoParsIndex[i]] < evoParMin[i]) || (parameter[evoParsIndex[i]] > evoParMax[i]))
        error("\nValue of evolutionary parameter #%d (%G) not in the appropriate range (%G - %G)!", evoParsIndex[i], parameter[evoParsIndex[i]],
              evoParMin[i], evoParMax[i]);
    }

  strcpy(pntstring, "c(");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(pntstring, ", ");
      snprintf(tmpstr, sizeof(tmpstr), "%.6G", REAL(initVals)[i]);
      strcat(pntstring, tmpstr);
    }
  strcat(pntstring, ")");

  //============================= Get the program name ===============================================================================
  // Get the name of the module file
  if (isString(moduleName) && (length(moduleName) == 1))
    strcpy(progname, CHAR(STRING_ELT(moduleName, 0)));
  else
    error("\nModel name argument must be a single string\n\n");

  // call the computational routine
  ComputeCurve(0, NULL);
  FreeHeapMemory();

  ResetCurve();

  if (errfile) fclose(errfile);
  if (outfile) fclose(outfile);

#if (RFUNCTIONS == 1)
  Rinterface_End();
#endif

  PROTECT(resfil = mkString(runname));

  UNPROTECT(1);

  return resfil;
}


/*==================================================================================================================================*/
#endif
