/***
  NAME
    PSPMdemo

  PURPOSE
    Generic, problem-independent specification for demographic analysis
    of models with POPULATION_NR structured population, whose individuals
    are characterized by I_STATE_DIM state variables. All problem-specific
    life-history functions are specified in an include file

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

    Last modification: AMdR - May 05, 2018
***/

#define PSPMDEMO                  1

#if (!defined(RFUNCTIONS) || (RFUNCTIONS != 1))
#define RFUNCTIONS                0
#endif
#if (!defined(MFUNCTIONS) || (MFUNCTIONS != 1))
#define MFUNCTIONS                0
#endif

#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
double                            parameter[PARAMETER_NR];
const char                        *parameternames[PARAMETER_NR];
#endif

#include "globals.h"

/*
 *====================================================================================================================================
 *  Import the population and environment dimension settings
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

#ifndef _MSC_VER                                                                    // Microsoft Visual C 8.0 can not handle the #warning preprocessor statement
#if defined(ENVIRON_DIM) && (ENVIRON_DIM > 0)
#warning ENVIRON_DIM is non-zero: You can use environment variables but their value is not set in demographic analysis!
#endif
#endif

#if !defined(PARAMETER_NR) || (PARAMETER_NR < 2)
#error PARAMETER_NR should be defined larger than 1
#endif


/*
 *====================================================================================================================================
 *  Definition of problem dimensions
 *====================================================================================================================================
 */

#undef PULSED
#if defined(REPRODUCTION_INTERVAL)
#define PULSED                    1
#else
#define PULSED                    0
#endif


/*
 *====================================================================================================================================
 *  Definition of global variables and parameters
 *====================================================================================================================================
 */

// Global dimension variables
static int                        BirthStateNr[POPULATION_NR];
static int                        MaxPntDim;

// These are the variables to solve for
#if defined(ENVIRON_DIM) && (ENVIRON_DIM > 0)
double                            Evar[ENVIRON_DIM];
#else
#define ENVIRON_DIM               0
static double                     *Evar = NULL;
#endif

static double                     R0[POPULATION_NR];
static double                     PGRvar[POPULATION_NR];
static double                     GenTime[POPULATION_NR];

// Global variables to hold variables shared among routines
static int                        LastMemAllocated = 0;

// Global pointers into the heap
static double                     *RightEigenvecMem = NULL;
static double                     *LeftEigenvecMem  = NULL;

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
static int                        FirstEstimate = 0;
static int                        ReportLevel   = 1;

// Global variables for other purposes

static char                       progname[MAXPATHLEN];
static double                     LogMinSurvival;
static double                     minbound1, maxbound1;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
static char                       parstring[MAX_STR_LEN];
static char                       optstring[MAX_STR_LEN];
static char                       curvestring[MAX_STR_LEN];
#endif


/*
 *====================================================================================================================================
 *  Include header files with generic routines
 *====================================================================================================================================
 */

#include "memory.h"
#include "lifehistory.h"
#include "dopri5.h"

/*==================================================================================================================================*/

int Equation(double *argument, double *result)

{
  int     indx = 0, i, j, b = 0, p, retval = SUCCES;
  int     MaxCohortDim = CohortDim + 1;
  double  norm;
  double  *NextGenMatrix = NULL, *FinalIstateMem = NULL;

  /*
   *===========================================================================
   *  Map current estimate of solution to global variables
   *===========================================================================
   */
  if (TestRun)
    memset(PGRvar, 0, PopulationNr*sizeof(double));
  else
    {
      if (Bifparone != -1)
        {
          parameter[Bifparone] = argument[indx]*pnt_scale[indx];
          indx++;
        }
      for (p = 0; p < PopulationNr; p++, indx++) PGRvar[p] = argument[indx]*pnt_scale[indx];
    }

  /*
   *===========================================================================
   *  Set the dimensions
   *===========================================================================
   */
  CurPopulationNr = PopulationNr;
  SetBirthStates(BirthStateNr, Evar);

  MaxStatesAtBirth = 0;
  for (p = 0; p < CurPopulationNr; p++)
    {
      BirthStateNr[p]  = max(BirthStateNr[p], 1);
      MaxStatesAtBirth = max(MaxStatesAtBirth, BirthStateNr[p]);
      MaxCohortDim     = CohortDim + MaxStatesAtBirth;
      PopDensCohortDim = IStateDim + 1 + MaxStatesAtBirth;
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
   *  Testing output
   *===========================================================================
   */

  if (TestRun)
    {
      if (Bifparone != -1) STDOUT("\n\nParameter #1:                %15.6G\n", parameter[Bifparone]);

      STDOUT("\n\n%40s", "");
      for (i = 0; i < IStateDim; i++) STDOUT("%12s%2d]", "Istate[", i);
      STDOUT("       Survival%15s%15s", "R0", "Gen. time");
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

  //==================================================================================================================================
  // Compute the final values of the fixed point equation F(y)=0,

  for (p = 0; p < CurPopulationNr; p++)
    {
      if (BirthStateNr[p] == 1)
        {
          R0[p] = FinalIstate(0, p, IStateDim + 1);
          RightEigenvec(p, 0) = 1.0;
          LeftEigenvec(p, 0)  = 1.0;

          if (R0[p] > 0)
            GenTime[p] = FinalIstate(0, p, IStateDim + 2)/R0[p];
          else
            GenTime[p] = INFINITY;

          result[p] = R0[p] - 1.0;
#if (FULLSTATEOUTPUT > 0)
          CohortLimit(0, p) = (FinalIstate(0, p, SortIndex) - BirthState(p, 0, SortIndex))/(COHORT_NR - 1);
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

          if (Eigenval(BirthStateNr[p], NextGenMatrix, 0, R0 + p, DOMINANT, &(RightEigenvec(p, 0)), &(LeftEigenvec(p, 0)), RHSTOL) != SUCCES)
            {
              ErrorMsg(__FILE__, __LINE__, "Computation of dominant eigenvalue failed for population %d!", p);
              if (FinalIstateMem) free(FinalIstateMem);
              if (NextGenMatrix) free(NextGenMatrix);
              return FAILURE;
            }
          result[p] = R0[p] - 1.0;

          if (TestRun || DoStateOutput || FirstEstimate)
            {
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

              // Scale the left eigenvector such that u^T v = 1
              // All elements of RightEigenvec are positive, while the elements of
              // LeftEigenvec are all negative or all positive. Scaling with the
              // inproduct should make all elements of  LeftEigenvec positive
              norm = DOT(BirthStateNr[p], &(RightEigenvec(p, 0)), 1, &(LeftEigenvec(p, 0)), 1);
              SCAL(BirthStateNr[p], 1.0/norm, &(LeftEigenvec(p, 0)), 1);

              // Check the validity of the left eigenvector
              for (b = 0; b < BirthStateNr[p]; b++)
                {
                  if (LeftEigenvec(p, b) < 0)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Negative and positive elements in left eigenvector of population %d!", p);
                      if (FinalIstateMem) free(FinalIstateMem);
                      if (NextGenMatrix) free(NextGenMatrix);
                      return FAILURE;
                    }
                }

              // NEEDS TO BE CHECKED!!: This might not be correct, but I do not know what the definition of the generation time is in the multidimensional context
              if (R0[p] > 0)
                {
                  GenTime[p] = 0;
                  for (b = 0; b < BirthStateNr[p]; b++) GenTime[p] += RightEigenvec(p, b)*FinalIstate(b, p, IStateDim + BirthStateNr[p] + 1);
                  GenTime[p] /= R0[p];
                }
              else GenTime[p] = INFINITY;
            }
#if (FULLSTATEOUTPUT == 1)
          CohortMin[p] = DBL_MAX;
          CohortMax[p] = -DBL_MAX;
          for (b = 0; b < BirthStateNr[p]; b++)
            {
              if (RightEigenvec(p, b))
                {
                  CohortMin[p] = min(CohortMin[p], BirthState(p, b, SortIndex));
                  CohortMax[p] = max(CohortMax[p], FinalIstate(b, p, SortIndex));
                }
            }
          for (b = 0; b < BirthStateNr[p]; b++) CohortLimit(b, p) = (CohortMax[p] - CohortMin[p])/(COHORT_NR - 1);
#elif (FULLSTATEOUTPUT == 2)
          for (b = 0; b < BirthStateNr[p]; b++) CohortLimit(b, p) = (FinalIstate(b, p, SortIndex) - BirthState(p, b, SortIndex))/(COHORT_NR - 1);
#endif
        }
    }

  if (TestRun)
    {
      for (p = 0; p < CurPopulationNr; p++)
        {
          STDOUT("\n\n\nPopulation #%2d:          ", p);
          for (i = 0; i < IStateDim; i++) STDOUT("%12s%2d]", "Istate[", i);
          STDOUT("       Survival%15s%15s%15s", "R0", "Gen. time", "PGR (approx.)");

          for (b = 0; b < BirthStateNr[p]; b++)
            {
              STDOUT("\nBstate %2d - (Final):     ", b);
              for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", FinalIstate(b, p, i));
              STDOUT("%15.6G", exp(FinalIstate(b, p, IStateDim)));
              STDOUT("%15.6G", R0[p]);
              STDOUT("%15.6G", GenTime[p]);
              if (R0[p] > 0)
                STDOUT("%15.6G", log(R0[p])/GenTime[p]);
              else
                STDOUT("%15s", "-");
            }

          if (BirthStateNr[p] > 1)
            {
              STDOUT("\n\n\ndet(M-I) : %12.6E\n", result[p]);
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

              STDOUT("\n\nStable birth distribution       Reproductive value\n");
              for (b = 0; b < BirthStateNr[p]; b++)
                {
                  STDOUT("%25.6G", RightEigenvec(p, b));
                  STDOUT("%25.6G", LeftEigenvec(p, b));
                  STDOUT("\n");
                }
            }
        }

      STDOUT("\n\n");
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  if (FinalIstateMem) free(FinalIstateMem);
  if (NextGenMatrix) free(NextGenMatrix);
  return SUCCES;
}


/*==================================================================================================================================*/

int DefineOutput(double *x, double *output)

{
  int     outnr = 0, i, b, j, k, p;
  double  reproval;
  double  rhsL[MaxPntDim];
  double  rhsH[MaxPntDim];
  double  drdp[ParameterNr][PopulationNr];
  double  old, dif, dR0dr, *parpntr;

  // Compute the sensitivities
  for (p = 0, k = (Bifparone != -1); p < PopulationNr; p++, k++)
    {
      old = x[k];
      dif = max(fabs(Jacobian_Step*x[k]), Jacobian_Min_Step);

      x[k] = old + dif;
      Equation(x, rhsH);

      x[k] = old - dif;
      Equation(x, rhsL);
      x[k]  = old;
      dR0dr = (rhsH[p] - rhsL[p])/(2*dif);
      dR0dr /= pnt_scale[k];

      for (j = 0; j < ParameterNr; j++)
        {
          if (j == Bifparone)
            parpntr = x;
          else
            parpntr = parameter + j;

          old =*parpntr;
          dif = max(fabs(Jacobian_Step*(*parpntr)), Jacobian_Min_Step);

          *parpntr = old + dif;
          dif      =*parpntr - old;
          Equation(x, rhsH);

          *parpntr = old - dif;
          if (*parpntr > 0)
            {
              Equation(x, rhsL);
              *parpntr = old;
              drdp[j][p] = (rhsH[p] - rhsL[p])/(2*dif);
            }
          else
            {
              *parpntr = old;
              Equation(x, rhsL);
              drdp[j][p] = (rhsH[p] - rhsL[p])/(dif);
            }

          drdp[j][p] /= -dR0dr;

          // Apply the scaling for the bifurcation parameter
          if (j == Bifparone) drdp[j][p] /= pnt_scale[0];
        }
    }

#if (FULLSTATEOUTPUT > 0)
  DoStateOutput = 1;
  Equation(x, rhsH);
  DoStateOutput = 0;

  // Compute the reproductive value using eq. (16) in Ecol. Letters (2008)
  for (p = 0; p < CurPopulationNr; p++)
    for (j = 0; j < CohortNr; j++)
      for (b = 0; b < BirthStateNr[p]; b++)
        {
          for (i = 0, reproval = 0; i < BirthStateNr[p]; i++)                       // Compute the sum of v_i(0)*(L_{ij}(Amax) - L_{ij}(a))
            reproval += LeftEigenvec(p, i)*(PopDens(p, b, IStateDim + 1 + i, COHORT_NR - 1) - PopDens(p, b, IStateDim + 1 + i, j));
          PopDens(p, b, IStateDim + 1, j) = reproval/PopDens(p, b, IStateDim, j);
        }

  WriteStateToFile(FULLSTATEOUTPUT);
#endif

  if (Bifparone != -1) output[outnr++] = parameter[Bifparone];

  for (p = 0; p < PopulationNr; p++)
    {
      output[outnr++] = PGRvar[p];
      output[outnr++] = GenTime[p];
      for (j = 0; j < ParameterNr; j++) output[outnr++] = drdp[j][p];
    }

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
  int           pntnr = 0, outmax, CurveEnd = 0, savedbifpar, pntdim;
  int           cycles, last = 1, retval = 0, hasjac = 0;
  double        point[MaxPntDim], oldpoint[MaxPntDim], *savedpntr;
  double        rhs[MaxPntDim];
  double        tanvec[MaxPntDim], Jacmat[MaxPntDim*MaxPntDim];
  char          csbname[MAX_STR_LEN], errname[MAX_STR_LEN], outname[MAX_STR_LEN];
  char          tmpstr[MAX_STR_LEN];
  struct stat   buffer;

  pntdim = PopulationNr + (Bifparone != -1);
  
  (void)SetScales(NULL, pntdim);

  memset((void *)point, 0, pntdim*sizeof(double));
  if (Bifparone != -1)
    {
      point[0] = parameter[Bifparone];
    }

  if (TestRun)
    {
      double rhsnorm;

#if (defined(R_PACKAGE))
      STDOUT("\n\nExecuting : ");
      STDOUT("PSPMdemo(\"%s\", %s, %s, %s)", progname, curvestring, parstring, optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      STDOUT("\n\nExecuting : ");
      STDOUT("PSPMdemo('%s', %s, %s, %s)", progname, curvestring, parstring, optstring);
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
      rhsnorm = NRM2(pntdim - 1, rhs, 1);
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

  // Calculate R0, generation time and initial estimate of population growth rate
  FirstEstimate = 1;
  Equation(point, rhs);
  FirstEstimate = 0;
  for (i = 0, j = pntdim - 1; i < PopulationNr; i++, j--) point[j] = log(R0[i])/GenTime[i];
  (void)SetScales(point, pntdim);

  if (strlen(runname))
    {
      sprintf(errname, "%s.err", runname);
      sprintf(outname, "%s.out", runname);
    }
  else
    {
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)
      strcpy(progname, argv[0]);
      progname[strlen(progname) - 4] = '\0';                                        // Cut off the 'demo' appendix
#endif
      i = 0;
      while (1)
        {
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
          sprintf(csbname, "%s-%s-%04d.mat", progname, "PGR", i);
#else
          sprintf(csbname, "%s-%s-%04d.csb", progname, "PGR", i);
#endif
          sprintf(errname, "%s-%s-%04d.err", progname, "PGR", i);
          sprintf(outname, "%s-%s-%04d.out", progname, "PGR", i);
          if (stat(csbname, &buffer) && stat(errname, &buffer) && stat(outname, &buffer)) break;
          i++;
        }
      sprintf(runname, "%s-%s-%04d", progname, "PGR", i);
    }

  errfile                      = fopen(errname, "w");
  if (Bifparone != -1) outfile = fopen(outname, "w");

  if ((Bifparone != -1) && outfile)
    {
      fprintf(outfile, "#\n# Executing : ");

#if (defined(R_PACKAGE))
      fprintf(outfile, "PSPMdemo(\"%s\", %s, %s, %s)", progname, curvestring, parstring, optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      fprintf(outfile, "PSPMdemo('%s', %s, %s, %s)", progname, curvestring, parstring, optstring);
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
      fprintf(outfile, "# Index and name of bifurcation parameter #1                   : %2d (%s)\n", Bifparone, parameternames[Bifparone]);
      fprintf(outfile, "#\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
      colnr = 0;
#else
      colnr = 1;
#endif
      sprintf(tmpstr, "%2d:%s", colnr++, parameternames[Bifparone]);
      fprintf(outfile, "#%15s", tmpstr);
      for (i = 0; i < PopulationNr; i++)
        {
          sprintf(tmpstr, "%d:PGR[%d]", colnr++, i);
          fprintf(outfile, "%16s", tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          sprintf(tmpstr, "%d:Tc[%d]", colnr++, i);
          fprintf(outfile, "%16s", tmpstr);
          for (j = 0; j < ParameterNr; j++)
            {
              sprintf(tmpstr, "%d:S[%d][%d]", colnr++, i, j);
              fprintf(outfile, "%16s", tmpstr);
            }
        }
      fprintf(outfile, "\n");
      fflush(outfile);
    }
  else
    {
      STDOUT("#\n# Executing : ");

#if (defined(R_PACKAGE))
      STDOUT("PSPMdemo(\"%s\", %s, %s, %s)", progname, curvestring, parstring, optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      STDOUT("PSPMdemo('%s', %s, %s, %s)", progname, curvestring, parstring, optstring);
#else
      for (i = 0; i < argc; i++) STDOUT("%s ", argv[i]);
#endif
      STDOUT("\n#\n");

      STDOUT("# Parameter values  : \n#");
      for (i = 0; i < ParameterNr; i++)
        {
          if (!(i % 3)) STDOUT("\n# ");
          STDOUT("\t%-10s:", parameternames[i]);
          STDOUT("  %-13G", parameter[i]);
        }
      STDOUT("\n#\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
      colnr = 0;
#else
      colnr = 1;
#endif
      STDOUT("#");
      sprintf(tmpstr, "%d:PGR[0]", colnr++);
      STDOUT("%15s", tmpstr);
      for (i = 1; i < PopulationNr; i++)
        {
          sprintf(tmpstr, "%d:PGR[%d]", colnr++, i);
          STDOUT("%16s", tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          sprintf(tmpstr, "%d:Tc[%d]", colnr++, i);
          STDOUT("%16s", tmpstr);
          for (j = 0; j < ParameterNr; j++)
            {
              sprintf(tmpstr, "%d:S[%d][%d]", colnr++, i, j);
              STDOUT("%16s", tmpstr);
            }
        }
      STDOUT("\n#\n");
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  ReportMsg("\n\nPrediction :\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
  ReportMsg("\n");

  memset((void *)tanvec, 0, MaxPntDim*sizeof(double));
  tanvec[0] = 1.0;

  // Continue the curve
  while (1)
    {
      // Compute fixed point with new varied parameter
      cycles = 0;
      while (Stepreduce <= MAX_STEPREDUCE)
        {
          if (Bifparone == -1)
            retval = FindPoint(pntdim, point, NULL, NULL, DYTOL, RHSTOL, MAXITER, Equation);
          else
            {
              parameter[Bifparone] = point[0]*pnt_scale[0];
              savedbifpar          = Bifparone;
              Bifparone            = -1;
              savedpntr            = pnt_scale;
              pnt_scale            = savedpntr + 1;
              if (hasjac)
                retval  = FindPoint(pntdim - 1, point + 1, Jacmat + (pntdim - 1), NULL, DYTOL, RHSTOL, MAXITER, Equation);
              else
                retval  = FindPoint(pntdim - 1, point + 1, NULL, NULL, DYTOL, RHSTOL, MAXITER, Equation);
              Bifparone = savedbifpar;
              pnt_scale = savedpntr;
            }
          hasjac = 0;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif

          if (retval == SUCCES) break;

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
          COPY(pntdim, oldpoint, 1, point, 1);

          AXPY(pntdim, (Maxcurvestep/(Stepreduce*tanvec[0]*pnt_scale[0])), tanvec, 1, point, 1);

          ReportMsg("\n\nPrediction :\t");
          for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
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

      // Prediction
      ReportMsg("\nNew point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
      ReportMsg("\n");

      if ((retval = SetScales(point, pntdim)))
        {
          ReportMsg("\n\nVariable %d rescaled!\n", retval);
          Stepchange = 0;
        }
      else                                                                          // Otherwise generate output and predict new point on the curve
        {
          // Signal curve stop if one of the components has become negative or parameter is out of bounds
          // and this is not the curve beginning
          if (pntnr > 20) CurveEnd = CurveEnd || (parameter[Bifparone] <= minbound1) || (parameter[Bifparone] >= maxbound1);

          // Generate output: invoked after setting CurveEnd to allow for output of last solution point on branch
          if ((Bifparone != -1) && (((pntnr + 1) % ReportLevel) == 0))
            {
              for (i = 0; i < pntdim; i++)
                {
#if (defined(R_PACKAGE))
                  if (i)
                    STDOUT(",%15.8E", point[i]*pnt_scale[i]);
                  else
#endif
                    STDOUT("%16.8E",  point[i]*pnt_scale[i]);
                }
              STDOUT("\n");
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
            }

          outmax = DefineOutput(point, Output);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif
          if (outmax) PrettyPrintArray(outfile, outmax, Output);

          // End the continuation at end of curve after generation of last output point
          if ((CurveEnd) || (Bifparone == -1)) break;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif

          pntnr++;

          // Store info on current solution point
          COPY(pntdim, point, 1, oldpoint, 1);

          // Compute the new tangent vector
          retval = TangentVec(pntdim, point, Jacmat, tanvec, Equation, NULL, DYTOL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) break;
#endif
          hasjac = 1;

          Stepchange = 1;

          AXPY(pntdim, (Maxcurvestep/(Stepreduce*tanvec[0]*pnt_scale[0])), tanvec, 1, point, 1);
        }
      // Prediction
      ReportMsg("\n\nPrediction :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
      ReportMsg("\n");

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
#if (defined(_MSC_VER) && (_MSC_VER < 1500)) || (defined(R_PACKAGE) && defined(_WIN32))
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
  epsMach = dlamch("Epsilon");

  MaxPntDim         = PopulationNr + 1;
  CurPopulationNr   = PopulationNr;
  CohortDim         = IStateDim + 2;
  PopDensCohortDim  = CohortDim;
  Bifparone         = -1;

  TestRun       = 0;
  DoStateOutput = 0;
  SortIndex     = 0;
  FirstEstimate = 0;
  ReportLevel   = 1;
  
  CurveType         = PGR;
  
  eVarPntr          = PGRvar;
  parPntr           = parameter;

  // The following variables can be modified by the user with optional #define statements 
  LogMinSurvival      = log(MIN_SURVIVAL);
  CohortNr            = COHORT_NR;
  
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

#if (PULSED == 1)
  if (REPRODUCTION_INTERVAL <= Odesolve_Min_Step)
    {
      ErrorMsg(__FILE__, __LINE__, "For pulsed reproduction REPRODUCTION_INTERVAL should be set to a non-negligible, positive value!");
      STDOUT("\n");
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
      return;
    }
#endif

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
  fprintf(stderr, "\nUsage:\t%s [-par1 <index> <initial value> <step> <min. par.1> <max. par.1>] [-isort <index>] [-test]", progname);
  fprintf(stderr, "\n\nAim:\tComputation of the population growth rate and its parameter sensitivities for (multiple) structured "
                  "populations\n\tpossibly as a function of a model parameter\n\n");
  fprintf(stderr, "Possible options are:\n\n");
  fprintf(stderr, "\t-par1   <index>: Index of the parameter, as a function of which to compute the population growth rate\n");
  fprintf(stderr, "\t                 If this option is used <initial value>, <step>, <min. par.1> and <max. par.1> specify\n");
  fprintf(stderr, "\t                 the initial value, stepsize value, minimum and maximum value of this parameter.\n");
  fprintf(stderr, "\t                 If not specified, a single computation of the population growth rate and parameter sensitivities\n");
  fprintf(stderr, "\t                 is carried out for the current parameter values.\n");
  fprintf(stderr, "\t-isort  <index>: Index of i-state variable to use as ruling variable for sorting the structured populations.\n");
  fprintf(stderr, "\t-report <value>: Interval of reporting computed output to console. Minimum value of 1 implies output of every point.\n");
  fprintf(stderr, "\t-test          : Perform only a single integration over the life history, reporting dynamics of survival, R0\n");
  fprintf(stderr, "\t                 and i-state variables.\n");

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
  int   my_argc, tmpint, cvfdone = 0;
  char  **argpnt1 = NULL, **argpnt2 = NULL, *my_argv[argc];
  char  *ch;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = POPULATION_NR;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

  //=========================== Process the command-line arguments ===================================================================
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
      else if (!strcmp(*argpnt1, "-par1"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo index of first bifurcation parameter specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= ParameterNr))
            {
              fprintf(stderr, "\nIndex of first parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, ParameterNr);
              Usage(argv[0]);
            }
          Bifparone = tmpint;
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
              fprintf(stderr, "CVF file %s not read\n", runname);
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

  if (Bifparone != -1)
    {
      if (my_argc == 5)
        {
          // Map all command-line variables into the argument vector
          parameter[Bifparone] = atof(my_argv[1]);

          Maxcurvestep = atof(my_argv[2]);

          minbound1 = atof(my_argv[3]);
          maxbound1 = atof(my_argv[4]);
        }
      else
        Usage(my_argv[0]);
    }

  ComputeCurve(argc, argv);
      
  return 0;
}


/*
 *====================================================================================================================================
 *  Matlab interface and exit/cleanup function
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

// The gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t          nrows, ncols;
  double          curvepars[5], tmpdouble;
  int             tmpint;
  const mxArray   *cell_element_ptr;
  char            optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN];
  mwIndex         i;
  size_t          total_num_of_cells, buflen;
  int             status;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = POPULATION_NR;
  ParameterNr  = PARAMETER_NR;

  InitialiseVars();

#if (MFUNCTIONS == 1)
  Minterface_Init();
#endif

  // check for proper number of arguments
  if (nrhs != 3)
    mexErrMsgIdAndTxt("MATLAB:PSPMdemo:nrhs", "\nIncorrect number of command-line arguments.\n\nUse: %s(%s)\n\n%12s: %s\n%12s: %s\n%12s: %s",
                      mexFunctionName(), "<curve settings>, <parameters>, <options>", "<curve settings>",
                      "Index, initial value, step size, minimum and maximum value of the independent parameter for the curve", "<parameters>",
                      "Array of parameter values to use (empty array or same length as parameter array)", "<options>", "Possible options: test or isort");

  // check for proper number of output variables
  if (nlhs != 1) mexErrMsgIdAndTxt("MATLAB:PSPMdemo:nlhs", "A single output argument is required.");

  //============================== Process the options argument ======================================================================
  // Extract the contents of MATLAB cell into the C array
  if (!mxIsCell(prhs[2])) mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options", "\nOptions should be specified as a cell array!\n");

  total_num_of_cells = mxGetNumberOfElements(prhs[2]);
  strcpy(optstring, "{");
  for (i = 0; i < total_num_of_cells; i++)
    {
      cell_element_ptr = mxGetCell(prhs[2], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optname, buflen);
      if (!((!strcmp(optname, "test")) || (!strcmp(optname, "isort")) || (!strcmp(optname, "report"))))
        mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options", "\nIllegal option %s!\n", optname);

      if (!strcmp(optname, "test"))
        {
          TestRun = 1;

          // optstring still equal to "{"
          if (strlen(optstring) == 1) strcat(optstring, "'");
          else strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "'");
        }
      else if (!strcmp(optname, "isort"))
        {
          if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options", "\nNo value specified for option %s!\n", optname);

          cell_element_ptr = mxGetCell(prhs[2], i);
          buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
          status           = mxGetString(cell_element_ptr, optval, buflen);
          if (status) mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options", "\nError in retrieving value for option %s!\n", optname);

          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= IStateDim))
            {
              mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options",
                                "\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!\n",
                                tmpint, IStateDim);
              continue;
            }
          SortIndex = tmpint;

          // optstring still equal to "{"
          if (strlen(optstring) == 1) strcat(optstring, "'");
          else strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "', '");
          strcat(optstring, optval);
          strcat(optstring, "'");
        }
      else if (!strcmp(optname, "report"))
        {
          if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options", "\nNo value specified for option %s!\n", optname);

          cell_element_ptr = mxGetCell(prhs[2], i);
          buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
          status           = mxGetString(cell_element_ptr, optval, buflen);
          if (status) mexErrMsgIdAndTxt("MATLAB:PSPMdemo:options", "\nError in retrieving value for option %s!\n", optname);

          tmpint = atoi(optval);
          ReportLevel = max(tmpint, 1);

          // optstring still equal to "{"
          if (strlen(optstring) == 1) strcat(optstring, "'");
          else strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "', '");
          strcat(optstring, optval);
          strcat(optstring, "'");
        }
    }
  strcat(optstring, "}");

  //============================== Process the parameters argument ===================================================================

  nrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if ((ncols == ParameterNr) && (nrows == 1))
    memcpy(parameter, mxGetPr(prhs[1]), ncols*mxGetElementSize(prhs[1]));
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:PSPMdemo:parameters", "\nParameter argument ignored as it is not a row vector of length %d\n", ParameterNr);

  total_num_of_cells = mxGetNumberOfElements(prhs[1]);
  strcpy(parstring, "[");
  for (i = 0; i < total_num_of_cells; i++)
    {
      if (i) strcat(parstring, " ");
      memcpy(&tmpdouble, mxGetPr(prhs[1]) + i, mxGetElementSize(prhs[1]));
      sprintf(tmpstr, "%.6G", tmpdouble);
      strcat(parstring, tmpstr);
    }
  strcat(parstring, "]");

  //============================== Process the curve parameters argument =============================================================

  nrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if ((ncols == 5) && (nrows == 1))
    {
      memcpy(curvepars, mxGetPr(prhs[0]), ncols*mxGetElementSize(prhs[0]));

      tmpint = (int)floor(curvepars[0] + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        mexErrMsgIdAndTxt("MATLAB:PSPMdemo:curvepars", "\nIndex of first parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint,
                          ParameterNr);
      Bifparone            = tmpint;
      parameter[Bifparone] = curvepars[1];
      Maxcurvestep         = curvepars[2];

      minbound1 = curvepars[3];
      maxbound1 = curvepars[4];
    }
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:PSPMdemo:curvepars",
                       "\nCurve parameters ignored as they are not in the form [index, initial value, step size, minimum value, maximum value].\n");

  total_num_of_cells = mxGetNumberOfElements(prhs[0]);
  strcpy(curvestring, "[");
  for (i = 0; i < total_num_of_cells; i++)
    {
      if (i) strcat(curvestring, " ");
      memcpy(&tmpdouble, mxGetPr(prhs[0]) + i, mxGetElementSize(prhs[0]));
      sprintf(tmpstr, "%.6G", tmpdouble);
      strcat(curvestring, tmpstr);
    }
  strcat(curvestring, "]");

  //============================= Get the program name ===============================================================================
  // Get the name of the mex file
  strcpy(progname, mexFunctionName());
  progname[strlen(progname) - 4] = '\0';                                            // Cut off the 'demo' appendix

  mexAtExit(CloseStreams);

  // call the computational routine
  ComputeCurve(0, NULL);
  FreeHeapMemory();

  if (TestRun) plhs[0] = mxCreateString("");
  else plhs[0] = mxCreateString(runname);

#if (MFUNCTIONS == 1)
  Minterface_End();
#endif

  return;
}

/*
 *====================================================================================================================================
 *  R interface function
 *====================================================================================================================================
 */

#elif defined(R_PACKAGE)

SEXP PSPMdemo(SEXP moduleName, SEXP curveVals, SEXP parVals, SEXP optVals)

{
  int     i, ncols, tmpint;
  double  curvepars[5];
  char    optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN];
  SEXP    resfil;

  PopulationNr = POPULATION_NR;
  Stages       = STAGES;
  IStateDim    = I_STATE_DIM;
  EnvironDim   = POPULATION_NR;
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
      if (!((!strcmp(optname, "test")) || (!strcmp(optname, "isort")) || (!strcmp(optname, "report")))) error("\nIllegal option %s!\n\n", optname);

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

      if (!strcmp(optname, "isort"))
        {
          if (!(++i < ncols)) error("\nNo value specified for option %s!\n\n", optname);

          strcpy(optval, CHAR(STRING_ELT(optVals, i)));
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= IStateDim))
            error("\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  IStateDim);

          SortIndex = tmpint;
          if (!strlen(optstring)) strcat(optstring, "c(\"");
          else strcat(optstring, ", \"");
          strcat(optstring, optname);
          strcat(optstring, "\", \"");
          strcat(optstring, optval);
          strcat(optstring, "\"");
          continue;
        }

      if (!strcmp(optname, "report"))
        {
          if (!(++i < ncols))
            error("\nInterval for reporting computed output to console (option \"%s\") not specified!\n\n", optname);

          strcpy(optval, CHAR(STRING_ELT(optVals, i)));
          tmpint = atoi(optval);

          ReportLevel = max(tmpint, 1);
          if (!strlen(optstring)) strcat(optstring, "c(\"");
          else strcat(optstring, ", \"");
          strcat(optstring, optname);
          strcat(optstring, "\", \"");
          strcat(optstring, optval);
          strcat(optstring, "\"");
          continue;
        }
    }
  if (strlen(optstring)) strcat(optstring, ")");
  else strcpy(optstring, "NULL");

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

  //============================== Process the curve parameters argument =============================================================

  ncols = length(curveVals);
  if (isReal(curveVals) && (ncols == 5))
    {
      memcpy(curvepars, REAL(curveVals), ncols*sizeof(double));

      tmpint = (int)floor(curvepars[0] + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        error("\nIndex of first parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint, ParameterNr);

      Bifparone            = tmpint;
      parameter[Bifparone] = curvepars[1];
      Maxcurvestep         = curvepars[2];

      minbound1 = curvepars[3];
      maxbound1 = curvepars[4];
    }
  else if (ncols)
    warning("\nCurve parameters ignored as they are not in the form [index, initial value, step size, minimum value, maximum value].\n\n");

  if (ncols) strcpy(curvestring, "c(");
  else strcpy(curvestring, "NULL");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(curvestring, ", ");
      sprintf(tmpstr, "%.6G", curvepars[i]);
      strcat(curvestring, tmpstr);
    }
  if (ncols) strcat(curvestring, ")");

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
