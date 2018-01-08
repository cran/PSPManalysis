/***
  NAME
    PSPMequi

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

    Last modification: AMdR - Oct 25, 2017
***/

#define PSPMEQUI                  1

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
 *  Import the population and environment dimension settings
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
#error Equilibrium analysis requires ENVIRON_DIM to be larger than 0
#endif

#if !defined(INTERACT_DIM) || (INTERACT_DIM < 1)
#error INTERACT_DIM should be defined larger than 0
#endif

#if !defined(PARAMETER_NR) || (PARAMETER_NR < 3)
#error PARAMETER_NR should be defined larger than 2
#endif


/*
 *====================================================================================================================================
 *  Definition of problem dimensions
 *====================================================================================================================================
 */

#undef PULSED

/*
 *====================================================================================================================================
 *  Definition of global variables and parameters
 *====================================================================================================================================
 */

// Global dimension variables
static int                        BirthStateNr[POPULATION_NR + 1];
static int                        MaxPntDim;

// These are the variables to solve for
static double                     Evar[ENVIRON_DIM];
static double                     EnvEquiCondition[ENVIRON_DIM];

static int                        EnvTrivEqui[ENVIRON_DIM];
static int                        EnvResIndex[ENVIRON_DIM];
static int                        EnvPntIndex[ENVIRON_DIM];

static double                     R0[POPULATION_NR + 1];
static double                     Beq[POPULATION_NR + 1];
static double                     InteractVars[POPULATION_NR + 1][INTERACT_DIM];

static int                        PopTrivEqui[POPULATION_NR + 1];
static int                        R0ResIndex[POPULATION_NR + 1];
static int                        PopPntIndex[POPULATION_NR + 1];

// Global variables to hold variables shared among routines
static int                        LastMemAllocated = 0;
static int                        pntdim;

// Global pointers into the heap
static double                     *RightEigenvecMem = NULL;

#if (FULLSTATEOUTPUT > 0)
static double                     *BirthStateMem = NULL;
static double                     *PopDensMem    = NULL;
static int                        *CohortsMem    = NULL;

static double *CohortLimitMem = NULL;
#if (FULLSTATEOUTPUT == 1)
static double                     CohortMin[POPULATION_NR + 1], CohortMax[POPULATION_NR + 1];
#endif
#endif

// Global flags to tailor execution
static int                        TestRun       = 0;
static int                        DoStateOutput = 0;
static int                        SortIndex     = 0;
static int                        evoParsIndex[PARAMETER_NR];
static int                        DoSingle    = 0;
static int                        BPdetection = 1;
static int                        LPdetection = 1;
static int                        ParEVOIndex = -1;

// Global variables for other purposes
static char                       progname[MAXPATHLEN];
static double                     LogMinSurvival;
static char                       ContinuationString[MAX_STR_LEN];
static double                     initpnt[ENVIRON_DIM + POPULATION_NR + 1 + PARAMETER_NR];
static double                     pntmin[ENVIRON_DIM + POPULATION_NR + 1 + PARAMETER_NR];
static double                     pntmax[ENVIRON_DIM + POPULATION_NR + 1 + PARAMETER_NR];

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
static char                       parstring[MAX_STR_LEN];
static char                       optstring[MAX_STR_LEN];
static char                       curvestring[MAX_STR_LEN];
static char                       pntstring[MAX_STR_LEN];
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
  int     indx = 0, e, i, j, b = 0, p, retval = SUCCES, openmpDone = 0;
  int     MaxCohortDim = CohortDim + 1 + InteractDim;
  double  norm;
  double  *NextGenMatrix = NULL, *FinalIstateMem = NULL;
  int     savedTestRun;
  double  savedparval, savedbirthstatenr, savedbirthrate;
  double  *respntr;

  /*
   *===========================================================================
   * Map current estimate of solution to global variables
   *===========================================================================
   */
  if (setBifParVal) parameter[Bifparone] = argument[indx]*pnt_scale[indx];
  indx++;

  for (e = 0; e < EnvironDim; e++)
    {
      if (EnvTrivEqui[e] || (e == EnvBPIndex))
        Evar[e] = 0.0;
      else
        {
          Evar[e] = argument[indx]*pnt_scale[indx];
          indx++;
        }
    }

  for (p = 0; p < CurPopulationNr; p++)
    {
      if (PopTrivEqui[p] || (p == PopBPIndex) || (p == PopulationNr))
        Beq[p] = 0.0;
      else
        {
          Beq[p] = argument[indx]*pnt_scale[indx];
          indx++;
        }
    }

  if (CurveType == PIP)
    MutantParVal = argument[indx]*pnt_scale[indx];
  else if (CurveType == ESS)
    {
      for (i = 0; i < evoParsDim; i++)
        {
          parameter[evoParsIndex[i]] = argument[indx]*pnt_scale[indx];
          indx++;
        }
    }
  else if (CurveType != EQ)
    parameter[Bifpartwo] = argument[indx]*pnt_scale[indx];

  /*
   *===========================================================================
   * Set the dimensions
   *===========================================================================
   */
  if (CurveType == PIP)
    {
      CurPopulationNr = PopulationNr + 1;
      savedparval     = parameter[Bifparone];

      parameter[Bifparone] = MutantParVal;
      // Set the mutant number of birth states
      SetBirthStates(BirthStateNr, Evar);
      BirthStateNr[PopulationNr] = BirthStateNr[PopEVOIndex];

      // Restore the resident parameter value
      parameter[Bifparone] = savedparval;
    }
  else
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
      STDOUT("\n\nParameter #1:                %15.6G", parameter[Bifparone]);
      for (e = 0; e < EnvironDim; e++) STDOUT("\nEnvironment variable #%d:     %15.6G", e, Evar[e]);
      for (p = 0; p < CurPopulationNr; p++) STDOUT("\nBirth rate of population #%d: %15.6G", p, Beq[p]);
      if (CurveType == PIP)
        STDOUT("\nMutant parameter value:      %15.6G", MutantParVal);
      else if ((CurveType != EQ) && (CurveType != ESS))
        STDOUT("\nParameter #2:                %15.6G", parameter[Bifpartwo]);

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
  if (CurveType == PIP)
    {
      // Save the resident number of birth states
      savedparval       = parameter[Bifparone];
      savedbirthstatenr = BirthStateNr[PopEVOIndex];
      savedbirthrate    = Beq[PopEVOIndex];
      Beq[PopEVOIndex]  = 0.0;

      // Set the mutant parameter value and number of birth states
      parameter[Bifparone]      = MutantParVal;
      BirthStateNr[PopEVOIndex] = BirthStateNr[PopulationNr];

      openmpDone = 0;
#if (defined(OPENMP) && (RFUNCTIONS != 1) && (MFUNCTIONS != 1))                     // If using R-defined or Matlab-defined model this parallelization is impossible
      if (!TestRun)
        {
#pragma omp parallel for private(i) if (MaxStatesAtBirth > 1)                       // Only use threading with multiple states at birth
          for (b = 0; b < BirthStateNr[PopEVOIndex]; b++)
            {
              if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;

              memcpy(FinalIstatePnt(b, PopulationNr, 0), FinalIstatePnt(b, PopEVOIndex, 0), (CohortDim + BirthStateNr[PopEVOIndex])*sizeof(double));
#if (FULLSTATEOUTPUT > 0)
              if (DoStateOutput)
                {
                  memcpy(BirthStatePnt(PopulationNr, b, 0), BirthStatePnt(PopEVOIndex, b, 0), IStateDim*sizeof(double));
                  for (i = 0; i < CohortDim; i++)                                  // Do not copy the zero population density
                    memcpy(&(PopDens(PopulationNr, b, i, 0)), &(PopDens(PopEVOIndex, b, i, 0)), Cohorts(PopEVOIndex, b)*sizeof(double));
                  Cohorts(PopulationNr, b) = Cohorts(PopEVOIndex, b);
                }
#endif
            }
          openmpDone = 1;
        }
#endif
      if (!openmpDone)
        {
          for (b = 0; b < BirthStateNr[PopEVOIndex]; b++)
            {
              if (!LifeHistory(BirthStateNr, b, MaxCohortDim, FinalIstatePnt(b, 0, 0))) retval = FAILURE;

              memcpy(FinalIstatePnt(b, PopulationNr, 0), FinalIstatePnt(b, PopEVOIndex, 0), (CohortDim + BirthStateNr[PopEVOIndex])*sizeof(double));
#if (FULLSTATEOUTPUT > 0)
              if (DoStateOutput)
                {
                  memcpy(BirthStatePnt(PopulationNr, b, 0), BirthStatePnt(PopEVOIndex, b, 0), IStateDim*sizeof(double));
                  for (i = 0; i < CohortDim; i++)                                  // Do not copy the zero population density
                    memcpy(&(PopDens(PopulationNr, b, i, 0)), &(PopDens(PopEVOIndex, b, i, 0)), Cohorts(PopEVOIndex, b)*sizeof(double));
                  Cohorts(PopulationNr, b) = Cohorts(PopEVOIndex, b);
                }
#endif
            }
        }

      // Restore the resident parameter value and number of birth states
      parameter[Bifparone]      = savedparval;
      BirthStateNr[PopEVOIndex] = savedbirthstatenr;
      Beq[PopEVOIndex]          = savedbirthrate;
    }

  openmpDone = 0;
#ifdef OPENMP
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
          R0[p] = FinalIstate(0, p, CohortDim);
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
            for (j = 0; j < BirthStateNr[p]; j++) NextGenMatrix[b*BirthStateNr[p] + j] = FinalIstate(j, p, CohortDim + b);

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
              STDOUT("%15.6G", ASUM(BirthStateNr[p], FinalIstatePnt(b, p, CohortDim), 1));
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
   * Remember parameter has index 0, hence the first environment variable has index 1 in pnt_scale
   */
  indx    = 1;
  respntr = result;
  for (e = 0; e < EnvironDim; e++)
    {
      if (EnvTrivEqui[e] && (e != EnvBPIndex)) continue;
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
      if (PopTrivEqui[p] && (p != PopBPIndex) && (p != PopulationNr)) continue;
      *respntr = R0[p] - 1.0;
      respntr++;
    }
  if ((CurveType == LP) && (!DoStateOutput))
    {
      savedTestRun = TestRun;
      TestRun      = 0;
      retval       = LPcondition(pntdim, argument, Equation, CENTRAL, 1, respntr, DYTOL);
      TestRun      = savedTestRun;
    }
  else if ((CurveType == ESS) && (!DoStateOutput))
    {
      savedTestRun = TestRun;
      TestRun      = 0;
      for (i = 0; i < evoParsDim; i++)
        {
          retval = SelectionGradient(pntdim, argument, Equation, (pntdim - evoParsDim) + i, R0ResIndex[PopEVOIndex], respntr);
          respntr++;
        }
      TestRun = savedTestRun;
    }

  // Add additional conditions in case of LP or ESS localization. Occurs only during EQ or ESS continuation, when called from LocateLP() or locateESS()
  savedTestRun = TestRun;
  TestRun      = 0;
  if (LocalizeType == LP)
    retval = LPcondition(pntdim, argument, Equation, CENTRAL, 0, respntr, DYTOL);
  else if (LocalizeType == ESS)
    retval = SelectionGradient(pntdim, argument, Equation, 0, R0ResIndex[PopEVOIndex], respntr);
  TestRun  = savedTestRun;

  if (FinalIstateMem) free(FinalIstateMem);
  if (NextGenMatrix) free(NextGenMatrix);
  return retval;
}


/*==================================================================================================================================*/

int DefineOutput(double *x, double *output)

{
  int     outnr = 0, i, j, retval;
  double  result[MaxPntDim];
  double  evJ, evH, zCz;

#if (FULLSTATEOUTPUT > 0)
  DoStateOutput = 1;
  Equation(x, result);
  DoStateOutput = 0;
  WriteStateToFile(FULLSTATEOUTPUT);
  if ((CurveType == LP) || (CurveType == ESS) || (CurveType == PIP)) Equation(x, result);
#else
  Equation(x, result);
#endif

  // There are maximally (EnvironDim+PopulationNr+2) values in the point vector
  // to solve for, which occurs when continuing an LP for this system. In all cases
  // the (EnvironDim+PopulationNr+2) values are written to the output file
  output[outnr++] = parameter[Bifparone];
  for (i = 0; i < EnvironDim; i++) output[outnr++] = Evar[i];
  for (i = 0; i < PopulationNr; i++) output[outnr++] = Beq[i];
  if (CurveType == PIP)
    output[outnr++] = MutantParVal;
  else if (CurveType == ESS)
    {
      for (i = 0; i < evoParsDim; i++) output[outnr++] = parameter[evoParsIndex[i]];
    }
  else if (CurveType != EQ)
    output[outnr++] = parameter[Bifpartwo];

  // Also add all interaction variables
  for (i = 0; i < PopulationNr; i++)
    for (j = 0; j < InteractDim; j++) output[outnr++] = InteractVars[i][j];

  for (i = 0; i < EnvironDim; i++)
    if (EnvironmentType[i] == PERCAPITARATE) output[outnr++] = EnvEquiCondition[i];

  for (i = 0; i < CurPopulationNr; i++) output[outnr++] = R0[i];

  if ((CurveType != PIP) && (PopEVOIndex != -1) && (R0ResIndex[PopEVOIndex] >= 0))
    {
      if (ParEVOIndex < 0)
        retval = SelectionGradient(pntdim, x, Equation, 0, R0ResIndex[PopEVOIndex], output + outnr);
      else
        retval = SelectionGradient(pntdim, x, Equation, pntdim + ParEVOIndex, R0ResIndex[PopEVOIndex], output + outnr);
      outnr++;
    }

  if (CurveType == ESS)
    {
      retval = ESSclassify(pntdim, x, Equation, DYTOL, RHSTOL, R0ResIndex[PopEVOIndex], &evJ, &evH, &zCz, 0);
      if (retval == SUCCES)
        {
          output[outnr++]                     = evJ;
          output[outnr++]                     = evH;
          if (evoParsDim > 1) output[outnr++] = zCz;
        }
      else
        {
          output[outnr++]                     = INFINITY;
          output[outnr++]                     = INFINITY;
          if (evoParsDim > 1) output[outnr++] = INFINITY;
        }
    }
  // Should we indeed add how accurate the solution was?
  output[outnr++] = NRM2(pntdim - 1, result, 1);

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
  int           pntnr = 0, outmax, CurveEnd = 0, hasjac = 0, TestBifs = 0, DoOutput = 1;
  int           cycles, last = 1, retval = 0;
  double        oldpoint[MaxPntDim];
  double        point[MaxPntDim], tanvec[MaxPntDim], oldvec[MaxPntDim], Jacmat[MaxPntDim*MaxPntDim];
  double        newEEC[EnvironDim], oldEEC[EnvironDim];
  double        newR0[PopulationNr + 1], oldR0[PopulationNr + 1];
  double        newdR0dp[PopulationNr + 1], olddR0dp[PopulationNr + 1];
  double        R0_xx, R0_yy, zCz;
  char          bifname[MAX_STR_LEN], csbname[MAX_STR_LEN], errname[MAX_STR_LEN], outname[MAX_STR_LEN];
  char          tmpstr[MAX_STR_LEN];
  struct stat   buffer;

#if defined(R_PACKAGE)
  STDOUT("\n");
#else
  fprintf(stderr, "\n");
#endif

  for (i = 0; i < (PopulationNr + 1); i++) oldR0[i] = 1.0;
  for (i = 0; i < (PopulationNr + 1); i++) newR0[i] = 1.0;
  memset(olddR0dp, 0, (PopulationNr + 1)*sizeof(double));
  memset(newdR0dp, 0, (PopulationNr + 1)*sizeof(double));

  COPY(pntdim, initpnt, 1, point, 1);
  (void)SetScales(point, pntdim);

  if (TestRun)
    {
      double rhs[MaxPntDim], rhsnorm;

#if defined(R_PACKAGE)
      STDOUT("\n\nExecuting : ");
      STDOUT("PSPMequi(\"%s\", '%s', %s, %s, %s, %s)", progname, ContinuationString, pntstring, curvestring, parstring, optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      STDOUT("\n\nExecuting : ");
      STDOUT("PSPMequi('%s', '%s', %s, %s, %s, %s)", progname, ContinuationString, pntstring, curvestring, parstring, optstring);
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

  if (strlen(runname))
    {
      sprintf(bifname, "%s.bif", runname);
      sprintf(errname, "%s.err", runname);
      sprintf(outname, "%s.out", runname);
    }
  else
    {
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)
      strcpy(progname, argv[0]);
      progname[strlen(progname) - 4] = '\0';                                        // Cut off the 'equi' appendix
#endif
      i = 0;
      while (1)
        {
          sprintf(bifname, "%s-%s-%04d.bif", progname, ContinuationString, i);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
          sprintf(csbname, "%s-%s-%04d.mat", progname, ContinuationString, i);
#else
          sprintf(csbname, "%s-%s-%04d.csb", progname, ContinuationString, i);
#endif
          sprintf(errname, "%s-%s-%04d.err", progname, ContinuationString, i);
          sprintf(outname, "%s-%s-%04d.out", progname, ContinuationString, i);
          if (stat(bifname, &buffer) && stat(csbname, &buffer) && stat(errname, &buffer) && stat(outname, &buffer)) break;
          i++;
        }
      sprintf(runname, "%s-%s-%04d", progname, ContinuationString, i);
    }

  if ((CurveType == EQ) || (CurveType == ESS)) biffile = fopen(bifname, "w");
  errfile                                              = fopen(errname, "w");
  outfile                                              = fopen(outname, "w");

  if (outfile)
    {
      fprintf(outfile, "#\n# Executing : ");

#if defined(R_PACKAGE)
      fprintf(outfile, "PSPMequi(\"%s\", \"%s\", %s, %G, %s, %s, %s)", progname, ContinuationString, pntstring, Maxcurvestep, curvestring, parstring,
              optstring);
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      fprintf(outfile, "PSPMequi('%s', '%s', %s, %G, %s, %s, %s)", progname, ContinuationString, pntstring, Maxcurvestep, curvestring, parstring,
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
      fprintf(outfile, "# Index and name of bifurcation parameter #1                   : %2d (%s)\n", Bifparone, parameternames[Bifparone]);
      if ((CurveType == BP) || (CurveType == LP) || (CurveType == BPE))
        fprintf(outfile, "# Index and name of bifurcation parameter #2                   : %2d (%s)\n", Bifpartwo, parameternames[Bifpartwo]);
      if (CurveType == ESS)
        {
          for (i = 0; i < evoParsDim; i++)
            fprintf(outfile, "# Index and name of parameter #%d at ESS value                  : %2d (%s)\n", i + 1, evoParsIndex[i],
                    parameternames[evoParsIndex[i]]);
        }
      if (CurveType == BPE) fprintf(outfile, "# Index of environment variable with transcritical bifurcation : %d\n", EnvBPIndex);
      if (CurveType == BP) fprintf(outfile, "# Index of structured population with transcritical bifurcation: %d\n", PopBPIndex);
      if ((CurveType == ESS) || (CurveType == PIP))
        fprintf(outfile, "# Index of structured population for evolutionary analysis     : %d\n", PopEVOIndex);
      if ((CurveType == EQ) && (PopEVOIndex != -1))
        {
          fprintf(outfile, "# Index of structured population for evolutionary analysis     : %d\n", PopEVOIndex);
          if (ParEVOIndex < 0)
            fprintf(outfile, "# Index of parameter for computation of selection gradient     : %d\n", Bifparone);
          else
            fprintf(outfile, "# Index of parameter for computation of selection gradient     : %d\n", ParEVOIndex);
        }
      fprintf(outfile, "#\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
      colnr = 0;
#else
      colnr = 1;
#endif
      sprintf(tmpstr, "%2d:%s", colnr++, parameternames[Bifparone]);
      fprintf(outfile, "#%15s", tmpstr);
      for (i = 0; i < EnvironDim; i++) fprintf(outfile, "        %2d:E[%2d]", colnr++, i);
      for (i = 0; i < PopulationNr; i++) fprintf(outfile, "        %2d:b[%2d]", colnr++, i);
      if (CurveType == ESS)
        {
          for (i = 0; i < evoParsDim; i++)
            {
              sprintf(tmpstr, "%2d:%s", colnr++, parameternames[evoParsIndex[i]]);
              fprintf(outfile, " %15s", tmpstr);
            }
        }
      else if (CurveType == PIP)
        {
          sprintf(tmpstr, "%2d:%s'", colnr++, parameternames[Bifpartwo]);
          fprintf(outfile, " %15s", tmpstr);
        }
      else if (CurveType != EQ)
        {
          sprintf(tmpstr, "%2d:%s", colnr++, parameternames[Bifpartwo]);
          fprintf(outfile, " %15s", tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        for (j = 0; j < InteractDim; j++) fprintf(outfile, "    %2d:I[%2d][%2d]", colnr++, i, j);
      for (i = 0; i < EnvironDim; i++)
        if (EnvironmentType[i] == PERCAPITARATE) fprintf(outfile, "     %2d:pcgE[%2d]", colnr++, i);
      for (i = 0; i < CurPopulationNr; i++) fprintf(outfile, "       %2d:R0[%2d]", colnr++, i);
      if ((CurveType != PIP) && (PopEVOIndex != -1) && (R0ResIndex[PopEVOIndex] >= 0)) fprintf(outfile, "     %2d:R0_x[%2d]", colnr++, Bifparone);
      if (CurveType == ESS)
        {
          if (evoParsDim == 1)
            {
              fprintf(outfile, "    %2d:R0_xx[%2d]", colnr++, PopEVOIndex);
              fprintf(outfile, "    %2d:R0_yy[%2d]", colnr++, PopEVOIndex);
            }
          else
            {
              fprintf(outfile, " %2d:Max. EigVal J", colnr++);
              fprintf(outfile, " %2d:Max. EigVal H", colnr++);
              fprintf(outfile, "     %2d:Z^T C01 Z", colnr++);
            }
        }
      fprintf(outfile, "     %2d:RHS norm\n", colnr++);
      fprintf(outfile, "#\n");
      fflush(outfile);
    }

  memset((void *)tanvec, 0, pntdim*sizeof(double));
  tanvec[0] = 1.0;

  // Continue the curve
  while (1)
    {
      // Compute fixed point with new varied parameter
      cycles = 0;
      while (Stepreduce <= MAX_STEPREDUCE)
        {
          if (hasjac)
            retval = FindPoint(pntdim, point, Jacmat, tanvec, DYTOL, RHSTOL, MAXITER, Equation);
          else
            retval = FindPoint(pntdim, point, NULL, tanvec, DYTOL, RHSTOL, MAXITER, Equation);
          hasjac   = 0;

          if (retval == SUCCES) break;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) return;
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
          COPY(pntdim, oldpoint, 1, point, 1);

          // AXPY(pntdim, (curvestep/(Stepreduce*pnt_scale[0])), tanvec, 1, point, 1);
          AXPY(pntdim, (curvestep/Stepreduce), tanvec, 1, point, 1);

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

      // Report on located point to stderr file
      ReportMsg("New point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
      ReportMsg("\n");

      if (!DoSingle)                                                                // Single point computation or curve continuation
        {
          // Increase step when both this and previous point were located with
          // the current step size
          if ((last && !cycles) && (Stepreduce > 1)) Stepreduce /= 2;
          last = (!cycles);

          // Signal curve stop if one of the components has become negative or parameter is out of bounds
          // and this is not the curve beginning
          if (pntnr > 20)
            for (i = 0; i < pntdim; i++)
              CurveEnd = CurveEnd || ((point[i]*pnt_scale[i]) < pntmin[i] - epsMach) || ((point[i]*pnt_scale[i]) > pntmax[i] + epsMach);

          // The R0 values and equilibrium conditions of the environmental variables are being recomputed at every call to Equation()
          // Therefore save these values for later use in bifurcation detection
          COPY(EnvironDim, EnvEquiCondition, 1, newEEC, 1);
          COPY(CurPopulationNr, R0, 1, newR0, 1);

          // Compute the new tangent vector
          retval = TangentVec(pntdim, point, Jacmat, tanvec, Equation, NULL, DYTOL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) return;
#endif
          hasjac = 1;
          /*
           * ===================================================================================================
           * The following code section is the only section that is specific to the problem of continuation
           * of structured population curves. All other parts of the main routine are generic.
           * The following lines detect bifurcation points
           * ===================================================================================================
           */
          if (TestBifs)
            {
              //******** LAST INDICES IN THE NEXT 3 FUNCTION CALLS HAVE TO BE CORRECTED FOR TOTAL VARS IN TRIVIAL
              for (i = 0; i < PopulationNr; i++)
                {
                  if (BPdetection && (PopTrivEqui[i] && ((oldR0[i] - 1)*(newR0[i] - 1) < -RHSTOL*RHSTOL)))
                    {
                      LPdetection = 0;
                      LocateBP(&pntdim, point, Equation, DYTOL, RHSTOL, i, -1);
                      newR0[i] = 1.0;
                    }
                  else if (PopTrivEqui[i])
                    continue;
                  else if (BPdetection && (oldpoint[PopPntIndex[i]]*point[PopPntIndex[i]] < -DYTOL*DYTOL))
                    {
                      LPdetection = 0;
                      LocateBP(&pntdim, point, Equation, DYTOL, RHSTOL, i, PopPntIndex[i]);
                    }
                  else if ((i == PopEVOIndex) && (R0ResIndex[PopEVOIndex] >= 0))
                    {
                      retval = SelectionGradient(pntdim, point, Equation, 0, R0ResIndex[PopEVOIndex], newdR0dp + PopEVOIndex);
                      if ((retval == SUCCES) && (olddR0dp[PopEVOIndex]*newdR0dp[PopEVOIndex] < -RHSTOL*RHSTOL))
                        LocateESS(pntdim, point, Equation, DYTOL, RHSTOL, R0ResIndex[PopEVOIndex]);
                    }
                }
              for (i = 0; BPdetection && (i < EnvironDim); i++)
                if (EnvironmentType[i] == PERCAPITARATE)
                  {
                    if (EnvTrivEqui[i] && (oldEEC[i]*newEEC[i] < -RHSTOL*RHSTOL))
                      {
                        LPdetection = 0;
                        LocateBPE(&pntdim, point, Equation, DYTOL, RHSTOL, i, -1);
                      }
                    else if (EnvTrivEqui[i])
                      continue;
                    else if (oldpoint[EnvPntIndex[i]]*point[EnvPntIndex[i]] < -DYTOL*DYTOL)
                      {
                        LPdetection = 0;
                        LocateBPE(&pntdim, point, Equation, DYTOL, RHSTOL, i, EnvPntIndex[i]);
                      }
                  }
              if (LPdetection && (tanvec[0]*oldvec[0] < -DYTOL*DYTOL)) LocateLP(pntdim, point, Equation, DYTOL, RHSTOL);
            }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) return;
#endif

          // When during the construction of a PIP we end up in the ESS on the line y=x,
          // make the tangent vector to the curve perpendicular to the line y=x to force the
          // solution away from the diagonal
          if ((CurveType == PIP) && (fabs(point[0] - point[pntdim - 1]) < DYTOL) && fabs(newdR0dp[PopEVOIndex]) < RHSTOL)
            {
              // By definition during a PIP continuation no other parameter are at their ESS value, hence evoParsDim = 0
              retval = ESSclassify(pntdim, point, Equation, DYTOL, RHSTOL, R0ResIndex[PopEVOIndex], &R0_xx, &R0_yy, &zCz, 1);
              /*
               * Notice that R0(E(x),y) with x and y the resident and mutant trait, respectively, can locally be approximated as:
               *
               * R0(E(x+dx),y+dy) = R0(E(x),y) + R0_x(E(x),y) dx + R0_y(E(x),y) dy + R0_xx(E(x),y)/2 dx^2 + R0_xy(E(x),y) dxdy + R0_yy(E(x),y)/2 dy^2 + h.o.t
               *
               * Because in equilibrium x=y=x* necessarily R0 = 1, both R0(E(x+dx),y+dy) and R0(E(x),y) equal 1
               * (we are always in equilibrium!!!) and for the above relation to hold  it is also necessarily true that:
               *
               *      R0_x(E(x*),x*) + R0_y(E(x*),x*) = 0                               (1)
               *
               *  and
               *
               *      R0_xx(E(x*),x*) + 2 R0_xy(E(x*),x*) + R0_yy(E(x*),x*) = 0         (2)
               *
               *  Because we are in an evolutionary singularity, it moreover holds that:
               *
               *      R0_x(E(x*),x*) = - R0_y(E(x*),x*) = 0                             (3)
               *
               *  and
               *
               *      R0_xy(E(x*),x*) = -(R0_xx(E(x*),x*) + R0_yy(E(x*),x*))/2
               *
               *  Therefore,
               *
               *      R0_xx(E(x),y)/2 dx^2 + R0_xy(E(x),y) dxdy + R0_yy(E(x),y)/2 dy^2 = (R0_xx(E(x),y) dx - R0_yy(E(x),y) dy)(dx - dy)/2 = 0
               *
               *  Therefore, dy/dx = 1 (which is the diagonal) or dy/dx = R0_xx(E(x),y)/R0_yy(E(x),y), i.e.
               *
               *      dy = (R0_xx(E(x),y)/R0_yy(E(x),y)) dx
               *
               */
              if ((retval == SUCCES) && R0_yy)
                {
                  tanvec[pntdim - 1] = zCz*tanvec[0];
                }
              else
                {
                  STDOUT("\nTangent vector element for mutant reversed to step away from the PIP diagonal\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
                  mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
                  R_FlushConsole();
                  R_ProcessEvents();
#endif
                  tanvec[pntdim - 1] = -tanvec[0];
                }
              hasjac = 0;
            }
        }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      if (checkInterrupt()) return;
#endif

      // Generate output: invoked after setting CurveEnd to allow for output of last solution point on branch
      if (DoOutput)
        {
          // Output located point to stdout and output file
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
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
          mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
          R_FlushConsole();
          R_ProcessEvents();
#endif

          outmax = DefineOutput(point, Output);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) return;
#endif
          if (outmax) PrettyPrintArray(outfile, outmax, Output);
        }

      // End the continuation at end of curve after generation of last output point
      if (CurveEnd || DoSingle) break;

      // Store all info on current solution point
      COPY(pntdim, point, 1, oldpoint, 1);
      COPY(pntdim, tanvec, 1, oldvec, 1);
      COPY(EnvironDim, newEEC, 1, oldEEC, 1);
      COPY(CurPopulationNr, newR0, 1, oldR0, 1);
      COPY(CurPopulationNr, newdR0dp, 1, olddR0dp, 1);

      // Scale the point vector anew if necessary and redo the current point
      if ((retval = SetScales(point, pntdim)))
        {
          ReportMsg("\n\nVariable %d rescaled!\n", retval);

          // Compute the new tangent vector
          // If the inner product with the previous tangent vector is negative,
          // the direction has been reversed and is hence not preserved. Therefore,
          // scale the tangent vector with a factor -1 to flip the direction.
          retval = TangentVec(pntdim, point, Jacmat, tanvec, Equation, NULL, DYTOL);
          hasjac = 1;
          if (DOT(pntdim, oldvec, 1, tanvec, 1) < 0) SCAL(pntdim, -1, tanvec, 1);

          Stepchange = DoOutput = 0;
          TestBifs              = 0;
        }
      // Otherwise generate output and predict new point on the curve
      else
        {
          Stepchange = DoOutput = 1;
          TestBifs              = (CurveType == EQ) || (CurveType == ESS);

          // Determine appropriate stepsize
          ReportMsg("\nTangent vector in component   0: %.8G\n", tanvec[0]);
          ReportMsg("Targeted  step in component   0: %.8G\n", curvestep*tanvec[0]/Stepreduce);

          if ((Stepreduce == 1) && (fabs(curvestep) > fabs(Maxcurvestep))) curvestep = sign(curvestep)*fabs(Maxcurvestep);

          AXPY(pntdim, (curvestep/Stepreduce), tanvec, 1, point, 1);

          ReportMsg("Realized  step in component   0: %.8G\n", curvestep*tanvec[0]/Stepreduce);
        }

      // Prediction
      ReportMsg("\n************************************************************************************************************\n");
      ReportMsg("\nPrediction :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", point[i]*pnt_scale[i]);
      ReportMsg("\n");

      pntnr++;

      fflush(NULL);
#if defined(R_PACKAGE)
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

#if (defined(_MSC_VER) && (_MSC_VER < 1500)) || (defined(R_PACKAGE) && defined(_WIN32))
  (void)_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  // Initialize some variables
  errfile    = NULL;
  outfile    = NULL;
  Stepchange = 0;
  Stepreduce = 1;
  strcpy(runname, "");

  // Get the machine precisions
  epsMach = dlamch("Epsilon");

  MaxPntDim         = EnvironDim + PopulationNr + 1 + ParameterNr;
  CurPopulationNr   = PopulationNr;
  CohortDim         = IStateDim + 1;
  PopDensCohortDim  = CohortDim;
  evoParsDim        = 0;

  PopBPIndex        = -1;
  EnvBPIndex        = -1;
  PopEVOIndex       = -1;
  ParEVOIndex       = -1;
  Bifparone         = -1;
  Bifpartwo         = -1;
  setBifParVal      = 1;
  
  CurveType         = UNDEFINED;
  LocalizeType      = UNDEFINED;
  
  eVarPntr          = Evar;
  birthRatePntr     = Beq;
  parPntr           = parameter;
  evoParsIndexPntr  = evoParsIndex;

  memset(PopTrivEqui, 0, (PopulationNr + 1)*sizeof(int));
  memset(EnvTrivEqui, 0, EnvironDim*sizeof(int));
#if (ALLOWNEGATIVE)
  for (i = 0; i < MaxPntDim; i++) pntmin[i] = -SAFETY*DBL_MAX;
#else
  for (i = 0; i < MaxPntDim; i++) pntmin[i] = 0.0;
#endif
  for (i = 0; i < MaxPntDim; i++) pntmax[i] = SAFETY*DBL_MAX;
  for (i = 0; i < ParameterNr; i++) evoParsIndex[i] = -1;


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
  char  bifstr[MAX_STR_LEN];
  char  desc[MAX_STR_LEN];
  char  varstr[MAX_STR_LEN];

  switch (CurveType)
    {
      case BP:
        strcpy(bifstr, "BP");
        strcpy(desc, "Aim:\tContinuation of a transcritical bifurcation of a structured population\n\tas a function of two parameters");
        break;
      case EQ:
        strcpy(bifstr, "EQ");
        strcpy(desc, "Aim:\tContinuation of a non-trivial equilibrium of a structured population\n\tas a function of a single parameter");
        break;
      case BPE:
        strcpy(bifstr, "BPE");
        strcpy(desc, "Aim:\tContinuation of a transcritical bifurcation in one of the environment variables\n\t");
        strcat(desc, "of a structured population as a function of two parameters.\n\tOnly works if the dynamics of the environment variable\n\tis "
                     "linear in the variable itself");
        break;
      case LP:
        strcpy(bifstr, "LP");
        strcpy(desc, "Aim:\tContinuation of a saddle-node bifurcation of a structured population\n\tas a function of two parameters");
        break;
      case ESS:
        strcpy(bifstr, "ESS");
        strcpy(desc, "Aim:\tContinuation of the ESS value(s) of life history parameter(s) for a structured population\n\tas a function of the first "
                     "bifurcation parameter");
        break;
      case PIP:
        strcpy(bifstr, "PIP");
        strcpy(desc, "Aim:\tContinuation of the boundary in resident-mutant values of the first parameter (PIP)\n\twhere R0=1 for both the resident "
                     "and mutant structured population");
        break;
      default:
        strcpy(bifstr, "<Type>");
        strcpy(desc, "Aim:\tContinuation of trivial or non-trivial equilibria, transcritical and saddle-node bifurcations\n\t");
        strcat(desc, "of structured populations as well as transcritical bifurcation in one of its environment variables\n\tand evolutionary "
                     "continuation as a function of one or two parameters");
    }
  if (CurveType == UNDEFINED)
    strcpy(varstr, "<Initial values>");
  else
    {
      strcpy(varstr, "Par.1");
      for (i = 0; i < EnvironDim; i++)
        {
          if (CurveType == BPE)
            {
              if ((i == EnvBPIndex) || (EnvTrivEqui[i])) continue;
            }
          else if (EnvTrivEqui[i])
            continue;
          sprintf(tmpstr, " E[%d]", i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          if (CurveType == BP)
            {
              if ((i == PopBPIndex) || (PopTrivEqui[i])) continue;
            }
          else if (PopTrivEqui[i])
            continue;
          sprintf(tmpstr, " b[%d]", i);
          strcat(varstr, tmpstr);
        }
      if (CurveType == ESS)
        strcat(varstr, " Par.2 ... Par.n");
      else if (!(CurveType == EQ))
        strcat(varstr, " Par.2");
    }

  fprintf(stderr, "Usage:\t%s [<options>] %s %s %s", progname, bifstr, varstr, "<max. parameter step> <index par.1> <min. par.1> <max. par.1>");
  if (CurveType == UNDEFINED)
    fprintf(stderr, " [<index par.2> <min. par.2> <max. par.2> ... <index par.n> <min. par.n> <max. par.n>]");
  else if (CurveType == ESS)
    fprintf(stderr, " <index par.2> <min. par.2> <max. par.2> ... <index par.n> <min. par.n> <max. par.n>");
  else if (!(CurveType == EQ))
    fprintf(stderr, " <index par.2> <min. par.2> <max. par.2>");
  fprintf(stderr, "\n\n%s\n\n", desc);
  if (CurveType == UNDEFINED) fprintf(stderr, "<Type>: Type of curve computation to be performed, either BP, EQ, LP, BPE, ESS or PIP\n\n");
  fprintf(stderr, "Possible options are:\n\n");
  fprintf(stderr, "\t-envBP   <index> : Index of environment variable, of which to continue the transcritical bifurcation\n");
  fprintf(stderr, "\t-popBP   <index> : Index of structured population, of which to continue the transcritical bifurcation\n");
  fprintf(stderr,
          "\t-popEVO  <index> : Index of structured population, for which to compute selection gradient or perform ESS or PIP continuation\n");
  fprintf(stderr, "\t-parEVO  <index> : Index of parameter to compute selection gradient for during EQ continuation\n");
  fprintf(stderr, "\t-envZE   <index> : Index of environment variable in trivial equilibrium (can be used multiple times)\n");
  fprintf(stderr, "\t-popZE   <index> : Index of structured population in trivial equilibrium (can be used multiple times)\n");
  fprintf(stderr, "\t-evoPars <number>: Number of life history parameters of structured population at their ESS value\n");
  fprintf(stderr, "\t-isort   <index> : Index of i-state variable to use as ruling variable for sorting the structured populations\n");
  fprintf(stderr, "\t-noBP            : Do not check for branching points while computing equilibrium curves\n");
  fprintf(stderr, "\t-noLP            : Do not check for limit points while computing equilibrium curves\n");
  fprintf(stderr, "\t-single          : Only compute the first point of the solution curve, do not continue the curve\n");
  fprintf(stderr, "\t-test            : Perform only a single integration over the life history, reporting dynamics of survival, R0,\n");
  fprintf(stderr, "\t                   i-state and interaction variables\n");
  fprintf(stderr, "\nThe value for -evoPars defaults to 0 unless set on the command-line\n");
  fprintf(stderr, "The values for -envBP, -popBP, -popEVO, -parEVO, -envZE and -popZE are undefined unless set on the command-line\n");

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
  register int  i, j, k, pind, rind;
  int           my_argc, tmpint, cvfdone = 0;
  double        minval, maxval;
  char          **argpnt1 = NULL, **argpnt2 = NULL, *my_argv[argc];
  char          *ch;

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
      else if (!strcmp(*argpnt1, "-envBP"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nIndex of environment variable for branching point continuation not specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            {
              fprintf(stderr, "\nIndex of environment variable for branching point continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, EnvironDim);
              Usage(argv[0]);
            }
          if (EnvironmentType[tmpint] != PERCAPITARATE)
            {
              fprintf(stderr, "\nDynamics of environment variable %d not specified as PERCAPITARATE!\n", tmpint);
              fprintf(stderr, "Branching point continuation not possible\n");
              Usage(argv[0]);
            }
          EnvBPIndex = tmpint;
        }
      else if (!strcmp(*argpnt1, "-popBP"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nIndex of structured population for branching point continuation not specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            {
              fprintf(stderr, "\nIndex of structured population for branching point continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, PopulationNr);
              Usage(argv[0]);
            }
          PopBPIndex = tmpint;
        }
      else if (!strcmp(*argpnt1, "-popEVO"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nIndex of structured population for ESS or PIP continuation not specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            {
              fprintf(stderr, "\nIndex of structured population for ESS or PIP continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, PopulationNr);
              Usage(argv[0]);
            }
          PopEVOIndex = tmpint;
        }
      else if (!strcmp(*argpnt1, "-parEVO"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nIndex of parameter to compute selection gradient for during EQ continuation not specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= ParameterNr))
            {
              fprintf(stderr, "\nIndex of parameter to compute selection gradient for (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint,
                      ParameterNr);
              Usage(argv[0]);
            }
          ParEVOIndex = tmpint;
        }
      else if (!strcmp(*argpnt1, "-evoPars"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo number specified for the number of life history parameters at their ESS values!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if (tmpint <= 0)
            {
              fprintf(stderr, "\nNumber of life history parameters at their ESS values should be larger than 0!\n");
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
      else if ((!strcmp(*argpnt1, "-noBP")))
        {
          BPdetection = 0;
        }
      else if ((!strcmp(*argpnt1, "-noLP")))
        {
          LPdetection = 0;
        }
      else if ((!strcmp(*argpnt1, "-test")))
        {
          TestRun = 1;
        }
      else if ((!strcmp(*argpnt1, "-single")))
        {
          DoSingle = 1;
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

  if (!strcmp(*(my_argv + 1), "BP"))
    {
      pntdim    = 2 + EnvironDim + (PopulationNr - 1);
      CurveType = BP;
      strcpy(ContinuationString, *(my_argv + 1));
      PopTrivEqui[PopBPIndex] = 0;
      EnvBPIndex              = -1;
    }
  else if (!strcmp(*(my_argv + 1), "BPE"))
    {
      pntdim    = 2 + (EnvironDim - 1) + PopulationNr;
      CurveType = BPE;
      strcpy(ContinuationString, *(my_argv + 1));
      EnvTrivEqui[EnvBPIndex] = 0;
      PopBPIndex              = -1;
    }
  else if (!strcmp(*(my_argv + 1), "EQ"))
    {
      pntdim    = 1 + EnvironDim + PopulationNr;
      CurveType = EQ;
      strcpy(ContinuationString, *(my_argv + 1));
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(*(my_argv + 1), "LP"))
    {
      pntdim    = 2 + EnvironDim + PopulationNr;
      CurveType = LP;
      strcpy(ContinuationString, *(my_argv + 1));
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(*(my_argv + 1), "ESS"))
    {
      if (evoParsDim <= 0)
        {
          fprintf(stderr, "\nNumber of life history parameters at their ESS values should be larger than 0! Define with option '-evoPars'.\n");
          Usage(argv[0]);
        }
      pntdim    = 1 + EnvironDim + PopulationNr + evoParsDim;
      CurveType = ESS;
      strcpy(ContinuationString, *(my_argv + 1));
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(*(my_argv + 1), "PIP"))
    {
      pntdim                     = 2 + EnvironDim + PopulationNr + 1;
      CurPopulationNr            = PopulationNr + 1;
      PopTrivEqui[PopulationNr]  = 1;
      CurveType                  = PIP;
      strcpy(ContinuationString, *(my_argv + 1));
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else
    Usage(my_argv[0]);

  if (((CurveType == ESS) || (CurveType == PIP)) && (PopEVOIndex == -1))
    {
      fprintf(stderr, "\nSet the structured population index for ESS or PIP continuation (0 <= i < %d) using the 'popEVO' option!\n\n",
              PopulationNr);
      Usage(argv[0]);
    }
  else if ((CurveType == BP) && (PopBPIndex == -1))
    {
      fprintf(stderr, "\nSet the structured population index for BP continuation (0 <= i < %d) using the 'popBP' option!\n\n", PopulationNr);
      Usage(argv[0]);
    }
  else if ((CurveType == BPE) && (EnvBPIndex == -1))
    {
      fprintf(stderr, "\nSet the environmental variable index for BPE continuation (0 <= i < %d) using the 'envBP' option!\n\n", EnvironDim);
      Usage(argv[0]);
    }

  // Determine the index of the environment variables and the birth rates of the structured populations in the
  // the point vector. Also determine the index of the equilibrium conditions for the environment variables and
  // the R0 values of the structured populations in the result vector
  for (i = 0, pind = 1, rind = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i])
        {
          pntdim--;
          EnvPntIndex[i] = -1;
        }
      else
        EnvPntIndex[i] = pind++;

      if (EnvTrivEqui[i] && (i != EnvBPIndex))
        EnvResIndex[i] = -1;
      else
        EnvResIndex[i] = rind++;
    }
  for (i = 0; i < CurPopulationNr; i++)
    {
      if (PopTrivEqui[i])
        {
          pntdim--;
          PopPntIndex[i] = -1;
        }
      else
        PopPntIndex[i] = pind++;

      if (PopTrivEqui[i] && (i != PopBPIndex) && (i != PopulationNr))
        R0ResIndex[i] = -1;
      else
        R0ResIndex[i] = rind++;
    }

  if ((CurveType == EQ) && (my_argc != pntdim + 6))
    Usage(argv[0]);
  else if ((CurveType == ESS) && (my_argc != pntdim + 6 + evoParsDim*3))
    Usage(argv[0]);
  else if ((CurveType != EQ) && (CurveType != ESS) && (my_argc != pntdim + 9))
    Usage(argv[0]);

  // Map all initial values of the variables into the argument vector
  memset((void *)initpnt, 0, pntdim*sizeof(double));
  for (i = 0, j = 2; i < pntdim; i++, j++) initpnt[i] = atof(my_argv[j]);

  // Map the next input value to the stepsize along the branch
  Maxcurvestep = atof(my_argv[j++]);
  curvestep    = Maxcurvestep;

  for (i = 0; j < my_argc; i++)
    {
      tmpint = (int)floor(atof(my_argv[j++]) + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        {
          fprintf(stderr, "\nIndex of parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint, ParameterNr);
          Usage(argv[0]);
        }
      minval = atof(my_argv[j++]);
      maxval = atof(my_argv[j++]);
      if (minval >= maxval)
        {
          fprintf(stderr, "\nMinimum parameter bound (%G not smaller than maximum (%G)!\n\n", minval, maxval);
          Usage(argv[0]);
        }
      if (!i)
        {
          Bifparone = tmpint;
          pntmin[0] = minval;
          pntmax[0] = maxval;
          if (CurveType == EQ) break;
        }
      else
        {
          pntmin[pntdim - 1] = minval;
          pntmax[pntdim - 1] = maxval;
          if (CurveType == ESS)
            {
              for (k = 0; k < (i - 1); k++)
                if (tmpint == evoParsIndex[k])
                  {
                    fprintf(stderr, "\nMultiple specifications for the same evolutionary parameter (%d) not allowed!\n\n", tmpint);
                    Usage(argv[0]);
                  }
              evoParsIndex[i - 1] = tmpint;
            }
          else
            {
              Bifpartwo = tmpint;
              if (CurveType == PIP)
                Bifpartwo = Bifparone;
              else if (Bifpartwo == Bifparone)
                {
                  fprintf(stderr, "\nIndex of first and second bifurcation parameter the same!\nTwo parameter continuation not possible!\n");
                  Usage(argv[0]);
                }
              break;
            }
        }
    }

  parameter[Bifparone] = initpnt[0];
  if ((CurveType == BP) || (CurveType == LP) || (CurveType == BPE))
    parameter[Bifpartwo] = initpnt[pntdim - 1];
  else if (CurveType == ESS)
    {
      for (j = 0; j < evoParsDim; j++) parameter[evoParsIndex[j]] = initpnt[pntdim - evoParsDim + j];
    }

  if ((CurveType != EQ) || (ParEVOIndex == Bifparone)) ParEVOIndex = -1;

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
  if (biffile) fclose(biffile);
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
  size_t        nrows, ncols;
  double        *curveVals, tmpdouble;
  int           tmpint, pind, rind, irhs;
  const mxArray *cell_element_ptr;
  char          optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN], varstr[MAX_STR_LEN];

  mwIndex       i, j;
  size_t        total_num_of_cells, buflen;
  int           status;

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
        "MATLAB:PSPMequi:nrhs",
        "\nIncorrect number of command-line arguments.\n\nUse: %s(%s)\n\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s",
        mexFunctionName(), "<type>, <point>, <step size>, <curve settings>, <parameters>, <options>", "<type>",
        "Type of bifurcation to perform (BP, BPE, EQ, LP, ESS or PIP)", "<point>", "Initial point of the bifurcation", "<step size>",
        "Step size along bifurcation curve", "<curve settings>", "Index, minimum and maximum value for the variable parameters (n*3 values)",
        "<parameters>", "Array of parameter values to use (empty array or of same length as parameter array)", "<options>",
        "Possible bifurcation options: envBP, popBP, popEVO, parEVO, envZE, popZE, isort, noBP, noLP, single or test");

  // check for proper number of output variables
  if (nlhs != 1) mexErrMsgIdAndTxt("MATLAB:PSPMequi:nlhs", "A single output argument is required.");

#if (MFUNCTIONS == 1)
  Minterface_Init();
#endif

  //============================== Process the options argument ======================================================================
  // Extract the contents of MATLAB cell into the C array
  irhs = 5;
  if (!mxIsCell(prhs[irhs])) mexErrMsgIdAndTxt("MATLAB:PSPMequi:options", "\nOptions should be specified as a cell array!\n");

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(optstring, "{");
  for (i = 0; i < total_num_of_cells; i++)
    {
      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optname, buflen);
      if (!((!strcmp(optname, "envBP")) || (!strcmp(optname, "popBP")) || (!strcmp(optname, "popEVO")) || (!strcmp(optname, "parEVO")) ||
            (!strcmp(optname, "envZE")) || (!strcmp(optname, "popZE")) || (!strcmp(optname, "isort")) || (!strcmp(optname, "noBP")) ||
            (!strcmp(optname, "noLP")) || (!strcmp(optname, "single")) || (!strcmp(optname, "test"))))
        mexErrMsgIdAndTxt("MATLAB:PSPMequi:options", "\nIllegal option %s!\n", optname);

      if (!strcmp(optname, "noBP") || !strcmp(optname, "noLP") || !strcmp(optname, "single") || !strcmp(optname, "test"))
        {
          if (!strcmp(optname, "noBP")) BPdetection = 0;
          if (!strcmp(optname, "noLP")) LPdetection = 0;
          if (!strcmp(optname, "single")) DoSingle  = 1;
          if (!strcmp(optname, "test")) TestRun     = 1;

          // optstring still equal to "{"
          if (strlen(optstring) == 1)
            strcat(optstring, "'");
          else
            strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "'");
          continue;
        }

      if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMequi:options", "\nNo value specified for option %s!\n", optname);

      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optval, buflen);
      if (status) mexErrMsgIdAndTxt("MATLAB:PSPMequi:options", "\nError in retrieving value for option %s!\n", optname);

      // optstring still equal to "{"
      if (strlen(optstring) == 1)
        strcat(optstring, "'");
      else
        strcat(optstring, ", '");
      strcat(optstring, optname);
      strcat(optstring, "', '");
      strcat(optstring, optval);
      strcat(optstring, "'");

      if (!strcmp(optname, "envBP"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of environment variable for branching point continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                              tmpint, EnvironDim);
          if (EnvironmentType[tmpint] != PERCAPITARATE)
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options", "\nDynamics of environment variable %d not specified as PERCAPITARATE!\n%s", tmpint,
                              "Branching point continuation not possible\n");
          EnvBPIndex = tmpint;
        }
      else if (!strcmp(optname, "popBP"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of structured population for branching point continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                              tmpint, PopulationNr);
          PopBPIndex = tmpint;
        }
      else if (!strcmp(optname, "popEVO"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of structured population for ESS or PIP continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                              tmpint, PopulationNr);
          PopEVOIndex = tmpint;
        }
      else if (!strcmp(optname, "parEVO"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= ParameterNr))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of parameter to compute selection gradient for (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint,
                              ParameterNr);
          ParEVOIndex = tmpint;
        }
      else if (!strcmp(optname, "envZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of environment variable in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n",
                              tmpint, EnvironDim);
          if (EnvironmentType[tmpint] == GENERALODE)
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nDynamics of environment variable %d not specified as PERCAPITARATE or POPULATIONINTEGRAL!\n%s", tmpint,
                              "Enforcing a boundary (zero-valued) equilibrium for it not possible\n");
          EnvTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(optname, "popZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of structured population in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n",
                              tmpint, PopulationNr);
          PopTrivEqui[tmpint] = 1;
        }
      else if (!strcmp(optname, "isort"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= IStateDim))
            mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                              "\nIndex of i-state variable for sorting structured populations (%d) not in the appropriate range (0 <= i < %d)!\n",
                              tmpint, IStateDim);
          SortIndex = tmpint;
        }
    }
  strcat(optstring, "}");

  //============================== Process the bifurcation type argument =============================================================
  // Input must be a string
  irhs = 0;
  if ((mxIsChar(prhs[irhs]) != 1) || (mxGetM(prhs[irhs]) != 1))
    mexErrMsgIdAndTxt("MATLAB:PSPMequi:inputNotString", "Bifurcation type must be a string (BP, BPE, EQ, LP, ESS or PIP).\n\n");

  // Copy the string data from prhs[irhs] into a C string input_ buf.
  strcpy(ContinuationString, mxArrayToString(prhs[irhs]));
  if (!((!strcmp(ContinuationString, "BP")) || (!strcmp(ContinuationString, "BPE")) || (!strcmp(ContinuationString, "EQ")) ||
        (!strcmp(ContinuationString, "LP")) || (!strcmp(ContinuationString, "ESS")) || (!strcmp(ContinuationString, "PIP"))))
    mexErrMsgIdAndTxt("MATLAB:PSPMequi:inputNotString", "\nBifurcation type must be a string (BP, BPE, EQ, LP, ESS or PIP).\n\n");

  if (!strcmp(ContinuationString, "BP"))
    {
      pntdim    = 2 + EnvironDim + (PopulationNr - 1);
      CurveType = BP;
      if (PopBPIndex == -1)
        mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                          "\nSet the structured population index for BP continuation (0 <= i < %d) using the 'popBP' option!\n\n", PopulationNr);
      PopTrivEqui[PopBPIndex] = 0;
      EnvBPIndex              = -1;
    }
  else if (!strcmp(ContinuationString, "BPE"))
    {
      pntdim    = 2 + (EnvironDim - 1) + PopulationNr;
      CurveType = BPE;
      if (EnvBPIndex == -1)
        mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                          "\nSet the environmental variable index for BPE continuation (0 <= i < %d) using the 'envBP' option!\n\n", EnvironDim);
      EnvTrivEqui[EnvBPIndex] = 0;
      PopBPIndex              = -1;
    }
  else if (!strcmp(ContinuationString, "EQ"))
    {
      pntdim     = 1 + EnvironDim + PopulationNr;
      CurveType  = EQ;
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(ContinuationString, "LP"))
    {
      pntdim     = 2 + EnvironDim + PopulationNr;
      CurveType  = LP;
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(ContinuationString, "ESS"))
    {
      pntdim     = 1 + EnvironDim + PopulationNr;                                 // Number of ESS parameters is added below
      CurveType  = ESS;
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(ContinuationString, "PIP"))
    {
      pntdim                     = 2 + EnvironDim + PopulationNr + 1;
      CurPopulationNr            = PopulationNr + 1;
      PopTrivEqui[PopulationNr]  = 1;
      CurveType                  = PIP;
      PopBPIndex                 = -1;
      EnvBPIndex                 = -1;
    }

  if (((CurveType == ESS) || (CurveType == PIP)) && (PopEVOIndex == -1))
    mexErrMsgIdAndTxt("MATLAB:PSPMequi:options",
                      "\nSet the structured population index for ESS or PIP continuation (0 <= i < %d) using the 'popEVO' option!\n\n",
                      PopulationNr);

  // Determine the index of the environment variables and the birth rates of the structured populations in the
  // the point vector. Also determine the index of the equilibrium conditions for the environment variables and
  // the R0 values of the structured populations in the result vector
  for (i = 0, pind = 1, rind = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i])
        {
          pntdim--;
          EnvPntIndex[i] = -1;
        }
      else
        EnvPntIndex[i] = pind++;

      if (EnvTrivEqui[i] && (i != EnvBPIndex))
        EnvResIndex[i] = -1;
      else
        EnvResIndex[i] = rind++;
    }
  for (i = 0; i < CurPopulationNr; i++)
    {
      if (PopTrivEqui[i])
        {
          pntdim--;
          PopPntIndex[i] = -1;
        }
      else
        PopPntIndex[i] = pind++;

      if (PopTrivEqui[i] && (i != PopBPIndex) && (i != PopulationNr))
        R0ResIndex[i] = -1;
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
    mexWarnMsgIdAndTxt("MATLAB:PSPMequi:parameters", "\nParameter argument ignored as it is not a row vector of length %d\n", ParameterNr);

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

  //============================== Process the step size argument ====================================================================

  irhs  = 2;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols != 1) || (nrows != 1))
    mexErrMsgIdAndTxt("MATLAB:PSPMequi:stepsize", "\nStep size along bifurcation curve must be a single double value.\n");

  memcpy(&Maxcurvestep, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  curvestep = Maxcurvestep;

  //============================== Process the curve parameters argument =============================================================

  irhs  = 3;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);

  if (CurveType == ESS)
    {
      if ((ncols < 6) || ((ncols - 6) % 3))
        mexErrMsgIdAndTxt("MATLAB:PSPMequi:bounds", "\nCurve argument must be a row vector of length (3*n) with n>=2 for ESS continuation, each "
                                                    "triplet consisting of index, minimum and maximum of the parameter.\n\n");
    }
  else if (CurveType == EQ)
    {
      if (ncols != 3) mexErrMsgIdAndTxt("MATLAB:PSPMequi:bounds", "\nCurve argument must be a row vector of length 3 for EQ continuation.\n\n");
    }
  else if (ncols != 6)
    mexErrMsgIdAndTxt("MATLAB:PSPMequi:bounds", "\nCurve argument must be a row vector of length 6 for BP, BPE, LP and PIP continuation.\n\n");

  curveVals = mxGetPr(prhs[irhs]);
  for (i = 0; i < ncols; i += 3)
    {
      tmpint = (int)floor(curveVals[i] + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        mexErrMsgIdAndTxt("MATLAB:PSPMequi:bounds", "\nIndex of parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                          ParameterNr);
      if (curveVals[i + 1] >= curveVals[i + 2])
        mexErrMsgIdAndTxt("MATLAB:PSPMequi:bounds", "\nMinimum parameter bound (%G not smaller than maximum (%G)!\n\n", curveVals[i + 1],
                          curveVals[i + 2]);
      if (!i)
        {
          Bifparone = tmpint;
          pntmin[0] = curveVals[i + 1];
          pntmax[0] = curveVals[i + 2];
        }
      else
        {
          if (CurveType == ESS)
            {
              evoParsIndex[evoParsDim] = tmpint;
              evoParsDim++;
              pntdim++;
            }
          else
            {
              Bifpartwo = tmpint;
              if (CurveType == PIP)
                Bifpartwo = Bifparone;
              else if (Bifpartwo == Bifparone)
                mexErrMsgIdAndTxt("MATLAB:PSPMequi:bounds",
                                  "\nIndex of first and second bifurcation parameter the same!\nTwo parameter continuation not possible!\n");
            }
          pntmin[pntdim - 1] = curveVals[i + 1];
          pntmax[pntdim - 1] = curveVals[i + 2];
        }
    }

  strcpy(curvestring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(curvestring, " ");
      sprintf(tmpstr, "%.6G", curveVals[i]);
      strcat(curvestring, tmpstr);
    }
  strcat(curvestring, "]");

  //============================== Process the initial point argument ================================================================

  irhs = 1;
  memset((void *)initpnt, 0, pntdim*sizeof(double));
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((nrows != 1) || (ncols != pntdim))
    {
      strcpy(varstr, "[ Par.1");
      for (i = 0; i < EnvironDim; i++)
        {
          if (CurveType == BPE)
            {
              if ((i == EnvBPIndex) || (EnvTrivEqui[i])) continue;
            }
          else if (EnvTrivEqui[i])
            continue;
          sprintf(tmpstr, " E[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          if (CurveType == BP)
            {
              if ((i == PopBPIndex) || (PopTrivEqui[i])) continue;
            }
          else if (PopTrivEqui[i])
            continue;
          sprintf(tmpstr, " b[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      if (!(CurveType == EQ)) strcat(varstr, " Par.2");
      strcat(varstr, " ]");
      mexErrMsgIdAndTxt("MATLAB:PSPMequi:point", "\nInitial point for bifurcation must be a row vector of length %d:  %s\n", pntdim, varstr);
    }
  memcpy(initpnt, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  parameter[Bifparone] = initpnt[0];

  if ((CurveType == BP) || (CurveType == LP) || (CurveType == BPE))
    parameter[Bifpartwo] = initpnt[pntdim - 1];
  else if (CurveType == ESS)
    {
      for (j = 0; j < evoParsDim; j++) parameter[evoParsIndex[j]] = initpnt[pntdim - evoParsDim + j];
    }

  strcpy(pntstring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(pntstring, " ");
      sprintf(tmpstr, "%.6G", initpnt[i]);
      strcat(pntstring, tmpstr);
    }
  strcat(pntstring, "]");

  //============================= Get the program name ===============================================================================
  // Get the name of the mex file
  strcpy(progname, mexFunctionName());
  progname[strlen(progname) - 4] = '\0';                                            // Cut off the 'equi' appendix

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
 *  R interface function
 *====================================================================================================================================
 */

#elif defined(R_PACKAGE)

SEXP PSPMequi(SEXP moduleName, SEXP bifType, SEXP initVals, SEXP stepsize, SEXP curveVals, SEXP parVals, SEXP optVals)

{
  int   i, j, ncols, tmpint, pind, rind;
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
      if (!((!strcmp(optname, "envBP")) || (!strcmp(optname, "popBP")) || (!strcmp(optname, "popEVO")) || (!strcmp(optname, "parEVO")) ||
            (!strcmp(optname, "envZE")) || (!strcmp(optname, "popZE")) || (!strcmp(optname, "isort")) || (!strcmp(optname, "noBP")) ||
            (!strcmp(optname, "noLP")) || (!strcmp(optname, "single")) || (!strcmp(optname, "test"))))
        error("\nIllegal option %s!\n\n", optname);

      if (!strcmp(optname, "noBP") || !strcmp(optname, "noLP") || !strcmp(optname, "single") || !strcmp(optname, "test"))
        {
          if (!strcmp(optname, "noBP")) BPdetection = 0;
          if (!strcmp(optname, "noLP")) LPdetection = 0;
          if (!strcmp(optname, "single")) DoSingle  = 1;
          if (!strcmp(optname, "test")) TestRun     = 1;

          if (!strlen(optstring))
            strcpy(optstring, "c(\"");
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
        strcpy(optstring, "c(\"");
      else
        strcat(optstring, ", \"");
      strcat(optstring, optname);
      strcat(optstring, "\", \"");
      strcat(optstring, optval);
      strcat(optstring, "\"");

      if (!strcmp(optname, "envBP"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            error("\nIndex of environment variable for branching point continuation (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  EnvironDim);
          if (EnvironmentType[tmpint] == GENERALODE)
            error("\nDynamics of environment variable %d not specified as PERCAPITARATE or POPULATIONINTEGRAL!\n%s", tmpint,
                  "Enforcing a boundary (zero-valued) equilibrium for it not possible\n");
          EnvBPIndex = tmpint;
        }
      else if (!strcmp(optname, "popBP"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            error("\nIndex of structured population for branching point continuation (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  PopulationNr);
          PopBPIndex = tmpint;
        }
      else if (!strcmp(optname, "popEVO"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PopulationNr))
            error("\nIndex of structured population for ESS or PIP continuation (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  PopulationNr);
          PopEVOIndex = tmpint;
        }
      else if (!strcmp(optname, "parEVO"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= ParameterNr))
            error("\nIndex of parameter to compute selection gradient for (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  ParameterNr);
          ParEVOIndex = tmpint;
        }
      else if (!strcmp(optname, "envZE"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EnvironDim))
            error("\nIndex of environment variable in boundary (trivial) equilibrium (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint,
                  EnvironDim);
          if (EnvironmentType[tmpint] == GENERALODE)
            error("\nDynamics of environment variable %d not specified as PERCAPITARATE or POPULATIONINTEGRAL!\n%s", tmpint,
                  "Enforcing a boundary (zero-valued) equilibrium for it not possible\n\n");
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
    }

  if (strlen(optstring))
    strcat(optstring, ")");
  else
    strcpy(optstring, "NULL");

  //============================== Process the bifurcation type argument =============================================================

  // Get the name of the module file
  if ((!isString(bifType)) || (length(moduleName) != 1))
    error("\nBifurcation type must be a single string\n\n");

  strcpy(ContinuationString, CHAR(STRING_ELT(bifType, 0)));
  
  if (!((!strcmp(ContinuationString, "BP")) || (!strcmp(ContinuationString, "BPE")) || (!strcmp(ContinuationString, "EQ")) ||
        (!strcmp(ContinuationString, "LP")) || (!strcmp(ContinuationString, "ESS")) || (!strcmp(ContinuationString, "PIP"))))
    error("\nBifurcation type must be a string (BP, BPE, EQ, LP, ESS or PIP).\n\n");

  if (!strcmp(ContinuationString, "BP"))
    {
      pntdim    = 2 + EnvironDim + (PopulationNr - 1);
      CurveType = BP;
      if (PopBPIndex == -1)
        error("\nSet the structured population index for BP continuation (0 <= i < %d) using the 'popBP' option!\n\n", PopulationNr);
      PopTrivEqui[PopBPIndex] = 0;
      EnvBPIndex              = -1;
    }
  else if (!strcmp(ContinuationString, "BPE"))
    {
      pntdim    = 2 + (EnvironDim - 1) + PopulationNr;
      CurveType = BPE;
      if (EnvBPIndex == -1)
        error("\nSet the environmental variable index for BPE continuation (0 <= i < %d) using the 'envBP' option!\n\n", EnvironDim);
      EnvTrivEqui[EnvBPIndex] = 0;
      PopBPIndex              = -1;
    }
  else if (!strcmp(ContinuationString, "EQ"))
    {
      pntdim     = 1 + EnvironDim + PopulationNr;
      CurveType  = EQ;
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(ContinuationString, "LP"))
    {
      pntdim     = 2 + EnvironDim + PopulationNr;
      CurveType  = LP;
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(ContinuationString, "ESS"))
    {
      pntdim     = 1 + EnvironDim + PopulationNr;                                 // Number of ESS parameters is added below
      CurveType  = ESS;
      PopBPIndex = -1;
      EnvBPIndex = -1;
    }
  else if (!strcmp(ContinuationString, "PIP"))
    {
      pntdim                     = 2 + EnvironDim + PopulationNr + 1;
      CurPopulationNr            = PopulationNr + 1;
      PopTrivEqui[PopulationNr]  = 1;
      CurveType                  = PIP;
      PopBPIndex                 = -1;
      EnvBPIndex                 = -1;
    }

  if (((CurveType == ESS) || (CurveType == PIP)) && (PopEVOIndex == -1))
    error("\nSet the structured population index for ESS or PIP continuation (0 <= i < %d) using the 'popEVO' option!\n\n", PopulationNr);

  // Determine the index of the environment variables and the birth rates of the structured populations in the
  // the point vector. Also determine the index of the equilibrium conditions for the environment variables and
  // the R0 values of the structured populations in the result vector
  for (i = 0, pind = 1, rind = 0; i < EnvironDim; i++)
    {
      if (EnvTrivEqui[i])
        {
          pntdim--;
          EnvPntIndex[i] = -1;
        }
      else
        EnvPntIndex[i] = pind++;

      if (EnvTrivEqui[i] && (i != EnvBPIndex))
        EnvResIndex[i] = -1;
      else
        EnvResIndex[i] = rind++;
    }
  for (i = 0; i < CurPopulationNr; i++)
    {
      if (PopTrivEqui[i])
        {
          pntdim--;
          PopPntIndex[i] = -1;
        }
      else
        PopPntIndex[i] = pind++;

      if (PopTrivEqui[i] && (i != PopBPIndex) && (i != PopulationNr))
        R0ResIndex[i] = -1;
      else
        R0ResIndex[i] = rind++;
    }

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

  //============================== Process the step size argument ====================================================================

  if (isReal(stepsize) && (length(stepsize) == 1))
    {
      Maxcurvestep = REAL(stepsize)[0];
      curvestep    = Maxcurvestep;
    }
  else
    error("\nStep size argument should be a single, real value\n\n");

  //============================== Process the curve parameters argument =============================================================

  ncols = length(curveVals);
  if (!isReal(curveVals) || (!length(curveVals)))
    error("\nThe curve settings argument should be a vector of real values\n\n");

  if (CurveType == ESS)
    {
      if ((ncols < 6) || ((ncols - 6) % 3))
        error("\nCurve argument must be a row vector of length (3*n) with n>=2 for ESS continuation, each triplet consisting of index, minimum and "
              "maximum of the parameter.\n\n");
    }
  else if (CurveType == EQ)
    {
      if (ncols != 3) error("\nCurve argument must be a row vector of length 3 for EQ continuation.\n\n");
    }
  else if (ncols != 6)
    error("\nCurve argument must be a row vector of length 6 for BP, BPE, LP and PIP continuation.\n\n");

  for (i = 0; i < ncols; i += 3)
    {
      tmpint = (int)floor(REAL(curveVals)[i] + MICRO);
      if ((tmpint < 0) || (tmpint >= ParameterNr))
        error("\nIndex of parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", tmpint, ParameterNr);
      if (REAL(curveVals)[i + 1] >= REAL(curveVals)[i + 2])
        error("\nMinimum parameter bound (%G not smaller than maximum (%G)!\n\n", REAL(curveVals)[i + 1], REAL(curveVals)[i + 2]);
      if (!i)
        {
          Bifparone = tmpint;
          pntmin[0] = REAL(curveVals)[i + 1];
          pntmax[0] = REAL(curveVals)[i + 2];
        }
      else
        {
          if (CurveType == ESS)
            {
              evoParsIndex[evoParsDim] = tmpint;
              evoParsDim++;
              pntdim++;
            }
          else
            {
              Bifpartwo = tmpint;
              if (CurveType == PIP)
                Bifpartwo = Bifparone;
              else if (Bifpartwo == Bifparone)
                error("\nIndex of first and second bifurcation parameter the same!\nTwo parameter continuation not possible!\n");
            }
          pntmin[pntdim - 1] = REAL(curveVals)[i + 1];
          pntmax[pntdim - 1] = REAL(curveVals)[i + 2];
        }
    }

  if (ncols)
    {
      strcpy(curvestring, "c(");
      for (i = 0; i < ncols; i++)
        {
          if (i) strcat(curvestring, ", ");
          sprintf(tmpstr, "%.6G", REAL(curveVals)[i]);
          strcat(curvestring, tmpstr);
        }
      strcat(curvestring, ")");
    }
  else
    strcpy(curvestring, "NULL");

  //============================== Process the initial point argument ================================================================

  memset((void *)initpnt, 0, pntdim*sizeof(double));
  ncols = length(initVals);
  if (!isReal(initVals) || (ncols != pntdim))
    {
      strcpy(varstr, "c( Par.1");
      for (i = 0; i < EnvironDim; i++)
        {
          if (CurveType == BPE)
            {
              if ((i == EnvBPIndex) || (EnvTrivEqui[i])) continue;
            }
          else if (EnvTrivEqui[i])
            continue;
          sprintf(tmpstr, " E[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      for (i = 0; i < PopulationNr; i++)
        {
          if (CurveType == BP)
            {
              if ((i == PopBPIndex) || (PopTrivEqui[i])) continue;
            }
          else if (PopTrivEqui[i])
            continue;
          sprintf(tmpstr, " b[%d]", (int)i);
          strcat(varstr, tmpstr);
        }
      if (!(CurveType == EQ)) strcat(varstr, " Par.2");
      strcat(varstr, " )");
      error("\nInitial point for bifurcation must be a row vector of length %d:  %s\n\n", pntdim, varstr);
    }
  memcpy(initpnt, REAL(initVals), ncols*sizeof(double));
  parameter[Bifparone] = initpnt[0];

  if ((CurveType == BP) || (CurveType == LP) || (CurveType == BPE))
    parameter[Bifpartwo] = initpnt[pntdim - 1];
  else if (CurveType == ESS)
    {
      for (j = 0; j < evoParsDim; j++) parameter[evoParsIndex[j]] = initpnt[pntdim - evoParsDim + j];
    }

  strcpy(pntstring, "c(");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(pntstring, ", ");
      sprintf(tmpstr, "%.6G", initpnt[i]);
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

  if (biffile) fclose(biffile);
  if (errfile) fclose(errfile);
  if (outfile) fclose(outfile);

  ResetCurve();
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
