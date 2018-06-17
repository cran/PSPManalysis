/***
  NAME
    lifehistory.h
  DESCRIPTION
    Header file with implementation of the routines describing the
    individual life history.
    This header file is included in the files PSPMdemo.c, PSPMequi.c
    PSPMevodyn.c and PSPMind.c.

    Copyright (C) 2015, Andre M. de Roos, University of Amsterdam

    This file is part of the PSPManalysis software package.

    This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - May 02, 2018
***/
/*
 *====================================================================================================================================
 *  Implementation of problem specification routines
 *====================================================================================================================================
 */

#if (RFUNCTIONS == 1)

#include "PSPMRinterface.h"

#elif (MFUNCTIONS == 1)

#include "Minterface.h"

#else
void SetDerivatives(const int totalOdeDim, const int *birthStateNr, const int curBirthState, int populationStart[PopulationNr],
                    int lifeStage[PopulationNr], double *bState[PopulationNr], double age, double *y, double *dyda, double *discard)
{
  int     p;
  double  *iStatePnt[PopulationNr];
  double  development[PopulationNr][IStateDim];
  double  mortality[PopulationNr];
  double  *derpnt;
#if (PULSED != 1)
  int     jj;
  double  *fecundPnt[PopulationNr];
  double  survival;
#if (PSPMDEMO == 1)
  double  cumrepro;
#else
  double  impact[PopulationNr][InteractDim];

  memset(impact, 0, PopulationNr*InteractDim*sizeof(double));
#endif
#endif

  memset(dyda, 0, totalOdeDim*sizeof(double));
  memset(development, 0, PopulationNr*IStateDim*sizeof(double));
  memset(mortality, 0, PopulationNr*sizeof(double));

  for (p = 0; p < PopulationNr; p++)
    {
      iStatePnt[p] = (populationStart[p] < 0) ? bState[p] : (y + populationStart[p]);
#if (PULSED != 1)
      fecundPnt[p] = (populationStart[p] < 0) ? discard : (dyda + populationStart[p] + (IStateDim + 1));
#endif
    }

  Development(lifeStage, iStatePnt, bState, curBirthState, Evar, development);
  Mortality(lifeStage, iStatePnt, bState, curBirthState, Evar, mortality);
#if (PULSED != 1)
  Fecundity(lifeStage, iStatePnt, bState, curBirthState, Evar, fecundPnt);
#endif
#if (PSPMDEMO != 1)
  Impact(lifeStage, iStatePnt, bState, curBirthState, Evar, impact);
#endif

  for (p = 0; p < PopulationNr; p++)
    {
      if (populationStart[p] < 0) continue;

      derpnt = (dyda + populationStart[p]);
      memcpy(derpnt, development[p], IStateDim*sizeof(double));
      derpnt += IStateDim;
#if (PSPMDEMO == 1)
      *derpnt = -(mortality[p] + PGRvar[p]);
      derpnt++;
#if (PULSED != 1)
      survival = exp(iStatePnt[p][IStateDim]);
      for (jj = 0, cumrepro = 0.0; jj < birthStateNr[p]; jj++, derpnt++)
        {
          *derpnt *= survival;                                                      // Cumulative reproduction
          cumrepro += age*(*derpnt);
        }
      *derpnt = cumrepro;                                                           // Generation time, last element
#endif
#else
      *derpnt = -(mortality[p]);
      derpnt++;
      survival = exp(iStatePnt[p][IStateDim]);
      for (jj = 0; jj < birthStateNr[p]; jj++, derpnt++) *derpnt *= survival;
      for (jj = 0; jj < InteractDim; jj++, derpnt++) *derpnt = impact[p][jj]*survival;

#if ((PSPMEQUI == 1) || (PSPMEVODYN == 1)) && (FULLSTATEOUTPUT > 0)
      if (DoStateOutput)
        {
          for (jj = 0; jj < IStateDim; jj++, derpnt++)
            *derpnt = iStatePnt[p][jj]*survival;                                    // Average i-state in cohort
          *derpnt   = survival;                                                     // Number of individuals in cohort is last
        }
#endif
#endif
    }

  return;
}
#endif // RFUNCTIONS == 0 or MFUNCTIONS == 0

/*==================================================================================================================================*/

double SetMaxLimit(const int *birthStateNr, const int curBirthState, int populationStart[PopulationNr], int lifeStage[PopulationNr],
                   double *bState[PopulationNr], double cohortBound[PopulationNr], double *y)
{
  int     p;
  double  maxlimit = -1.0;
  double  *iStatePnt[PopulationNr];
  double  Limits[PopulationNr];

  for (p = 0; p < PopulationNr; p++)
    {
      iStatePnt[p] = (populationStart[p] < 0) ? bState[p] : y + populationStart[p];
      Limits[p]    = -1.0;
    }

  IntervalLimit(lifeStage, iStatePnt, bState, curBirthState, Evar, Limits);

  for (p = 0; p < PopulationNr; p++)
    {
      if (populationStart[p] < 0) continue;
      Limits[p] = max(Limits[p], LogMinSurvival - iStatePnt[p][IStateDim]);
      maxlimit  = max(maxlimit, Limits[p]);
#if (FULLSTATEOUTPUT > 0)
      if (DoStateOutput) maxlimit = max(maxlimit, (iStatePnt[p][SortIndex] - cohortBound[p]));
#endif
    }

  return maxlimit;
}


/*==================================================================================================================================*/

int SwitchLifeStage(const int *birthStateNr, const int curBirthState, int populationStart[PopulationNr], int lifeStage[PopulationNr],
                    double *bState[PopulationNr], double *y, double cohortBound[PopulationNr], double cohortLims[PopulationNr], const int initial)
{
  int     i, j, p, s;
  int     newstage[PopulationNr];
  double  newiStates[PopulationNr][IStateDim + 1];
  double *iStatePnt[PopulationNr], *newiStatePnt[PopulationNr];
  double  Limits[PopulationNr];

  for (p = 0; p < PopulationNr; p++)
    {
      iStatePnt[p]    = (populationStart[p] < 0) ? bState[p] : y + populationStart[p];
      newiStatePnt[p] = newiStates[p];
      memcpy(newiStatePnt[p], iStatePnt[p], (IStateDim + 1)*sizeof(double));
      newstage[p] = -1;
    }

  for (p = 0; p < PopulationNr; p++)
    {
      if (populationStart[p] < 0) continue;

      if ((LogMinSurvival - iStatePnt[p][IStateDim]) > -DYTOL)                    // Individual is dead
        {
          if (TestRun)
            {
              STDOUT("\nPop. #%2d - Bstate %2d - Stage %2d   (End):", p, curBirthState, lifeStage[p]);
              for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", iStatePnt[p][i]);
              STDOUT("%15.6G", exp(iStatePnt[p][IStateDim]));
              STDOUT("%15.6G", SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
#if (PSPMDEMO == 1)
              if (fabs(SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1)) > 0)
                STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p]]/SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
              else
                STDOUT("%15s", "Undefined");
#else
              for (i = 0; i < InteractDim; i++) STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p] + i]);
#endif
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
            }
          lifeStage[p]       = Stages;
          populationStart[p] = -1;
          continue;
        }

      for (s = lifeStage[p]; s < Stages; s++)
        {
          for (j = 0; j < PopulationNr; j++) Limits[j] = -1.0;
          IntervalLimit(lifeStage, newiStatePnt, bState, curBirthState, Evar, Limits);
          Limits[p] = max(Limits[p], LogMinSurvival - newiStates[p][IStateDim]);

          if (initial && (Limits[p] < 0))
            break;
          else if (!initial && (fabs(Limits[p]) > DYTOL))
            break;

          if (TestRun)
            {
              STDOUT("\nPop. #%2d - Bstate %2d - Stage %2d   (End):", p, curBirthState, lifeStage[p]);
              for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", newiStates[p][i]);
              STDOUT("%15.6G", exp(newiStates[p][IStateDim]));
              // Unchanged values
              STDOUT("%15.6G", SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
#if (PSPMDEMO == 1)
              if (fabs(SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1)) > 0)
                STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p]]/SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
              else
                STDOUT("%15s", "Undefined");
#else
              for (i = 0; i < InteractDim; i++) STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p] + i]);
#endif
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
            }

          // Carry out the discrete changes for population p that are required on transitioning into the new life history stage
          lifeStage[p]++;
          newstage[p] = lifeStage[p];
          DiscreteChanges(newstage, newiStatePnt, bState, curBirthState, Evar);
          newstage[p] = -1;

          if (TestRun)
            {
              for (i = 0; i < (IStateDim + 1); i++)
                if (fabs(iStatePnt[p][i] - newiStates[p][i]) > DYTOL) break;
              if (i < (IStateDim + 1))
                {
                  STDOUT("\nPop. #%2d - Bstate %2d - Stage %2d (Start):", p, curBirthState, lifeStage[p]);
                  for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", newiStates[p][i]);
                  STDOUT("%15.6G", exp(newiStates[p][IStateDim]));
                  // Unchanged values
                  STDOUT("%15.6G", SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
#if (PSPMDEMO == 1)
                  if (fabs(SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1)) > 0)
                    STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p]]/SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
                  else
                    STDOUT("%15s", "Undefined");
#else
                  for (i = 0; i < InteractDim; i++) STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p] + i]);
#endif
#if (defined(R_PACKAGE))
                  R_FlushConsole();
                  R_ProcessEvents();
#endif
                }
            }
          memcpy(iStatePnt[p], newiStatePnt[p], (IStateDim + 1)*sizeof(double));

#if (PSPMIND == 1)
          if (DoStateOutput && !initial)
            {
              for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = iStatePnt[p][i];

              PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(iStatePnt[p][IStateDim]);

              // Store the cumulative contribution to the impact
              for (i = 0; i < InteractDim; i++)
                PopDens(p, curBirthState, (IStateDim + 1) + i, Cohorts(p, curBirthState)) = iStatePnt[p][IStateDim + 1 + birthStateNr[p] + i];

              // Store the cumulative number of offspring in all birth states
              for (i = 0; i < birthStateNr[p]; i++)
                PopDens(p, curBirthState, CohortDim + i, Cohorts(p, curBirthState)) = iStatePnt[p][IStateDim + 1 + i];
              Cohorts(p, curBirthState)++;
              cohortBound[p] += cohortLims[p];
            }
#endif
        }
      if (lifeStage[p] >= Stages)
        populationStart[p] = -1;
      else if (Limits[p] > 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Integration failed to stop at end of life stage %d, starting from birth state %d in population %d",
                   lifeStage[p], curBirthState, p);

          if (TestRun)
            {
              for (j = 0; j < PopulationNr; j++)
                {
                  if (populationStart[j] < 0) continue;
                  STDOUT("\nPop. #%2d - Bstate %2d - Stage %2d   (End):", j, curBirthState, lifeStage[j]);
                  for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", iStatePnt[j][i]);
                  STDOUT("%15.6G", exp(iStatePnt[j][IStateDim]));
                  STDOUT("%15.6G", SUM(birthStateNr[j], iStatePnt[j] + (IStateDim + 1), 1));
#if (PSPMDEMO == 1)
                  if (fabs(SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1)) > 0)
                    STDOUT("%15.6G", iStatePnt[p][(IStateDim + 1) + birthStateNr[p]]/SUM(birthStateNr[p], iStatePnt[p] + (IStateDim + 1), 1));
                  else
                    STDOUT("%15s", "Undefined");
#else
                  for (i = 0; i < InteractDim; i++) STDOUT("%15.6G", iStatePnt[j][(IStateDim + 1) + birthStateNr[j] + i]);
#endif
#if (defined(R_PACKAGE))
                  R_FlushConsole();
                  R_ProcessEvents();
#endif
                }
            }
          return FAILURE;
        }
    }

  return SUCCES;
}


/*==================================================================================================================================*/
