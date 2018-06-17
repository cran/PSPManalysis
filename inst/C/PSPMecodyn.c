/*
  NAME
     PSPMecodyn.c

     This file provides the backbone program file for the numerical integration with the Escalator Boxcar
     Train software of a physiologically structured population model that is specificied in the type of
     header file used for the PSPManalysis package.

  Last modification: AMdR - May 05, 2018
*/

#define PSPMECODYN                1                                                 // File identification
#define EBTLIB					                                                            // and file grouping

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

#include "escbox.h"
#include "ebtcohrt.h"
#include "ebtmain.h"
#include "ebtutils.h"

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)
#include "ebtinit.h"
#endif


/*================================ PREPARE IMPORTING THE MODEL IMPLEMENTATION ======================================================*/

#define PERCAPITARATE             2001
#define GENERALODE                2002
#define POPULATIONINTEGRAL        2003

/*================================ IMPORT THE MODEL DIMENSIONS AND IMPLEMENTATION ==================================================*/

#define Survival(p)               (exp(-(istate[p][-1] - 1.0)))  
#define SetSurvival(p, s)         istate[p][-1] = (((s) >= 0) && ((s) < exp(-(istate[p][-1] - 1.0)))) ? (1.0 - log(max((s), DBL_EPSILON))) : (istate[p][-1])

#if ((RFUNCTIONS != 1) && (MFUNCTIONS != 1))
#undef POPULATION_NR
#undef STAGES
#undef I_STATE_DIM
#undef ENVIRON_DIM
#undef INTERACT_DIM
#undef PARAMETER_NR

#if defined(PROBLEMHEADER)                                                          // Include header file
#define HEADERNAME <PROBLEMHEADER>
#include HEADERNAME
#else
#error No header file defined!
#endif
#endif

#ifndef EBTMETHOD
#define EBTMETHOD                 1
#endif


/*=================================== END IMPORTED DIMENSION SETTINGS  =============================================================*/

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
  I-constant variables include:
    Index:
    0                   : Initial cohort density
    1
    :                   : All i-state variable values at birth
    I_STATE_DIM
    I_STATE_DIM + 1     : Current life stage
    I_STATE_DIM + 2     : Number of the state at birth
    I_STATE_DIM + 3     : Time that the cohort was initiated
*/
#ifndef I_CONST_DIM
#define I_CONST_DIM	              (I_STATE_DIM + 4)
#endif

/* 
  Output variables include:
    Index:
    0
    :                                                                               : All environmental variables
    ENVIRON_DIM_PSPM - 1
    ENVIRON_DIM_PSPM
    :                                                                               : Total birth rates of all populations
    ENVIRON_DIM_PSPM + POPULATION_NR - 1
    ENVIRON_DIM_PSPM + POPULATION_NR
    :                                                                               : Interaction variables of all populations
    ENVIRON_DIM_PSPM + POPULATION_NR + POPULATION_NR*INTERACT_DIM - 1
*/
#undef  OUTPUT_VAR_NR
#define OUTPUT_VAR_NR             (ENVIRON_DIM + POPULATION_NR + POPULATION_NR*INTERACT_DIM)

// For each population the end of each life stage
#ifndef EVENT_NR
#define EVENT_NR                  (POPULATION_NR * STAGES)
#endif

/*==================================================================================================================================*/

#define popIDlifestage(m, n)      popIDcard[m][n][I_STATE_DIM + 1]
#define popIDbstatenr(m, n)       popIDcard[m][n][I_STATE_DIM + 2]
#define popIDbirthtime(m, n)      popIDcard[m][n][I_STATE_DIM + 3]

#define ofsIDlifestage(m, n)      ofsIDcard[m][n][I_STATE_DIM + 1]
#define ofsIDbstatenr(m, n)       ofsIDcard[m][n][I_STATE_DIM + 2]
#define ofsIDbirthtime(m, n)      ofsIDcard[m][n][I_STATE_DIM + 3]

#define getPopIDlifestage(m, n)   (lrint(popIDlifestage(m, n)))
#define getPopIDbstatenr(m, n)    (lrint(popIDbstatenr(m, n)))

#define getOfsIDlifestage(m, n)   (lrint(ofsIDlifestage(m, n)))
#define getOfsIDbstatenr(m, n)    (lrint(ofsIDbstatenr(m, n)))

#define isdead(m, n)              ((pop[m][n][number] > logMinSurvival) || (popIDcard[m][n][I_STATE_DIM + 3] < (env[0] - MAX_AGE)))

static int                        MaxCohortNr = -1;
static int                        MaxStatesAtBirth = -1;
static double                     curPSPMEnv[ENVIRON_DIM];
static double                     logMinSurvival;
int                               EBTMethod;
static double                     Time = 0;

static void setCurrentEnvironmentValues(double *env, population *pop, population *ofs, double popImpacts[POPULATION_NR][INTERACT_DIM]);


#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
static char                       parstring[MAXFILENAMELEN];
static char                       optstring[MAXFILENAMELEN];
static char                       timestring[MAXFILENAMELEN];
static char                       bifstring[MAXFILENAMELEN];
#endif

#if (RFUNCTIONS == 1)

#include "PSPMRinterface.h"

#elif (MFUNCTIONS == 1)

#include "Minterface.h"

#endif

/*
 *====================================================================================================================================
 *
 * USER INITIALIZATION ROUTINE ALLOWS OPERATIONS ON INITIAL POPULATIONS
 *
 *====================================================================================================================================
 */

void	UserInit( int argc, char **argv, double *env,  population *pop)
  
{
  int     ii, j, jj, pp;
  int     lifeStage[POPULATION_NR], bstatenr;
  double *istatePnt[POPULATION_NR];
  double *bstatePnt[POPULATION_NR];
  double  limitVals[POPULATION_NR];
  int     bstateNr[POPULATION_NR];

  Time = env[0];
  logMinSurvival = 1.0 - log(MIN_SURVIVAL);

  // Set the maximum number of cohorts
  for (pp = 0, MaxCohortNr = 0; pp < POPULATION_NR; pp++)
    MaxCohortNr = max(MaxCohortNr, cohort_no[pp]);
  
  memcpy(curPSPMEnv, env + 1, ENVIRON_DIM*sizeof(double));

  // Dummy call to set the names attributes to the i-state variable and the value of MaxStatesAtBirth
  memset(bstateNr, 0, POPULATION_NR*sizeof(int));
  SetBirthStates(bstateNr, curPSPMEnv);
  for (pp = 0; pp < POPULATION_NR; pp++)
    MaxStatesAtBirth = max(MaxStatesAtBirth, bstateNr[pp]);

  // Set the number entry of the cohort data to 1 - log(Survival) if popIDcard[pp][jj][number] is missing
  for (pp = 0; pp < POPULATION_NR; pp++)
    {
      for (jj = 0; jj < cohort_no[pp]; jj++)
        {
          if (ismissing(popIDcard[pp][jj][number]))
            {
              for (ii = cohort_no[pp] - 1; ii > 0; ii--)
                if (getPopIDbstatenr(pp, jj) == getPopIDbstatenr(pp, ii)) break;
              popIDcard[pp][jj][number] = pop[pp][ii][number];
              pop[pp][jj][number]       = 1 - log(pop[pp][jj][number]/popIDcard[pp][jj][number]);    
            }
        }
    }

  // Set the appropriate stage index
  for (j = 0; j < MaxCohortNr; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;

          istatePnt[pp] = pop[pp][jj] + 1;                                          // Index 0 is cohort density
          bstatePnt[pp] = popIDcard[pp][jj] + 1;

          popIDlifestage(pp, jj) = 0;
          lifeStage[pp]          = 0;

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getPopIDbstatenr(pp, jj);

          limitVals[pp] = -1.0;
        }

      // Call the life history routine
      IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);

      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if ((j >= cohort_no[pp]) || (isdead(pp, jj))) continue;

          while ((limitVals[pp] > 0) || iszero(limitVals[pp]))
            {
              popIDlifestage(pp, jj) += 1;
              lifeStage[pp] = getPopIDlifestage(pp, jj);
              for (ii = 0; ii < POPULATION_NR; ii++) limitVals[ii] = -1.0;
              IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);
            }
        }
    }

  return;
}



/*
 *====================================================================================================================================
 *
 *	SPECIFICATION OF THE NUMBER AND VALUES OF BOUNDARY POINTS
 *
 *====================================================================================================================================
 */

void SetBpointNo(double *env, population *pop, int *bpoint_no)

{
  int pp;

  Time = env[0];
  for (pp = 0, MaxCohortNr = 0; pp < POPULATION_NR; pp++)
    MaxCohortNr = max(MaxCohortNr, cohort_no[pp]);

  // Set the current values of the environmental variables
  setCurrentEnvironmentValues(env, pop, NULL, NULL);

  SetBirthStates(bpoint_no, curPSPMEnv);

  for (pp = 0, MaxStatesAtBirth = 0; pp < POPULATION_NR; pp++)
    MaxStatesAtBirth = max(MaxStatesAtBirth, bpoint_no[pp]);

#if (POPULATION_NR > 1)
  if (MaxStatesAtBirth > 1)
    ErrorExit(0, "Current version of PSPMecodyn() does not allow for multiple populations with multiple states at birth!");
#endif

  return;
}


/*==================================================================================================================================*/

void	SetBpoints(double *env, population *pop, population *bpoints)

{
  int     bb, pp;
  double *bstatePnt[POPULATION_NR];
  double  dummy[I_STATE_DIM];

  Time = env[0];
  // No need to update the POPULATIONINTEGRAL variables

  for (bb = 0; bb < MaxStatesAtBirth; bb++)
    {
      for (pp = 0; pp < POPULATION_NR; pp++)
        bstatePnt[pp] = (bb < bpoint_no[pp]) ? bpoints[pp][bb] + 1 : dummy;

      StateAtBirth(bstatePnt, bb, curPSPMEnv);

      for (pp = 0; pp < POPULATION_NR; pp++)
        if (bb < bpoint_no[pp])
          {
            memcpy(ofsIDcard[pp][bb] + 1, bpoints[pp][bb] + 1, I_STATE_DIM*sizeof(double));
#if (EBTMETHOD == 0)
            memcpy(ofs[pp][bb] + 1, bpoints[pp][bb] + 1, I_STATE_DIM*sizeof(double));
#endif
            ofsIDbstatenr(pp, bb)  = bb;
            ofsIDbirthtime(pp, bb) = env[0];
          }
    }
  
  return;
}


/*
 *====================================================================================================================================
 *
 *			SPECIFICATION OF DERIVATIVES
 *
 *====================================================================================================================================
 */

void Gradient(double *env, population *pop, population *ofs, double *envgrad, population *popgrad, population *ofsgrad, population *bpoints)

{
  int     bb, ii, j, jj, pp;
  double *istatePnt[POPULATION_NR];
  double *bstatePnt[POPULATION_NR];
  int     lifeStage[POPULATION_NR], bstatenr;
  double  istateVals[POPULATION_NR][I_STATE_DIM];
  double  development[POPULATION_NR][I_STATE_DIM];
  double  mortality[POPULATION_NR];
#if (EBTMETHOD != 0)
  double  istateValsdI[POPULATION_NR][I_STATE_DIM], istateDiff;
  double  developmentdI[POPULATION_NR][I_STATE_DIM];
  double  mortalitydI[POPULATION_NR];
#endif
  double *fecundity[POPULATION_NR], *fecundityMem;
  double *curBirthRate[POPULATION_NR], *curBirthRateMem;
  double  curPopIntegrals[POPULATION_NR][INTERACT_DIM];
  double  curEnvCondition[ENVIRON_DIM];
  double  cohortDensity;

  Time = env[0];
  fecundityMem = malloc(POPULATION_NR*MaxStatesAtBirth*sizeof(double));
  if (!fecundityMem) ErrorExit(0, "Memory allocation failure in Gradient()!");

  curBirthRateMem = malloc(POPULATION_NR*MaxStatesAtBirth*sizeof(double));
  if (!curBirthRateMem) ErrorExit(0, "Memory allocation failure in Gradient()!");
  memset(curBirthRateMem, 0, POPULATION_NR*MaxStatesAtBirth*sizeof(double));

  for (pp = 0; pp < POPULATION_NR; pp++)
    {
      fecundity[pp]    = fecundityMem + pp*MaxStatesAtBirth;
      curBirthRate[pp] = curBirthRateMem + pp*MaxStatesAtBirth;
    }

  // Set the current values of the environmental variables
  setCurrentEnvironmentValues(env, pop, ofs, curPopIntegrals);

  for (j = 0; j < MaxCohortNr; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;

          istatePnt[pp] = pop[pp][jj] + 1;                                                // Index 0 is cohort density
          bstatePnt[pp] = popIDcard[pp][jj] + 1;
          lifeStage[pp] = getPopIDlifestage(pp, jj);

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getPopIDbstatenr(pp, jj);
        }

      // Call the life history routines
      memset(development, 0, POPULATION_NR *I_STATE_DIM*sizeof(double));
      memset(fecundityMem, 0, POPULATION_NR*MaxStatesAtBirth*sizeof(double));
      memset(mortality, 0, POPULATION_NR*sizeof(double));

#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
      SetLifeHistoryRates(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv,
                          development, fecundity, mortality, NULL);   
#else
      Development(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, development);      
      Fecundity(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, fecundity);
      Mortality(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, mortality);
#endif

      // Assign the derivatives
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if (j < cohort_no[pp])
            {
              popgrad[pp][jj][number] = mortality[pp];
              memcpy(popgrad[pp][jj] + 1, development[pp], I_STATE_DIM*sizeof(double));

              if  (isdead(pp, jj)) continue;

              cohortDensity = popIDcard[pp][jj][number]*exp(-(pop[pp][jj][number] - 1.0));
              for (bb = 0; bb < bpoint_no[pp]; bb++)
                curBirthRate[pp][bb] += fecundity[pp][bb]*cohortDensity;
            }
        }
    }

  // Process the offspring cohorts
  for (j = 0; j < MaxStatesAtBirth; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < bpoint_no[pp]) ? j : bpoint_no[pp] - 1;

#if (EBTMETHOD != 0)
          memcpy(istateVals[pp], ofsIDcard[pp][jj] + 1, I_STATE_DIM*sizeof(double));
#else
          memcpy(istateVals[pp], ofs[pp][jj] + 1, I_STATE_DIM*sizeof(double));
#endif
          istatePnt[pp] = istateVals[pp];
          bstatePnt[pp] = ofsIDcard[pp][jj] + 1;
          lifeStage[pp] = getOfsIDlifestage(pp, jj);

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getOfsIDbstatenr(pp, jj);
        }

      // Call the life history routines
      memset(development, 0, POPULATION_NR *I_STATE_DIM*sizeof(double));
      memset(mortality, 0, POPULATION_NR*sizeof(double));

#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
      SetLifeHistoryRates(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv,
                          development, NULL, mortality, NULL);   
#else
      Development(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, development);
      Mortality(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, mortality);
#endif

      // Assign the derivatives
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if (j < bpoint_no[pp])
            {
              ofsgrad[pp][jj][number]  = -mortality[pp]*ofs[pp][jj][number];
              ofsgrad[pp][jj][number] += curBirthRate[pp][jj];
              for (ii = 0; ii < I_STATE_DIM; ii++)
                {
                  ofsgrad[pp][jj][i_state(ii)] = development[pp][ii];
#if (EBTMETHOD != 0)
                  ofsgrad[pp][jj][i_state(ii)] *= ofs[pp][jj][number];
                  ofsgrad[pp][jj][i_state(ii)] -= mortality[pp]*ofs[pp][jj][i_state(ii)];
#endif
                }
            }
        }

#if (EBTMETHOD != 0)
      // Numerically compute the forward derivative of the development and mortality function
      for (ii = 0; ii < I_STATE_DIM; ii++)
        {
          memcpy(istateValsdI,  istateVals,  POPULATION_NR*I_STATE_DIM*sizeof(double));
          for (pp = 0; pp < POPULATION_NR; pp++)
            {
              if (j < bpoint_no[pp])
                istateValsdI[pp][ii] += max(JACOBIAN_STEP*istateVals[pp][ii], JACOBIAN_MIN_STEP);

              // Reset the pointer for the derivative computation
              istatePnt[pp] = istateValsdI[pp];
            }
            
          // Call the life history routines
          memset(developmentdI, 0, POPULATION_NR *I_STATE_DIM*sizeof(double));
          memset(mortalitydI, 0, POPULATION_NR*sizeof(double));

#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
          SetLifeHistoryRates(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv,
                              developmentdI, NULL, mortalitydI, NULL);   
#else
          Development(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, developmentdI);
          Mortality(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, mortalitydI);
#endif

          // Add the partial derivative terms to the derivatives
          for (pp = 0; pp < POPULATION_NR; pp++)
            {
              if (j < bpoint_no[pp])
                {
                  istateDiff = istateValsdI[pp][ii] - istateVals[pp][ii];
                  for (bb = 0; bb < I_STATE_DIM; bb++)
                    {
                      developmentdI[pp][bb] -= development[pp][bb];
                      developmentdI[pp][bb] /= istateDiff;
                    }
                  mortalitydI[pp] -= mortality[pp];
                  mortalitydI[pp] /= istateDiff;

                  ofsgrad[pp][jj][number] += -mortalitydI[pp]*ofs[pp][jj][i_state(ii)];
                  for (bb = 0; bb < I_STATE_DIM; bb++)
                    ofsgrad[pp][jj][bb + 1] += developmentdI[pp][bb]*ofs[pp][jj][i_state(ii)];
                }
            }
        }
#endif
    }

  memset(curEnvCondition, 0, ENVIRON_DIM*sizeof(double));
  EnvEqui(curPSPMEnv, curPopIntegrals, curEnvCondition);

  memset(envgrad, 0, ENVIRON_DIM_EBT*sizeof(double));

  envgrad[0] = 1;
  for (ii = 0; ii < ENVIRON_DIM; ii++)
    {
      switch (EnvironmentType[ii])
        {
          case GENERALODE:
            envgrad[ii + 1] = curEnvCondition[ii];
            break;
          case PERCAPITARATE:
            envgrad[ii + 1] = curEnvCondition[ii] * env[ii + 1];
            break;
          default:
            break;
        }
    }

  free(fecundityMem);
  free(curBirthRateMem);

  return;
}


/*
 *====================================================================================================================================
 *
 *	SPECIFICATION OF EVENT LOCATION AND DYNAMIC COHORT CLOSURE
 *
 *====================================================================================================================================
 */

void	EventLocation(double *env, population *pop, population *ofs, population *bpoints, double *events)

{
  int     j, jj, pp;
  int     lifeStage[POPULATION_NR], bstatenr;
  double *istatePnt[POPULATION_NR];
  double *bstatePnt[POPULATION_NR];
  double  limitVals[POPULATION_NR];

  Time = env[0];
  // Set the current values of the environmental variables
  setCurrentEnvironmentValues(env, pop, ofs, NULL);

  for (jj = 0; jj < EVENT_NR; jj++) events[jj] = -1;

  for (j = 0; j < MaxCohortNr; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;

          istatePnt[pp] = pop[pp][jj] + 1;                                                // Index 0 is cohort density
          bstatePnt[pp] = popIDcard[pp][jj] + 1;
          lifeStage[pp] = getPopIDlifestage(pp, jj);

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getPopIDbstatenr(pp, jj);
          limitVals[pp] = -1.0;
        }

      // Call the life history routine
      IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);

      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if ((j < cohort_no[pp]) && (!isdead(pp, jj)))
            events[pp*STAGES + lifeStage[pp]] = max(limitVals[pp], events[pp*STAGES + lifeStage[pp]]);
        }
    }

  return;
}


/*==================================================================================================================================*/

int	ForceCohortEnd(double *env, population *pop, population *ofs, population *bpoints)
  
{
  int     ii, j, jj, pp;
  int     lifeStage[POPULATION_NR], bstatenr;
  int     newStage[POPULATION_NR];
  double *istatePnt[POPULATION_NR];
  double *bstatePnt[POPULATION_NR];
  double  limitVals[POPULATION_NR];

  Time = env[0];
  // No need to set the current values of the environmental variables, have been set in EventLocation()

  for (pp = 0; pp < POPULATION_NR; pp++) newStage[pp] = -1;

  for (j = 0; j < MaxCohortNr; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;

          istatePnt[pp] = pop[pp][jj] + 1;                                                // Index 0 is cohort density
          bstatePnt[pp] = popIDcard[pp][jj] + 1;
          lifeStage[pp] = getPopIDlifestage(pp, jj);

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getPopIDbstatenr(pp, jj);

          limitVals[pp] = -1.0;
        }

      // Call the life history routine
      IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);

      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if ((j >= cohort_no[pp]) || (isdead(pp, jj))) continue;

          while  ((limitVals[pp] > 0) || iszero(limitVals[pp]))
            {
              popIDlifestage(pp, jj) += 1;

              newStage[pp] = getPopIDlifestage(pp, jj);
              if (newStage[pp] >= STAGES) break;

              DiscreteChanges(newStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv);
              newStage[pp] = -1;

              lifeStage[pp] = getPopIDlifestage(pp, jj);
              for (ii = 0; ii < POPULATION_NR; ii++) limitVals[ii] = -1.0;
              IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);
            }
        }
    }

  return NO_COHORT_END;
}


/*
 *====================================================================================================================================
 *
 *		SPECIFICATION OF BETWEEN COHORT CYCLE DYNAMICS
 *
 *====================================================================================================================================
 */

void	InstantDynamics(double *env, population *pop, population *ofs)
  
{
  int     ii, j, jj, pp;
  int     lifeStage[POPULATION_NR], bstatenr;
  int     newStage[POPULATION_NR];
  double *istatePnt[POPULATION_NR];
  double *bstatePnt[POPULATION_NR];
  double  limitVals[POPULATION_NR];

  Time = env[0];
  for (pp = 0; pp < POPULATION_NR; pp++) newStage[pp] = -1;
  for (j = 0; j < MaxCohortNr; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;

          istatePnt[pp] = pop[pp][jj] + 1;                                                // Index 0 is cohort density
          bstatePnt[pp] = popIDcard[pp][jj] + 1;
          lifeStage[pp] = getPopIDlifestage(pp, jj);

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getPopIDbstatenr(pp, jj);
          limitVals[pp] = -1.0;
        }

      // Call the life history routine
      IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);

      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if ((j >= cohort_no[pp]) || (isdead(pp, jj))) continue;

          while ((limitVals[pp] > 0) || iszero(limitVals[pp]))
            {
              popIDlifestage(pp, jj) += 1;

              newStage[pp] = getPopIDlifestage(pp, jj);
              if (newStage[pp] >= STAGES) break;

              DiscreteChanges(newStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv);
              newStage[pp] = -1;

              lifeStage[pp] = getPopIDlifestage(pp, jj);
              for (ii = 0; ii < POPULATION_NR; ii++) limitVals[ii] = -1.0;
              IntervalLimit(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, limitVals);
            }
        }
    }

  for (pp = 0; pp < POPULATION_NR; pp++)
    {
      for (jj = 0; jj < cohort_no[pp]; jj++)
        {
          if (isdead(pp, jj)) pop[pp][jj][number] = 0;
        }

      for (jj = 0; jj < bpoint_no[pp]; jj++)
        {
          ofsIDcard[pp][jj][number] = ofs[pp][jj][number];
          ofs[pp][jj][number]       = 1.0;
#if (EBTMETHOD == 0)
          // In addition to the method proposed by Brännström et al., who suggest to use dx_b/dt=g(E,x_b) as ODE for the boundary
          // cohort, here we at the end of the cohort interval the value of x_b(t) is averaged with the value x_b(0), the 
          // initial value of the state at birth.
          for (ii = 0; ii < I_STATE_DIM; ii++)
            {
              ofs[pp][jj][1 + ii] += ofsIDcard[pp][jj][1 + ii];
              ofs[pp][jj][1 + ii] /= 2.0;
            }
#endif
        }
    }


  return;
}


/*
 *====================================================================================================================================
 *
 *			SPECIFICATION OF OUTPUT VARIABLES
 *
 *====================================================================================================================================
 */

void	DefineOutput(double *env, population *pop, double *output)
  
{
  int            bb, j, jj, pp, outnr;
  int            lifeStage[POPULATION_NR], bstatenr;
  double        *istatePnt[POPULATION_NR];
  double        *bstatePnt[POPULATION_NR];
  double        *fecundity[POPULATION_NR], *fecundityMem;
  double         curPopIntegrals[POPULATION_NR][INTERACT_DIM];
  double         cohortDensity;
  static double  nextreport = 0.0;

  Time = env[0];
  // Compute MaxCohortNr anew after SievePop()
  for (pp = 0, MaxCohortNr = 0; pp < POPULATION_NR; pp++)
    MaxCohortNr = max(MaxCohortNr, cohort_no[pp]);

  memset(output, 0, OUTPUT_VAR_NR*sizeof(double));

  // Set the current values of the environmental variables
  setCurrentEnvironmentValues(env, pop, NULL, curPopIntegrals);

  // Output all environmental variables
  for (jj = 0, outnr = 0; jj < ENVIRON_DIM; jj++) output[outnr++] = curPSPMEnv[jj];

  fecundityMem = malloc(POPULATION_NR*MaxStatesAtBirth*sizeof(double));
  if (!fecundityMem) ErrorExit(0, "Memory allocation failure in DefineOutput()!");

  for (pp = 0; pp < POPULATION_NR; pp++)
    fecundity[pp] = fecundityMem + pp*MaxStatesAtBirth;

  // Compute and output the total population birth rates
  for (j = 0; j < MaxCohortNr; j++)
    {
      // Set up the pointers
      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;

          istatePnt[pp] = pop[pp][jj] + 1;                                                // Index 0 is cohort density
          bstatePnt[pp] = popIDcard[pp][jj] + 1;
          lifeStage[pp] = getPopIDlifestage(pp, jj);

          // For multiple populations and multiple birth states the next statement is incorrect
          bstatenr      = getPopIDbstatenr(pp, jj);
        }

      memset(fecundityMem, 0, POPULATION_NR*MaxStatesAtBirth*sizeof(double));
        
      // Call the life history routines
#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
      SetLifeHistoryRates(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv,
                          NULL, fecundity, NULL, NULL);   
#else
      Fecundity(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, fecundity);
#endif

      for (pp = 0; pp < POPULATION_NR; pp++)
        {
          if (j < cohort_no[pp])
            {
              cohortDensity = popIDcard[pp][jj][number]*exp(-(pop[pp][jj][number] - 1.0));              
              for (bb = 0; bb < bpoint_no[pp]; bb++)
                output[outnr + pp] += fecundity[pp][bb]*cohortDensity;
            }
        }
    }
  free(fecundityMem);

  if (isequal(env[0], nextreport) || (env[0] >= nextreport))
    {
      // Output environmental variables and population birth rates to stdout
      if (env[0] < 1.0E7)
        STDOUT("%12.2f", env[0]);
      else
        STDOUT("%12.6f", env[0]);

      for (j = 0; j < outnr; j++)
        {
#if (defined(R_PACKAGE))
          STDOUT(",%15.8E", output[j]);
#else
          STDOUT("%16.8E",  output[j]);
#endif
        }
      STDOUT("\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
      nextreport = env[0] + report_level;
    }

  // Output the values of all population level impact variables
  for (pp = 0, outnr += POPULATION_NR; pp < POPULATION_NR; pp++)
    for (jj = 0; jj < INTERACT_DIM; jj++)
      output[outnr++] = curPopIntegrals[pp][jj];

  return;
}


/*==================================================================================================================================*/
// Computes the current value of all POPULATIONINTEGRAL environmental variables

/*
  The macro constant DECEPENDENCYLEVEL indicates whether the POPULATIONINTEGRAL environmental variables are independent of each
  other (DECEPENDENCYLEVEL == 1; default) or have a uni-directional dependence (DECEPENDENCYLEVEL == 2), in which variable
  I_i may depend on other variables I_j, but none of these variables I_j in turn depends on I_i
*/

#ifndef DECEPENDENCYLEVEL
#define DECEPENDENCYLEVEL         1
#endif

static void setCurrentEnvironmentValues(double *env, population *pop, population *ofs, double  popImpacts[POPULATION_NR][INTERACT_DIM])

{
  int     dd, ii, jj, j, pp;
  int     lifeStage[POPULATION_NR], bstatenr;
  double  istateVals[POPULATION_NR][I_STATE_DIM];
  double *istatePnt[POPULATION_NR];
  double *bstatePnt[POPULATION_NR];
  double  curImpacts[POPULATION_NR][INTERACT_DIM];
  double  curPopIntegrals[POPULATION_NR][INTERACT_DIM];
  double  curEnvCondition[ENVIRON_DIM];
  double  cohortDensity;

  Time = env[0];
  memcpy(curPSPMEnv, env + 1, ENVIRON_DIM*sizeof(double));
  for (dd = 0; dd < DECEPENDENCYLEVEL; dd++)
    {
      memset(curPopIntegrals, 0, POPULATION_NR*INTERACT_DIM*sizeof(double));
      for (j = 0; j < MaxCohortNr; j++)
        {
          for (pp = 0; pp < POPULATION_NR; pp++)
            {
              jj = (j < cohort_no[pp]) ? j : cohort_no[pp] - 1;
              istatePnt[pp] = pop[pp][jj] + 1;                                      // Index 0 is cohort density
              bstatePnt[pp] = popIDcard[pp][jj] + 1;
              lifeStage[pp] = getPopIDlifestage(pp, jj);

              // For multiple populations and multiple birth states the next statement is incorrect
              bstatenr      = getPopIDbstatenr(pp, jj);    
            } 

          memset(curImpacts, 0, POPULATION_NR*INTERACT_DIM*sizeof(double));
#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
          SetLifeHistoryRates(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv,
                              NULL, NULL, NULL, curImpacts);   
#else
          Impact(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, curImpacts);
#endif

          for (pp = 0; pp < POPULATION_NR; pp++)
            if ((j < cohort_no[pp]) && (!isdead(pp, jj)))
              {
                cohortDensity = popIDcard[pp][jj][number]*exp(-(pop[pp][jj][number] - 1.0));                              
                for (ii = 0; ii < INTERACT_DIM; ii++)
                  curPopIntegrals[pp][ii] +=  curImpacts[pp][ii]*cohortDensity;
              }
        }

      // Add offspring cohort contributions
      if (ofs)
        {
          for (j = 0; j < MaxStatesAtBirth; j++)
            {
              for (pp = 0; pp < POPULATION_NR; pp++)
                {
                  jj = (j < bpoint_no[pp]) ? j : bpoint_no[pp] - 1;
#if (EBTMETHOD == 0)
                  memcpy(istateVals[pp], ofs[pp][jj] + 1, I_STATE_DIM*sizeof(double));
#else
                  memcpy(istateVals[pp], ofsIDcard[pp][jj] + 1, I_STATE_DIM*sizeof(double));

                  if (ofs[pp][jj][number] > DYTOL)
                    {
                      for (ii = 0; ii < I_STATE_DIM; ii++)
                        istateVals[pp][ii] += ofs[pp][jj][i_state(ii)]/ofs[pp][jj][number];
                    }
#endif
                  istatePnt[pp] = istateVals[pp];
                  bstatePnt[pp] = ofsIDcard[pp][jj] + 1;
                  lifeStage[pp] = getOfsIDlifestage(pp, jj);

                  // For multiple populations and multiple birth states the next statement is incorrect
                  bstatenr      = getOfsIDbstatenr(pp, jj);
                }

              memset(curImpacts, 0, POPULATION_NR*INTERACT_DIM*sizeof(double));
#if ((RFUNCTIONS == 1) || (MFUNCTIONS == 1))
              SetLifeHistoryRates(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv,
                                  NULL, NULL, NULL, curImpacts);
#else  
              Impact(lifeStage, istatePnt, bstatePnt, bstatenr, curPSPMEnv, curImpacts);
#endif

              for (pp = 0; pp < POPULATION_NR; pp++)
                if (j < bpoint_no[pp])
                  {
                    for (ii = 0; ii < INTERACT_DIM; ii++)
                      curPopIntegrals[pp][ii] +=  curImpacts[pp][ii]*ofs[pp][jj][number];
                  }
            }
        }

      memset(curEnvCondition, 0, ENVIRON_DIM*sizeof(double));
      EnvEqui(curPSPMEnv, curPopIntegrals, curEnvCondition);
      for (ii = 0; ii < ENVIRON_DIM; ii++)
        if (EnvironmentType[ii] == POPULATIONINTEGRAL) curPSPMEnv[ii] = curEnvCondition[ii];
    }

  if (popImpacts)
    memcpy(popImpacts, curPopIntegrals, POPULATION_NR*INTERACT_DIM*sizeof(double));
  
  return;
}


/*==================================================================================================================================*/
/*
 * Start of main routine implementations.
 */
/*==================================================================================================================================*/
 
/*
 * The error messages that occur in the routines in the present file.
 */

#define CSB  "Unable to open CSB file for complete state output!"
#define DBG  "Unable to open DBG file for event location information!"
#define MAFC "Memory allocation failure for cohort variables!"
#define MAFI "Memory allocation failure for cohort constants!"
#define OUT "Failure in opening OUT file for writing!"
#define SGE "Error in installing the signal handlers!"
 
 static void CatchSig(int signl)
 
 /*
  * CatchSig - Routine catches and handles the signals received from the
  *            ebttool application running this program as a child.
  */
 
 {
   char fpe_mes[80];
 
   if (signl == SIGFPE)
     {
       if (initState && currentState)
         {
           // In the middle of integration, save most recent data
           if (!dbgfile)                                                            // Open file if necessary
             {
               char filename[MAXFILENAMELEN];
 
               (void)strcpy(filename, runname);
               (void)strcat(filename, "dbg");
               dbgfile = fopen(filename, "a");
               if (!dbgfile) Warning("Unable to open DBG file!");
             }
           if (dbgfile)(void)fprintf(dbgfile, "==============================================================\n");
           if (dbgfile)(void)fprintf(dbgfile, "\nFPE EXCEPTION AT T = %.10f, STEP SIZE = %.6E\n", env[0], step_size);
           if (dbgfile && initState)
             {
               (void)fprintf(dbgfile, "\nInitial state of current integration step\n\n");
               WriteStateToFile(dbgfile, initState);
             }
           if (dbgfile && currentState)
             {
               (void)fprintf(dbgfile, "\nCurrent state in integration step (rk_level = %d)\n\n", rk_level);
               WriteStateToFile(dbgfile, currentState);
             }
           if (dbgfile && currentDers[rk_level - 1])
             {
               (void)fprintf(dbgfile, "\nLast derivative in integration step (rk_level = %d)\n\n", rk_level - 1);
               WriteStateToFile(dbgfile, currentDers[rk_level - 1]);
             }
           if (dbgfile && currentDers[rk_level])
             {
               (void)fprintf(dbgfile, "\nCurrent derivative in integration step (rk_level = %d)\n\n", rk_level);
               WriteStateToFile(dbgfile, currentDers[rk_level]);
             }
           if (dbgfile)(void)fprintf(dbgfile, "==============================================================\n");
         }
       if (EBTMethod) TransBcohorts();
       sprintf(fpe_mes, "Floating point error at time %.4f", env[0]);
       ErrorExit(1, fpe_mes);
     }
 
   if (signal(SIGFPE, CatchSig) == SIG_ERR) Warning(SGE);
 
   return;
 }
 
/*==================================================================================================================================*/

static void InitVars()

/* 
   * InitVars - Routine initializes all global variables to 0.
   */

{
  register int i;

  PopulationNr  = POPULATION_NR;
  Stages        = STAGES;
  IStateDim     = I_STATE_DIM;
  EnvironDim    = ENVIRON_DIM;
  InteractDim   = INTERACT_DIM;
  ParameterNr   = PARAMETER_NR;

  debug_level   = 0;
  report_level  = 1;
  EBTEnvironDim = ENVIRON_DIM_EBT;
  IConstDim     = I_CONST_DIM;
  output_var_nr = OUTPUT_VAR_NR;

  EBTMethod     = EBTMETHOD;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  CtrlCPressed  = 0;
#endif

  Odesolve_Init_Step  = ODESOLVE_INIT_STEP;
  Odesolve_Fixed_Step = ODESOLVE_FIXED_STEP;
  Odesolve_Min_Step   = ODESOLVE_MIN_STEP;
  Odesolve_Max_Step   = ODESOLVE_MAX_STEP;
  Odesolve_Abs_Err    = ODESOLVE_ABS_ERR;
  Odesolve_Rel_Err    = ODESOLVE_REL_ERR;
  Odesolve_Func_Tol   = ODESOLVE_FUNC_TOL;
  
  step_size = Odesolve_Init_Step;

  for (i = 0; i < ENVIRON_DIM_EBT; i++) env[i] = 0.0;                               // Create zero environment

  initState    = NULL;
  currentState = NULL;
  for (i = 0; i < MAXDERS; i++) currentDers[i] = NULL;
  ForcedRunEnd = 0;
  BifurcationRun = 0;
  csbnew = 1;
  identical_zero = DYTOL;

  for (i = 0; i < POPULATION_NR; i++)
    {
      CohortNo[i]         = 0;
      BpointNo[i]         = 0;
      pop_extinct[i]      = 0;
      pop[i]              = NULL;
      DataMemAllocated[i] = 0;
      popIDcard[i]        = NULL;
      IDMemAllocated[i]   = 0;
    }

  csbfile   = NULL;
  outfile   = NULL;
  averages  = NULL;
  gaverages = NULL;
  variances = NULL;
  extrema   = NULL;
  dbgfile   = NULL;

  return;
}

 
/*==================================================================================================================================*/
 
 int ComputeCurve(int argc, char **argv)
{
  int             i, j, colnr;
  char            csbname[MAXFILENAMELEN], dbgname[MAXFILENAMELEN], outname[MAXFILENAMELEN];
  double          oldBifParVal = 0;
  struct stat     buffer;
  char            tmpstr[MAX_STR_LEN];
  
  i = 0;
  while (1)
    {
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      sprintf(csbname, "%s-%s-%04d.mat", progname, "ECODYN", i);
#else
      sprintf(csbname, "%s-%s-%04d.csb", progname, "ECODYN", i);
#endif
      sprintf(dbgname, "%s-%s-%04d.err", progname, "ECODYN", i);
      sprintf(outname, "%s-%s-%04d.out", progname, "ECODYN", i);
      if (stat(csbname, &buffer) && stat(dbgname, &buffer) && stat(outname, &buffer)) break;
      i++;
    }
  sprintf(runname, "%s-%s-%04d", progname, "ECODYN", i);

  outfile = fopen(outname, "w");
  if (!outfile) ErrorAbort(OUT);

  // Open DBG file with lower case extension
  if (debug_level)
    {
      dbgfile = fopen(dbgname, "w");
      if (!dbgfile) Warning(DBG);                                                   // On error warn only
    }
  else dbgfile = NULL;

  csbnew  = 1;
  csbfile = fopen(csbname, "wb");
  if (!csbfile)
    {
      Warning(CSB);
      csbnew = 0;
    }

  fprintf(outfile, "#\n# Executing : ");

#if defined(R_PACKAGE)
#if (RFUNCTIONS == 1)
  fprintf(outfile, "PSPMecodyn(\"%s.R\", %s, %s, %s, %s, %s)", progname, "<Initial state>", timestring, bifstring, parstring, optstring);
#else
  fprintf(outfile, "PSPMecodyn(\"%s\", %s, %s, %s, %s, %s)", progname, "<Initial state>", timestring, bifstring, parstring, optstring);
#endif
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#if (MFUNCTIONS == 1)
  fprintf(outfile, "PSPMecodyn('%s.m', '%s', %s, %s, %s, %s)", progname, "<Initial state>", timestring, bifstring, parstring, optstring);
#else
  fprintf(outfile, "PSPMecodyn('%s', '%s', %s, %s, %s, %s)", progname, "<Initial state>", timestring, bifstring, parstring, optstring);
#endif
#else
  for (i = 0; i < argc; i++) fprintf(outfile, "%s ", argv[i]);
#endif
  fprintf(outfile, "\n#\n");

  fprintf(outfile, "# Parameter values  : \n#");
  for (i = 0; i < PARAMETER_NR; i++)
    {
      if (!(i % 3)) fprintf(outfile, "\n# ");
      fprintf(outfile, "\t%-10s:", parameternames[i]);
      fprintf(outfile, "  %-13G", parameter[i]);
    }
  fprintf(outfile, "\n#\n");
  fprintf(outfile, "#\tCohort cycle time interval                                         : %.1f\n", cohort_limit);
  fprintf(outfile, "#\tOutput time interval                                               : %.1f\n", delt_out);
  fprintf(outfile, "#\tComplete state output interval                                     : %.1f\n", state_out);
  fprintf(outfile, "#\tMaximum integration time                                           : %.1f\n", max_time);
  fprintf(outfile, "#\n");
  if (BifurcationRun)
    {
      fprintf(outfile, "#\tIndex of the bifurcation parameter                                 : %d\n", BifParIndex);
      fprintf(outfile, "#\tStarting value of the bifurcation parameter                        : %-13G\n", parameter[BifParIndex]);
      fprintf(outfile, "#\tStep size in  the bifurcation parameter                            : %-13G\n", BifParStep);
      fprintf(outfile, "#\tFinal value of the bifurcation parameter                           : %-13G\n", BifParLastVal);
      fprintf(outfile, "#\tPeriod of producing data output during each bifurcation interval   : %-13G\n", BifOutput);
      fprintf(outfile, "#\tPeriod of producing state output during each bifurcation interval  : %-13G\n", BifStateOutput);
    }
  fprintf(outfile, "#\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
  colnr = 0;
#else
  colnr = 1;
#endif
  sprintf(tmpstr, "%d:%s  ", colnr++, "Time");
  fprintf(outfile, "#%11s", tmpstr);
  for (i = 0; i < ENVIRON_DIM; i++)
    {
      sprintf(tmpstr, "%d:E[%d]", colnr++, i);
      fprintf(outfile, "%16s", tmpstr);
    }
  for (i = 0; i < POPULATION_NR; i++)
    {
      sprintf(tmpstr, "%d:b[%d]", colnr++, i);
      fprintf(outfile, "%16s", tmpstr);
    }
  for (i = 0; i < POPULATION_NR; i++)
    for (j = 0; j < INTERACT_DIM; j++)
      {
        sprintf(tmpstr, "%d:I[%d][%d]", colnr++, i, j);
        fprintf(outfile, "%16s", tmpstr);
      }
  if (BifurcationRun)
    {
      sprintf(tmpstr, "%d:Bif. par.", colnr++);
      fprintf(outfile, "%16s", tmpstr);
    }
  fprintf(outfile, "\n");
  fflush(outfile);

  if (signal(SIGFPE, CatchSig) == SIG_ERR) Warning(SGE);

  next_cohort_end = (floor(env[0]/cohort_limit) + 1)*cohort_limit;

  // By default (state) output at start up
  next_output = (floor((env[0] - identical_zero)/delt_out) + 1)*delt_out;
  if (state_out > 0.0)
    next_state_output = (floor((env[0] - identical_zero)/state_out) + 1)*state_out;
  else
    next_state_output = 0.0;

  if (BifurcationRun)
    {
      // Initialise the bifurcation and set the current parameter value
      SetBifOutputTimes(env);

      BifParBase = parameter[BifParIndex];

      parameter[BifParIndex] = (BifParBase + floor((env[0] + BIFTINY)/BifPeriod)*BifParStep);

      oldBifParVal = parameter[BifParIndex];
    }

  // Allow the user to carry out some initialization
  UserInit(argc, argv, env, pop);

  SievePop();

  if (BifurcationRun)
    {
      // Bifurcation parameter as additional output
      output_var_nr++;

      // Redo the initialisation in case the user has changed some of the bifurcation parameters
      SetBifOutputTimes(env);

      // Change the BifParBase only when user has changed the parameter value in UserInit()
      if (oldBifParVal != parameter[BifParIndex]) BifParBase = parameter[BifParIndex];

      parameter[BifParIndex] = (BifParBase + floor((env[0] + BIFTINY)/BifPeriod)*BifParStep);
    }

  if ((next_output - env[0]) < identical_zero)
    {
      // Increment next time. Produce output if required at this time
      next_output += delt_out;
      FileOut();
    }

  if ((state_out > 0.0) && ((next_state_output - env[0]) < identical_zero))
    {
      // Increment next time. Output complete state if required at this time
      next_state_output += state_out;
      FileState();
    }

  while (!((env[0] >= max_time) || ForcedRunEnd))
    {
      CohortCycle(next_cohort_end);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      if (checkInterrupt()) break;
#endif
    }

  // Program shut down procedure
  ShutDown(0);

  return 0;
}
 
  
/*
 *====================================================================================================================================
 *  UNIX shell interface function main() and supporting functions.
 *====================================================================================================================================
 */

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)

int main(int argc, char **argv)
{
  /*
   * Possible exception types to set are
   *   (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW | FE_INEXACT);
   * The last two are, however, frequently occurring in mathematical functions
   * without any consequences. Only enable testing the first 3.
   */
  feclearexcept(FE_ALL_EXCEPT);
#if defined(__APPLE__)
  static fenv_t fenv;
  unsigned int  new_excepts;                                                        // previous masks

  new_excepts = (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW) & FE_ALL_EXCEPT,

  fegetenv(&fenv);

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr &= ~(new_excepts << 7);
  fesetenv(&fenv);

#else
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif                                                                              // defined(__APPLE__)

  if (argc)
    {
      strcpy(progname, argv[0]);
      progname[strlen(progname) - 6] = '\0';                                        // Cut off the 'ecodyn' appendix    
    }

  // Initialization of the environment, population and output devices
  InitVars();
  Initialize(argc, argv);

  // call the computational routine
  ComputeCurve(0, NULL);

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
  // Program shut down procedure
  ShutDown(0);

#if defined(MATLAB_MEX_FILE) && (FULLSTATEOUTPUT > 0)
   if (pmat) matClose(pmat);
   pmat = NULL;
#endif

  return;
 }
 
 
 void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
  size_t          nrows, ncols;
  double          tmpdouble;
  int             irhs, field_num, status;
  const mxArray   *element_ptr;
  int             bb, ii, jj, pp, num, maxBstates, maxCohorts;
  int             parsdefined, bstatesdefined;
  int             bStateNr[POPULATION_NR];
  double          bstateMem[POPULATION_NR*I_STATE_DIM], *bstatePnt[POPULATION_NR];
  double          *dblpnt, *bStates = NULL;
  double          val_tmp[COHORT_SIZE + I_CONST_DIM];
  long            mem_req;
  char            optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN];
  unsigned char   pmag, bmag;

  mwIndex         i;
  size_t          total_num_of_cells, buflen;

  InitVars();
  for (pp = 0; pp < POPULATION_NR; pp++) bstatePnt[pp] = bstateMem + pp*I_STATE_DIM;

  // check for proper number of arguments
  if (nrhs != 5)
    mexErrMsgIdAndTxt(
        "MATLAB:PSPMecodyn:nrhs",
        "\nIncorrect number of command-line arguments.\n\nUse: %s(%s)\n\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s",
        mexFunctionName(), "<state>, <time settings>, <bifurcation settings>, <parameters>, <options>",
        "<state>", "Initial state of the environment and populations",
        "<time settings>", "Array of 4 values: cohort cyle time, interval of data output, interval of state output, maximum integration time",
        "<bifurcation settings>", "Empty array or array of 6 values: index, starting value, step size and final value of bifurcation parameter followed by period of data and state output"
        "<parameters>", "Array of parameter values to use (empty array or PARAMETER_NR in length)",
        "<options>", "Cell array of possible options: info (followed by info level 1, 2, 3 or 4)");

  // check for proper number of output variables
  if (nlhs != 1) mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:nlhs", "A single output argument is required.");

#if (MFUNCTIONS == 1)
  Minterface_Init();
#endif

  //============================== Process the options argument ======================================================================
  // Extract the contents of MATLAB cell into the C array
  irhs = 4;
  if (!mxIsCell(prhs[irhs])) mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nOptions should be specified as a cell array!\n");

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(optstring, "{");
  for (i = 0; i < total_num_of_cells; i++)
    {
      element_ptr = mxGetCell(prhs[irhs], i);
      buflen      = mxGetN(element_ptr)*sizeof(mxChar) + 1;
      status      = mxGetString(element_ptr, optname, buflen);
      // optstring still equal to "{"
      if (strlen(optstring) == 1)
        { 
          strcat(optstring, "'"); strcat(optstring, optname); strcat(optstring, "'");
        }
      else
        {
          strcat(optstring, ", '"); strcat(optstring, optname); strcat(optstring, "'");
        }
      if (!strcmp(optname, "info"))
        {
          if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nNo value specified for option %s!\n", optname);

          element_ptr = mxGetCell(prhs[irhs], i);
          buflen      = mxGetN(element_ptr)*sizeof(mxChar) + 1;
          status      = mxGetString(element_ptr, optval, buflen);
          if (status) mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nError in retrieving value for option %s!\n", optname);

          strcat(optstring, ", '"); strcat(optstring, optval); strcat(optstring, "'");

          debug_level = atoi(optval);
          if ((debug_level < 1) || (debug_level > 4))
            mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nOutput level of the DOPRI5 time integrator (%d) not in the appropriate range (0 < i < 5)!\n\n", debug_level);
        }
      else if (!strcmp(optname, "report"))
        {
          if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nNo value specified for option %s!\n", optname);

          element_ptr = mxGetCell(prhs[irhs], i);
          buflen      = mxGetN(element_ptr)*sizeof(mxChar) + 1;
          status      = mxGetString(element_ptr, optval, buflen);
          if (status) mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nError in retrieving value for option %s!\n", optname);

          strcat(optstring, ", '"); strcat(optstring, optval); strcat(optstring, "'");

          report_level = atoi(optval);
          if (report_level < 1)
            mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nInterval of reporting data output to console (%d) should be at least 1!\n\n", report_level);
        }
      else
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:options", "\nIllegal option %s!\n", optname);
    }
  strcat(optstring, "}");

  //============================== Process the parameters argument ===================================================================

  parsdefined = 0;
  irhs        = 3;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols == PARAMETER_NR) && (nrows == 1))
    {
      memcpy(parameter, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));   // Overrides parameters from initial state
      parsdefined = 1;
    }
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:PSPMecodyn:parameters", "\nParameter argument ignored as it is not a row vector of length PARAMETER_NR\n");

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

  //============================== Process the time settings argument ================================================================

  irhs  = 1;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols == 4) && (nrows == 1))
    {
      dblpnt = mxGetPr(prhs[irhs]);

      cohort_limit  = dblpnt[0];
      delt_out      = dblpnt[1];
      state_out     = dblpnt[2];
      max_time      = dblpnt[3];
      if (cohort_limit <= 0.0)
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:times", "\nCohort cylce time should be positive!\n\n");
      if (delt_out < cohort_limit)
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:times", "\nInterval for data output to .out file should be larger than or equal to cohort cylce time!\n\n");
      if ((state_out != 0) && (state_out < cohort_limit))
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:times", "\nInterval for complete state output to .csb file should either be 0 or larger than or equal to cohort cylce time!\n\n");
      if (max_time < cohort_limit)
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:times", "\nMaximum integration time should be larger than cohort or equal to cylce time!\n\n");
    }
  else if (ncols)
    mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:times", "\n%s\n%s\n"
                      "Time settings argument ignored as it is not a row vector of length 4 consisting of",
                      "cohort cycle time, interval between data output, interval between state output and maximum integration time");

  strcpy(timestring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(timestring, " ");
      sprintf(tmpstr, "%.6G", dblpnt[i]);
      strcat(timestring, tmpstr);
    }
  strcat(timestring, "]");

  //============================== Process the bifurcation settings argument =========================================================
  // max_time has to be set before this point

  irhs  = 2;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  BifurcationRun = 0;
  if ((ncols == 6) && (nrows == 1))
    {
      dblpnt = mxGetPr(prhs[irhs]);

      BifParIndex             = lrint(dblpnt[0]);
      if ((BifParIndex < 0) || (BifParIndex > PARAMETER_NR))
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:bifurcation", "\nIndex of bifurcation parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", BifParIndex, PARAMETER_NR);

      parameter[BifParIndex]  = dblpnt[1];

      BifParStep              = dblpnt[2];
      // Sanitize the bifurcation control variable: Make it absolute
      BifParStep = fabs(BifParStep);
      if (BifParStep < Odesolve_Min_Step)
        mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:bifurcation", "\nStep size in bifurcation parameter (%G) should be larger than %G!\n\n", BifParStep, Odesolve_Min_Step);

      BifParLastVal           = dblpnt[3];
      // Set the sign of bifurcation step size
      if (parameter[BifParIndex] > BifParLastVal) BifParStep *= -1;

      BifOutput               = dblpnt[4];
      // Sanitize the bifurcation control variable
      BifOutput = max(BifOutput, 0.0);
      BifOutput = min(BifOutput, max_time);

      BifStateOutput          = dblpnt[5];
      // Sanitize the bifurcation control variable
      BifStateOutput = max(BifStateOutput, 0.0);
      BifStateOutput = min(BifStateOutput, max_time);

      BifurcationRun = 1;
    }
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:PSPMecodyn:bifurcation", "\n%s\n%s\n"
                       "Bifurcation settings argument ignored as it is not a row vector of length 6 consisting of",
                       "index, starting value, step size and final value of the bifurcation parameter and the period of data and state output");

  strcpy(bifstring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(bifstring, " ");
      sprintf(tmpstr, "%.6G", dblpnt[i]);
      strcat(bifstring, tmpstr);
    }
  strcat(bifstring, "]");

  //============================== Process the initial point argument ================================================================

  irhs = 0;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);

  if ((!mxIsStruct(prhs[irhs])) || (nrows != 1) || (ncols != 1))
    mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:init", "%s\n%s\n",
                      "Initial state should be a single structure with elements 'Parameters' (optional), 'Environment'",
                      "'PopNN_Birthstates' (optional) and 'PopNN' (or 'PopNN_BstateMM' in case of multiple states at birth) for all populations");

  // Environment is required in the initial state
  field_num = mxGetFieldNumber(prhs[irhs], "Environment");
  if (field_num > -1)
    {
      element_ptr = mxGetFieldByNumber(prhs[irhs], 0, field_num);
      if ((element_ptr != NULL) && (mxGetClassID(element_ptr) == mxDOUBLE_CLASS)  && (!mxIsComplex(element_ptr)) &&
          (mxGetM(element_ptr) == 1) && (mxGetN(element_ptr) == ENVIRON_DIM))
        {
          memcpy(env + 1, mxGetPr(element_ptr), ENVIRON_DIM*sizeof(double));
        }
      else
        field_num = -1;
    }
  if (field_num < 0)
    mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:init", "%s\n",
                      "Initial state structure does not contain a valid element 'Environment' of length ENVIRON_DIM");
                     
  // Parameters are optional. Command-line parameter overrule parameters stored in the initial state
  if (parsdefined != 1)
    {
      field_num = mxGetFieldNumber(prhs[irhs], "Parameters");
      if (field_num > -1)
        {
          element_ptr = mxGetFieldByNumber(prhs[irhs], 0, field_num);
          if ((element_ptr != NULL) && (mxGetClassID(element_ptr) == mxDOUBLE_CLASS)  && (!mxIsComplex(element_ptr)) &&
              (mxGetM(element_ptr) == 1) && (mxGetN(element_ptr) == PARAMETER_NR))
            {
              memcpy(parameter, mxGetPr(element_ptr), PARAMETER_NR*sizeof(double));
              parsdefined = 1;
            }
        }
      if (parsdefined != 1)
        mexWarnMsgIdAndTxt("MATLAB:PSPMecodyn:init", "%s\n%s\n",
                           "Initial state structure does not contain a valid element 'Parameters' of length PARAMETER_NR",
                           "Parameter values taken from model file");
    }

  pmag = 1; num  = POPULATION_NR; while (num > 0) { pmag++; num = num/10; }
  for (pp = 0; pp < POPULATION_NR; pp++)
    {
      bstatesdefined = 0;
      sprintf(tmpstr, "Pop%0*d_BirthStates", pmag, pp);
      field_num = mxGetFieldNumber(prhs[irhs], tmpstr);
      if (field_num > -1)
        {
          element_ptr = mxGetFieldByNumber(prhs[irhs], 0, field_num);
          if ((element_ptr != NULL) && (mxGetClassID(element_ptr) == mxDOUBLE_CLASS)  && (!mxIsComplex(element_ptr)) &&
              (mxGetM(element_ptr) > 0) && (mxGetN(element_ptr) == I_STATE_DIM))
            {
              maxBstates = mxGetM(element_ptr);
              bStates = mxGetPr(element_ptr);
              bstatesdefined = 1;
            }
        }
      if (bstatesdefined != 1)
        {
          mexWarnMsgIdAndTxt("MATLAB:PSPMecodyn:init",
                             "\nBirth states of population Pop%0*d not specified, will be obtained through calls to SetBirthStates() and StateAtBirth()\n\n",
                             pmag, pp);
          SetBirthStates(bStateNr, env + 1);
          maxBstates = bStateNr[pp];
        }

      bmag = 1; num  = maxBstates; while (num > 0) { bmag++; num = num/10; }
      for (bb = 0; bb < maxBstates; bb++)
        {
          if (maxBstates == 1) sprintf(tmpstr, "Pop%0*d", pmag, pp);
          else sprintf(tmpstr, "Pop%0*d_Bstate%0*d", pmag, pp, bmag, bb);

          field_num = mxGetFieldNumber(prhs[irhs], tmpstr);
          if (field_num > -1)
            {
              element_ptr = mxGetFieldByNumber(prhs[irhs], 0, field_num);
              if ((element_ptr != NULL) && (mxGetClassID(element_ptr) == mxDOUBLE_CLASS)  && (!mxIsComplex(element_ptr)) &&
                  (mxGetM(element_ptr) > 0) && (mxGetN(element_ptr) >= (I_STATE_DIM + 1)))
                {
                  maxCohorts = mxGetM(element_ptr);
                  dblpnt = mxGetPr(element_ptr);
                  for (jj = 0; jj < maxCohorts; jj++)
                    {
                      for (ii = 0; ii < COHORT_SIZE; ii++)
                        val_tmp[ii] = dblpnt[ii*maxCohorts + jj];
                      val_tmp[I_STATE_DIM + 1] = MISSING_VALUE;
                      if (bstatesdefined == 1)
                        {
                          for (ii = 0; ii < I_STATE_DIM; ii++)
                            val_tmp[I_STATE_DIM + 2 + ii] = bStates[ii*maxBstates + bb];
                        }
                      else
                        {
                          StateAtBirth(bstatePnt, bb, env + 1);
                          for (ii = 0; ii < I_STATE_DIM; ii++)
                            val_tmp[I_STATE_DIM + 2 + ii] = bstateMem[pp*I_STATE_DIM + ii];
                        }
                      val_tmp[I_STATE_DIM + 2 + I_STATE_DIM    ] = 0;
                      val_tmp[I_STATE_DIM + 2 + I_STATE_DIM + 1] = bb;
                      val_tmp[I_STATE_DIM + 2 + I_STATE_DIM + 2] = 0.0;

                      mem_req = (CohortNo[pp] + 1)*COHORT_SIZE;
                      if (!(mem_req < DataMemAllocated[pp]))
                        {
                          DataMemAllocated[pp] = MemBlocks(mem_req);
                          pop[pp]              = (population)Myalloc((void *)pop[pp], (size_t)DataMemAllocated[pp], sizeof(double));
                          if (!(pop[pp])) ErrorAbort(MAFC);
                        }
                      mem_req = (CohortNo[pp] + 1)*I_CONST_DIM;
                      if (!(mem_req < IDMemAllocated[pp]))
                        {
                          IDMemAllocated[pp] = MemBlocks(mem_req);
                          popIDcard[pp]      = (popID)Myalloc((void *)popIDcard[pp], (size_t)IDMemAllocated[pp], sizeof(double));
                          if (!(popIDcard[pp])) ErrorAbort(MAFI);
                        }
                      if (CohortNo[pp])
                        CohortNo[pp] = InsCohort(val_tmp, val_tmp + COHORT_SIZE, pp);
                      else
                        {
                          for (ii = 0; ii < COHORT_SIZE; ii++) pop[pp][CohortNo[pp]][ii] = val_tmp[ii];
                          for (ii = 0; ii < I_CONST_DIM; ii++) popIDcard[pp][CohortNo[pp]][ii] = val_tmp[ii + COHORT_SIZE];
                          CohortNo[pp]++;
                        }
                  }
                }
              else
                field_num = -1;
            }
          if (field_num < 0)
            mexErrMsgIdAndTxt("MATLAB:PSPMecodyn:init",
                              "Initial state structure does not contain a valid matrix element %s with %d columns\n",
                              tmpstr, I_STATE_DIM + 1);
        }
      cohort_no[pp] = CohortNo[pp];
    }

  //============================= Get the program name ===============================================================================
  // Get the name of the mex file
  strcpy(progname, mexFunctionName());
  progname[strlen(progname) - 6] = '\0';                                            // Cut off the 'ecodyn' appendix

  mexAtExit(CloseStreams);

  // call the computational routine
  ComputeCurve(0, NULL);

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

SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int  i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
      {
        elmt = VECTOR_ELT(list, i);
        break;
      }

  if (elmt == R_NilValue) return R_NilValue;

  return elmt;
}


SEXP PSPMecodyn(SEXP moduleName, SEXP initState, SEXP timePars, SEXP bifPars, SEXP parVals, SEXP optVals)

{
  char            tmpstr[MAX_STR_LEN], optname[MAX_STR_LEN], optval[MAX_STR_LEN];
  int             bb, ii, jj, ncols, pp, num, maxBstates, maxCohorts;
  int             parsdefined, bstatesdefined;
  int             bStateNr[POPULATION_NR];
  double          bstateMem[POPULATION_NR*I_STATE_DIM], *bstatePnt[POPULATION_NR];
  double          *dblpnt, *bStates = NULL;
  double          val_tmp[COHORT_SIZE + I_CONST_DIM];
  long            mem_req;
  unsigned char   pmag = 1, bmag = 1;
  SEXP            resfil;
  SEXP            R_listel, R_pnt;
  SEXP            dim;


  InitVars();
  for (pp = 0; pp < POPULATION_NR; pp++) bstatePnt[pp] = bstateMem + pp*I_STATE_DIM;

#if (RFUNCTIONS == 1)
  Rinterface_Init();
#endif

  // init <- csbread("Checks/Initstate", 1); init5 <- csbread("Checks/Initstate5.csb", 1); PSPMecodyn("PNAS2002", init, c(1.0, 1.0, 0.0, 2000.0), NULL, NULL, c("report", "10"), clean = TRUE, force = TRUE)

  //============================== Process the options argument ======================================================================

  strcpy(optstring, "");

  ncols =length(optVals);

  for (ii = 0; ii < ncols; ii++)
    {
      strcpy(optname, CHAR(STRING_ELT(optVals, ii)));

      // optstring still equal to ""
      if (!strlen(optstring))
        {
          strcpy(optstring, "c(\""); strcat(optstring, optname); strcat(optstring, "\"");
        }
      else
        {
          strcat(optstring, ", \""); strcat(optstring, optname); strcat(optstring, "\"");
        }

      if (!strcmp(optname, "info"))
        {
          if (!(++ii < ncols))
            error("\nOutput level of the DOPRI5 time integrator (option \"%s\") not specified!\n\n", optname);

          strcpy(optval, CHAR(STRING_ELT(optVals, ii)));
          strcat(optstring, ", \""); strcat(optstring, optval); strcat(optstring, "\"");

          debug_level = atoi(optval);
          if ((debug_level < 1) || (debug_level > 4))
            error("\nOutput level of the DOPRI5 time integrator (%d) not in the appropriate range (0 < i < 4)!\n\n", debug_level);
        }
      else if (!strcmp(optname, "report"))
        {
          if (!(++ii < ncols))
            error("\nInterval for reporting data output to console (option \"%s\") not specified!\n\n", optname);

          strcpy(optval, CHAR(STRING_ELT(optVals, ii)));
          strcat(optstring, ", \""); strcat(optstring, optval); strcat(optstring, "\"");

          report_level = atoi(optval);
          if (report_level < 1)
            error("\nInterval of reporting data output to console (%d) should be at least 1!\n\n", report_level);
        }
      else
        error("\nIllegal option %s!\n", optname);
    }

  if (strlen(optstring))
    strcat(optstring, ")");
  else
    strcpy(optstring, "NULL");

  //============================== Process the parameters argument ===================================================================

  strcpy(parstring, "NULL");
  parsdefined = 0;
  ncols       = length(parVals);
  if (isReal(parVals))
    {
      dblpnt = REAL(parVals);
      if (ncols == PARAMETER_NR)
        {
          memcpy(parameter, dblpnt, ncols*sizeof(double));
          parsdefined = 1;
        }

      if (ncols)
        {
          strcpy(parstring, "c(");
          for (ii = 0; ii < ncols; ii++)
            {
              if (ii) strcat(parstring, ", ");
              sprintf(tmpstr, "%.6G", dblpnt[ii]);
              strcat(parstring, tmpstr);
            }
          strcat(parstring, ")");
        }
    }

  if (ncols && (!parsdefined))
    warning("\nParameter argument ignored as it is not a real-valued vector of length PARAMETER_NR\n\n");

  //============================== Process the time parameter argument ===============================================================

  strcpy(timestring, "NULL");
  ncols  = length(timePars);

  if ((!isReal(timePars)) || (ncols != 4))
    error("\nTime settings argument is not a real-valued vector of length 4\n\n");

  dblpnt = REAL(timePars);
  cohort_limit  = dblpnt[0];
  delt_out      = dblpnt[1];
  state_out     = dblpnt[2];
  max_time      = dblpnt[3];
  if (cohort_limit <= 0.0)
    error("\nCohort cylce time should be positive!\n\n");
  if (delt_out < cohort_limit)
    error("\nInterval for data output to .out file should be larger than or equal to cohort cylce time!\n\n");
  if ((state_out != 0) && (state_out < cohort_limit))
    error("\nInterval for complete state output to .csb file should either be 0 or larger than or equal to cohort cylce time!\n\n");
  if (max_time < cohort_limit)
    error("\nMaximum integration time should be larger than cohort or equal to cylce time!\n\n");
  strcpy(timestring, "c(");
  for (ii = 0; ii < ncols; ii++)
    {
      if (ii) strcat(timestring, ", ");
      sprintf(tmpstr, "%.6G", dblpnt[ii]);
      strcat(timestring, tmpstr);
    }
  strcat(timestring, ")");

  //============================== Process the bifurcation parameter argument ========================================================
  // max_time has to be set before this point

  strcpy(bifstring, "NULL");
  BifurcationRun = 0;
  ncols  = length(bifPars);

  if (ncols && ((!isReal(bifPars)) || (ncols != 6)))
    error("\nBifurcation settings argument is not a real-valued vector of length 6\n\n");
  else if (ncols)
    {
      dblpnt = REAL(bifPars);
      BifParIndex             = lrint(dblpnt[0]);
      if ((BifParIndex < 0) || (BifParIndex > PARAMETER_NR))
        error("\nIndex of bifurcation parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", BifParIndex, PARAMETER_NR);

      parameter[BifParIndex]  = dblpnt[1];

      BifParStep              = dblpnt[2];
      // Sanitize the bifurcation control variable: Make it absolute
      BifParStep = fabs(BifParStep);
      if (BifParStep < Odesolve_Min_Step)
        error("\nStep size in bifurcation parameter (%G) should be larger than %G!\n\n", BifParStep, Odesolve_Min_Step);

      BifParLastVal           = dblpnt[3];
      // Set the sign of bifurcation step size
      if (parameter[BifParIndex] > BifParLastVal) BifParStep *= -1;

      BifOutput               = dblpnt[4];
      // Sanitize the bifurcation control variable
      BifOutput = max(BifOutput, 0.0);
      BifOutput = min(BifOutput, max_time);

      BifStateOutput          = dblpnt[5];
      // Sanitize the bifurcation control variable
      BifStateOutput = max(BifStateOutput, 0.0);
      BifStateOutput = min(BifStateOutput, max_time);

      BifurcationRun = 1;

      strcpy(bifstring, "c(");
      for (ii = 0; ii < ncols; ii++)
        {
          if (ii) strcat(bifstring, ", ");
          sprintf(tmpstr, "%.6G", dblpnt[ii]);
          strcat(bifstring, tmpstr);
        }
      strcat(bifstring, ")");
    }

  //============================== Process the initial state argument ================================================================

  // Environment is required in the initial state
  PROTECT(R_listel = getListElement(initState, "Environment"));
  if (R_listel == R_NilValue)
    error("\nInitial state list should contain an element \"Environment\" with initial values for the environmental variables\n\n");
  PROTECT(R_pnt = coerceVector(R_listel, REALSXP));
  if (length(R_pnt) == ENVIRON_DIM)
    memcpy(env + 1, REAL(R_pnt), ENVIRON_DIM*sizeof(double));
  else
    error("\nInitial state for the environmental variables should be a vector of length %d\n\n", ENVIRON_DIM);
  UNPROTECT(2);

  // Parameters are optional. Command-line parameter overrule parameters stored in the initial state 
  if (parsdefined != 1)
    {
      PROTECT(R_listel = getListElement(initState, "Parameters"));
      if (R_listel != R_NilValue)
        {
          PROTECT(R_pnt = coerceVector(R_listel, REALSXP));
          if (length(R_pnt) == PARAMETER_NR)
            memcpy(parameter, REAL(R_pnt), PARAMETER_NR*sizeof(double));
          else
            warning("\nInitial state element \"Parameters\" with parameter values is not of length PARAMETER_NR\n%s\n\n",
                    "Parameter values taken from model file");
          UNPROTECT(1);
        }
      else
        warning("\nInitial state list does not contain an element \"Parameters\" with parameter values\n%s\n\n",
                "Parameter values taken from model file");
      UNPROTECT(1);
    }

  pmag = 1; num  = POPULATION_NR; while (num > 0) { pmag++; num = num/10; }
  for (pp = 0; pp < POPULATION_NR; pp++)
    {
      // If birth states not specified, generate them via a call to StateAtBirth()
      bstatesdefined = 0;
      sprintf(tmpstr, "Pop%0*d_BirthStates", pmag, pp);
      PROTECT(R_listel = getListElement(initState, tmpstr));
      if (R_listel != R_NilValue)
        {
          PROTECT(R_pnt = coerceVector(R_listel, REALSXP));
          maxBstates = length(R_pnt)/I_STATE_DIM;
          bmag = 1; num  = maxBstates; while (num > 0) { bmag++; num = num/10; }

          if (length(R_pnt) == maxBstates*I_STATE_DIM)
            {
              bStates = REAL(R_pnt);
              bstatesdefined = 1;
            }
          UNPROTECT(1);
        }
      UNPROTECT(1);
      if (bstatesdefined != 1)
        {
          warning("\nBirth states of population Pop%0*d not specified, will be obtained through calls to SetBirthStates() and StateAtBirth()\n\n", pmag, pp);
          SetBirthStates(bStateNr, env + 1);
          maxBstates = bStateNr[pp];
        }

      for (bb = 0; bb < maxBstates; bb++)
        {
          if (maxBstates == 1)
            sprintf(tmpstr, "Pop%0*d", pmag, pp);
          else 
            sprintf(tmpstr, "Pop%0*d_Bstate%0*d", pmag, pp, bmag, bb);

          PROTECT(R_listel = getListElement(initState, tmpstr));
          if (R_listel == R_NilValue)
            error("\nInitial state list should contain an element \"%s\" with the initial state values for population %d and birth state %d\n\n", tmpstr, pp, bb);

          dim = getAttrib(R_listel, R_DimSymbol);
          maxCohorts = INTEGER(dim)[0];
          ncols      = INTEGER(dim)[1];
          if ((ncols < (I_STATE_DIM + 1)) || (maxCohorts < 1))
            error("\nThe initial state values for population %d and birth state %d should be specified as a matrix with at least %d columns\n\n", pp, bb, I_STATE_DIM + 1);

          PROTECT(R_pnt = coerceVector(R_listel, REALSXP));

          dblpnt = REAL(R_pnt);
          for (jj = 0; jj < maxCohorts; jj++)
            {
              for (ii = 0; ii < (I_STATE_DIM + 1); ii++)
                val_tmp[ii] = dblpnt[ii*maxCohorts + jj];
              val_tmp[I_STATE_DIM + 1] = MISSING_VALUE;
              if (bstatesdefined == 1)
                {
                  for (ii = 0; ii < I_STATE_DIM; ii++)
                    val_tmp[I_STATE_DIM + 2 + ii] = bStates[ii*maxBstates + bb];
                }
              else
                {
                  StateAtBirth(bstatePnt, bb, env + 1);
                  for (ii = 0; ii < I_STATE_DIM; ii++)
                    val_tmp[I_STATE_DIM + 2 + ii] = bstateMem[pp*I_STATE_DIM + ii];
                }
              val_tmp[I_STATE_DIM + 2 + I_STATE_DIM    ] = 0;
              val_tmp[I_STATE_DIM + 2 + I_STATE_DIM + 1] = bb;
              val_tmp[I_STATE_DIM + 2 + I_STATE_DIM + 2] = 0.0;

              mem_req = (CohortNo[pp] + 1)*COHORT_SIZE;
              if (!(mem_req < DataMemAllocated[pp]))
                {
                  DataMemAllocated[pp] = MemBlocks(mem_req);
                  pop[pp]              = (population)Myalloc((void *)pop[pp], (size_t)DataMemAllocated[pp], sizeof(double));
                  if (!(pop[pp])) ErrorAbort(MAFC);
                }
              mem_req = (CohortNo[pp] + 1)*I_CONST_DIM;
              if (!(mem_req < IDMemAllocated[pp]))
                {
                  IDMemAllocated[pp] = MemBlocks(mem_req);
                  popIDcard[pp]      = (popID)Myalloc((void *)popIDcard[pp], (size_t)IDMemAllocated[pp], sizeof(double));
                  if (!(popIDcard[pp])) ErrorAbort(MAFI);
                }
              if (CohortNo[pp])
                CohortNo[pp] = InsCohort(val_tmp, val_tmp + COHORT_SIZE, pp);
              else
                {
                  for (ii = 0; ii < COHORT_SIZE; ii++) pop[pp][CohortNo[pp]][ii] = val_tmp[ii];
                  for (ii = 0; ii < I_CONST_DIM; ii++) popIDcard[pp][CohortNo[pp]][ii] = val_tmp[ii + COHORT_SIZE];
                  CohortNo[pp]++;
                }
            }
          UNPROTECT(2);
        }
      cohort_no[pp] = CohortNo[pp];
    }

  //============================= Get the program name ===============================================================================
  // Get the name of the module file
  if (isString(moduleName) && (length(moduleName) == 1))
    strcpy(progname, CHAR(STRING_ELT(moduleName, 0)));
  else
    error("\nModel name argument must be a single string\n\n");

  // call the computational routine
  ComputeCurve(0, NULL);

#if (RFUNCTIONS == 1)
  Rinterface_End();
#endif

  PROTECT(resfil = mkString(runname));

  UNPROTECT(1);

  return resfil;
}


/*==================================================================================================================================*/
#endif
