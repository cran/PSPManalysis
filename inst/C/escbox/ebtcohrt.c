/*
  NAME
     ebtcohrt.c

     This file contains NewCohort() routine and all the functions used by
     it to keep the set of population cohorts organized and to deal with
     their increasing and decreasing number. 

  Last modification: AMdR - Dec 15, 2017
*/

#define EBTCOHRT_C				                                                          // Identification of file
#define EBTLIB					                                                            // and file grouping

#include "escbox.h"
#include "ebtmain.h"
#include "ebtcohrt.h"
#include "ebtdopri5.h"
#include "ebtutils.h"





/*==================================================================================================================================*/
/*
 * The error messages that occur in the routines in the present file.
 */

#define AEXT "All populations are extinct!"
#define MAFB "Memory allocation failure for boundary points!"
#define MAFC "Memory allocation failure for i-state variables!"
#define MAFI "Memory allocation failure for cohort constants!"
#define MAFT "Memory allocation failure while inserting boundary cohorts!"



/*==================================================================================================================================*/
/*
 * Definitions of static variables, restricted to this file.
 */

static population	                BPData = NULL;		                                // Pointer to bpoints data
static long		                    BPAllocated = 0L;	                                // # of doubles allocated
static cohort_pnt                 BcohortsPopPnt = NULL;
static cohortID_pnt               BcohortsIDPnt = NULL;
static int                        MaxBcohorts = -1;
static double		                  maxss, minss;




/*==================================================================================================================================*/
/*
 * Start of function implementations.
 */
/*==================================================================================================================================*/

int	  InsCohort(cohort newcoh, cohortID id, int pop_nr)

  /* 
   * InsCohort - Routine inserts a new cohort on its appropriate, ordered 
   *		 position in the population with number "pop_nr". The 
   *		 pointer "new" points to an array of i-state variables of
   *		 length COHORT_SIZE, the pointer "id" points to an array
   *		 of i-state constants of length I_CONST_DIM.
   */

{
  register int      found, equal;
  register int      i, pos;
  cohort_pnt        p;
  cohortID_pnt      pid;

  // Top-down search for appropriate cohort place
  p = pop[pop_nr];
  for (pos = CohortNo[pop_nr] - 1, found = 0; pos >= 0; pos--)
    {
      for (i = 1, equal = 1; (i < COHORT_SIZE) && equal; i++)
        {
          equal = (p[pos][i] == newcoh[i]);
          found = (p[pos][i] > newcoh[i]);
        }
      if (found) break;
    }

  // Shift all entries lying above position upwards
  pos++;
  p = pop[pop_nr] + pos;
  if (pos < CohortNo[pop_nr]) (void)memmove(*(p + 1), *p, (CohortNo[pop_nr] - pos)*COHORT_SIZE*sizeof(double));

  // Place new cohort
  (void)memcpy((DEF_TYPE *)*p, (DEF_TYPE *)newcoh, COHORT_SIZE*sizeof(double));
  pid = popIDcard[pop_nr] + pos;
  if (pos < CohortNo[pop_nr]) (void)memmove(*(pid + 1), *pid, (CohortNo[pop_nr] - pos)*I_CONST_DIM*sizeof(double));

  // Place new cohort id
  (void)memcpy((DEF_TYPE *)*pid, (DEF_TYPE *)id, I_CONST_DIM*sizeof(double));

  return (++CohortNo[pop_nr]);
}


/*==================================================================================================================================*/

static void	CreateBcohorts()

  /* 
   * CreateBcohorts - Routine sets up the boundary cohorts and the fixed
   *		      boundary points.
   */

{
  register int i, base;
  long         mem_req, bp_req;

  for (i = 0; i < POPULATION_NR; i++)
    {
      if (BpointNo[i])
        {
          mem_req = (CohortNo[i] + BpointNo[i])*COHORT_SIZE;

          // Create memory for new cohorts if necessary
          if (!(mem_req < DataMemAllocated[i]))
            {
              DataMemAllocated[i] = MemBlocks(mem_req);
              pop[i]              = (population)Myalloc((void *)pop[i], (size_t)DataMemAllocated[i], sizeof(double));
              if (!(pop[i])) ErrorAbort(MAFC);
            }

          // set up pointer to offspring cohort
          ofs[i] = pop[i] + CohortNo[i];
          (void)memset((DEF_TYPE *)(ofs[i][0]), 0, BpointNo[i]*COHORT_SIZE*sizeof(double));
          mem_req = (CohortNo[i] + BpointNo[i])*I_CONST_DIM;
          if (!(mem_req < IDMemAllocated[i]))
            {
              IDMemAllocated[i] = MemBlocks(mem_req);
              popIDcard[i]      = (popID)Myalloc((void *)popIDcard[i], (size_t)IDMemAllocated[i], sizeof(double));
              if (!(popIDcard[i])) ErrorAbort(MAFI);
            }

          ofsIDcard[i] = popIDcard[i] + CohortNo[i];
          (void)memset((DEF_TYPE *)(ofsIDcard[i][0]), 0, BpointNo[i]*I_CONST_DIM*sizeof(double));
        }
      else
        {
          ofs[i]       = NULL;
          ofsIDcard[i] = NULL;
        }
    }

  for (i = 0, bp_req = 0; i < POPULATION_NR; i++) bp_req += BpointNo[i];
  bp_req *= COHORT_SIZE;

  if ((bp_req > 0) && !(bp_req < BPAllocated))
    {
      BPAllocated = MemBlocks(bp_req);
      BPData      = (population)Myalloc((void *)BPData, (size_t)BPAllocated, sizeof(double));
      if (!(BPData)) ErrorAbort(MAFB);
    }

  for (i = 0, base = 0; i < POPULATION_NR; i++)
    {
      if (BpointNo[i] > 0)
        bpoints[i] = BPData + base;
      else
        bpoints[i] = NULL;
      bpoint_no[i] = BpointNo[i];
      base += BpointNo[i];
    }

  if (BPData && BPAllocated) (void)memset((DEF_TYPE *)(BPData), 0, BPAllocated*sizeof(double));

  return;
}


/*==================================================================================================================================*/

void		TransBcohorts()

  /* 
   * TransBcohorts - Routine transforms the boundary cohorts into internal 
   *		     cohorts and clears the current boundary points.
   */


{
  register int        i;
  register int        j, k;
  register cohort_pnt p, bp;

  for (i = 0; i < POPULATION_NR; i++)
    {
      if (!bpoints[i]) continue;
      p  = ofs[i];
      bp = bpoints[i];
      for (j = 0; j < BpointNo[i]; j++)                                             // If the new cohorts are not empty, transform i-state
        {
          if (p[j][number] > 0)
            {
              for (k = 1; k < COHORT_SIZE; k++) p[j][k] = bp[j][k] + p[j][k]/p[j][number];
            }
          else                                                                      // Reset cohort to empty
            (void)memset((DEF_TYPE *)p[j], 0, COHORT_SIZE*sizeof(double));
        }
      bpoints[i] = NULL;                                                            // Set the pointers back to NULL pointers
    }

  return;
}


/*==================================================================================================================================*/

static void	InsertBcohorts()

  /* 
   * InsertBcohorts - Routine inserts the boundary cohorts into the 
   *		      populations, ignoring empty ones.
   */

{
  register int        i;
  register int        j;
  int                 cbc;

  for (i = 0, cbc = 0; i < POPULATION_NR; i++) cbc = imax(cbc, BpointNo[i]);
  if (!cbc) return;

  if (cbc > MaxBcohorts)
    {
      MaxBcohorts = cbc;
      BcohortsPopPnt = (cohort_pnt)Myalloc((void *)BcohortsPopPnt,  (size_t)(MaxBcohorts*COHORT_SIZE), sizeof(double));
      if (!BcohortsPopPnt) ErrorAbort(MAFT);

      BcohortsIDPnt  = (cohortID_pnt)Myalloc((void *)BcohortsIDPnt, (size_t)(MaxBcohorts*I_CONST_DIM), sizeof(double));
      if (!BcohortsIDPnt) ErrorAbort(MAFT);
    }

  for (i = 0; i < POPULATION_NR; i++)
    {
      (void)memcpy((DEF_TYPE *)BcohortsPopPnt, (DEF_TYPE *)(ofs[i]), BpointNo[i]*COHORT_SIZE*sizeof(double));
      (void)memset((DEF_TYPE *)(ofs[i]), 0, BpointNo[i]*COHORT_SIZE*sizeof(double));
      (void)memcpy((DEF_TYPE *)BcohortsIDPnt, (DEF_TYPE *)(ofsIDcard[i]), BpointNo[i]*I_CONST_DIM*sizeof(double));
      (void)memset((DEF_TYPE *)(ofsIDcard[i]), 0, BpointNo[i]*I_CONST_DIM*sizeof(double));

      // If the new cohorts are not empty, insert in population
      for (j = 0; j < BpointNo[i]; j++)
        {
          if (BcohortsPopPnt[j][number] > MIN_ACCURACY) CohortNo[i] = InsCohort(BcohortsPopPnt[j], BcohortsIDPnt[j], i);
        }

      // Set the pointers back to NULL pointers
      ofs[i]      = NULL;
      BpointNo[i] = 0;
    }

  return;
}


/*==================================================================================================================================*/

void	SievePop()

  /* 
   * SievePop - Routine that scans all the cohorts, removing the ones that 
   *		are below the minimum size and conjugating cohorts that 
   *		have become too similar.
   */

{
  register int          i;
  register int          j, ind1, ind2;
  register cohort_pnt   p;
  register cohortID_pnt pid;

  for (i = 0; i < POPULATION_NR; i++)
    {
      // Test the size of the cohorts
      j   = 0;
      p   = pop[i];
      pid = popIDcard[i];
      while ((j < CohortNo[i]) && (p[j][0] > MIN_ACCURACY)) j++;
      ind1 = j;
      while (j < CohortNo[i])
        {
          while ((j < CohortNo[i]) && !(p[j][0] > MIN_ACCURACY)) j++;
          ind2 = j;
          while ((j < CohortNo[i]) && (p[j][0] > MIN_ACCURACY)) j++;
          if (j > ind2)
            {
              (void)memmove((DEF_TYPE *)p[ind1], (DEF_TYPE *)p[ind2], (j - ind2)*COHORT_SIZE*sizeof(double));
              (void)memmove((DEF_TYPE *)pid[ind1], (DEF_TYPE *)pid[ind2], (j - ind2)*I_CONST_DIM*sizeof(double));
              ind1 += (j - ind2);
            }
        }
      CohortNo[i] = ind1;
    }

  for (i = 0; i < POPULATION_NR; i++) cohort_no[i] = CohortNo[i];

  return;
}


/*==================================================================================================================================*/

void	    CohortCycle(double next)

  /* 
   * CohortCycle - Routine loops through one cycle of cohort creation and output
   *               generation. First new boundary cohorts are created, then
   *               integration of the continuous dynamics takes place, after
   *               which instantaneous dynamics can occur at the end of a cohort
   *               cycle. Finally the boundary cohorts are merged with  the
   *               existing population and output is produced if necessary.
   */

{
  register int i;
  int          all_zero = 1;
  char         pext[80];

  for (i = 0; i < POPULATION_NR; i++) /* If all populations are   */
    {                                 /* extinct exit	            */
      if ((!CohortNo[i]) && (!pop_extinct[i]))
        {
          sprintf(pext, "Population %d is extinct at time %.2f", i, env[0]);
          Warning(pext);
          pop_extinct[i] = 1;
        }
      else if (CohortNo[i])
        pop_extinct[i] = 0;
      all_zero *= pop_extinct[i];
    }
  if (all_zero) ErrorExit(0, AEXT);

  for (i = 0; i < POPULATION_NR; i++) cohort_no[i] = CohortNo[i];

  if (BifurcationRun)                                                               // Set the current bifurcation parameter value
    parameter[BifParIndex] = (BifParBase + floor((env[0] + BIFTINY)/BifPeriod)*BifParStep);

  // User specifies no. of birth points at the start of each cohort cycle
  SetBpointNo(env, pop, BpointNo);

  // Create boundary cohorts and bpoints
  CreateBcohorts();

  // User defines bpoints
  SetBpoints(env, pop, bpoints);

  if (BifurcationRun) measureBifstats(env, pop, 1, 0);

  if (step_size <= Odesolve_Min_Step) step_size = cohort_limit;
  minss = maxss = step_size;

  // Do as many integration steps as possible with adaptable stepsize and end with an optional rest step
  PrepareCycle();
  ForcedCohortEnd = cohort_end = 0;

  while ((!((next - env[0]) < Odesolve_Min_Step)) && (!cohort_end))
    {
      (void)IntegrationStep(step_size, (next - (env[0])), 0);

      initState    = NULL;
      currentState = NULL;
      for (i = 0; i < MAXDERS; i++) currentDers[i] = NULL;

      minss = min(minss, step_size);
      maxss = max(maxss, step_size);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      if (checkInterrupt()) return;
#endif
    }
  ForcedCohortEnd = cohort_end;

  if (EBTDEBUG(2))
    {
      (void)fprintf(dbgfile, "Cohort end: T = %15.8f   min. dt: %12.7E  max. dt: %12.7E\n", env[0], minss, maxss);
      (void)fflush(dbgfile);
    }

  // Transform boundary	cohorts
  TransBcohorts();

  // Instantaneous dynamics at the end of cycle
  InstantDynamics(env, pop, ofs);

  // Insert boundary cohorts
  InsertBcohorts();

  // Delete all cohorts that are too small
  SievePop();

  if (BifurcationRun) SetBifOutputTimes(env);

  outputDefined = 0;
  if (env[0] >= (next_output - identical_zero))
    { 
      // Increment next time and produce output if required at this time
      next_output = (floor((env[0] + identical_zero)/delt_out))*delt_out;
      next_output += delt_out;
      FileOut();
      outputDefined = 1;
    }

  if (BifurcationRun) outputMeasureBifstats(env);

  if (state_out > 0.0)
    {
      if (env[0] >= (next_state_output - identical_zero))
        {
          // Increment next time and produce complete state output if required at this time
          next_state_output = (floor((env[0] + identical_zero)/state_out))*state_out;
          next_state_output += state_out;
          FileState();
        }
    }

  if (fabs(next_cohort_end - env[0]) < Odesolve_Min_Step)
    next_cohort_end += cohort_limit;
  else if ((next_cohort_end - env[0]) < Odesolve_Min_Step)
    next_cohort_end = env[0] + cohort_limit;

  return;
}


/*==================================================================================================================================*/

void ResetCohorts(void)
{
  int             pp;

  for (pp = 0; pp < POPULATION_NR; pp++)
    {
      if (pop[pp])        free(pop[pp]); 
      if (popIDcard[pp])  free(popIDcard[pp]);
      pop[pp]              = NULL;
      popIDcard[pp]        = NULL;
      DataMemAllocated[pp] = 0;
      IDMemAllocated[pp]   = 0;
      CohortNo[pp]         = 0;
      BpointNo[pp]         = 0;
    }
  if (BPData) free(BPData);
  BPData = NULL;
  BPAllocated = 0L;

  if (BcohortsPopPnt) free(BcohortsPopPnt);
  BcohortsPopPnt = NULL;
  if (BcohortsIDPnt)  free(BcohortsIDPnt);
  BcohortsIDPnt  = NULL;
  MaxBcohorts = -1;

  return;
}


/*==================================================================================================================================*/

