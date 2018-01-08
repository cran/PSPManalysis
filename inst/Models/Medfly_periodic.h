/*
  Medfly_periodic.h - Header file specifying the elementary life-history functions of
                      the Medfly example model with periodic mortality as described in

                      A.M. de Roos, 2008. Demographic analysis of continuous-time
                      life-history models. Ecol. Lett. 11(1): 1-15

  Layout for the i-state variables:

    istate[0][0]  : Age
    istate[0][1]  : Phase

  Last modification: AMdR - Jun 07, 2017
*/


/*
 *====================================================================================================================================
 *  SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
 *====================================================================================================================================
 */
// Dimension settings: Required
#define POPULATION_NR             1                                                 // Structured consumer population
#define STAGES                    2                                                 // Juvenile & adult
#define I_STATE_DIM               2                                                 // See below
#define PARAMETER_NR              8

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL              1.0E-9                                            // Survival at which individual is considered dead
#define MAX_AGE                   1000                                              // Give some absolute maximum for individual age

#define DYTOL                     1.0E-5                                            // Variable tolerance
#define RHSTOL                    1.0E-6                                            // Function tolerance

#define ODESOLVE_ABS_ERR          1.0E-7
#define ODESOLVE_REL_ERR          1.0E-8


// Descriptive names of parameters in parameter array (at least two parameters are required)
char *parameternames[PARAMETER_NR] = {"Beta0", "Beta1", "AJ", "Mu0", "Mu1", "Chi0", "Chi1", "Period"};

// Default values of all parameters
double parameter[PARAMETER_NR] = {47.0, 0.04, 11.0, 0.00095, 0.0581, 20.0, 2.0, 11.6};


// Aliases definitions for all istate variables
#define AGE                       istate[0][0]
#define PHASE                     istate[0][1]

// Aliases definitions for all parameters
#define BETA0                     parameter[0]                                      // Default: 47.0
#define BETA1                     parameter[1]                                      // Default:  0.04
#define AJ                        parameter[2]                                      // Default: 11.0
#define MU0                       parameter[3]                                      // Default: 0.00095
#define MU1                       parameter[4]                                      // Default: 0.0581

#define CHI0                      parameter[5]                                      // Default: 20
#define CHI1                      parameter[6]                                      // Default: 2
#define PERIOD                    parameter[7]                                      // Default: 40


/*
 *====================================================================================================================================
 *  SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
 *====================================================================================================================================
 */

/*
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 */

int States;

void SetBirthStates(int BirthStates[POPULATION_NR], double E[])
{
  States = (int)floor((PERIOD + 1.0E-6)/0.2);
  BirthStates[0] = States;

  // Adjust the variable Odesolve_Fixed_Step, which represents the fixed step size in time.
  // As a result, the time integrator will stop exactly on all the intermediate nodal
  // time points that cover the period of the oscillation in mortality
  Odesolve_Fixed_Step = (PERIOD/((double)States));

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

void StateAtBirth(double *istate[POPULATION_NR], int BirthStateNr, double E[])
{
  AGE   = 0.0;
  PHASE = BirthStateNr*PERIOD/((double)States);

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

void IntervalLimit(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                   double limit[POPULATION_NR])
{
  if (lifestage[0] == 0) limit[0] = AGE - AJ;

  return;
}


/*
 * Specify the development of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * Notice that the first index of the variables 'istate[][]' and 'development[][]'
 * refers to the number of the structured population, the second index refers
 * to the number of the individual state variable. The interpretation of the
 * latter is up to the user.
 */

void Development(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                 double development[POPULATION_NR][I_STATE_DIM])
{
  development[0][0] = 1.0;
  development[0][1] = 0.0;

  return;
}


/*
 * Specify the possible discrete changes (jumps) in the individual state
 * variables when ENTERING the stage specified by 'lifestage[]'.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 */

void DiscreteChanges(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[])
{
  return;
}


/*
 * Specify the fecundity of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * The number of offspring produced has to be specified for every possible
 * state at birth in the variable 'fecundity[][]'. The first index of this
 * variable refers to the number of the structured population, the second
 * index refers to the number of the birth state.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 */

void Fecundity(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double *fecundity[POPULATION_NR])
{
  int    i, index;
  double tau, dindex, frac;

  for (i = 0; i < States; i++) fecundity[0][i] = 0;
  if (lifestage[0] == 1)                                                            // Only for adults
    {
      // tau is the current timing in the mortality cycle given the age and the phase at birth
      // 0 <= tau < PERIOD. If tau equals 0, offspring is produced at phase 0 (i.e. index 0 in
      // the array of birth points).
      tau = fmod(AGE + birthstate[0][1] + 1.0E-9, PERIOD);

      // dindex converts tau to the appropriate index in double format between 0...States-1
      // Linear subdivision over period
      dindex = (tau/PERIOD)*States;

      // dindex is between index and index+1. If dindex = index all offspring is assigned to
      // birth point index, if dindex = index+1 zero offspring is assigned to birth point index
      index                              = (int)floor(dindex + 1.0E-9);
      frac                               = (dindex - index);
      fecundity[0][index % States]       = (1.0 - frac)*BETA0*exp(-BETA1*(AGE - AJ));
      fecundity[0][(index + 1) % States] = frac*BETA0*exp(-BETA1*(AGE - AJ));
    }

  return;
}


/*
 * Specify the mortality of individuals as a function of the i-state
 * variables and the individual's state at birth for all populations in every
 * life stage.
 *
 * Notice that the first index of the variable 'istate[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the individual state variable. The interpretation of the latter
 * is up to the user.
 */

void Mortality(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double mortality[POPULATION_NR])
{
  double tau;

  tau = fmod(AGE + birthstate[0][1], PERIOD);
  if (lifestage[0] == 0)                                                            // Only for juveniles
    {
      mortality[0] = MU0*exp(MU1*AGE) + CHI0*exp(-CHI1*tau);
    }
  else
    mortality[0] = MU0*exp(MU1*AGE);

  return;
}


/*==================================================================================================================================*/
