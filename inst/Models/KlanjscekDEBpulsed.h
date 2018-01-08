/*
  KlanjscekDEBpulsed.h -  Header file specifying the elementary life-history functions
                          of a structured consumer-unstructured resource model, in which the
                          consumer follows the Kooijman DEB model as specified in:

                          Klanjscek, T., Caswell, H., Neubert, M.G. & Nisbet, R.M. (2006).
                          Integrating dynamic energy budgets into matrix population models.
                          Ecol. Modell., 196, 407-420.

  Layout for the i-state variables:

    istate[0][0]  : Age
    istate[0][1]  : Volume
    istate[0][2]  : Q
    istate[0][3]  : H
    istate[0][4]  : E: Eggs in the egg buffer

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
#define I_STATE_DIM               5                                                 // See below
#define PARAMETER_NR              10

// The following definition will force the program to consider reproduction pulses
#define REPRODUCTION_INTERVAL     1.0

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL              1.0E-9                                            // Survival at which individual is considered dead
#define MAX_AGE                   100000                                            // Give some absolute maximum for individual age

#define DYTOL                     1.0E-7                                            // Variable tolerance
#define RHSTOL                    1.0E-8                                            // Function tolerance


// Descriptive names of parameters in parameter array (at least two parameters are required)
char *parameternames[PARAMETER_NR] = {"Food", "Kappa", "Kappa_R", "Nu", "m", "g", "Vb", "Vp", "[Em]", "ha"};

// Default values of all parameters
double parameter[PARAMETER_NR] = {1.0, 0.8, 0.001, 0.075, 0.583, 1.286, 1.0E-9, 1.73E-6, 0.7, 0.15};


// Aliases definitions for all istate variables
#define AGE                       istate[0][0]
#define VOLUME                    istate[0][1]
#define Q                         istate[0][2]
#define H                         istate[0][3]
#define EGGS                      istate[0][4]

// Aliases definitions for all parameters
#define FOOD                      parameter[0]                                      // Default: 0.5

// DEB parameters
#define KAPPA                     parameter[1]                                      // Default: 0.8
#define KAPPA_R                   parameter[2]                                      // Default: 0.001
#define NU                        parameter[3]                                      // Default: 0.075
#define M                         parameter[4]                                      // Default: 0.583
#define G                         parameter[5]                                      // Default: 1.286
#define VB                        parameter[6]                                      // Default: 1.0E-9
#define VP                        parameter[7]                                      // Default: 1.73E-6

// Aging parameters
#define EM                        parameter[8]                                      // Default: 0.7
#define HA                        parameter[9]                                      // Default: 0.15


/*
 *====================================================================================================================================
 *  SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
 *====================================================================================================================================
 */

/*
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 */

void SetBirthStates(int BirthStates[POPULATION_NR], double E[])
{
  BirthStates[0] = 1;

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
  AGE    = 0.0;
  VOLUME = VB;                                                                      // Tiny length at birth to start with
  Q      = 0.0;
  H      = 0.0;

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
  switch (lifestage[0])
    {
      // Juvenile stage ends when V reaches Vp
      case 0:
        limit[0] = VOLUME/VP - 1.0;
        break;
    }

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
  double dVda, dQda, dHda;
  double Er;

  //  Assume growth always occurs
  dVda = max((FOOD*NU*pow(VOLUME, 2.0/3.0) - M*G*VOLUME)/(FOOD + G), 0);
  dQda = G*EM*(dVda + M*VOLUME);
  dHda = HA*Q/VOLUME;

  development[0][0] = 1.0;
  development[0][1] = dVda;                                                         // dV/da
  development[0][2] = dQda;                                                         // dQ/da
  development[0][3] = dHda;                                                         // dH/da

  if (lifestage[0] == 1)                                                            // Only for adults
    {
      Er                = (1 - KAPPA)*EM*FOOD*G*(NU*pow(VOLUME, 2.0/3.0) + M*VOLUME)/(FOOD + G);
      development[0][4] = max(Er - (1 - KAPPA)*EM*M*G*VP, 0);
      development[0][4] /= EM*(KAPPA*G + FOOD)*VB/KAPPA_R;
    }
  else
    development[0][4] = 0;

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
  if (lifestage[0] == 1)                                                            // Only for adults
    {
      fecundity[0][0] = EGGS;
    }
  else
    fecundity[0][0] = 0;

  EGGS = 0.0;                                                                       // Empty the egg buffer

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
  mortality[0] = H;

  return;
}


/*==================================================================================================================================*/
