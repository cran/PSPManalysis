/*
  Indet_growth.h -  Header file specifying the elementary life-history functions of
                    a structured consumer-unstructured resource model, in which the
                    consumer follows a net-production DEB model with arbitrary, allometric
                    scaling functions of ingestion and maintenance as a function of body size.

                    The variable that is solved for is the biomass of the resource.

  Layout for the i-state variables:

    istate[i][0]  : Age(i)
    istate[i][1]  : Size(i)

  Layout for the environment variables:

    E[0]          : Resource

  Layout for the interaction (and output) variables:

    I[0][0]       : Total ingestion

    I[0][1]       : Total biomass of juveniles
    I[0][2]       : Total biomass of growing adults
    I[0][3]       : Total biomass of non-growing adults

  Last modification: AMdR - Jun 07, 2017
*/

/*
 *====================================================================================================================================
 *  SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
 *====================================================================================================================================
 */
// Dimension settings: Required
#define POPULATION_NR             1
#define STAGES                    2
#define I_STATE_DIM               2
#define ENVIRON_DIM               1
#define INTERACT_DIM              3
#define PARAMETER_NR              11

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL              1.0E-9                                            // Survival at which individual is considered dead
#define MAX_AGE                   100000                                            // Give some absolute maximum for individual age

#define DYTOL                     1.0E-7                                            // Variable tolerance
#define RHSTOL                    1.0E-6                                            // Function tolerance

#define ODESOLVE_REL_ERR          1.0E-6


// Descriptive names of parameters in parameter array (at least two parameters are required)
char *parameternames[PARAMETER_NR] = {"Delta", "Rmax", "Sb", "Sj", "Sm", "Imax", "q", "Sigma", "T", "p", "Mu"};

// Default values of all parameters
double parameter[PARAMETER_NR] = {0.1, 2.0, 0.05, 1.0, 2.0, 1.0, 1.0, 0.5, 0.1, 1.0, 0.01};


// Aliases definitions for all istate variables
#define AGE(n)                    istate[n][0]
#define SIZE(n)                   istate[n][1]

// Aliases definitions for all environment variables
#define R                         E[0]

// Aliases definitions for all parameters
#define DELTA                     parameter[0]                                      // Default: 0.1
#define RMAX                      parameter[1]                                      // Default: 2

#define SB                        parameter[2]                                      // Default: 0.05  Mass of newborn individual
#define SJ                        parameter[3]                                      // Default: 1.0   Mass at maturation
#define SM                        parameter[4]                                      // Default: 2.0   Maximum mass

#define IMAX                      parameter[5]                                      // Default: 1     Attack rate constant
#define Q                         parameter[6]                                      // Default: 1     Attack rate exponent

#define SIGMA                     parameter[7]                                      // Default: 0.5   Food conversion efficiency

#define T                         parameter[8]                                      // Default: 0.1   Maintenance rate constant
#define P                         parameter[9]                                      // Default: 1     Maintenance rate exponent

#define MU                        parameter[10]                                     // Default: 0.01  Mortality rate


/*
 *====================================================================================================================================
 *  SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
 *====================================================================================================================================
 */

/*
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 */

#define BIRTHSTATES               5
#define BIRTHSPREAD               0.02

void SetBirthStates(int BirthStates[POPULATION_NR], double E[])
{
  BirthStates[0] = BIRTHSTATES;

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

  AGE(0)  = 0.0;
  SIZE(0) = SB + (((double)BirthStateNr)/((double)(BIRTHSTATES - 1)) - 0.5)*BIRTHSPREAD;

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
      case 0:
        limit[0] = SIZE(0) - SJ;
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
  double Ingest, netproduction;
  double L, kappa;

  development[0][0] = 1.0;

  Ingest        = IMAX*pow(SIZE(0), Q)*R;
  netproduction = SIGMA*Ingest - T*pow(SIZE(0), P);

  if (lifestage[0] == 0)
    {
      development[0][1] = netproduction;
    }
  else
    {
      L                 = (SIZE(0) - SJ)/(SM - SJ);
      kappa             = 1 - 3*L*L + 2*L*L*L;
      development[0][1] = kappa*netproduction;
    }

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
  double Ingest, netproduction;
  double L, kappa, fec;

  Ingest        = IMAX*pow(SIZE(0), Q)*R;
  netproduction = SIGMA*Ingest - T*pow(SIZE(0), P);
  L             = (SIZE(0) - SJ)/(SM - SJ);
  kappa         = 1 - 3*L*L + 2*L*L*L;

  if (lifestage[0] == 0)
    {
      fecundity[0][0] = 0.0;
      fecundity[0][1] = 0.0;
      fecundity[0][2] = 0.0;
      fecundity[0][3] = 0.0;
      fecundity[0][4] = 0.0;
    }
  else
    {
      fec             = (1 - kappa)*netproduction/SB;
      fecundity[0][0] = 0.1*fec;
      fecundity[0][1] = 0.2*fec;
      fecundity[0][2] = 0.5*fec;
      fecundity[0][3] = 0.2*fec;
      fecundity[0][4] = 0.1*fec;
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
  mortality[0] = MU;

  return;
}


/*
 *====================================================================================================================================
 *  SECTION 3: FEEDBACK ON THE ENVIRONMENT
 *====================================================================================================================================
 */

/*
 * For all the integrals (measures) that occur in interactions of the
 * structured populations with their environments and for all the integrals
 * that should be computed for output purposes (e.g. total juvenile or adult
 * biomass), specify appropriate weighing function dependent on the i-state
 * variables, the individual's state at birth, the environment variables and
 * the current life stage of the individuals. These weighing functions should
 * be specified for all structured populations in the problem. The number of
 * weighing functions is the same for all of them.
 *
 * Notice that the first index of the variables 'istate[][]' and 'impact[][]'
 * refers to the number of the structured population, the second index of the
 * variable 'istate[][]' refers to the number of the individual state variable,
 * while the second index of the variable 'impact[][]' refers to the number of
 * the interaction variable. The interpretation of these second indices is up
 * to the user.
 */

void Impact(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
            double impact[POPULATION_NR][INTERACT_DIM])
{
  impact[0][0] = IMAX*pow(SIZE(0), Q)*R;
  switch (lifestage[0])
    {
      case 0:
        impact[0][1] = SIZE(0);
        impact[0][2] = 0;
        break;
      case 1:
        impact[0][1] = 0;
        impact[0][2] = SIZE(0);
        break;
    }

  return;
}


/*
 * Specify the type of each of the environment variables by setting
 * the entries in EnvironmentType[ENVIRON_DIM] to PERCAPITARATE, GENERALODE
 * or POPULATIONINTEGRAL based on the classification below:
 *
 * Set an entry to PERCAPITARATE if the dynamics of E[j] follow an ODE and 0
 * is a possible equilibrium state of E[j]. The ODE is then of the form
 * dE[j]/dt = P(E,I)*E[j], with P(E,I) the per capita growth rate of E[j].
 * Specify the equilibrium condition as condition[j] = P(E,I), do not include
 * the multiplication with E[j] to allow for detecting and continuing the
 * transcritical bifurcation between the trivial and non-trivial equilibrium.
 *
 * Set an entry to GENERALODE if the dynamics of E[j] follow an ODE and 0 is
 * NOT an equilibrium state of E. The ODE then has a form dE[j]/dt = G(E,I).
 * Specify the equilibrium condition as condition[j] = G(E,I).
 *
 * Set an entry to POPULATIONINTEGRAL if E[j] is a (weighted) integral of the
 * population distribution, representing for example the total population
 * biomass. E[j] then can be expressed as E[j] = I[p][i]. Specify the
 * equilibrium condition in this case as condition[j] = I[p][i].
 *
 * Notice that the first index of the variable 'I[][]' refers to the
 * number of the structured population, the second index refers to the
 * number of the interaction variable. The interpretation of the latter
 * is up to the user. Also notice that the variable 'condition[j]' should
 * specify the equilibrium condition of environment variable 'E[j]'.
 */

const int EnvironmentType[ENVIRON_DIM] = {GENERALODE};

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM], double condition[ENVIRON_DIM])
{
  condition[0] = DELTA*(RMAX - R) - I[0][0];

  return;
}


/*==================================================================================================================================*/
