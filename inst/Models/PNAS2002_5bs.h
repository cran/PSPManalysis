/*
  PNAS2002_5bs.h -  Header file specifying the elementary life-history functions of
                    the tri-trophic model analyzed in:

                    A.M. de Roos & L. Persson, 2002. Size-dependent life-history
                    traits promote catastrophic collapses of top predators.
                    Proc. Natl. Acad. Sciences 99(20): 12907-12912.

                    The model includes a basic resource, a size-structured consumer
                    and an unstructured predator population. The difference with the
                    original model is that there are 5 possible states at birth

  Layout for the i-state variables:

    istate[0][0]  : Age
    istate[0][1]  : Length

  Layout for the environment variables:

    E[0]          : Resource
    E[1]          : Predators
    E[2]          : Vulnerable consumer biomass

  Layout for the interaction (and output) variables:

    I[0][0]       : Total ingestion

    I[0][1]       : Total biomass of small juveniles
    I[0][2]       : Total biomass of non-vulnerable juveniles
    I[0][3]       : Total adult biomass

  Last modification: AMdR - Apr 18, 2018
*/

/*
 *====================================================================================================================================
 *  SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
 *====================================================================================================================================
 */
// Dimension settings: Required
#define POPULATION_NR             1
#define STAGES                    3
#define I_STATE_DIM               2
#define ENVIRON_DIM               3
#define INTERACT_DIM              4
#define PARAMETER_NR              16

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL              1.0E-9                                            // Survival at which individual is considered dead
#define MAX_AGE                   100000                                            // Give some absolute maximum for individual age

#define DYTOL                     1.0E-7                                            // Variable tolerance
#define RHSTOL                    1.0E-6                                            // Function tolerance

#define ALLOWNEGATIVE             0                                                 // Negative solution values allowed?
#define COHORT_NR                 100                                               // Number of cohorts in state output

#define ODESOLVE_REL_ERR          1.0E-7

#define FULLSTATEOUTPUT           2

// Descriptive names of parameters in parameter array (at least two parameters are required)
char *parameternames[PARAMETER_NR] = {"Rho", "Rmax",  "Lb", "Lv",  "Lj", "Lm", "Omega",   "Imax",
                                      "Rh",  "Gamma", "Rm", "Mub", "A",  "Th", "Epsilon", "Delta"};

// Default values of all parameters
double parameter[PARAMETER_NR] = {0.1, 3.0E-4, 7.0, 27.0, 110.0, 300.0, 9.0E-6, 1.0E-4, 1.5E-5, 0.006, 0.003, 0.01, 5000.0, 0.1, 0.5, 0.01};


// Aliases definitions for all istate variables
#define AGE                       istate[0][0]
#define LENGTH                    istate[0][1]

// Aliases definitions for all environment variables
#define R                         E[0]
#define P                         E[1]
#define B                         E[2]

// Aliases definitions for all parameters
#define RHO                       parameter[0]                                      // Default: 0.1
#define RMAX                      parameter[1]                                      // Default: 3.0E-4

#define LB                        parameter[2]                                      // Default: 7
#define LV                        parameter[3]                                      // Default: 27
#define LJ                        parameter[4]                                      // Default: 110
#define LM                        parameter[5]                                      // Default: 300

#define OMEGA                     parameter[6]                                      // Default: 9.0E-6

#define IMAX                      parameter[7]                                      // Default: 1.0E-4
#define RH                        parameter[8]                                      // Default: 1.5E-5

#define GAMMA                     parameter[9]                                      // Default: 0.006
#define RM                        parameter[10]                                     // Default: 0.003

#define MUB                       parameter[11]                                     // Default: 0.01

#define A                         parameter[12]                                     // Default: 5000.0
#define TH                        parameter[13]                                     // Default: 0.1
#define EPSILON                   parameter[14]                                     // Default: 0.5
#define DELTA                     parameter[15]                                     // Default: 0.01


/*
 *====================================================================================================================================
 *  SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
 *====================================================================================================================================
 */

/*
 * Specify the number of states at birth for the individuals in all structured
 * populations in the problem in the vector BirthStates[].
 */

#define BIRTHSTATES 5
#define BIRTHSPREAD 2.0

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
  AGE    = 0.0;
  LENGTH = LB + (((double)BirthStateNr)/((double)(BIRTHSTATES - 1)) - 0.5)*BIRTHSPREAD;

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
        limit[0] = LENGTH - LV;
        break;
      case 1:
        limit[0] = LENGTH - LJ;
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
  development[0][0] = 1.0;
  development[0][1] = GAMMA*(LM*R/(R + RH) - LENGTH);

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
  double fec;

  fecundity[0][0] = 0.0;
  fecundity[0][1] = 0.0;
  fecundity[0][2] = 0.0;
  fecundity[0][3] = 0.0;
  fecundity[0][4] = 0.0;
  if (lifestage[0] == 2)
    {
      fec             = RM*R/(R + RH)*LENGTH*LENGTH;
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
  if (lifestage[0] == 0)
    mortality[0] = MUB + A*P/(1 + A*TH*B);
  else
    mortality[0] = MUB;

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
  impact[0][0] = IMAX*R/(R + RH)*LENGTH*LENGTH;

  switch (lifestage[0])
    {
      case 0:
        impact[0][1] = OMEGA*LENGTH*LENGTH*LENGTH;
        impact[0][2] = 0;
        impact[0][3] = 0;
        break;
      case 1:
        impact[0][1] = 0;
        impact[0][2] = OMEGA*LENGTH*LENGTH*LENGTH;
        impact[0][3] = 0;
        break;
      case 2:
        impact[0][1] = 0;
        impact[0][2] = 0;
        impact[0][3] = OMEGA*LENGTH*LENGTH*LENGTH;
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

const int EnvironmentType[ENVIRON_DIM] = {GENERALODE, PERCAPITARATE, POPULATIONINTEGRAL};

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM], double condition[ENVIRON_DIM])
{
  condition[0] = RHO*(RMAX - R) - I[0][0];
  condition[1] = EPSILON*A*I[0][1]/(1 + A*TH*I[0][1]) - DELTA;
  condition[2] = I[0][1];

  return;
}


/*==================================================================================================================================*/
