/*
  Salmon.h -  Header file specifying the elementary life-history functions of
              the model analyzed in:

                P. Catalina Chaparro-Pedraza & Andre M. de Roos.
                Density-dependent effects of mortality on the optimal body 
                size to shift habitat: Why smaller is better despite increased 
                mortality risk. 
                Evolution 74-5: 831-841 (2020)

              The model includes a basic resource and a size-structured 
              consumer population that switches habitat from a nursery
              habitat to a growth habitat.               
              In addition, a size-selective predator forages on consumers
              in the growth habitat. Vulnerability of consumers to predation
              in the growth habitat scales with L^(-D).
 
  Layout for the i-state variables:

    istate[0][0]  : Age
    istate[0][1]  : Length

  Layout for the environment variables:

    E[0]          : Resource
    E[1]          : Predators
 
  Layout for the interaction (and output) variables:

    I[0][0]       : Total resource ingestion by the consumer population
    I[0][1]       : Consumer biomass availability for predators
 
    I[0][2]       : Total biomass of small juveniles in the nursery habitat
    I[0][3]       : Total biomass of juveniles in the growth habitat
    I[0][4]       : Total adult biomass
 
  Last modification: AMdR - Jul 16, 2020
*/

/*
 *==================================================================================
 *  SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
 *==================================================================================
 */
// Dimension settings: Required
#define POPULATION_NR             1
#define STAGES                    3
#define I_STATE_DIM               2
#define ENVIRON_DIM               2
#define INTERACT_DIM              5
#define PARAMETER_NR              15

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL              1.0E-9      // Minimum individual survival
#define MAX_AGE                   200000      // Absolute maximum individual age

#define DYTOL                     1.0E-7      // Variable tolerance
#define RHSTOL                    1.0E-6      // Function tolerance

#define ALLOWNEGATIVE             0           // Negative solution values allowed?
#define COHORT_NR                 200         // Number of cohorts in state output


// Descriptive names of parameters in parameter array (at least two are required)
char *parameternames[PARAMETER_NR] = {"Rho", "Xmax",   "K", "Imax",  "Bmax", 
                                      "L0",  "Ls",    "Lm", "Linf",  "Xi", 
                                      "Mu1", "Mu2",   "D", "Alpha", "Mup"};

// Default values of all parameters
double parameter[PARAMETER_NR] = {0.01, 0.5, 1.0, 0.0025, 0.002725, 
                                  2.0, 20.0, 30.0, 115.0, 0.00051, 
                                  0.002, 0.006, 0.75, 0.001, 0.006};

// Aliases definitions for all istate variables
#define AGE       istate[0][0]
#define LENGTH    istate[0][1]

// Aliases definitions for all environment variables
#define X         E[0]
#define P         E[1]

// Aliases definitions for all parameters

// Resource in the nursery habitat
#define RHO       parameter[ 0]       // Resource growth rate
#define XMAX      parameter[ 1]       // Maximum resource density

// Population with habitat shift
#define K         parameter[ 2]       // Half saturation resource density
#define IMAX      parameter[ 3]       // Maximum ingestion proportionality constant

#define BMAX      parameter[ 4]       // Maximum fecundity proportionality constant

#define L0        parameter[ 5]       // Body size of a newborn
#define LS        parameter[ 6]       // Body size at the habitat shift
#define LM        parameter[ 7]       // Body size at maturation
#define LINF      parameter[ 8]       // Maximum body size at maximum feeding rate

#define XI        parameter[ 9]       // Von Bertalanffy growth rate parameter

#define MU1       parameter[10]       // Background mortality rate in habitat 1 
#define MU2       parameter[11]       // Background mortality rate in habitat 2

#define D         parameter[12]       // Exponent of size-dependent mortality

#define PHI       parameter[13]       // Scaled predator attack rate
#define MUP       parameter[14]       // Scaled predator mortality rate

/*
 *==================================================================================
 *  SECTION 2: DEFINITION OF THE INDIVIDUAL LIFE HISTORY
 *==================================================================================
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
  LENGTH = L0;

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

void IntervalLimit(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                   double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                   double limit[POPULATION_NR])
{
  switch (lifestage[0])
    {
      case 0:
        limit[0] = LENGTH - LS;
        break;
      case 1:
        limit[0] = LENGTH - LM;
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

void Development(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                 double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
                 double development[POPULATION_NR][I_STATE_DIM])
{
  development[0][0] = 1.0;
  if (lifestage[0] == 0)
    development[0][1] = XI * (LINF * X / (K + X) - LENGTH);
  else
    development[0][1] = XI * (LINF - LENGTH);

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

void DiscreteChanges(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
                     double *birthstate[POPULATION_NR], int BirthStateNr,
                     double E[])
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

void Fecundity(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], 
               double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double *fecundity[POPULATION_NR])
{
  fecundity[0][0] = 0.0;

  if (lifestage[0] > 1) fecundity[0][0] = BMAX * LENGTH * LENGTH;

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

void Mortality(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], 
               double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
               double mortality[POPULATION_NR])
{
  if (lifestage[0] == 0)
    mortality[0] = MU1;
  else
    mortality[0] = MU2 + PHI * P * pow(LENGTH, -D);

  return;
}


/*
 *==================================================================================
 *  SECTION 3: FEEDBACK ON THE ENVIRONMENT
 *==================================================================================
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

void Impact(int lifestage[POPULATION_NR], double *istate[POPULATION_NR],
            double *birthstate[POPULATION_NR], int BirthStateNr, double E[],
            double impact[POPULATION_NR][INTERACT_DIM])
{
  impact[0][0] = IMAX * X / (K + X) * LENGTH * LENGTH;

  switch (lifestage[0])
    {
      case 0:
        impact[0][1] = 0;
        impact[0][2] = LENGTH*LENGTH*LENGTH;
        impact[0][3] = 0;
        impact[0][4] = 0;
        break;
      case 1:
        impact[0][1] = PHI * pow(LENGTH, -D) * LENGTH * LENGTH * LENGTH;
        impact[0][2] = 0;
        impact[0][3] = LENGTH*LENGTH*LENGTH;
        impact[0][4] = 0;
        break;
      case 2:
        impact[0][1] = PHI * pow(LENGTH, -D) * LENGTH * LENGTH * LENGTH;
        impact[0][2] = 0;
        impact[0][3] = 0;
        impact[0][4] = LENGTH*LENGTH*LENGTH;
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

const int EnvironmentType[ENVIRON_DIM] = {GENERALODE, PERCAPITARATE};

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM], 
             double condition[ENVIRON_DIM])
{
  condition[0] = RHO * (XMAX - X) - I[0][0];
  condition[1] = I[0][1] - MUP;

  return;
}


/*================================================================================*/
