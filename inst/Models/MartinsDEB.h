/*
  MartinsDEB.h -  Header file specifying the elementary life-history functions of
                  a structured consumer-unstructured resource model, in which the
                  consumer follows the variant of the Kooijman DEB model as specified in:

                  Martins et al., 2013, The American Naturalist, http://www.jstor.org/stable/10.1086/669904

                  The variable that are solved for are the biomass of the resource
                  value of E at birth.

  Layout for the i-state variables:

    istate[0][0]  : Age
    istate[0][1]  : Length
    istate[0][2]  : Ue
    istate[0][3]  : Uh
    istate[0][4]  : Q
    istate[0][5]  : H

  Layout for the environment variables:

    E[0]          : Resource
    E[1]          : Ue0

  Layout for the interaction (and output) variables:

    I[0][0]       : Total ingestion

    I[0][1]       : Total number of embryos
    I[0][2]       : Total number of juveniles
    I[0][3]       : Total adult number

    I[0][4]       : Total biomass of embryos
    I[0][5]       : Total biomass of juveniles
    I[0][6]       : Total adult biomass

    I[0][7]       : Stage-specific net-biomass balance embryos
    I[0][8]       : Stage-specific net-biomass balance juveniles
    I[0][9]       : Stage-specific net-biomass balance adults

  Last modification: AMdR - Oct 30, 2017
*/

#define EATBIRTH                  0                                                 // Value should be 0 to use f (=X/(X+K)) or 1 to use 1, as in Martin et al. (2013)

#ifndef _MSC_VER                                                                    // Microsoft Visual C 8.0 can not handle the #warning preprocessor statement
#if (EATBIRTH == 0)
#warning
#warning E at birth assumed to equal X/(X+K) following Kooijman (2009)
#warning
#elif (EATBIRTH == 1)
#warning
#warning E at birth assumed to equal 1 following Martin et al. (2013)
#warning
#endif
#endif


/*
 *====================================================================================================================================
 *  SECTION 1: PROBLEM DIMENSIONS, NUMERICAL SETTINGS AND MODEL PARAMETERS
 *====================================================================================================================================
 */
// Dimension settings: Required
#define POPULATION_NR             1                                                 // Structured consumer population
#define STAGES                    3                                                 // Embryo, juvenile & adult
#define I_STATE_DIM               6                                                 // See below
#define ENVIRON_DIM               2                                                 // Unstructured resource, value of Ue0
#define INTERACT_DIM              10
#define PARAMETER_NR              22

// Numerical settings: Optional (default values adopted otherwise)
#define MIN_SURVIVAL              1.0E-9                                            // Survival at which individual is considered dead
#define MAX_AGE                   100000                                            // Give some absolute maximum for individual age

#define DYTOL                     1.0E-7                                            // Variable tolerance
#define RHSTOL                    1.0E-6                                            // Function tolerance


// Descriptive names of parameters in parameter array (at least two parameters are required)
char *parameternames[PARAMETER_NR] = {"Delta", "Xmax",    "Kappa", "KappaR", "KMdot", "KJdot",  "Uhb",   "Uhp",   "Vdot",  "G",     "HAdotdot",
                                      "Sg",    "JXamdot", "K",     "M",      "Muc",   "Volume", "I-Eff", "M-Eff", "G-Eff", "R-Eff", "ES-Eff"};

// Default values of all parameters
double parameter[PARAMETER_NR] = {0.1,   10.0,   0.678,  0.95,  0.3314, 0.1921, 0.1108, 2.555, 18.1, 10.0, 3.04E-6,
                                  0.019, 3.8E+5, 1.5850, 0.090, 0.0,    1.0E+9, 0.0,    0.0,   0.0,  0.0,  0.0};


// Aliases definitions for all istate variables
#define AGE                       istate[0][0]
#define LENGTH                    istate[0][1]
#define UE                        istate[0][2]
#define UH                        istate[0][3]
#define Q                         istate[0][4]
#define H                         istate[0][5]

// Aliases definitions for all environment variables
#define X                         E[0]
#define UE0                       E[1]

// Aliases definitions for all parameters
// These are their problem-specific meaning (I assume semi-chemostat resource dynamics))
#define DELTA                     parameter[0]                                      // Default: 0.1
#define XMAX                      parameter[1]                                      // Default: 1.0E4
// DEB parameters
#define KAPPA                     parameter[2]                                      // Default: 0.678
#define KAPPAR                    parameter[3]                                      // Default: 0.95
#define KMDOT                     parameter[4]                                      // Default: 0.3314
#define KJDOT                     parameter[5]                                      // Default: 0.1921
#define UHB                       parameter[6]                                      // Default: 0.1108
#define UHP                       parameter[7]                                      // Default: 2.555
#define VDOT                      parameter[8]                                      // Default: 18.1
#define G                         parameter[9]                                      // Default: 10.0

// Aging parameters
#define HADOTDOT                  parameter[10]                                     // Default: 3.04E-6
#define SG                        parameter[11]                                     // Default: 0.019

// Foraging parameters
#define JXAMDOT                   parameter[12]                                     // Default: 3.8E+5
#define K                         parameter[13]                                     // Default: 1585

// Daphnia-specific parameters
#define M                         parameter[14]                                     // Default: 0.090

#define MUC                       parameter[15]                                     // Default: 0.0

#define VOLUME                    parameter[16]                                     // Scaling values for all densities (in cm^3 because of dimension of X and K)

// Toxic effect parameters
#define INTAKEEFFECT              parameter[17]                                     // Varies between 0 and 1
#define MAINTEFFECT               parameter[18]                                     // Larger than 0
#define GROWTHEFFECT              parameter[19]                                     // Larger than 0
#define REPROEFFECT               parameter[20]                                     // Larger than 0
#define EGGSURVEFFECT             parameter[21]                                     // Larger than 0


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
  LENGTH = 1.0E-4;                                                                  // Tiny length at birth to start with
  UE     = UE0;                                                                     // Ue at birth is variable to solve for
  UH     = 0.0;
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
      // Embryo stage ends when Uh reaches Uh at birth (UHB)
      case 0:
        limit[0] = UH - UHB;
        break;
      // Juvenile stage ends when Uh reaches Uhp
      case 1:
        limit[0] = UH - UHP;
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
  double KMDOT_s, KJDOT_s, G_s;
  double f, Lm, L2;
  double Sa, Sc;
  double e, rdot;
  double growthL, maturation, aging, hazardchange;

  KMDOT_s = KMDOT*(1 + MAINTEFFECT)/(1 + GROWTHEFFECT);
  KJDOT_s = KJDOT*(1 + MAINTEFFECT);
  G_s     = G*(1 + GROWTHEFFECT);

  f = X/(X + K);
  f *= (1 - INTAKEEFFECT);
  Lm = VDOT/(G_s*KMDOT_s);
  L2 = LENGTH*LENGTH;

  if (lifestage[0] == 0)
    {
      // In the embryo we always have to track Ue dynamics irrespective of whether
      // e = 1 or e = f at birth
      Sa = 0;

      // Rewrite the original equations
      //      e  = VDOT*Ue/(LENGTH*LENGTH*LENGTH);
      //      Sc = L2*G*e*(1 + L*KMDOT/VDOT)/(G + e);
      // for numerical reasons
      Sc = L2*G_s*UE*(VDOT + LENGTH*KMDOT_s)/((LENGTH*LENGTH*LENGTH)*G_s + VDOT*UE);

      //  growthL = (VDOT*Sc/(G*L2) - KMDOT*L)/3;
      growthL = (VDOT*UE*(VDOT + LENGTH*KMDOT_s)/((LENGTH*LENGTH*LENGTH)*G_s + VDOT*UE) - KMDOT_s*LENGTH)/3;

      maturation = (1 - KAPPA)*Sc - KJDOT_s*UH;
      // I assume that embryos do not experience any mortality
      aging        = 0;
      hazardchange = 0;
    }
  else
    {
      Sa = f*L2;
#if (EATBIRTH == 0)
      e = f;
#elif (EATBIRTH == 1)
      e        = VDOT*UE/(LENGTH*LENGTH*LENGTH);
#endif
      Sc = L2*G_s*e*(1 + LENGTH*KMDOT_s/VDOT)/(G_s + e);

      //  growthL = (VDOT*Sc/(G*L2) - KMDOT*L)/3;
      //  growthL = (e*(VDOT + KMDOT*L)/(G + e) - KMDOT*L)/3;
      growthL = (VDOT*e/(G_s + e) - KMDOT_s*LENGTH*G_s/(G_s + e))/3;

      rdot         = 3*growthL/LENGTH;
      aging        = (Q*(LENGTH*LENGTH*LENGTH)*SG/(Lm*Lm*Lm) + HADOTDOT)*e*(VDOT/LENGTH - rdot) - rdot*Q;
      hazardchange = Q - rdot*H;

      if (lifestage[0] == 1)                                                        // Juveniles
        maturation = (1 - KAPPA)*Sc - KJDOT_s*UH;
      else
        maturation = 0;
    }

  development[0][0] = 1.0;
  development[0][1] = growthL;                                                      // dL/da
  development[0][2] = Sa - Sc;                                                      // dUe/da
  development[0][3] = maturation;                                                   // dUh/da
  development[0][4] = aging;                                                        // dQ/da
  development[0][5] = hazardchange;                                                 // dH/da

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

// The following global variable is needed as it is part of the equilibrium condition
static double EBirth;

void DiscreteChanges(int lifestage[POPULATION_NR], double *istate[POPULATION_NR], double *birthstate[POPULATION_NR], int BirthStateNr, double E[])
{
#if (EATBIRTH == 0)
  double f;
#endif

  if (lifestage[0] == 1)                                                            // End of the embryo stage
    {
      // Store the computed life history value (scaled energy density at birth) used in the iteration
      // The realized scaled energy density at birth equals either 1 or f
      EBirth = VDOT*UE/(LENGTH*LENGTH*LENGTH);

// Reset Ue now to the value it is assumed to have at birth
// e = VDOT*Ue/L3 = 1 or = f
#if (EATBIRTH == 0)
      f = X/(X + K);
      f *= (1 - INTAKEEFFECT);
      UE = f*LENGTH*LENGTH*LENGTH/VDOT;
#elif (EATBIRTH == 1)
      UE       = LENGTH*LENGTH*LENGTH/VDOT;
#endif
    }

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
  double KMDOT_s, KJDOT_s, G_s, KAPPAR_s;
  double Sc;
  double e;

  KMDOT_s  = KMDOT*(1 + MAINTEFFECT)/(1 + GROWTHEFFECT);
  KJDOT_s  = KJDOT*(1 + MAINTEFFECT);
  G_s      = G*(1 + GROWTHEFFECT);
  KAPPAR_s = KAPPAR*exp(-EGGSURVEFFECT)/(1 + REPROEFFECT);

  if (lifestage[0] == 2)                                                            // Only for adults
    {
#if (EATBIRTH == 0)
      e = (1 - INTAKEEFFECT)*X/(X + K);
#elif (EATBIRTH == 1)
      e        = VDOT*UE/(LENGTH*LENGTH*LENGTH);
#endif
      Sc = LENGTH*LENGTH*G_s*e*(1 + LENGTH*KMDOT_s/VDOT)/(G_s + e);

      fecundity[0][0] = KAPPAR_s*((1 - KAPPA)*Sc - KJDOT_s*UHP)/UE0;                // In number of offspring
    }
  else
    fecundity[0][0] = 0;

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
  double e;

  if (lifestage[0] == 0)
    // I assume that embryos do not experience any mortality
    mortality[0] = 0;
  else
    {
#if (EATBIRTH == 0)
      e = (1 - INTAKEEFFECT)*X/(X + K);
#elif (EATBIRTH == 1)
      e        = VDOT*UE/(LENGTH*LENGTH*LENGTH);
#endif
      if (lifestage[0] == 1)                                                        // Juveniles
        mortality[0] = Q + M*(1 - e) + MUC;
      else                                                                          // Adults
        mortality[0] = Q + MUC;
    }

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
  double KMDOT_s, KJDOT_s, G_s, KAPPAR_s;
  double f, e, Sa, Sc, growthL;
  double fecundity;

  KMDOT_s  = KMDOT*(1 + MAINTEFFECT)/(1 + GROWTHEFFECT);
  KJDOT_s  = KJDOT*(1 + MAINTEFFECT);
  G_s      = G*(1 + GROWTHEFFECT);
  KAPPAR_s = KAPPAR*exp(-EGGSURVEFFECT)/(1 + REPROEFFECT);

  switch (lifestage[0])
    {
      case 0:
        impact[0][0] = 0.0;                                                         // No embryo resource intake

        // Number in size class
        impact[0][1] = 1.0;
        impact[0][2] = 0;
        impact[0][3] = 0;

        // Biomass in size class
        impact[0][4] = LENGTH*LENGTH*LENGTH;
        impact[0][5] = 0;
        impact[0][6] = 0;

        /*
       * Stage-specific net-biomass balance
       * Rewrite the original equations
       *
       * Sc = L^2*G_s*e*(1 + L*KMDOT/VDOT)/(G_s + e)
       *
       * with e  = VDOT*Ue/L3 (only valid for embryos)
       *
       * to obtain
       *
       * Sc = L^2*G_s*Ue*(VDOT + L*KMDOT_s)/(L3*G_s + VDOT*Ue)
       *
       * Substitute in
       *
       * growthL = (VDOT*Sc/(G_s*L^2) - KMDOT*L)/3
       *
       * to obtain:
       */
        growthL = (VDOT*UE*(VDOT + LENGTH*KMDOT_s)/((LENGTH*LENGTH*LENGTH)*G_s + VDOT*UE) - KMDOT_s*LENGTH)/3;

        impact[0][7] = 3.0*LENGTH*LENGTH*growthL;
        impact[0][8] = 0;
        impact[0][9] = 0;
        break;
      case 1:
        f = X/(X + K);
        f *= (1 - INTAKEEFFECT);
        Sa           = f*LENGTH*LENGTH;
        impact[0][0] = JXAMDOT*Sa;                                                  // Juvenile feeding

        // Number in size class
        impact[0][1] = 0;
        impact[0][2] = 1.0;
        impact[0][3] = 0;

        // Biomass in size class
        impact[0][4] = 0;
        impact[0][5] = LENGTH*LENGTH*LENGTH;
        impact[0][6] = 0;

// Stage-specific net-biomass balance
#if (EATBIRTH == 0)
        e = f;
#elif (EATBIRTH == 1)
        e      = VDOT*UE/(LENGTH*LENGTH*LENGTH);
#endif
        /*
       * Substitute
       *
       * Sc = L^2*G_s*e*(1 + L*KMDOT_s/VDOT)/(G_s + e)
       *
       * into
       *
       * growthL = (VDOT*Sc/(G_s*L^2) - KMDOT*L)/3
       *
       * to obtain
       *
       * growthL = (e*(VDOT + KMDOT*L)/(G + e) - KMDOT*L)/3
       *
       * and rewrite:
       */
        growthL = (VDOT*e/(G_s + e) - KMDOT_s*LENGTH*G_s/(G_s + e))/3;

        impact[0][7] = 0;
        impact[0][8] = 3.0*LENGTH*LENGTH*growthL + (Q + M*(1 - e) + MUC)*(LENGTH*LENGTH*LENGTH);
        impact[0][9] = 0;
        break;
      case 2:
        f = X/(X + K);
        f *= (1 - INTAKEEFFECT);
        Sa           = f*LENGTH*LENGTH;
        impact[0][0] = JXAMDOT*Sa;                                                  // Adult feeding

        // Number in size class
        impact[0][1] = 0;
        impact[0][2] = 0;
        impact[0][3] = 1.0;

        // Biomass in size class
        impact[0][4] = 0;
        impact[0][5] = 0;
        impact[0][6] = LENGTH*LENGTH*LENGTH;

// Stage-specific net-biomass balance
#if (EATBIRTH == 0)
        e = f;
#elif (EATBIRTH == 1)
        e      = VDOT*UE/(LENGTH*LENGTH*LENGTH);
#endif

        Sc = LENGTH*LENGTH*G_s*e*(1 + LENGTH*KMDOT_s/VDOT)/(G_s + e);

        /*
       * Substitute the above expression for Sc into
       *
       * growthL = (VDOT*Sc/(G_s*L^2) - KMDOT*L)/3
       *
       * to obtain
       *
       * growthL = (e*(VDOT + KMDOT*L)/(G + e) - KMDOT*L)/3
       *
       * and rewrite:
       */
        growthL   = (VDOT*e/(G_s + e) - KMDOT_s*LENGTH*G_s/(G_s + e))/3;
        fecundity = KAPPAR_s*((1 - KAPPA)*Sc - KJDOT_s*UHP);

        fecundity    = KAPPAR_s*((1 - KAPPA)*Sc - KJDOT_s*UHP);
        impact[0][7] = 0;
        impact[0][8] = 0;
        impact[0][9] = 3.0*LENGTH*LENGTH*growthL + fecundity - (Q + MUC)*(LENGTH*LENGTH*LENGTH);
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

const int EnvironmentType[ENVIRON_DIM] = {GENERALODE, GENERALODE};

void EnvEqui(double E[], double I[POPULATION_NR][INTERACT_DIM], double condition[ENVIRON_DIM])
{
#if (EATBIRTH == 0)
  double f;

  f = X/(X + K);
  f *= (1 - INTAKEEFFECT);

  condition[0] = DELTA*(XMAX - E[0]) - I[0][0]/VOLUME;
  condition[1] = EBirth - f;
#elif (EATBIRTH == 1)
  condition[0] = DELTA*(XMAX - E[0]) - I[0][0]/VOLUME;
  condition[1] = EBirth - 1;
#endif

  return;
}


/*==================================================================================================================================*/
