#
# StageStructuredBiomass.R -  R file specifying the elementary life-history functions
#                   of the size-structured consumer-unstructured resource model from
#
#                   A.M. de Roos et al., 2008. Simplifying a physiologically structured
#                   population model to a stage-structured biomass model.
#                   Theor. Pop. Biol. 73: 47–62.
#
#                   This model forms the basis for the stage-structured biomass model
#                   analyzed in:
#
#                   A.M. de Roos et al., 2007. Food-dependent growth leads to
#                   overcompensation in stage- specific biomass when mortality
#                   increases: the influence of maturation versus reproduction
#                   regulation. American Naturalist 170: E59-E76.
#
#                   Consumers hence follow a net-production DEB model with scaling
#                   functions of ingestion and maintenance proportional to body size.
#                   Juveniles and adults use all their net-production for growth and
#                   reproduction, respectively.
#
#                   The variable that is solved for is the biomass of the resource.
#
#   Model i-state variables:
#     1 : Age
#     2 : Size
#
#   Model of environmental state variables:
#     1 : Resource
#
#   Model interaction (and output) variables:
#     1 : Total resource ingestion by the consumer population
#     2 : Total biomass of juveniles
#     3 : Total biomass of adults
#
#   Last modification: AMdR - Nov 22, 2017
#

#
# Model dimension variables: PSPMdimensions (required)
#
# Define a numerical (integer) vector called 'PSPMdimensions' that specifies the
# dimensions of the model. The vector should include the following named vecctor
# elements:
#
# PopulationNr:       The number of populations in the model that are structured,
#                     that is, of which the life history of individuals is explicitly
#                     modelled
# IStateDimension:    The number of individual state variables that characterise the
#                     individuals in the structured populations. This number should be
#                     the same for all structured populations in the model
# LifeHistoryStages:  The number of distinct and discrete stages in the life history
#                     of the individuals. A part of the individual life history is
#                     considered a stage, when at the boundary of this part one of the
#                     life history processes (development, fecundity, mortality or
#                     impact) changes discontinuously
# ImpactDimension:    The number of functions that represent the impact of an individual
#                     on its environment. At the population level these impacts
#                     affect the dynamics of the environemtal state variables and hence
#                     determine their equilibrium values. Impact functions can, however,
#                     also be used to extract certain output information on the
#                     structured populations (such as number of biomass of juveniles and
#                     adults) as all population-level values of the impact functions
#                     are reported as output.
#
# An error will occur when one of the above named elements is missing.
#
PSPMdimensions <- c(PopulationNr = 1, IStateDimension = 2, LifeHistoryStages = 2, ImpactDimension = 3)


#
# Variable name: NumericalOptions (optional)
#
# Define a numerical vector called 'NumericalOptions' the elements of which specify
# the optional numerical settings to tweak the computations. The elements of the
# vector should have names corresponding to one of the possible numerical settings
# (see the vignette). Examples of such numerical settings are 'MIN_SURVIVAL = 1.0E-9'
# and 'RHSTOL = 1.0E-6', which set the survivial at which the individual is considered
# dead and the tolerance value determining when a solution has been found, respectively.
#

NumericalOptions <- c(MIN_SURVIVAL  = 1.0E-9,                                       # Survival at which individual is considered dead
                      MAX_AGE       = 100000,                                       # Give some absolute maximum for individual age
                      DYTOL         = 1.0E-7,                                       # Variable tolerance
                      RHSTOL        = 1.0E-6)                                       # Function tolerance

#
# Variable name: EnvironmentState (required)
#
# Define a vector called 'EnvironmentState' with a length equal to the number of
# environmental state variables in the problem. Each element of this vector should be
# one of the strings "PERCAPITARATE", "GENERALODE" or "POPULATIONINTEGRAL", defining
# the nature or type of environmental state variable. The 3 different types are defined
# as follows:
#
# Set an entry to "PERCAPITARATE" if the dynamics of E[j] follow an ODE and 0
# is a possible equilibrium state of E[j]. The ODE is then of the form
# dE[j]/dt = P(E,I)*E[j], with P(E,I) the per capita growth rate of E[j].
# Specify the equilibrium condition as condition[j] = P(E,I), do not include
# the multiplication with E[j] to allow for detecting and continuing the
# transcritical bifurcation between the trivial and non-trivial equilibrium.
#
# Set an entry to "GENERALODE" if the dynamics of E[j] follow an ODE and 0 is
# NOT an equilibrium state of E. The ODE then has a form dE[j]/dt = G(E,I).
# Specify the equilibrium condition as condition[j] = G(E,I).
#
# Set an entry to "POPULATIONINTEGRAL" if E[j] is a (weighted) integral of the
# population distribution, representing for example the total population
# biomass. E[j] then can be expressed as E[j] = I[p][i]. Specify the
# equilibrium condition in this case as condition[j] = I[p][i].
#
# If the members of the vector 'EnvironmentState' are given names, these names can be
# used in the functions below that define the life history processes.

EnvironmentState <- c(R = "GENERALODE")

#
# Variable name: DefaultParameters (required)
#
# Define a vector called 'DefaultParameters' with a length equal to the number of
# parameters in the model. Each element of this vector should be given the default
# for the particualr parameter. If the members of the vector 'DefaultParameters' are
# given names, these names can be used conveniently in the functions below that define
# the life history processes.

DefaultParameters <- c(Rho = 0.1, Rmax = 30.0, z = 0.01, Mc = 1.0, Hc = 3.0, q = 0.5, Sigma = 0.5, Tc = 0.1, Muc = 0.015, Muj = 0.0, Mua = 0.0)

#
# Function name: StateAtBirth  (required)
#
# Specify for all structured populations in the problem all possible values of the individual state
# at birth.
#
# Function arguments:
#
#   E:    Vector with the current values of the environmental state variables.
#   pars: Vector with the model parameters
#
# Required return:
#
# A vector of length equal to the number of i-state variables. The biological interpretation of
# each of the i-state variables is completely up to the user. Each element should specify the
# numeric value of the particular i-state variable with which the individual is born.
# at birth for the particular i-state variable. If the members of the vector are given
# meaningful names, these names can be used conveniently in the functions below that define
# the life history processes.
#
# If individuals can differ in their individual state at birth this function should return a
# matrix with the number of rows equal to the number of possible states at birth and the number
# of columns equal to the number of i-state variables. Each row then specifies the value of the
# individual state variable of the particualr state at birth.
#
# In case the model accounts for multiple, structured populations this function should return a
# a matrix with the number of rows equal to the number of structured populations in the problem
# and the number of columns equal to the number of i-state variables.
#
# In case the model accounts for multiple, structured populations and individuals can differ in
# their individual state at birth this function should return a 3-dimensional array with the
# first dimension having a length equal to the number of structured populations in the problem,
# the second dimension equal to the number of possible states at birth and the third dimension
# equal to the number of i-state variables.

StateAtBirth <- function(E, pars)
{
 with(as.list(c(E, pars)),{
    # We model a single structured population with two i-state variables:
    # 1: age (initial value 0); 2: length (initial value equal to parameter z)
    c(Age = 0.0, Size = z)
   })
}

#
# Function name: LifeStageEndings  (required)
#
# Specify for all structured populations in the problem the threshold value at which the current life
# stage of the individual ends and the individual matures to the next life history stage. The threshold
# value may depend on the current i-state variables of the individual, its state at birth and the life
# stage that it currently is in.
#
# Function arguments:
#
#  lifestage:     Integer value specifying the life stage that the individual is currently in.
#                 These stages are numbered 1 (youngest) and up.
#                 In case the model accounts for multiple, structured populations this argument
#                 is a vector of integer values.
#  istate:        Vector of length equal to the number of i-state variables that charaterize
#                 the state of an individual. The biological interpretation of the i-state
#                 variables is up to the user. Each element specifies the current value of the
#                 particular i-state variable of the individual.
#                 In case the model accounts for multiple, structured populations this argument
#                 is a matrix with the number of rows equal to the number of structured populations
#                 in the problem and the number of columns equal to the number of i-state variables.
#  birthstate:    Vector of length equal to the number of i-state variables that charaterize
#                 the state of an individual. Each element specifies the value of the particular
#                 i-state variables at which the individual was born.
#                 In case the model accounts for multiple, structured populations this argument
#                 is a matrix with the number of rows equal to the number of structured populations
#                 in the problem and the number of columns equal to the number of i-state variables.
#  BirthStateNr:  The integer index of the state of birth to be specified, ranging from 1 and up.
#  E:             Vector with the current values of the environmental state variables.
#  pars:          Vector with the model parameters
#
# Required return:
#
# maturation:   A single value specifying when the current life stage of the individual ends and
#               the individual matures to the next life history stage. The end of the current life
#               history stage occurs when this threshold value becomes 0 and switches sign from
#               negative to  positive. For the final life stage (which never ends) return a constant
#               negative value (for example, -1)
#               In case the model accounts for multiple, structured populations this argument
#               is a vector with the number of elements equal to the number of structured populations
#               in the problem.

LifeStageEndings <- function(lifestage, istate, birthstate, BirthStateNr, E, pars) {
  with(as.list(c(E, pars, istate)),{
        maturation  = switch(lifestage, Size - 1.0, -1)
      })
}

#
# Function name: LifeHistoryRates  (required)
#
# Specify for all structured populations in the problem the rates of the various life history
# processes (development, fecundity, mortality and impact on the environment)  of an individual
# as a function of its i-state variables, the individual's state at birth and the life
# stage that the individual is currently in.
#
# Function arguments:
#
#  lifestage:     Integer value specifying the life stage that the individual is currently in.
#                 These stages are numbered 1 (youngest) and up.
#                 In case the model accounts for multiple, structured populations this argument
#                 is a vector of integer values.
#  istate:        Vector of length equal to the number of i-state variables that charaterize
#                 the state of an individual. The biological interpretation of the i-state
#                 variables is up to the user. Each element specifies the current value of the
#                 particular i-state variable of the individual.
#                 In case the model accounts for multiple, structured populations this argument
#                 is a matrix with the number of rows equal to the number of structured populations
#                 in the problem and the number of columns equal to the number of i-state variables.
#  birthstate:    Vector of length equal to the number of i-state variables that charaterize
#                 the state of an individual. Each element specifies the value of the particular
#                 i-state variables at which the individual was born.
#                 In case the model accounts for multiple, structured populations this argument
#                 is a matrix with the number of rows equal to the number of structured populations
#                 in the problem and the number of columns equal to the number of i-state variables.
#  BirthStateNr:  The integer index of the state of birth to be specified, ranging from 1 and up.
#  E:             Vector with the current values of the environmental state variables.
#  pars:          Vector with the model parameters
#
# Required return:
#
# A list with 4 components, named "development", "fecundity", "mortality" and "impact". The
# components should have the following structure:
#
# development:  A vector of length equal to the number of i-state variables. Each element
#               specifies the rate of development for the particular i-state variable.
#               In case the model accounts for multiple, structured populations this component
#               should be a matrix with the number of rows equal to the number of structured
#               populations in the problem and the number of columns equal to the number of
#               i-state variables.
# fecundity:    The value of the current fecundity of the individual.
#               In case the model accounts for multiple, structured populations this argument
#               is a matrix of fecundities with the number of rows equal to the number
#               of structured populations in the problem and a single column.
#               In case individuals can be born with different states at birth the component
#               should have a number of columns equal to the number of states at birth. Each
#               column should specify the number of offspring produced with the particular
#               state at birth.
# mortality:    A single value specifying the current mortality rate that the individual experiences.
#               In case the model accounts for multiple, structured populations this argument
#               is a vector of mortality rates with the number of elements equal to the number
#               of structured populations in the problem.
# impact:       A single value or a vector of a length equal to the number of impact functions that
#               need to be monitored for the individual. The value (or the values of the vector)
#               should specify the current contribution of the individual to this population-level
#               impact.
#               In case the model accounts for multiple, structured populations this component
#               should be a matrix with the number of rows equal to the number of structured
#               populations in the problem and the number of columns equal to the number of impact
#               functions.

LifeHistoryRates <- function(lifestage, istate, birthstate, BirthStateNr, E, pars) {
        with(as.list(c(E, pars, istate)),{
                Ingest        = switch(lifestage, Mc*Size*R/(Hc + R), q*Mc*Size*R/(Hc + R))
                netproduction = Sigma*Ingest - Tc*Size

                list(
                        # We model a single structured population (nrow=1) with two i-state variables:
                        development = switch(lifestage, c(1.0, netproduction), c(1.0, 0.0)),

                        fecundity   = switch(lifestage, 0, netproduction/z),

                        mortality   = switch(lifestage, Muc + Muj, Muc + Mua),

                        impact      = switch(lifestage, c(Ingest, Size, 0), c(Ingest, 0, Size))
                )
        })
}

#
# Function name: DiscreteChanges  (optional)
#
# If discrete changes in the individual state occur at the moment individuals mature to the next
# life history stage, specify for all structured populations in the problem these discrete changes
# (jumps) in the individual state variables on ENTERING a particular life stage. The life stage that is
# entered is specified by the vector 'lifestage'.
#
# NOTE: LEAVE THIS FUNCTION UNSPECIFIED IF THERE ARE NO DISCRETE CHANGES, AS THIS SAVES COMPUTATION TIME
#
# Function arguments:
#
# lifestage:    Integer value specifying the life stage that the individual is currently in.
#               These stages are numbered 1 (youngest) and up.
#               In case the model accounts for multiple, structured populations this argument
#               is a vector of integer values.
# istate:       Vector of length equal to the number of i-state variables that charaterize
#               the state of an individual. The biological interpretation of the i-state
#               variables is up to the user. Each element specifies the current value of the
#               particular i-state variable of the individual.
#               In case the model accounts for multiple, structured populations this argument
#               is a matrix with the number of rows equal to the number of structured populations
#               in the problem and the number of columns equal to the number of i-state variables.
# birthstate:   Vector of length equal to the number of i-state variables that charaterize
#               the state of an individual. Each element specifies the value of the particular
#               i-state variables at which the individual was born.
#               In case the model accounts for multiple, structured populations this argument
#               is a matrix with the number of rows equal to the number of structured populations
#               in the problem and the number of columns equal to the number of i-state variables.
# BirthStateNr: The integer index of the state of birth to be specified, ranging from 1 and up.
# E:            Vector with the current values of the environmental state variables.
# pars:         Vector with the model parameters
#
# Required return:
#
# A vector of length equal to the number of i-state variables. Each element should specify the value
# of the particular i-state variable after the transition to the current state, as given in
# the vector 'lifestage'
# In case the model accounts for multiple, structured populations this function should return
# a matrix with the number of rows equal to the number of structured populations in the problem
# and the number of columns equal to the number of i-state variables.

# DiscreteChanges <- function(lifestage, istate, birthstate, BirthStateNr, E, pars) {
#  with(as.list(c(E, pars)),{
#    # No discrete changes in this problem, function is commented out, which
#    # would be equivalent to returning a copy of the input argument 'istate'
#    istate
#   })
# }

# Specify for each of the environmental state variables the condition that determines its
# equilibrium value, dependent on the values of the environmental state variables themselves,
# the population impact values and the paramenters.
#
# I:    Matrix of population impact values. The number of rows of this matrix equals
#       the number of structured populations in the problem, while the number of columns
#       equals the number of impact factors of an individual on its environment, as specified
#       by the function Impat(). The interpretation of the latter is up to the user.
# E:    List with named current values of the environmental state variables
# pars: List with named model parameters
#
# Required return:
#
# A vector with a length equal to the number of environmental state variables in the problem.
# Each entry i of this vector should specify the equilibrium condition for the particular environmental
# state variable (i).

EnvEqui <- function(I, E, pars) {
 with(as.list(c(E, pars)),{
    Rho*(Rmax - R) - I[1]
   })
}

