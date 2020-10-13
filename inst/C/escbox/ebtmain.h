/*
  NAME
     ebtmain.h

     Interface header file to the functions defined in ebtmain.c

  Last modification: AMdR - Oct 13, 2020
*/

#ifndef EBTMAIN_H
#define EBTMAIN_H

/*==================================================================================================================================*/
/*
 * Definitions of global variables exported to the linker. All the global 
 * variables are defined in the ebtmain.h file and exported to the other 
 * modules. This facilitates location of their definitions.
 */

#ifdef 	PSPMECODYN
#undef  EXTERN
#define EXTERN
#else
#undef  EXTERN
#define EXTERN	                  extern
#endif

EXTERN int                        PopulationNr;                                     // Number of populations
EXTERN int                        Stages;                                           // Number of life history stages
EXTERN int                        IStateDim;                                        // I-state dimension
EXTERN int                        EnvironDim;                                       // Environment dimension
EXTERN int                        InteractDim;                                      // Interaction Dimension
EXTERN int                        ParameterNr;                                      // Parameter number

EXTERN int	                      IConstDim;                                        // I-constant dimension
EXTERN int	                      EBTEnvironDim;                                    // Environment dimension for EBT

EXTERN int	                      output_var_nr;                                    // Output variable number

EXTERN double                     *initState;                                       // Pointer to initial system state
EXTERN double                     *currentState;                                    // Pointer to most current system state
EXTERN double                     *currentDers[MAXDERS];                            // Pointers to computed derivatives
EXTERN double                     env[ENVIRON_DIM_EBT];                             // Becomes vector of environmental values
EXTERN population                 pop[POPULATION_NR];                               // Array of pointers to the population data
EXTERN population                 ofs[POPULATION_NR];                               // Array of pointers to the offspring data
EXTERN population                 bpoints[POPULATION_NR];                           // Array of pointers to the boundary points
EXTERN popID                      popIDcard[POPULATION_NR];                         // Population ID structure
EXTERN popID                      ofsIDcard[POPULATION_NR];                         // & offspring ID structure

EXTERN int                        CohortNo[POPULATION_NR];                          // Number of cohorts in the various populations
EXTERN int                        cohort_no[POPULATION_NR];                         // Copy of CohortNo vector accessible to user
EXTERN int                        BpointNo[POPULATION_NR];                          // Number of boundary points for boundary cohorts
EXTERN int                        bpoint_no[POPULATION_NR];                         // Copy of BpointNo vector accessible to user

EXTERN double                     cohort_limit;                                     // Time limit new cohorts

EXTERN double                     max_time;                                         // Maximum integration time

EXTERN double                     Odesolve_Init_Step;                               // Dopri5 initial step size
EXTERN double                     Odesolve_Fixed_Step;                              // Dopri5 fixed step size
EXTERN double                     Odesolve_Min_Step;                                // Dopri5 minimum step size
EXTERN double                     Odesolve_Max_Step;                                // Dopri5 maximum step size
EXTERN double                     Odesolve_Abs_Err;                                 // Dopri5 absolute error in step
EXTERN double                     Odesolve_Rel_Err;                                 // Dopri5 relative error in size
EXTERN double                     Odesolve_Func_Tol;                                // Dopri5 function tolerance

EXTERN int                        cohort_end;                                       // Flag indicating dynamic ending of cohort cycle
EXTERN double                     step_size;                                        // The step size for the time integration
EXTERN double                     output[OUTPUT_VAR_NR + 2];                        // Array with output values

EXTERN int                        outputDefined;                                    // Output is defined flag

EXTERN double                     delt_out;                                         // Output time interval

EXTERN double                     state_out;                                        // The time interval for complete state output
EXTERN double                     next_cohort_end;                                  // Time of cohort closure

EXTERN double                     next_output, next_state_output;                   // Time of next (state) output
EXTERN FILE                       *outfile;                                         // Pointer to outputfile

EXTERN FILE                       *csbfile;                                         // Pointer binary state file

EXTERN FILE                       *averages;                                        // File pointers for bifurcation output
EXTERN FILE                       *gaverages;
EXTERN FILE                       *variances;
EXTERN FILE                       *extrema;

EXTERN int                        csbnew;                                           // Flag new state file

EXTERN FILE                       *dbgfile;                                         // Pointer debug report file

EXTERN double                     identical_zero;                                   // Tolerance value, determining identity with 0.0
EXTERN int                        pop_extinct[POPULATION_NR];                       // Flag indicating whether population extinction has been signalled already

EXTERN int                        BifurcationRun;                                   // Flag indicating whether or not this is a bifurcation run
EXTERN int                        BifParIndex;                                      // Index of bifurcation parameter
EXTERN double                     BifParStep;                                       // Step size in bifurcation parameter
EXTERN double                     BifParLastVal;                                    // Last value of bifurcation parameter
EXTERN double                     BifOutput;                                        // Output period during bifurcation run
EXTERN double                     BifStateOutput;                                   // State output period during bifurcation run
EXTERN double                     BifParBase;                                       // First value of bifurcation parameter
EXTERN double                     BifPeriod;                                        // Parameter change period during bifurcation run
EXTERN char                       progname[MAXFILENAMELEN];                         // Name of the program

EXTERN char                       runname[2*MAXFILENAMELEN];                        // Name of the current run

EXTERN int                        debug_level;                                      // Level of debug info
EXTERN int                        report_level;                                     // Interval of reporting data to console

EXTERN long                       DataMemAllocated[POPULATION_NR];                  // Total number of doubles currently allocated
EXTERN long                       IDMemAllocated[POPULATION_NR];                    // Total number of doubles currently allocated

EXTERN int                        rk_level;                                         // Number of current Runge- Kutta evaluation (1-6)

EXTERN int                        ForcedCohortEnd;                                  // Flag indicating forced dynamic cohort closure

EXTERN int                        ForcedRunEnd;                                     // Flag indicating forced ending of entire run

EXTERN int                        LocatedEvent;                                     // Index of located event


/*==================================================================================================================================*/
#endif // PSPMECODYN
