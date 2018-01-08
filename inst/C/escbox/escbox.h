/*
  NAME
     escbox.h

     This file contains all include statements for the program, some
     constant and type definitions that are used in the population data
     representations and the prototypes of various functions.
     This file is meant to be included in all escalator boxcar train
     modules, as it contains the global definitions that are accessible to
     all modules. 

  Last modification: AMdR - Oct 25, 2017
*/

#ifndef ESCBOX_H
#define ESCBOX_H

#define MAXDERS                   20                                                // Max. stages in ODE solver

#include  "ctype.h"
#include  "math.h"
#include  "stdio.h"
#include  "string.h"
#include  "float.h"
#include  "signal.h"
#include  "fenv.h"
#include  "limits.h"
#include  "stddef.h"
#include  "stdlib.h"
#include  "stdarg.h"
#include  "sys/stat.h"

#ifndef MY_OMP_H_
#define MY_OMP_H_

#ifdef _OPENMP                                                                       // Invocation with -fopenmp or -openmp argument
#undef OPENMP
#define OPENMP                    1
#include <omp.h>
#if _OPENMP >= 201307
#define OMP_VER_4
#endif
#endif

// Insert SIMD pragma if supported
#ifdef OMP_VER_4
#define SAFE_SIMD _Pragma("omp simd")
#define SAFE_FOR_SIMD _Pragma("omp for simd")
#else
#define SAFE_SIMD
#define SAFE_FOR_SIMD
#endif

#endif // MY_OMP_H_


/*==================================================================================================================================*/
/*
 * Definition of some global constants. Most of these definitions only take
 * effect if the flag EBTLIB is defined, which occurs in the
 * Escalator Boxcar Train module files.
 */

#define  COHORT_SIZE              (1+I_STATE_DIM)	                                  // size of cohort
#define  ENVIRON_DIM_EBT          (1+ENVIRON_DIM)	                                  // size of EBT environment vector

#ifndef I_CONST_DIM
#define I_CONST_DIM	              (I_STATE_DIM + 4)
#endif

#ifndef OUTPUT_VAR_NR
#define OUTPUT_VAR_NR             (ENVIRON_DIM + POPULATION_NR + POPULATION_NR*INTERACT_DIM)
#endif

#ifndef EVENT_NR
#define EVENT_NR                  (POPULATION_NR * STAGES)
#endif

#if defined(DBL_MAX)				                                                        // Stub value for missing data point
#define   MISSING_VALUE		        (DBL_MAX)
#elif defined(MAXDOUBLE)
#define   MISSING_VALUE		        (MAXDOUBLE)
#else
#define   MISSING_VALUE		        (1.23456789e+307)
#endif

#define   BIFTINY	 	              1.0E-4		                                        // Tiny value used in bifurcation computations

#ifdef EBTLIB
#if defined(PATH_MAX)
#define   MAXFILENAMELEN	 	      (PATH_MAX)                                        // Maximum length of the full path to a file
#elif defined(MAXPATHLEN)
#define   MAXFILENAMELEN	 	      (MAXPATHLEN)
#elif defined(_MAX_PATH)
#define   MAXFILENAMELEN	 	      (_MAX_PATH)
#else
#define   MAXFILENAMELEN	 	      (512)
#endif

#define   MAX_STR_LEN             1024
#define	  MEM_BLOCK_SIZE          256			                                          // Number of doubles in allocated memory block
#define	  MIN_ACCURACY	          1.0E-16		                                        // The minimum accuracy allowed in integration
#endif // EBTLIB

#define   NO_EVENT                -1
#define   NO_COHORT_END           0
#define   COHORT_END              1

#define   EBTDEBUG(a)             (dbgfile && (debug_level >= (a)))

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#include "mex.h"
#define STDOUT                    mexPrintf
#elif defined(R_PACKAGE)
#define R_USE_C99_IN_CXX
#include "R.h"
#include <Rinternals.h>
#include <Rmath.h>

#include "R_ext/Print.h"
#include <R_ext/Utils.h>

#define STDOUT                    Rprintf
#else
#define STDOUT                    printf
#endif

/*==================================================================================================================================*/
/*
 * Definition of the cohort structure that will hold the data of one 
 * population cohort. Also defined is a pointer type to a single cohort and
 * the population data type.
 */

typedef double                    cohort[COHORT_SIZE];
typedef	cohort                    *cohort_pnt;
typedef	cohort_pnt	              population;
typedef double                    cohortID[I_CONST_DIM];
typedef	cohortID	                *cohortID_pnt;
typedef	cohortID_pnt	            popID;
typedef int		                    cohort_ind[COHORT_SIZE];
typedef	cohort_ind	              *cohort_ind_pnt;
typedef	cohort_ind_pnt	          population_index;



/*==================================================================================================================================*/
/* 
 * Defining macro's for global use in all files
 */

#define number			              (0)
#define i_state(a)		            ((a+1))
#define i_const(a)		            ((a))

#define	DEF_TYPE		              void		                                          // Basic pointer type
#define CONST			const
#define SIZE_TYPE		size_t

#define MemBlocks(a)              (((a/MEM_BLOCK_SIZE)+1)*MEM_BLOCK_SIZE)


/*==================================================================================================================================*/
/*
 * Definition of global variables and functions accessible in user
 * specified problem file.
 */

#ifndef EBTLIB
extern double		                  cohort_limit;
extern double		                  next_output, next_state_output;
extern int		                    cohort_no[POPULATION_NR];
extern int		                    bpoint_no[POPULATION_NR];
extern int		                    rk_level;
extern int		                    ForcedCohortEnd;
extern int		                    ForcedRunEnd;
extern int		                    LocatedEvent;
extern int		                    ParamerNr;
extern int		                    BifParIndex;
extern double		                  BifOutput;
extern double		                  BifStateOutput;
extern double		                  BifPeriod;
extern popID		                  popIDcard[POPULATION_NR];
extern popID		                  ofsIDcard[POPULATION_NR];

extern int		                    imax(int, int);
extern int		                    imin(int, int);
extern double		                  max(double, double);
extern double		                  min(double, double);

extern int		                    iszero(double);
extern int		                    isequal(double, double);
extern void		                    SievePop(void);
extern void		                    ErrorAbort(const char *);
extern void		                    ErrorExit(const int, const char *);
extern void		                    Warning(const char *);
extern void		                    measureBifstats(double *env, population *pop, int Nstores, ...);
#else
extern double		                  parameter[PARAMETER_NR];
#endif // EBTLIB




/*==================================================================================================================================*/
/*
 * Prototyping all functions that are defined in the problem-specific
 * program file
 */

#ifdef	EBTLIB
#undef	EXTERN
#define EXTERN	                  extern
#else
#undef	EXTERN
#define EXTERN
#endif

EXTERN void	                      UserInit(int, char **, double *, population *);
EXTERN void	                      SetBpointNo(double *, population *, int *);
EXTERN void	                      SetBpoints(double *, population *, population *);
EXTERN void	                      Gradient(double *, population *, population *, double *, population *, population *, population *);
EXTERN void	                      EventLocation(double *, population *, population *, population *, double *);
EXTERN int	                      ForceCohortEnd(double *, population *, population *, population *);
EXTERN void	                      InstantDynamics(double *, population *, population *);
EXTERN void	                      DefineOutput(double *, population *, double *);

/*==================================================================================================================================*/
#endif // ESCBOX_H
