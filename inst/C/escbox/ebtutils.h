/*
  NAME
     ebtutils.h

     Interface header file to the functions defined in ebtutils.c

  Last modification: AMdR - Jan 19, 2023
*/

#ifndef EBTUTILS_H
#define EBTUTILS_H

#ifdef	EBTUTILS_C
#undef	EXTERN
#define EXTERN
#else
#undef	EXTERN
#define EXTERN	                  extern
#endif

EXTERN int	                      imin(int, int);
EXTERN int	                      imax(int, int);
EXTERN double	                    min(double, double);
EXTERN double	                    max(double, double);
EXTERN int	                      ismissing(double);
EXTERN int	                      isequal2zero(double);
EXTERN int	                      isequal(double, double);
EXTERN void                       ShutDown(int exitcode);
EXTERN void	                      ErrorAbort(const char *);
EXTERN void	                      ErrorExit(const int, const char *);
EXTERN void	                      Warning(const char *);
EXTERN void                       FileOut(void);
EXTERN void                       FileState(void);
EXTERN void	                      *Myalloc(void *, size_t, size_t);
EXTERN void	                      PrettyPrint(FILE *fp, double output);
EXTERN void	                      WriteStateToFile(FILE *fp, double *data);
EXTERN void	                      SetBifOutputTimes(double *);
EXTERN void	                      outputMeasureBifstats(double *env);
EXTERN void	                      measureBifstats(double *env, population *pop, int Nstores, ...);
EXTERN void                       ResetBifStats(void);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
// For Ctrl-C detection
int                               checkInterrupt(void);

EXTERN int                        CtrlCPressed;
#endif


/*==================================================================================================================================*/
#endif // EBTUTILS_H
