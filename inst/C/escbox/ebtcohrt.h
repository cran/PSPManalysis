/*
   NAME
     ebtcohrt.h

     Interface header file to the functions defined in ebtcohrt.c

  Last modification: AMdR - Aug 30, 2017
*/

/*==================================================================================================================================*/
#ifndef EBTCOHRT_H
#define EBTCOHRT_H

#ifdef	EBTCOHRT_C
#undef	EXTERN
#define EXTERN
#else
#undef	EXTERN
#define EXTERN	                  extern
#endif

EXTERN void     	                CohortCycle(double);
EXTERN int		                    InsCohort(cohort, cohortID, int);
EXTERN void		                    TransBcohorts(void);
EXTERN void		                    SievePop(void);
EXTERN void                       ResetCohorts(void);


/*==================================================================================================================================*/
#endif // EBTCOHRT_H
