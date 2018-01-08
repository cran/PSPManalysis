/*
   NAME
     ebtdopri5.h

     Interface header file to the functions defined in ebtdopri5.c

  Last modification: AMdR - Sep 06, 2017
*/

#ifndef EBTDOPRI5_H
#define EBTDOPRI5_H

#ifdef	EBTDOPRI5_C
#undef	EXTERN
#define EXTERN
#else
#undef	EXTERN
#define EXTERN	                  extern
#endif

EXTERN void	                      PrepareCycle(void);
EXTERN double	                    IntegrationStep(double, double, int);
EXTERN void                       ResetDopri5(void);


/*==================================================================================================================================*/
#endif // EBTDOPRI5_H
