/*
  NAME
     ebtinit.h

     Interface header file to the functions defined in ebtinit.c
     
   Last modification: AMdR - Aug 30, 2017
*/

/*==================================================================================================================================*/
#ifndef EBTINIT_H
#define EBTINIT_H

#ifdef	EBTINIT_C
#undef	EXTERN
#define EXTERN
#else
#undef	EXTERN
#define EXTERN	                  extern
#endif
EXTERN void	                      Initialize(int, char **);



/*==================================================================================================================================*/
#endif // EBTINIT_H
