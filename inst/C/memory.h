/***
  NAME
    memory.h
  DESCRIPTION
    Header file with implementation of the memory handling routines.
    This header file is included in the files PSPMdemo.c, PSPMequi.c
    PSPMevodyn.c and PSPMind.c.

    Copyright (C) 2015, Andre M. de Roos, University of Amsterdam

    This file is part of the PSPManalysis software package.

    This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Dec 05, 2017
***/

/*
 *====================================================================================================================================
 *  Definition of addressing macros
 *====================================================================================================================================
 */

#define FinalIstate(b, p, i)      (*(FinalIstateMem + ((b)*CurPopulationNr + (p))*MaxCohortDim + (i)))  // FinalIstate[MaxStatesAtBirth][CurPopulationNr][MaxCohortDim]
#define FinalIstatePnt(b, p, i)   (FinalIstateMem + ((b)*CurPopulationNr + (p))*MaxCohortDim + (i))     // FinalIstate[MaxStatesAtBirth][CurPopulationNr][MaxCohortDim]
#if (PSPMDEMO == 1)
#define LeftEigenvec(p, b)        (*(LeftEigenvecMem + (p)*MaxStatesAtBirth + (b)))                     // LeftEigenvec[CurPopulationNr][MaxStatesAtBirth]
#endif

#if (FULLSTATEOUTPUT > 0)
#define CohortLimit(b, p)         (*(CohortLimitMem + (b)*CurPopulationNr + (p)))                       // CohortLimit[MaxStatesAtBirth][CurPopulationNr]
#endif

/*
 *====================================================================================================================================
 *  Allocation and deallocation of global memory
 *====================================================================================================================================
 */

int AllocateHeapMemory(void)

{
  if ((CurPopulationNr*MaxStatesAtBirth) > LastMemAllocated)
    {
      if (!(RightEigenvecMem = realloc(RightEigenvecMem, CurPopulationNr*MaxStatesAtBirth*sizeof(double)))) return ReportMemError("AllocateHeapMemory");
      rightEigenvecMem       = RightEigenvecMem;

#if (PSPMDEMO == 1)
      if (!(LeftEigenvecMem = realloc(LeftEigenvecMem, CurPopulationNr*MaxStatesAtBirth*sizeof(double)))) return ReportMemError("AllocateHeapMemory");
#endif

#if (FULLSTATEOUTPUT > 0)
      if (!(PopDensMem = realloc(PopDensMem, CurPopulationNr*MaxStatesAtBirth*PopDensCohortDim*CohortNr*sizeof(double))))
        return ReportMemError("AllocateHeapMemory");
      if (!(CohortLimitMem = realloc(CohortLimitMem, CurPopulationNr*MaxStatesAtBirth*sizeof(double)))) return ReportMemError("AllocateHeapMemory");

      if (!(CohortsMem = realloc(CohortsMem, CurPopulationNr*MaxStatesAtBirth*sizeof(int)))) return ReportMemError("AllocateHeapMemory");

      popDensMem   = PopDensMem;
      birthStateNr = BirthStateNr;
      cohortsMem   = CohortsMem;

      if (!(BirthStateMem = realloc(BirthStateMem, CurPopulationNr*MaxStatesAtBirth*IStateDim*sizeof(double))))
        return ReportMemError("AllocateHeapMemory");
      birthStateMem = BirthStateMem;
#endif

      LastMemAllocated = CurPopulationNr*MaxStatesAtBirth;
    }

  return SUCCES;
}


int FreeHeapMemory(void)

{
  if (RightEigenvecMem) free(RightEigenvecMem);
  RightEigenvecMem  = NULL;
#if (PSPMDEMO == 1)
  if (LeftEigenvecMem) free(LeftEigenvecMem);
  LeftEigenvecMem   = NULL;
#endif
  rightEigenvecMem  = NULL;

#if (FULLSTATEOUTPUT > 0)
  if (PopDensMem) free(PopDensMem);
  PopDensMem        = NULL;
  if (CohortLimitMem) free(CohortLimitMem);
  CohortLimitMem    = NULL;

  if (CohortsMem) free(CohortsMem);
  CohortsMem        = NULL;

  popDensMem        = NULL;
  birthStateNr      = NULL;
  cohortsMem        = NULL;

#if (PSPMIND != 1)
  if (BirthStateMem) free(BirthStateMem);
  BirthStateMem     = NULL;
  birthStateMem     = NULL;
#endif
#endif

  LastMemAllocated  = 0;

  return SUCCES;
}


/*==================================================================================================================================*/
