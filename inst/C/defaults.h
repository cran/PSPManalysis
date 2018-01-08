/***
  NAME
    defaults
  DESCRIPTION
    Header file with global #define statements

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

    Last modification: AMdR - Jun 07, 2017
***/
#ifndef DEFAULTS_H_
#define DEFAULTS_H_

#if !defined(POPULATION_NR) || (POPULATION_NR < 1)
#error POPULATION_NR should be defined larger than 0
#endif

#if !defined(I_STATE_DIM) || (I_STATE_DIM < 1)
#error I_STATE_DIM should be defined larger than 0
#endif

#if !defined(FULLSTATEOUTPUT) || (FULLSTATEOUTPUT < 0) || (FULLSTATEOUTPUT > 2)
#undef FULLSTATEOUTPUT
// 0: No state output at all, 1: Condensed state output, 2: State output for each birth state
#define FULLSTATEOUTPUT           2
#endif

#ifndef MAX_AGE
#define MAX_AGE                   1.0E6
#endif

#ifndef MIN_SURVIVAL
#define MIN_SURVIVAL              1.0E-9
#endif

#ifndef COHORT_NR
#define COHORT_NR                 100
#endif
#ifndef DYTOL
#define DYTOL                     1.0E-7
#endif
#ifndef RHSTOL
#define RHSTOL                    1.0E-8
#endif
#ifndef ESSTOL
#define ESSTOL                    1.0E-8
#endif

#ifndef ALLOWNEGATIVE
#define ALLOWNEGATIVE             0                                                 // Allow negative solutions?
#endif
#define MAX_STEPREDUCE            2048
#define MICRO                     1.0E-6


#ifndef ODESOLVE_INIT_STEP
#define ODESOLVE_INIT_STEP        0.1
#endif
#ifndef ODESOLVE_FIXED_STEP
#define ODESOLVE_FIXED_STEP       MAX_AGE
#endif
#ifndef ODESOLVE_MIN_STEP
#define ODESOLVE_MIN_STEP         1.0E-8
#endif
#ifndef ODESOLVE_MAX_STEP
#define ODESOLVE_MAX_STEP         10.0
#endif
#ifndef ODESOLVE_ABS_ERR
#define ODESOLVE_ABS_ERR          1.0E-10
#endif
#ifndef ODESOLVE_REL_ERR
#define ODESOLVE_REL_ERR          1.0E-8
#endif
#ifndef ODESOLVE_FUNC_TOL
#define ODESOLVE_FUNC_TOL         1.0E-8
#endif


#ifndef JACOBIAN_MIN_STEP
#define JACOBIAN_MIN_STEP         1.0E-6
#endif
#ifndef JACOBIAN_STEP
#define JACOBIAN_STEP             1.0E-3                                            // % change for jacobian
#endif
#if !defined(JACOBIAN_UPDATES) || (JACOBIAN_UPDATES < 1)
#define JACOBIAN_UPDATES          5
#endif

#ifndef FASTNUMERICS
#define FASTNUMERICS              1
#endif

#endif /* DEFAULTS_H_ */
