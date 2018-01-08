/***
  NAME
    dopri5.h
  DESCRIPTION
    Header file with implementation of DOPRI5 integration method.
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

    Last modification: AMdR - Dec 15, 2017
***/

static double LocateSwitch(const int totalOdeDim, const int *birthStateNr, const int curBirthState, int populationStart[PopulationNr],
                           int lifeStage[PopulationNr], double *bState[PopulationNr], double cohortBound[PopulationNr], double *ystart,
                           double *yend, double *ymid, double *rc)

{
  int    jj, iter;
  double a, b, c = 0.0, d = 0.0, e = 0.0, min1, min2;
  double fa, fb, fc, p, q, r, s, tol1, xm;

  a  = 0.0;
  fa = SetMaxLimit(birthStateNr, curBirthState, populationStart, lifeStage, bState, cohortBound, ystart);
  b  = 1.0;
  fb = SetMaxLimit(birthStateNr, curBirthState, populationStart, lifeStage, bState, cohortBound, yend);

  if (fb*fa > 0.0) return -1.0;

  fc = fb;
  for (iter = 0; iter < MAXITER; iter++)
    {
      if (fb*fc > 0.0)
        {
          c  = a;
          fc = fa;
          e = d = b - a;
        }
      if (fabs(fc) < fabs(fb))
        {
          a  = b;
          fa = fb;
          b  = c;
          fb = fc;
          c  = a;
          fc = fa;
        }
      tol1 = 2.0*epsMach*fabs(b) + 0.5*Odesolve_Func_Tol;
      xm   = 0.5*(c - b);

      if ((fabs(xm) <= tol1 && fb*fc <= 0.0) || fb == 0.0) return b;

      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
          s = fb/fa;
          if (a == c)
            {
              p = 2.0*xm*s;
              q = 1.0 - s;
            }
          else
            {
              q = fa/fc;
              r = fb/fc;
              p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
              q = (q - 1.0)*(r - 1.0)*(s - 1.0);
            }
          if (p > 0.0) q = -q;
          p              = fabs(p);
          min1           = 3.0*xm*q - fabs(tol1*q);
          min2           = fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2))
            {
              e = d;
              d = p/q;
            }
          else
            {
              d = xm;
              e = d;
            }
        }
      else
        {
          d = xm;
          e = d;
        }
      a  = b;
      fa = fb;
      if (fabs(d) > tol1)
        b += d;
      else
        b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
      for (jj = 0; jj < totalOdeDim; jj++)
        ymid[jj] =
            rc[jj] + b*(rc[totalOdeDim + jj] + (1.0 - b)*(rc[2*totalOdeDim + jj] + b*(rc[3*totalOdeDim + jj] + (1.0 - b)*rc[4*totalOdeDim + jj])));
      fb = SetMaxLimit(birthStateNr, curBirthState, populationStart, lifeStage, bState, cohortBound, ymid);
    }

  return -1.0;
}


/*
 *====================================================================================================================================
 * Macro definitions needed in integration routine (system independent)
 *====================================================================================================================================
 */

#ifndef BETAFAC
#define BETAFAC                   0.04
#endif
#ifndef FAC1
#define FAC1                      0.2
#endif
#ifndef FAC2
#define FAC2                      10.0
#endif
#ifndef SAFETY
#define SAFETY                    0.9
#endif
#ifndef FACOLD
#define FACOLD                    1.0E-4
#endif


#define Istate(p, i)              (*(ynew + PopulationIndex[p] + (i)))
#define IstatePnt(p, i)           (ynew + PopulationIndex[p] + (i))
#define rcont(i, j)               (rcontpnt[(i)*TotalOdeDim + (j)])


int LifeHistory(const int *birthStateNr, const int curBirthState, const int maxCohortDim, double *iStateOut)
{
  int     i, jj, p, TotalOdeDim, retval;
  double  limitval;

  int     LifeStage[PopulationNr];
  int     PopulationStart[PopulationNr], PopulationIndex[PopulationNr];

  double  Bstates[PopulationNr][IStateDim];
  double  *BstatePnt[PopulationNr];
  double  CohortBound[PopulationNr], CohortLims[PopulationNr];

  int     keepgoing = 1, adjust_step = 1;
  double  tval, tnext, dt, dtfrac, fac, facold = FACOLD, fac11;
  double  *yold, *ytmp, *ynew;
  double  *k0, *k1, *k2, *k3, *k4, *k5, *k6, *scratch;
  double  *rcontpnt;

  double  sqr, err;

  static const double a10 = 0.2, a20 = 3.0/40.0, a21 = 9.0/40.0, a30 = 44.0/45.0, a31 = -56.0/15.0, a32 = 32.0/9.0, a40 = 19372.0/6561.0,
                      a41 = -25360.0/2187.0, a42 = 64448.0/6561.0, a43 = -212.0/729.0, a50 = 9017.0/3168.0, a51 = -355.0/33.0, a52 = 46732.0/5247.0,
                      a53 = 49.0/176.0, a54 = -5103.0/18656.0, a60 = 35.0/384.0, a62 = 500.0/1113.0, a63 = 125.0/192.0, a64 = -2187.0/6784.0,
                      a65 = 11.0/84.0, c1 = 0.2, c2 = 0.3, c3 = 0.8, c4 = 8.0/9.0, d0 = -12715105075.0/11282082432.0,
                      d2  = 87487479700.0/32700410799.0, d3 = -10690763975.0/1880347072.0, d4 = 701980252875.0/199316789632.0,
                      d5  = -1453857185.0/822651844.0, d6 = 69997945.0/29380423.0, e0 = 71.0/57600.0, e2 = -71.0/16695.0, e3 = 71.0/1920.0,
                      e4  = -17253.0/339200.0, e5 = 22.0/525.0, e6 = -1.0/40.0;
#if ((PSPMDEMO != 1) && (PSPMIND != 1))
  double  tmpval;
#endif
#if (PULSED == 1)
  double  t_next_repro = REPRODUCTION_INTERVAL;
  double *iStatePnt[PopulationNr];
  double *fecundPnt[PopulationNr];
  double  survival, cumrepro;
#elif (PSPMIND == 1)
  int CohortStored = 0;
#endif

  memset(LifeStage, 0, PopulationNr*sizeof(int));
  memset(CohortLims, 0, PopulationNr*sizeof(double));

  for (p = 0, TotalOdeDim = 0; p < PopulationNr; p++)
    {
      if (curBirthState < birthStateNr[p])
        {
          PopulationIndex[p] = PopulationStart[p] = TotalOdeDim;
          TotalOdeDim += CohortDim + birthStateNr[p];

#if ((PSPMEQUI == 1) || (PSPMEVODYN == 1))
          TotalOdeDim += InteractDim;
#if (FULLSTATEOUTPUT > 0)
          if (DoStateOutput) TotalOdeDim += CohortDim;
#endif
#endif
        }
      else
        {
          PopulationIndex[p] = PopulationStart[p] = -1;
          LifeStage[p]                            = Stages;
        }
    }
  if (!(ynew = calloc(15*TotalOdeDim + maxCohortDim, sizeof(double))))
    {
      ReportMemError("AllocateHeapMemory");
      return FAILURE;
    }
  yold     = ynew + TotalOdeDim;
  ytmp     = yold + TotalOdeDim;
  k0       = ytmp + TotalOdeDim;
  k1       = k0 + TotalOdeDim;
  k2       = k1 + TotalOdeDim;
  k3       = k2 + TotalOdeDim;
  k4       = k3 + TotalOdeDim;
  k5       = k4 + TotalOdeDim;
  k6       = k5 + TotalOdeDim;
  rcontpnt = k6 + TotalOdeDim;
  scratch  = rcontpnt + 5*TotalOdeDim;

  for (p = 0; p < PopulationNr; p++)
    {
      BstatePnt[p] = Bstates[p];
      memset(BstatePnt[p], 0, IStateDim*sizeof(double));
    }

  StateAtBirth(BstatePnt, curBirthState, Evar);

  for (p = 0; p < PopulationNr; p++)
    {
      if (PopulationStart[p] < 0) continue;
      memcpy(IstatePnt(p, 0), BstatePnt[p], IStateDim*sizeof(double));

#if (FULLSTATEOUTPUT > 0)
      memcpy(BirthStatePnt(p, curBirthState, 0), BstatePnt[p], IStateDim*sizeof(double));
      if (DoStateOutput)
        {
          CohortLims[p] = CohortLimit(curBirthState, p);
          Cohorts(p, curBirthState) = 0;

#if (FULLSTATEOUTPUT == 1)
          CohortBound[p] = CohortMin[p];
#if ((PSPMEQUI == 1) || (PSPMEVODYN == 1))
          CohortBound[p] += CohortLims[p];
#endif
          while (CohortBound[p] < Bstates[p][SortIndex] - DYTOL)
            {
              Cohorts(p, curBirthState)++;
              Cohorts(p, curBirthState) = min(Cohorts(p, curBirthState), COHORT_NR);
              CohortBound[p] += CohortLims[p];
            }
#else
          CohortBound[p] = Bstates[p][SortIndex];
#if ((PSPMEQUI == 1) || (PSPMEVODYN == 1))
          CohortBound[p] += CohortLims[p];
#endif
#endif

#if (PSPMDEMO == 1)
          if (fabs(Istate(p, SortIndex) - CohortBound[p]) <= DYTOL)
            {
              for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, i);

              PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(Istate(p, IStateDim));

              // Store the cumulative number of offspring in all birth states
              for (jj = 0; jj < birthStateNr[p]; jj++)
                PopDens(p, curBirthState, (IStateDim + 1) + jj, Cohorts(p, curBirthState)) = Istate(p, IStateDim + 1 + jj);
              Cohorts(p, curBirthState)++;
              CohortBound[p] = Bstates[p][SortIndex] + CohortLims[p];
            }
#elif (PSPMIND == 1)
          for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, i);

          PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(Istate(p, IStateDim));

          // Store the cumulative contribution to the impact
          for (i = 0, jj = IStateDim + 1 + birthStateNr[p]; i < InteractDim; i++, jj++)
            PopDens(p, curBirthState, (IStateDim + 1) + i, Cohorts(p, curBirthState)) = Istate(p, jj);

          // Store the cumulative number of offspring in all birth states
          for (i = 0, jj = IStateDim + 1; i < birthStateNr[p]; i++, jj++)
            PopDens(p, curBirthState, CohortDim + i, Cohorts(p, curBirthState)) = Istate(p, jj);
          Cohorts(p, curBirthState)++;
          CohortBound[p] += CohortLims[p];
          CohortStored = 1;
#endif
        }
#endif

      if (TestRun)
        {
          STDOUT("\nPop. #%2d - Bstate %2d - Stage %2d (Start):", p, curBirthState, LifeStage[p]);
          for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", Istate(p, i));
          STDOUT("%15.6G", exp(Istate(p, IStateDim)));
          STDOUT("%15.6G", ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1));

#if (PSPMDEMO == 1)
          if (ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1) > 0)
            STDOUT("%15.6G", Istate(p, (IStateDim + 1) + birthStateNr[p])/ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1));
          else
            STDOUT("%15s", "Undefined");
#else
          for (i   = 0; i < InteractDim; i++) STDOUT("%15.6G", Istate(p, (IStateDim + 1) + birthStateNr[p] + i));
#endif

#if (defined(R_PACKAGE))
          R_FlushConsole();
          R_ProcessEvents();
#endif
        }
    }

  // Set the initial life stage values
  SwitchLifeStage(birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, ynew, CohortBound, CohortLims, 1);

  tval = 0.0;
  memcpy(yold, ynew, TotalOdeDim*sizeof(double));
  dt = Odesolve_Init_Step;
  for (p = 0, keepgoing = 0; p < PopulationNr; p++) keepgoing = keepgoing || (PopulationStart[p] >= 0);
  while (keepgoing && (tval < MAX_AGE))
    {
#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt())
        {
          free(ynew);
          return FAILURE;
        }
#endif
      tnext = Odesolve_Fixed_Step*(1 + floor((1 + Odesolve_Rel_Err)*tval/Odesolve_Fixed_Step));
      if (tnext < tval + dt)
        {
          dt          = tnext - tval;
          adjust_step = 0;
        }

#if (PULSED == 1)
      if (t_next_repro < tval + dt)
        {
          dt          = t_next_repro - tval;
          adjust_step = 0;
        }
#elif (PSPMIND == 1)
      CohortStored = 0;
#endif

      /*==========================    Stage 1-5 of DOPRI5   ====================================================*/

      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval, ynew, k0, scratch);

      for (jj = 0; jj < TotalOdeDim; jj++) ynew[jj] = yold[jj] + dt*a10*k0[jj];

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif
      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval + c1*dt, ynew, k1, scratch);

      for (jj = 0; jj < TotalOdeDim; jj++) ynew[jj] = yold[jj] + dt*(a20*k0[jj] + a21*k1[jj]);

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif
      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval + c2*dt, ynew, k2, scratch);

      for (jj = 0; jj < TotalOdeDim; jj++) ynew[jj] = yold[jj] + dt*(a30*k0[jj] + a31*k1[jj] + a32*k2[jj]);

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif
      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval + c3*dt, ynew, k3, scratch);

      for (jj = 0; jj < TotalOdeDim; jj++) ynew[jj] = yold[jj] + dt*(a40*k0[jj] + a41*k1[jj] + a42*k2[jj] + a43*k3[jj]);

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif
      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval + c4*dt, ynew, k4, scratch);

      for (jj = 0; jj < TotalOdeDim; jj++) ynew[jj] = yold[jj] + dt*(a50*k0[jj] + a51*k1[jj] + a52*k2[jj] + a53*k3[jj] + a54*k4[jj]);

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif
      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval + dt, ynew, k5, scratch);

      for (jj = 0; jj < TotalOdeDim; jj++) ynew[jj] = yold[jj] + dt*(a60*k0[jj] + a62*k2[jj] + a63*k3[jj] + a64*k4[jj] + a65*k5[jj]);

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif
      SetDerivatives(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, tval + dt, ynew, k6, scratch);

      err = 0.0;
      for (jj = 0; jj < TotalOdeDim; jj++)
        {
          sqr = dt*(e0*k0[jj] + e2*k2[jj] + e3*k3[jj] + e4*k4[jj] + e5*k5[jj] + e6*k6[jj]);
          sqr /= Odesolve_Abs_Err + Odesolve_Rel_Err*max(fabs(yold[jj]), fabs(ynew[jj]));
          err += sqr*sqr;
        }
      err = sqrt(err/(double)TotalOdeDim);

#if ((defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)) && (OPENMP != 1))
      if (checkInterrupt()) return FAILURE;
#endif

      if (err > 1.0)                                                                // Step rejected
        {
          if (dt <= Odesolve_Min_Step)                                              // Minimum step has failed
            {
              ErrorMsg(__FILE__, __LINE__, "At t = %f the error test failed repeatedly with the minimum step size.", tval);
              free(ynew);
              return FAILURE;
            }
          memcpy(ynew, yold, TotalOdeDim*sizeof(double));
          // If bigger than accuracy take smaller step and restart
          fac11 = pow(err, 0.2 - BETAFAC*0.75);
          dt /= min(1.0/FAC1, fac11/SAFETY);
          dt          = max(dt, Odesolve_Min_Step);
          adjust_step = 0;
          continue;                                                                 // Skip the remainder of the while (keepgoing) {} loop
        }

      // Integration step has been successful. Check whether some event has occurred and process them
      limitval = SetMaxLimit(birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, CohortBound, ynew);
      if (limitval > 0)                                                             // Switching has occurred
        {
          // Update variables for event location and continuous output
          for (jj = 0; jj < TotalOdeDim; jj++)
            {
              rcont(0, jj) = yold[jj];
              rcont(1, jj) = ynew[jj] - yold[jj];
              rcont(2, jj) = dt*k0[jj] - rcont(1, jj);
              rcont(3, jj) = -dt*k6[jj] + rcont(1, jj) - rcont(2, jj);
              rcont(4, jj) = dt*(d0*k0[jj] + d2*k2[jj] + d3*k3[jj] + d4*k4[jj] + d5*k5[jj] + d6*k6[jj]);
            }

          dtfrac =
              LocateSwitch(TotalOdeDim, birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, CohortBound, yold, ynew, ytmp, rcontpnt);
          if (dtfrac < 0.0)
            {
              ErrorMsg(__FILE__, __LINE__, "Root location in LocateSwitch() failed. Integration aborted.");
              free(ynew);
              return FAILURE;
            }

          tval += dtfrac*dt;

          for (jj    = 0; jj < TotalOdeDim; jj++)
            ynew[jj] = rcont(0, jj) + dtfrac*(rcont(1, jj) + (1.0 - dtfrac)*(rcont(2, jj) + dtfrac*(rcont(3, jj) + (1.0 - dtfrac)*rcont(4, jj))));

          retval = SwitchLifeStage(birthStateNr, curBirthState, PopulationStart, LifeStage, BstatePnt, ynew, CohortBound, CohortLims, 0);
          if (retval != SUCCES)
            {
              free(ynew);
              return FAILURE;
            }

#if (FULLSTATEOUTPUT > 0)
          if (DoStateOutput)
            {
              for (p = 0; p < PopulationNr; p++)
                {
                  if (PopulationIndex[p] < 0) continue;
                  if (fabs(Istate(p, SortIndex) - CohortBound[p]) < DYTOL)
                    {
#if (PSPMDEMO == 1)
                      for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, i);

                      PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(Istate(p, IStateDim));

                      // Store the cumulative number of offspring in all birth states
                      for (jj = 0; jj < birthStateNr[p]; jj++)
                        PopDens(p, curBirthState, (IStateDim + 1) + jj, Cohorts(p, curBirthState)) = Istate(p, IStateDim + 1 + jj);
                      Cohorts(p, curBirthState)++;
                      Cohorts(p, curBirthState) = min(Cohorts(p, curBirthState), COHORT_NR - 1);
                      CohortBound[p] += CohortLims[p];
#elif (PSPMIND == 1)
                      for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, i);

                      PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(Istate(p, IStateDim));

                      // Store the cumulative contribution to the impact
                      for (i = 0, jj = IStateDim + 1 + birthStateNr[p]; i < InteractDim; i++, jj++)
                        PopDens(p, curBirthState, (IStateDim + 1) + i, Cohorts(p, curBirthState)) = Istate(p, jj);

                      // Store the cumulative number of offspring in all birth states
                      for (i = 0, jj = IStateDim + 1; i < birthStateNr[p]; i++, jj++)
                        PopDens(p, curBirthState, CohortDim + i, Cohorts(p, curBirthState)) = Istate(p, jj);
                      Cohorts(p, curBirthState)++;
                      CohortBound[p] += CohortLims[p];
                      CohortStored = 1;
#else
                      jj = CohortDim + birthStateNr[p] + InteractDim;
                      if (Istate(p, jj + IStateDim) > MIN_SURVIVAL*CohortLims[p])
                        {
                          if (Cohorts(p, curBirthState) == COHORT_NR)
                            Istate(p, jj + IStateDim) += PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState));
                          for (i = 0; i < IStateDim; i++)
                            {
                              if (Cohorts(p, curBirthState) == COHORT_NR)
                                Istate(p, jj + i) += PopDens(p, curBirthState, i, Cohorts(p, curBirthState))*
                                                     PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState));
                              PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, jj + i)/Istate(p, jj + IStateDim);
                              Istate(p, jj + i) = 0.0;
                            }
                          PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = Beq[p]*Istate(p, jj + IStateDim);
                          Istate(p, jj + IStateDim) = 0.0;

                          Cohorts(p, curBirthState)++;
                          Cohorts(p, curBirthState) = min(Cohorts(p, curBirthState), COHORT_NR);
                        }
                      CohortBound[p] += CohortLims[p];
#endif
                    }
                  else if (Istate(p, SortIndex) > CohortBound[p])
                    {
                      ErrorMsg(__FILE__, __LINE__, "Integration failed to stop at I-state = %.2f! (I-state = %.2f)", CohortBound[p],
                               Istate(p, SortIndex));
                      free(ynew);
                      return FAILURE;
                    }
                }
            }
#endif
        }
      else                                                                          // No switching has occurred
        {
          tval += dt;

          // Step size adjustment using Lund-stabilization: we require fac1 <=  hnew/h <= fac2. No increase if failed just before
          if (adjust_step)
            {
              fac11 = pow(err, 0.2 - BETAFAC*0.75);
              fac   = fac11/pow(facold, BETAFAC);
              dt /= max(1.0/FAC2, min(1.0/FAC1, fac/SAFETY));
              facold = max(err, FACOLD);
              dt     = min(dt, Odesolve_Max_Step);
              dt     = max(dt, Odesolve_Min_Step);
            }
        }

#if (PULSED == 1)
      if ((fabs(tval - t_next_repro) <= Odesolve_Min_Step))
        {
          (void)memcpy(ytmp, ynew, TotalOdeDim*sizeof(double));                     // Make a copy to protect the actual values
          for (p = 0; p < PopulationNr; p++)
            {
              iStatePnt[p] = (PopulationStart[p] < 0) ? Bstates[p] : (ytmp + PopulationStart[p]);
              fecundPnt[p] = (PopulationStart[p] < 0) ? scratch : (ytmp + PopulationStart[p] + (IStateDim + 1));
            }
          Fecundity(LifeStage, iStatePnt, BstatePnt, curBirthState, Evar, fecundPnt);
          (void)memcpy(ynew, ytmp, IStateDim*sizeof(double));                     // Copy the i-state variables back in case they have changed
          for (p = 0; p < PopulationNr; p++)
            {
              if (PopulationStart[p] < 0) continue;

              survival = exp(Istate(p, IStateDim));
              for (jj = 0, cumrepro = 0.0; jj < birthStateNr[p]; jj++)
                {
                  Istate(p, (IStateDim + 1) + jj) += fecundPnt[p][jj]*survival;   // Cumulative reproduction
                  cumrepro += tval*(fecundPnt[p][jj]*survival);
                }
              Istate(p, (IStateDim + 1) + birthStateNr[p]) += cumrepro;           // Generation time, last element

              if (TestRun)
                {
                  STDOUT("\nPop. #%2d - Bstate %2d - Repro. t = %5.1f:", p, curBirthState, tval);
                  for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", Istate(p, i));
                  STDOUT("%15.6G", exp(Istate(p, IStateDim)));
                  STDOUT("%15.6G", ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1));
                  if (ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1) > 0)
                    STDOUT("%15.6G", Istate(p, (IStateDim + 1) + birthStateNr[p])/ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1));
                  else
                    STDOUT("%15s", "Undefined");
#if (defined(R_PACKAGE))
                  R_FlushConsole();
                  R_ProcessEvents();
#endif
                }
            }

          t_next_repro += REPRODUCTION_INTERVAL;
          t_next_repro = min(MAX_AGE, t_next_repro);
        }
#endif

      (void)memcpy(yold, ynew, TotalOdeDim*sizeof(double));                         // Update all variables
      adjust_step = 1;

      for (p = 0, keepgoing = 0; p < PopulationNr; p++) keepgoing = keepgoing || (PopulationStart[p] >= 0);
    }                                                                               // while (keepgoing) {

  if (TestRun)
    {
      for (p = 0; p < PopulationNr; p++)
        {
          if (PopulationIndex[p] < 0) continue;
          STDOUT("\n===> Pop. #%2d - Bstate %2d (End of life):", p, curBirthState);
          for (i = 0; i < IStateDim; i++) STDOUT("%15.6G", Istate(p, i));
          STDOUT("%15.6G", exp(Istate(p, IStateDim)));
          STDOUT("%15.6G", ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1));

#if (PSPMDEMO == 1)
          if (ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1) > 0)
            STDOUT("%15.6G", Istate(p, (IStateDim + 1) + birthStateNr[p])/ASUM(birthStateNr[p], IstatePnt(p, (IStateDim + 1)), 1));
          else
            STDOUT("%15s", "Undefined");
#else
          for (i = 0; i < InteractDim; i++) STDOUT("%15.6G", Istate(p, (IStateDim + 1) + birthStateNr[p] + i));
#endif
        }
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

// Always add the final cohorts, as long as the final density is non-zero
#if (FULLSTATEOUTPUT > 0)
  if (DoStateOutput)
    {
      for (p = 0; p < PopulationNr; p++)
        {
#if (PSPMDEMO == 1)
          if (PopulationIndex[p] < 0) continue;

          for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, i);

          PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(Istate(p, IStateDim));

          // Store the cumulative number of offspring in all birth states
          for (jj = 0; jj < birthStateNr[p]; jj++)
            PopDens(p, curBirthState, (IStateDim + 1) + jj, Cohorts(p, curBirthState)) = Istate(p, IStateDim + 1 + jj);
          Cohorts(p, curBirthState)++;
#elif (PSPMIND == 1)
          if (CohortStored) continue;
          if ((PopulationIndex[p] < 0) || (Cohorts(p, curBirthState) >= CohortNr)) continue;

          for (i = 0; i < IStateDim; i++) PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, i);

          PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = exp(Istate(p, IStateDim));

          // Store the cumulative contribution to the impact
          for (i = 0, jj = IStateDim + 1 + birthStateNr[p]; i < InteractDim; i++, jj++)
            PopDens(p, curBirthState, (IStateDim + 1) + i, Cohorts(p, curBirthState)) = Istate(p, jj);

          // Store the cumulative number of offspring in all birth states
          for (i = 0, jj = IStateDim + 1; i < birthStateNr[p]; i++, jj++)
            PopDens(p, curBirthState, CohortDim + i, Cohorts(p, curBirthState)) = Istate(p, jj);
          Cohorts(p, curBirthState)++;
          CohortBound[p] += CohortLims[p];
#else
          if (PopulationIndex[p] < 0) continue;

          jj = CohortDim + birthStateNr[p] + InteractDim;
          tmpval = Beq[p]*Istate(p, jj + IStateDim);
          if (tmpval)
            {
              if (Cohorts(p, curBirthState) == COHORT_NR)
                Istate(p, jj + IStateDim) += PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState));
              for (i = 0; i < IStateDim; i++)
                {
                  if (Cohorts(p, curBirthState) == COHORT_NR)
                    Istate(p, jj + i) +=
                        PopDens(p, curBirthState, i, Cohorts(p, curBirthState))*PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState));
                  PopDens(p, curBirthState, i, Cohorts(p, curBirthState)) = Istate(p, jj + i)/Istate(p, jj + IStateDim);
                }
              PopDens(p, curBirthState, IStateDim, Cohorts(p, curBirthState)) = Beq[p]*Istate(p, jj + IStateDim);
              Cohorts(p, curBirthState)++;
            }
#endif
        }
    }
#endif

  // Now copy the locally stored values back to the global memory
  for (p = 0; p < PopulationNr; p++)
    {
      if (PopulationIndex[p] < 0) continue;
#if (PSPMDEMO == 1)
      memcpy(iStateOut + p*maxCohortDim, ynew + PopulationIndex[p], (CohortDim + birthStateNr[p])*sizeof(double));
#elif (PSPMIND == 1)
      memcpy(iStateOut + p*maxCohortDim, ynew + PopulationIndex[p], (CohortDim + birthStateNr[p])*sizeof(double));
      CohortLimit(curBirthState, p) = (Istate(p, SortIndex) - Bstates[p][SortIndex])/COHORT_NR;
#else
      memcpy(iStateOut + p*maxCohortDim, ynew + PopulationIndex[p], (CohortDim + birthStateNr[p] + InteractDim)*sizeof(double));
#endif
    }

  free(ynew);
  return SUCCES;
}


/*==================================================================================================================================*/
