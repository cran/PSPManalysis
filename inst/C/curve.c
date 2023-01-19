/***
  NAME
    curve
  DESCRIPTION
    This module implements routines that are specifically used in locating
    points on the equilibrium branch of a structured population model.

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

    Last modification: AMdR - Jan 19, 2023
***/
#ifndef CURVE
#define CURVE
#endif
#include "globals.h"


/*
 *====================================================================================================================================
 *      Some numerical settings
 *====================================================================================================================================
 */

#ifndef RHSMAX
#define RHSMAX              0.99
#endif
#ifndef CENTRALDIFF
#define CENTRALDIFF         1                                                       // Compute derivative by central difference
#endif
#ifndef FUNCTOL
#define FUNCTOL             1.0E-8
#endif

#define STEP_DOUBLE         4
#define STEP_HALF           (MAXITER/2)

#define FEMTO               1.0E-15

#define SAFETY              0.8
#define MIN_SCALE           1.0E-6
#define MAX_SCALE           1.0E6

#ifdef _MSC_VER
#include <float.h>
#define issane(a)           ((_fpclass(a) == _FPCLASS_NN) || (_fpclass(a) == _FPCLASS_NZ) || (_fpclass(a) == _FPCLASS_PZ) || (_fpclass(a) == _FPCLASS_PN))
#else
#define issane(a)           ((fpclassify(a) == FP_ZERO) || (fpclassify(a) == FP_NORMAL))
#endif

static int                  *oldscale;
static int                  FirstTangent = 1;
static int                  fast_iters = 0, slow_iters = 0;
static int                  LPImmediateReturn  = 0;
static int                  ESSImmediateReturn = 0;
const double                upper = ((2.0 - SAFETY)*10.0), lower = SAFETY;
static int                  pnt_scale_dim = -1;

static double               *dDETBaseMem    = NULL;
static LAPACK_SIZE_T        *iDETBaseMem    = NULL;
static int                  dDETBaseMemDim = 0, iDETBaseMemDim = 0;

static double               *dEVBaseMem    = NULL;
static LAPACK_SIZE_T        *iEVBaseMem    = NULL;
static int                  dEVBaseMemDim = 0, iEVBaseMemDim = 0;

static double               *dSLSBaseMem    = NULL;
static LAPACK_SIZE_T        *iSLSBaseMem    = NULL;
static int                  dSLSBaseMemDim = 0, iSLSBaseMemDim = 0;


/*==================================================================================================================================*/

double anorm(int rows, int cols, double *a)
{
  register int i;
  double       tmp, maxval = 0.0;

  for (i = 0; i < rows; i++)
    {
      tmp    = ASUM(cols, a + i*cols, 1);
      maxval = max(maxval, tmp);
    }

  return maxval;
}


/*==================================================================================================================================*/

int SetScales(double *point, int pntdim)

/*
 * This routine scales the vector of problem variables to within reasonable
 * bounds. This makes them more comparable, which is advantageous for the
 * computations.
 */

{
  register int i, scaleset = 0, newscale;
  double       tmp, val;

  if (!pnt_scale)
    {
      pnt_scale = calloc(pntdim, sizeof(double));
      for (i = 0; i < pntdim; i++) pnt_scale[i] = 1.0;
      pnt_scale_dim                             = pntdim;
    }
  if (!oldscale)
    {
      oldscale = calloc(pntdim, sizeof(int));
      for (i = 0; i < pntdim; i++) oldscale[i] = 0;
    }
  if (!point) return scaleset;

  for (i = 0; i < pntdim; i++)
    {
      val = fabs(point[i]);
      if ((val < lower) || (val > upper))
        {
          tmp      = max(val*pnt_scale[i], MIN_SCALE);
          tmp      = min(tmp, MAX_SCALE);
          tmp      = floor(log10(tmp) + FUNCTOL);
          newscale = (int)tmp;

          if (newscale != oldscale[i])
            {
              tmp = pow(10.0, tmp);
              point[i] *= pnt_scale[i]/tmp;
              pnt_scale[i] = tmp;
              oldscale[i]  = newscale;
              scaleset     = i + 1;
            }
        }
    }

  return scaleset;
}


/*==================================================================================================================================*/

int FindPoint(const int pntdim, double *guess, double *JacImport, double *tanvec, double ytol, double rhstol, const int max_iter,
              int (*fnc)(double *, double *))

/*
 * FindPoint -  Routine locates a point on a curve determined by a
 *              system of non-linear, algebraic equations.
 *              The iteration adjusts the vector-elements following a simple
 *              Newton-Chord method with Broyden update (see Kuznetsov pg. 418).
 *              Pseudo-arclength continuation is used to continue past curve folds.
 *
 * Arguments -  pntdim    : The dimension of the solution point on the curve.
 *                          The dimension of the system of equations is
 *                          assumed to be exactly 1 less.
 *              guess     : Pointer to an array containing the initial point
 *                          to start the iteration from. The first element of
 *                          the vector is assumed to be non-adjustable parameter.
 *              ytol      : Tolerance determining when change in y equals zero.
 *              rhstol    : Tolerance determining when RHS equals zero.
 *              max_iter  : Maximum stepnumber allowed in iteration.
 *              fnc       : Pointer to function specifying the system of
 *                          equations. The function must have a (double)
 *                          pointer as first argument, containing the point
 *                          in which to evaluate the system and a (double)
 *                          pointer as second argument, containing the
 *                          results after evaluation of the equations.
 */

{
  register int  iter, i, j;
  int           rhsdim;
  int           pntdim2 = pntdim*pntdim;
  int           retcode = NO_CONVERGENCE;
  double        ynorm, dynorm, rhsnorm;
  double        *basemem, *y, *tv, *dy, *rhs;
  double        *Jac, *JacCopy;

  y = basemem = calloc((4*pntdim + 2*pntdim2), sizeof(double));
  if (!basemem) return ReportMemError("FindPoint");

  tv      = y + pntdim;
  dy      = tv + pntdim;
  rhs     = dy + pntdim;
  Jac     = rhs + pntdim;
  JacCopy = Jac + pntdim2;

  // If tangent vector is given, we are doing pseudo-arc length continuation
  rhsdim = pntdim - (tanvec != NULL);

  ReportMsg("\nLocating new solution point\n");
  COPY(pntdim, guess, 1, y, 1);
  memset((void *)dy, 0, pntdim*sizeof(double));

  // The iteration loop
  for (iter = 0; iter < max_iter; iter++)
    {
      // Compute norm of Y and of dY
      ynorm  = NRM2(pntdim, y, 1);
      dynorm = NRM2(pntdim, dy, 1);
      if (!issane(ynorm) || !issane(dynorm))
        {
          ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
          retcode = NORM_OVERFLOW;
          break;
        }

      // Compute RHS and its norm
      memset((void *)rhs, 0, pntdim*sizeof(double));
      if ((*fnc)(y, rhs) == FAILURE)
        {
          ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
          retcode = FAILED_EVALUATION;
          break;
        }
      // The dimension of rhs is pntdim, for example in case of BP localisation or PGR calculations
      rhsnorm = NRM2(pntdim, rhs, 1);
      if (!iter) ReportMsg("Start:");
      for (i = 0; i < pntdim; i++) ReportMsg("\t%12.5E", y[i]*pnt_scale[i]);
      ReportMsg("\tdY  norm: %12.5E\tRHS norm: %12.5E\n", dynorm, rhsnorm);

      // Return if converged or diverged
      if ((!issane(rhsnorm)) || (rhsnorm/(1.0 + rhsnorm) > RHSMAX))
        {
          ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
          retcode = NORM_OVERFLOW;
          break;
        }
      // The dimension of rhs is pntdim, for example in case of BP localisation or PGR calculations
      else if ((rhsnorm < pntdim*rhstol) && (dynorm < pntdim*ytol))
        {
          COPY(pntdim, y, 1, guess, 1);

          if ((Stepchange) && (Stepreduce == 1))
            {
              if (iter < STEP_DOUBLE)
                {
                  fast_iters++;
                  slow_iters = 0;
                  if (fast_iters == 2)
                    {
                      curvestep *= 1.5;
                      fast_iters = 0;
                    }
                }
              else if (iter > STEP_HALF)
                {
                  slow_iters++;
                  fast_iters = 0;
                  if (slow_iters == 2)
                    {
                      curvestep *= 0.5;
                      slow_iters = 0;
                    }
                }
              else
                {
                  fast_iters = 0;
                  slow_iters = 0;
                }
            }

          retcode = SUCCES;
          break;
        }

      // Compute Jacobian every Jacobian_Updates steps, otherwise the Jacobian is updated
      // via a Broyden update (see below)
      if ((iter == 0) && JacImport)
        memcpy(Jac, JacImport, pntdim*rhsdim*sizeof(double));
      else if (!(iter % Jacobian_Updates))
        {
          ReportMsg("%-s", "Computing jacobian");
          memset((void *)Jac, 0, (pntdim*pntdim)*sizeof(double));
          /*
        * Notice that the Jacobian is stored as
        *
        *                    |dF1/dy1 ... dFn/dy1|
        *                    |dF1/dy2 ... dFn/dy2|
        *               J =     |   .          .   |
        *                    |   .          .   |
        *                    |dF1/dyn ... dFn/dyn|
        *
        * The matrix is hence stored in column-wise (fortran) style.
        * From a C perspective this means that all coefficients pertaining to yi are to be found
        * in ROW i (as opposed to column i).
        *
        * Solving J.dy = -F(y) with dy = (dy1 ... dyn) and F(y) = (F1(y) ... Fn(y)) requires the variable
        * trans[1] to be defined as {"N"} (see programs/various/testlapack.c for details).
        */
          Jacobian(pntdim, y, rhsdim, Jac, fnc, FORWARD);                           // Compute J = F_x(X^k)
          ReportMsg(".....Ok!\n");
        }
      else                                                                          // Broyden update of Jacobian
        {
          dynorm = DOT(pntdim, dy, 1, dy, 1);
          // See 10.7 on pg. 419 in Kuznetsov. Notice though that Jac is the transposed jacobian
          for (i = 0; i < pntdim; i++)
            for (j = 0; j < rhsdim; j++) Jac[j + i*rhsdim] += rhs[j]*dy[i]/dynorm;
        }

      // Extend the Jacobian matrix to include an additional row for the tangent
      // vector if we are doing pseudo-archlength continuation
      // If tangent is not given rhsdim == pntdim and we follow simple Newton
      memset((void *)JacCopy, 0, (pntdim*pntdim)*sizeof(double));
      for (i = 0; i < pntdim; i++)
        COPY(rhsdim, Jac + i*rhsdim, 1, JacCopy + i*pntdim, 1);                     // Extract dF/dx
      memset((void *)dy, 0, pntdim*sizeof(double));
      AXPY(rhsdim, -1.0, rhs, 1, dy, 1);

      if (tanvec != NULL)
        {
          // When tangent is present, find new point via pseudo-arclength continuation
          for (i = 0; i < pntdim; i++) *(JacCopy + i*pntdim + rhsdim) = tanvec[i];
          COPY(pntdim, y, 1, tv, 1);
          AXPY(pntdim, -1.0, guess, 1, tv, 1);
          dy[rhsdim] = DOT(pntdim, tv, 1, tanvec, 1);
        }

      // Solve the linear system
      retcode = SolveLinearSystem(pntdim, JacCopy, dy, ytol);
      if (retcode != SUCCES) break;

      // Adjust point
      AXPY(pntdim, 1.0, dy, 1, y, 1);
      retcode = NO_CONVERGENCE;
    }

  free(basemem);

  if (retcode == SUCCES)
    ReportMsg("New solution point found\n\n");
  else
    ReportMsg("Locating new solution point failed\n\n");

  return retcode;
}


/*==================================================================================================================================*/

int TangentVec(const int pntdim, double *sol, double *JacExport, double *tanvec, int (*fnc)(double *, double *), double *det, const double ytol)

/*
 * TangentVec - routine determines the direction of the curve defined by the
 *              system of equations
 *
 *                    F(y) = 0
 *
 *              The point y is considered to have a dimension of exactly 1
 *              larger than the number of equations (i.e. the dimension of
 *              F(y)).
 *
 * Arguments -  pntdim  : The dimension of the solution point on the curve.
 *              y       : Pointer to an array containing the fixed point
 *              tanvec  : Pointer to return tangent vector
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 */

{
  register int  j;
  int           rhsdim = pntdim - 1, pntdim2 = pntdim*pntdim, retcode;
  double        norm;
  double        *basemem, *y, *Jac, *JacCopy;

  y = basemem = calloc((pntdim + 2*pntdim2), sizeof(double));
  if (!basemem) return ReportMemError("TangentVec");

  Jac     = y + pntdim;
  JacCopy = Jac + pntdim2;

  // Initialize
  COPY(pntdim, sol, 1, y, 1);
  norm = NRM2(pntdim, y, 1);
  if (!issane(norm))
    {
      ErrorMsg(__FILE__, __LINE__, "Norm overflow in curvedir");
      free(basemem);
      return NORM_OVERFLOW;
    }

  // Determine the Jacobian of the extended system (variable plus parameter
  // dependence).
  ReportMsg("\nComputing curve direction     ");
  Jacobian(pntdim, y, rhsdim, JacCopy, fnc, CENTRALDIFF);
  if (JacExport) memcpy(JacExport, JacCopy, pntdim*rhsdim*sizeof(double));
  ReportMsg(".....Ok!\n");

  // Append the current tangent vector as the last row to the jacobian to
  // preserve direction. See the matcont manual at
  // http://www.matcont.ugent.be/manual.pdf, page 10 & 11
  // Notice, however, it is here added as the last COLUMN because of the
  // Fortran column-wise storage!

  for (j = 0; j < pntdim; j++)
    {
      COPY(rhsdim, JacCopy + j*rhsdim, 1, Jac + j*pntdim, 1);                       // Extract dF/dy
      *(Jac + j*pntdim + rhsdim) = tanvec[j];
    }

  memset((void *)JacCopy, 0, (pntdim*pntdim)*sizeof(double));
  COPY(pntdim2, Jac, 1, JacCopy, 1);
  memset((void *)tanvec, 0, pntdim*sizeof(double));
  tanvec[rhsdim] = 1.0;
  Stepchange     = 0;

  // Solve the linear system
  retcode = SolveLinearSystem(pntdim, JacCopy, tanvec, ytol);
  if (retcode != SUCCES)
    {
      ErrorMsg(__FILE__, __LINE__, "Failed to solve for tangent vector in TangentVec()");
      memset((void *)tanvec, 0, pntdim*sizeof(double));
      tanvec[0] = 1.0;
      free(basemem);
      return retcode;
    }

  if (det)
    {
      // Replace the last row of the (saved) Jacobian with the newly computed tangent vector
      // to compute the determinant for BP detection
      for (j = 0; j < pntdim; j++)
        {
          COPY(rhsdim, Jac + j*pntdim, 1, JacCopy + j*pntdim, 1);
          *(JacCopy + j*pntdim + rhsdim) = tanvec[j];
        }
      Determinant(pntdim, JacCopy, det, NULL);
    }
  norm = NRM2(pntdim, tanvec, 1); /* Normalize and store      */
  SCAL(pntdim, 1.0/norm, tanvec, 1);

  if (FirstTangent && (tanvec[0] < 0.0)) SCAL(pntdim, -1.0, tanvec, 1);
  FirstTangent = 0;

  free(basemem);

  return SUCCES;
}


/*==================================================================================================================================*/

int CentralDerivative(int fncdim, int (*fnc)(double *, double *), double *farg, double *frhs, double *x, double h0, double *feq, double *result,
                      int fast)
{
  // Compute the derivative using the 5-point rule (x-h, x-h/2, x, x+h/2, x+h). Note that the central point is not used.
  // Compute the error using the difference between the 5-point and the 3-point rule (x-h,x,x+h). Again the central point is not used.
  int    i;
  double old, h, r[2], trunc[2], round[2], error[2];
  double fm1, fp1, fmh, fph, r3, r5, e3, e5, dy;

  old =*x;
  h   = h0;
  for (i = 0; i < 2; i++)
    {
      *x = old + h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

      // Trick to reduce precision errors if non-optimized computation is used. See Num. Recipes 9.7, pg. 388.
      if (fast)
        {
          COPY(fncdim, feq, 1, result, 1);
          h =*x - old;
        }
      else
        fp1 =*feq;

      *x = old - h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

      if (fast)
        {
          AXPY(fncdim, -1.0, feq, 1, result, 1);
          SCAL(fncdim, 0.5/h, result, 1);
          *x = old;
          return SUCCES;
        }
      else
        fm1 =*feq;

      *x = old - h/2;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      fmh =*feq;

      *x = old + h/2;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      fph =*feq;

      r3 = 0.5*(fp1 - fm1);
      r5 = (4.0/3.0)*(fph - fmh) - (1.0/3.0)*r3;

      e3 = (fabs(fp1) + fabs(fm1))*DBL_EPSILON;
      e5 = 2.0*(fabs(fph) + fabs(fmh))*DBL_EPSILON + e3;

      // The next term is due to finite precision in x+h = O (eps * x)
      dy = max(fabs(r3/h), fabs(r5/h))*(fabs(*x)/h)*DBL_EPSILON;

      // The truncation error in the r5 approximation itself is O(h^4). However, for safety, we estimate the error from r5-r3,
      // which is O(h^2).  By scaling h we will minimise this estimated error, not the actual truncation error in r5.
      r[i]     = r5/h;
      trunc[i] = fabs((r5 - r3)/h);                                                 // Estimated truncation error O(h^2)
      round[i] = fabs(e5/h) + dy;                                                   // Rounding error (cancellations)
      error[i] = round[i] + trunc[i];

      // Compute an optimised stepsize to minimize the total error, using the scaling of the truncation error (O(h^2)) and
      // rounding error (O(1/h)).
      if ((i == 0) && ((round[0] < trunc[0]) && (round[0] > 0 && trunc[0] > 0)))
        h = (h0)*pow(round[0]/(2.0*trunc[0]), 1.0/3.0);
      else
        break;
    }

  // Check that the new error is smaller, and that the new derivative is consistent with the error bounds of the original estimate.
  if ((i == 1) && (error[1] < error[0]) && (fabs(r[1] - r[0]) < 4.0*error[0]))
    *result = r[1];
  else
    *result = r[0];

  *x = old;
  return SUCCES;
}

static int ForwardDerivative(int fncdim, int (*fnc)(double *, double *), double *farg, double *frhs, double *x, double h0, double *feq,
                             double *result, int fast)
{
  int    i;
  double old, h, r[2], trunc[2], round[2], error[2];
  double f1, f2, f3, f4, r2, r4, e4, dy;

  old =*x;
  h   = h0;
  for (i = 0; i < 2; i++)
    {
      *x = old + h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

      // Trick to reduce precision errors if non-optimized computation is used. See Num. Recipes 9.7, pg. 388.
      if (fast)
        {
          h =*x - old;
          COPY(fncdim, feq, 1, result, 1);
          *x = old;
          if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;

          AXPY(fncdim, -1.0, feq, 1, result, 1);
          SCAL(fncdim, 1.0/h, result, 1);
          return SUCCES;
        }
      else
        f4 =*feq;

      *x = old + (3.0/4.0)*h;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      f3 =*feq;

      *x = old + h/2.0;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      f2 =*feq;

      *x = old + h/4.0;
      if ((*fnc)(farg, frhs) == FAILURE) return FAILURE;
      f1 =*feq;

      r2 = 2.0*(f4 - f2);
      r4 = (22.0/3.0)*(f4 - f3) - (62.0/3.0)*(f3 - f2) + (52.0/3.0)*(f2 - f1);

      // Estimate the rounding error for r4
      e4 = 2*20.67*(fabs(f4) + fabs(f3) + fabs(f2) + fabs(f1))*DBL_EPSILON;

      // The next term is due to finite precision in x+h = O (eps * x)
      dy = max(fabs(r2/h), fabs(r4/h))*fabs(*x/h)*DBL_EPSILON;

      // The truncation error in the r4 approximation itself is O(h^3). However, for safety, we estimate the error from r4-r2,
      // which is O(h).  By scaling h we will minimise this estimated error, not the actual truncation error in r4.
      r[i]     = r4/h;
      trunc[i] = fabs((r4 - r2)/h);                                                 // Estimated truncation error O(h)
      round[i] = fabs(e4/h) + dy;
      error[i] = round[i] + trunc[i];

      // Compute an optimised stepsize to minimize the total error, using the scaling of the estimated truncation error
      // (O(h)) and rounding error (O(1/h)).
      if ((i == 0) && ((round[0] < trunc[0]) && (round[0] > 0 && trunc[0] > 0)))
        h = (h0)*pow(round[0]/(trunc[0]), 1.0/2.0);
      else
        break;
    }
  *x = old;

  // Check that the new error is smaller, and that the new derivative is consistent with the error bounds of the original estimate.
  if ((i == 1) && (error[1] < error[0]) && (fabs(r[1] - r[0]) < 4.0*error[0]))
    *result = r[1];
  else
    *result = r[0];

  return SUCCES;
}


int Jacobian(const int pntdim, double *pnt, const int fncdim, double *jac, int (*fnc)(double *, double *), int method)
/*
 * Routine determines the Jacobian of the n-dimensional function F(y) w.r.t. the m-dimensional
 * variable y at the current point given by 'pnt'. The routine hence returns in 'jac' the
 * following matrix of partial derivatives:
 *
 *                    |dF1/dy1 ... dFn/dy1|
 *                    |   .          .    |
 *               Df = |   .          .    |
 *                    |   .          .    |
 *                    |dF1/dym ... dFn/dym|
 *
 * Notice that all coefficients pertaining to yi are to be found in ROW i (as opposed to column i).
 * The matrix is hence stored in column-wise (fortran) style.
 */


{
  register int  j, k;
  double        *basemem, *y, *rhs, ydif;

  y = basemem = calloc(2*pntdim, sizeof(double));
  if (!basemem) return ReportMemError("Jacobian");

  rhs = y + pntdim;

  // Initialize
  COPY(pntdim, pnt, 1, y, 1);
  memset((void *)jac, 0, (pntdim*fncdim)*sizeof(double));

  for (j = 0; j < pntdim; j++)
    {
      ydif = fabs(Jacobian_Step*y[j]);
      ydif = (j < pnt_scale_dim) ? max(ydif, Jacobian_Min_Step/pnt_scale[j]) : max(ydif, Jacobian_Min_Step);

      if (FastNumerics == 1)
        {
          if (method == FORWARD)
            {
              if (ForwardDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs, jac + j*fncdim, 1) == FAILURE)
                {
                  ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                  free(basemem);
                  return FAILURE;
                }
            }
          else
            {
              if (CentralDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs, jac + j*fncdim, 1) == FAILURE)
                {
                  ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                  free(basemem);
                  return FAILURE;
                }
            }
        }
      else
        {
          for (k = 0; k < fncdim; k++)
            {
              if (method == FORWARD)
                {
                  if (ForwardDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs + k, jac + j*fncdim + k, 0) == FAILURE)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                      free(basemem);
                      return FAILURE;
                    }
                }
              else
                {
                  if (CentralDerivative(fncdim, fnc, y, rhs, y + j, ydif, rhs + k, jac + j*fncdim + k, 0) == FAILURE)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                      free(basemem);
                      return FAILURE;
                    }
                }
            }
        }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      if (checkInterrupt())
        {
          free(basemem);
          return FAILURE;
        }
#endif
    }

  free(basemem);
  return SUCCES;
}


/*==================================================================================================================================*/

int LPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), const int method, const int lpcurve, double *curvedir, const double ytol)

/*
 * LPcondition -  Routine computes the factor determining the location of a limit point, i.e.
 *                the parameter component of the tangent vector. This component always has
 *                index 0 in the vector of the solution point.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the parameters) plus the dimension of the vector of
 *                        state variables in case of lpcurve = 1, otherwise the dimension
 *                        of y equals 1 plus the vector of state variables.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables y. The last element is assumed to be the
 *                        second parameter in the LP continuation
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *              lpcurve : Routine called for detection of LP in EQ or ESS curve (0) or
 *                        during computation of LP curve (1)
 */
{
  register int  j;
  const int     lppntdim = pntdim - lpcurve;
  int           rhsdim = lppntdim - 1, lppntdim2 = lppntdim*lppntdim, retcode;
  int           maxind;
  double        *basemem, *rhs, *Jac, *JacCopy, *tanvec;
  double        norm, maxcond, cond;

  // Prevent recurrence
  *curvedir = 0.0;
  if (LPImmediateReturn) return SUCCES;

  rhs = basemem = calloc((pntdim + pntdim*pntdim + lppntdim*lppntdim + lppntdim), sizeof(double));
  if (!basemem) return ReportMemError("LPcondition");

  LPImmediateReturn = 1;

  Jac     = rhs + pntdim;
  JacCopy = Jac + pntdim*pntdim;
  tanvec  = JacCopy + lppntdim*lppntdim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  // Notice that when continuing a LP curve (lpcurve = 1) we have to call Jacobian() with the
  // full dimension (pntdim) of y to pass the entire argument vector to Equation().
  // As a result, the Jacobian matrix will have one row more than we really need,
  // representing the derivatives w.r.t. to the 2nd parameter of the LP continuation.
  // This additional row, with index lppntdim = pntdim-1 will be ignored.
  // Also notice that the last column will be unassigned, as the systems of equations returned
  // only has size rhsdim = pntdim-lpcurve-1. This last column is hence explicitly set to 0
  Jacobian(pntdim, y, lppntdim, Jac, fnc, method);
  for (j = 0; j < lppntdim; j++) Jac[j*lppntdim + rhsdim] = 0.0;

  /*
    * When lpcurve = 1 the Jacobian equals the following (n+2)x(n+1) matrix of partial
    * derivatives:
    *
    *                    |dF1/dp1 ... dFn/dp1  0|
    *                    |dF1/dy1 ... dFn/dy1  0|
    *                    |   .          .      0|
    *               Df = |   .          .      0|
    *                    |   .          .      0|
    *                    |dF1/dyn ... dFn/dyn  0|
    *                    |dF1/dp2 ... dFn/dp2  0|
    *
    * Otherwise, when lpcurve = 0 the (n+1)x(n+1) matrix
    *
    *
    *                    |dF1/dp1 ... dFn/dp1  0|
    *                    |dF1/dy1 ... dFn/dy1  0|
    *                    |   .          .      0|
    *               Df = |   .          .      0|
    *                    |   .          .      0|
    *                    |dF1/dyn ... dFn/dyn  0|
    *
    * In which n = pntdim-lpcurve-1 (i.e. equal to the number of state variables).
    * Notice that all coefficients pertaining to yi are to be found in ROW i (as 
    * opposed to column i). The matrix is hence stored in column-wise (fortran) style.
    */

  // Additional call to reset global variables
  (*fnc)(y, rhs);

  // Find the most non-singular matrix (largest inverse condition)
  for (j = 0, maxcond = 0.0, maxind = -1; j < lppntdim; j++)
    {
      // LU decompose the matrix and compute determinant
      COPY(lppntdim2, Jac, 1, JacCopy, 1);
      JacCopy[j*lppntdim + rhsdim] = 1.0;

      if (Determinant(lppntdim, JacCopy, NULL, &cond) == SUCCES)
        {
          if (cond > maxcond + FEMTO)
            {
              maxcond = cond;
              maxind  = j;
            }
        }
    }

  if (maxind == -1)
    {
      ErrorMsg(__FILE__, __LINE__, "No non-singular matrix found in LPcondition()");
      free(basemem);
      LPImmediateReturn = 0;
      return FAILURE;
    }

  Jac[maxind*lppntdim + rhsdim] = 1.0;
  COPY(lppntdim2, Jac, 1, JacCopy, 1);
  memset((void *)tanvec, 0, lppntdim*sizeof(double));
  tanvec[rhsdim] = 1.0;

  // Solve the linear system
  retcode = SolveLinearSystem(lppntdim, JacCopy, tanvec, ytol);
  if (retcode != SUCCES)
    {
      ErrorMsg(__FILE__, __LINE__, "Failed to solve for tangent vector in LPcondition()");
      free(basemem);
      LPImmediateReturn = 0;
      return retcode;
    }

  norm = NRM2(lppntdim, tanvec, 1); /* Normalize and store      */
  SCAL(lppntdim, 1.0/norm, tanvec, 1);

  *curvedir = tanvec[0];
  free(basemem);
  LPImmediateReturn = 0;

  return SUCCES;
}


/*==================================================================================================================================*/

int Determinant(const int N, double *M, double *det, double *cond)

{
  char          whichnorm;
  int           j;
  int           retval = SUCCES, memneeded;
  double        *A, *work, norm;
  LAPACK_SIZE_T nc = N, lwork = 4*N, *ipiv, *iwork, liwork = N, info;

  // Allocate temporarily minimally allowed size for workspace arrays
  memneeded = N*N + lwork;
  if (memneeded > dDETBaseMemDim)
    {
      dDETBaseMem = realloc(dDETBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dDETBaseMem == NULL)
        {
          if (dDETBaseMem) free(dDETBaseMem);
          dDETBaseMem    = NULL;
          dDETBaseMemDim = 0;
          if (iDETBaseMem) free(iDETBaseMem);
          iDETBaseMem    = NULL;
          iDETBaseMemDim = 0;
          return ReportMemError("Determinant");
        }
      dDETBaseMemDim = memneeded;
    }

  memneeded = N + liwork;
  if (memneeded > iDETBaseMemDim)
    {
      iDETBaseMem = realloc(iDETBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iDETBaseMem == NULL)
        {
          if (dDETBaseMem) free(dDETBaseMem);
          dDETBaseMem    = NULL;
          dDETBaseMemDim = 0;
          if (iDETBaseMem) free(iDETBaseMem);
          iDETBaseMem    = NULL;
          iDETBaseMemDim = 0;
          return ReportMemError("Determinant");
        }
      iDETBaseMemDim = memneeded;
    }

  A    = dDETBaseMem;
  work = A + N*N;

  // Copy the matrix
  COPY(N*N, M, 1, A, 1);

  memset((void *)iDETBaseMem, 0, iDETBaseMemDim*sizeof(LAPACK_SIZE_T));
  ipiv  = iDETBaseMem;
  iwork = ipiv + N;

  dgetrf(&nc, &nc, A, &nc, ipiv, &info);
  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dgetrf", abs((int)info));
      return ILLEGAL_INPUT;
    }

  if (det)
    {
      *det = 1.0;
      if (!info)
        {
          for (j = 0; j < N; j++)
            {
              if (ipiv[j] != (j + 1))
                *det *= -A[j*N + j];
              else
                *det *= A[j*N + j];
            }
        }
    }

  if (info > 0) return SINGULARITY;

  if (cond)
    {
      norm      = anorm(N, N, M);
      whichnorm = '1';
      dgecon(&whichnorm, &nc, A, &nc, &norm, cond, work, iwork, &info FCONE);
      if (info < 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in DGECON", abs((int)info));
          return ILLEGAL_INPUT;
        }
    }

  return retval;
}


/*==================================================================================================================================*/

int Eigenval(const int N, double *X, const int symmetric, double *eigval, const int eigenvaltype, double *righteigvec, double *lefteigvec,
             const double tol)

{
/*
 * This function calculates the dominant, real eigenvalue of the N*N general
 * matrix X and the corresponding right and left eigenvector.
 *
 * The matrix has to be in normal C-wise row format and is transferred in this
 * routine to Format vector format.
 *
 * The eigenvalue is returned in *eigval, whereas the vector eigevec (length N)
 * contains the eigenvector. The content of X is not changed.
 *
 * This function first queries the Lapack routines for optimal workspace sizes.
 * These memoryblocks are then allocated and the decomposition is calculated
 * using the Lapack function "dgeevx". The allocated memory is preserved for
 * subsequent calls.
 */

  // balanc = 'B'           : Do scaling and permutation to make computation more accurate
  // jobvl  = jobvr  = 'V'  : Both eigenvectors are needed to compute error estimates of eigenvalues
  // sense  = 'E'           : Compute error estimates on eigenvalues only

  char          balanc = 'B', jobvl = 'V', jobvr = 'V', sense = 'E';
  char          jobz = 'V', range = 'I', uplo = 'U';
  int           i, j, dbasemem, ibasemem, memneeded;
  double        *Xc = NULL, *wr = NULL, *wi = NULL, *vr = NULL, *vl = NULL, *scale = NULL, *rconde = NULL, *rcondv = NULL, *work = NULL;
  double        abnrm, abstol, ddummy;
  int           retval = SUCCES, eigvalindx = 0;
  LAPACK_SIZE_T nc = N, ilo, ihi, lwork, info, *iwork, liwork, nfound, *isuppz;

  // Allocate temporarily minimally allowed size for workspace arrays
  if (symmetric)
    {
      dbasemem = 2*N*N + N;
      lwork    = 26*N;
    }
  else
    {
      dbasemem = 3*N*N + 5*N;
      lwork    = N*(N + 6);
    }

  memneeded = dbasemem + lwork;
  if (memneeded > dEVBaseMemDim)
    {
      dEVBaseMem = realloc(dEVBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      dEVBaseMemDim = memneeded;
    }

  if (symmetric)
    {
      ibasemem = 2*N;
      liwork   = 10*N;
    }
  else
    {
      ibasemem = 0;
      liwork   = (2*N - 2);
    }

  memneeded = ibasemem + liwork;
  if (memneeded > iEVBaseMemDim)
    {
      iEVBaseMem = realloc(iEVBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      iEVBaseMemDim = memneeded;
    }

  memset((void *)dEVBaseMem, 0, dEVBaseMemDim*sizeof(double));
  Xc = dEVBaseMem;
  vr = Xc + N*N;
  wr = vr + N*N;
  if (symmetric)
    work = wr + N;
  else
    {
      vl     = wr + N;
      wi     = vl + N*N;
      scale  = wi + N;
      rconde = scale + N;
      rcondv = rconde + N;
      work   = rcondv + N;
    }

  memset((void *)iEVBaseMem, 0, iEVBaseMemDim*sizeof(LAPACK_SIZE_T));
  iwork  = iEVBaseMem;
  isuppz = iwork + liwork;                                                          // Only used in case of symmetric matrix

  // Get the machine precisions
  abstol = dlamch("Safe minimum" FCONE);

  // Rewrite C-style row vector format to Fortran-style column vector format
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) Xc[j*N + i] = X[i*N + j];

  ilo = ihi = N;
  // Query the Lapack routine for optimal sizes for workspace arrays
  if (symmetric)
    {
      lwork = liwork = -1;
      dsyevr(&jobz, &range, &uplo, &nc, Xc, &nc, &ddummy, &ddummy, &ilo, &ihi, &abstol, &nfound, wr, vr, &nc, isuppz, work, &lwork, iwork, &liwork,
             &info FCONE FCONE FCONE);
      lwork  = (int)work[0];
      liwork = (int)iwork[0];
    }
  else
    {
      lwork = -1;
      dgeevx(&balanc, &jobvl, &jobvr, &sense, &nc, Xc, &nc, wr, wi, vl, &nc, vr, &nc, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork,
             &info FCONE FCONE FCONE FCONE);
      lwork = (int)work[0];
    }

  // Free previous allocation and reallocate preferable workspaces, Check result
  memneeded = dbasemem + lwork;
  if (memneeded > dEVBaseMemDim)
    {
      dEVBaseMem = realloc(dEVBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      dEVBaseMemDim = memneeded;
    }

  memneeded = ibasemem + liwork;
  if (memneeded > iEVBaseMemDim)
    {
      iEVBaseMem = realloc(iEVBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      iEVBaseMemDim = memneeded;
    }

  memset((void *)dEVBaseMem, 0, dEVBaseMemDim*sizeof(double));
  Xc = dEVBaseMem;
  vr = Xc + N*N;
  wr = vr + N*N;
  if (symmetric)
    work = wr + N;
  else
    {
      vl     = wr + N;
      wi     = vl + N*N;
      scale  = wi + N;
      rconde = scale + N;
      rcondv = rconde + N;
      work   = rcondv + N;
    }

  memset((void *)iEVBaseMem, 0, iEVBaseMemDim*sizeof(LAPACK_SIZE_T));
  iwork  = iEVBaseMem;
  isuppz = iwork + liwork;                                                          // Only used in case of symmetric matrix

  // Rewrite C-style row vector format to Fortran-style column vector format
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) Xc[j*N + i] = X[i*N + j];

  // Now calculate the eigenvalues and vectors using optimal workspaces
  if (symmetric)
    dsyevr(&jobz, &range, &uplo, &nc, Xc, &nc, &ddummy, &ddummy, &ilo, &ihi, &abstol, &nfound, wr, vr, &nc, isuppz, work, &lwork, iwork, &liwork,
           &info FCONE FCONE FCONE);
  else
    dgeevx(&balanc, &jobvl, &jobvr, &sense, &nc, Xc, &nc, wr, wi, vl, &nc, vr, &nc, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork,
           &info FCONE FCONE FCONE FCONE);

  // Check for convergence
  if (info < 0)
    {
      if (symmetric)
        ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dsyevr()", abs((int)info));
      else
        ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dgeevx()", abs((int)info));
      retval = ILLEGAL_INPUT;
    }
  else if (info > 0)
    {
      ErrorMsg(__FILE__, __LINE__, "The algorithm failed to compute eigenvalues!");
      retval = NO_CONVERGENCE;
    }
  else
    {
      switch (eigenvaltype)
        {
          case DOMINANT:
            // Find the real eigenvalue with the largest real part
            if (symmetric)
              eigvalindx = 0;
            else
              {
                eigvalindx = -1;
                for (i = 0; i < N; i++)
                  if ((eigvalindx == -1) || (wr[i] > wr[eigvalindx])) eigvalindx = i;
                if ((epsMach*abnrm/rconde[eigvalindx]) > tol)
                  {
                    ErrorMsg(__FILE__, __LINE__, "The estimated error bound in the largest eigenvalue (%G) exceeds the tolerance level %G!",
                             (epsMach*abnrm/rconde[eigvalindx]), tol);
                    retval = NO_CONVERGENCE;
                  }
              }
            break;
          case MINIMUMNORM:
            // Find the real eigenvalue with the smallest norm
            eigvalindx = -1;
            for (i = 0; i < N; i++)
              if (wi[i] == (double)0.0)
                if ((eigvalindx == -1) || (fabs(wr[i]) < fabs(wr[eigvalindx]))) eigvalindx = i;

            if (eigvalindx == -1)
              {
                ErrorMsg(__FILE__, __LINE__, "Did not find any real eigenvalue with smallest norm in Eigenval!");
                retval = FAILURE;
              }
            else if (!symmetric && ((epsMach*abnrm/rconde[eigvalindx]) > tol))
              {
                ErrorMsg(__FILE__, __LINE__, "The estimated error bound in the minimum eigenvalue (%G) exceeds the tolerance level %G!",
                         (epsMach*abnrm/rconde[eigvalindx]), tol);
                retval = NO_CONVERGENCE;
              }

            break;
        }

      if (retval == SUCCES)
        {
          *eigval = wr[eigvalindx];
          COPY(nc, vr + eigvalindx*nc, 1, righteigvec, 1);
          // Left and right eigenvectors are the same
          if ((!symmetric) && lefteigvec) COPY(nc, vl + eigvalindx*nc, 1, lefteigvec, 1);
        }
    }

  return retval;
}


/*==================================================================================================================================*/

int SolveLinearSystem(const int N, double *A, double *B, double tol)

{
/*
 * This function solves the linear equation system A*x = B, where A is a NxN
 * matrix and B a N-dimensional vector.
 *
 * The matrix A has to be in Fortran column-wise format
 *
 * The solution is returned in *B, whereas the content of A is not changed.
 */

  char          fact, trans, equed, errstr[MAX_STR_LEN];
  int           i, j, memneeded;
  double        *Ac, *Af, *r, *c, *Bc, *x, *work, rcond, ferr = 0, berr;
  int           retval = SUCCES;
  LAPACK_SIZE_T nc = N, nrhs = 1, lwork = 4*N, *ipiv, *iwork, info;

  // Allocate temporarily minimally allowed size for workspace arrays
  memneeded = 2*N*N + 4*N + lwork;
  if (memneeded > dSLSBaseMemDim)
    {
      dSLSBaseMem = realloc(dSLSBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dSLSBaseMem == NULL)
        {
          if (dSLSBaseMem) free(dSLSBaseMem);
          dSLSBaseMem    = NULL;
          dSLSBaseMemDim = 0;
          if (iSLSBaseMem) free(iSLSBaseMem);
          iSLSBaseMem    = NULL;
          iSLSBaseMemDim = 0;
          return ReportMemError("SolveLinearSystem");
        }
      dSLSBaseMemDim = memneeded;
    }

  memneeded = 2*N;
  if (memneeded > iSLSBaseMemDim)
    {
      iSLSBaseMem = realloc(iSLSBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iSLSBaseMem == NULL)
        {
          if (dSLSBaseMem) free(dSLSBaseMem);
          dSLSBaseMem    = NULL;
          dSLSBaseMemDim = 0;
          if (iSLSBaseMem) free(iSLSBaseMem);
          iSLSBaseMem    = NULL;
          iSLSBaseMemDim = 0;
          return ReportMemError("SolveLinearSystem");
        }
      iSLSBaseMemDim = memneeded;
    }

  memset((void *)dSLSBaseMem, 0, dSLSBaseMemDim*sizeof(double));
  Ac   = dSLSBaseMem;
  Af   = Ac + N*N;
  r    = Af + N*N;
  c    = r + N;
  Bc   = c + N;
  x    = Bc + N;
  work = x + N;

  memset((void *)iSLSBaseMem, 0, iSLSBaseMemDim*sizeof(LAPACK_SIZE_T));
  ipiv  = iSLSBaseMem;
  iwork = ipiv + N;

  fact  = 'E';
  trans = 'N';

  // Fill the matrix and the right-hand side vector
  COPY(N*N, A, 1, Ac, 1);
  COPY(N, B, 1, Bc, 1);

  dgesvx(&fact, &trans, &nc, &nrhs, Ac, &nc, Af, &nc, ipiv, &equed, r, c, Bc, &nc, x, &nc, &rcond, &ferr, &berr, work, iwork, &info FCONE FCONE FCONE);

  // Check for singularity of the matrix
  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dgesvx()", abs((int)info));
      retval = ILLEGAL_INPUT;
    }
  else if (info > 0)
    {
      ErrorMsg(__FILE__, __LINE__, "(Nearly) Singular matrix in SolveLinearSystem(), solving the linear system A*x = B:\n");
      for (i = 0; i < N; i++)
        {
          if ((2*i == (N - 1)) || (2*i == N))
            STDOUT(" A = |");
          else
            STDOUT("     |");
          for (j = 0; j < N; j++) STDOUT("%16.8E", A[j*N + i]);
          if ((2*i == (N - 1)) || (2*i == N))
            STDOUT("|,     B = | %16.8E|\n", B[i]);
          else
            STDOUT("|          | %16.8E|\n", B[i]);
        }
      retval = SINGULARITY;
    }
  else
    {
      if (ferr > tol)
        {
          // Matlab does not handle correctly the direct printing of the double number via ErrorMsg
          snprintf(errstr, sizeof(errstr), "Warning: The estimated error bound in the solution of the linear system A*x = B (%G) exceeds the tolerance level %G!\n",
                  ferr, tol);
          ErrorMsg(__FILE__, __LINE__, errstr);
        }
      COPY(N, x, 1, B, 1);
    }

  return retval;
}


/*==================================================================================================================================*/

void ResetCurve(void)
{
  if (pnt_scale) free(pnt_scale);
  pnt_scale           = NULL;
  if (oldscale) free(oldscale);
  oldscale            = NULL;

  if (dDETBaseMem) free(dDETBaseMem);
  dDETBaseMem         = NULL;
  dDETBaseMemDim      = 0;
  if (iDETBaseMem) free(iDETBaseMem);
  iDETBaseMem         = NULL;
  iDETBaseMemDim      = 0;

  if (dEVBaseMem) free(dEVBaseMem);
  dEVBaseMem          = NULL;
  dEVBaseMemDim       = 0;
  if (iEVBaseMem) free(iEVBaseMem);
  iEVBaseMem          = NULL;
  iEVBaseMemDim       = 0;

  if (dSLSBaseMem) free(dSLSBaseMem);
  dSLSBaseMem         = NULL;
  dSLSBaseMemDim      = 0;
  if (iSLSBaseMem) free(iSLSBaseMem);
  iSLSBaseMem         = NULL;
  iSLSBaseMemDim      = 0;

  fast_iters          = 0;
  slow_iters          = 0;
  FirstTangent        = 1;

  LPImmediateReturn   = 0;
  ESSImmediateReturn  = 0;

  return;
}


/*==================================================================================================================================*/

int SelectionGradient(const int pntdim, double *pnt, int (*fnc)(double *, double *), const int parindex, const int R0index, double *dR0dp)

/*
 * SelectionGradient -  Routine computes the factor determining the location of an ESS, i.e.
 *                      the derivative of R0 w.r.t. a parameter.
 *
 * Arguments -  pntdim    : The dimension of the argument vector 'y'. Notice that this
 *                          equals 2 (for the parameters) plus the dimension of the vector of
 *                          state variables.
 *              y         : Pointer to an array containing as first element the value
 *                          of the parameter p and as subsequent elements the values of
 *                          the state variables y. The last element is assumed to be the
 *                          second parameter in the ESS continuation
 *              fnc       : Pointer to function specifying the system of
 *                          equations. The function must have a (double)
 *                          pointer as first argument, containing the point
 *                          in which to evaluate the system and a (double)
 *                          pointer as second argument, containing the
 *                          results after evaluation of the equations.
 *              parindex  : Index of the parameter in the solution point 'pnt' w.r.t. which the derivative has to be calculated
 *                          If parindex >= pntdim, the difference (parindex-pntdim) is taken as an index into the parameter
 *                          vector and the selection gradient w.r.t. this parameter is calculated.
 *              R0index   : Index of the equilibrium condition (R0-1 = 0) in the result vector
 *                          returned by (*fnc)().
 */
{
  double  pardif;
  double  *basemem, *y, *rhs, *parpnt;

  // Prevent recurrence
  *dR0dp = 0.0;
  if (ESSImmediateReturn) return SUCCES;
  ESSImmediateReturn = 1;

  y = basemem = calloc(2*pntdim, sizeof(double));
  if (!basemem)
    {
      ESSImmediateReturn = 0;
      return ReportMemError("SelectionGradient");
    }

  rhs = y + pntdim;

  // Initialize
  COPY(pntdim, pnt, 1, y, 1);

  if (parindex < pntdim)
    parpnt = y + parindex;
  else
    parpnt = parPntr + (parindex - pntdim);

  pardif = fabs(Jacobian_Step*(*parpnt));
  pardif = (parindex < pntdim) ? max(pardif, Jacobian_Min_Step/pnt_scale[parindex]) : max(pardif, Jacobian_Min_Step);

  if (CentralDerivative(1, fnc, y, rhs, parpnt, pardif, rhs + R0index, dR0dp, FastNumerics) == FAILURE)
    {
      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
      free(basemem);
      ESSImmediateReturn = 0;
      return FAILURE;
    }
  // If part of the solution point compute the derivative w.r.t. **unscaled** parameter
  if (parindex < pntdim) *dR0dp /= pnt_scale[parindex];

  // Additional call to reset global variables
  (*fnc)(y, rhs);

  free(basemem);
  ESSImmediateReturn = 0;

  // Return the derivative w.r.t. the scaled or unscaled parameter?
  // par=scaledpar*scale => dR0/dscaledpar = dR0/d(par/scale) = scale*DR0/dpar

  return SUCCES;
}


/*==================================================================================================================================*/
