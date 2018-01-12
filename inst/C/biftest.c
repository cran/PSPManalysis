/***
  NAME
    biftest
  DESCRIPTION
    This module implements routines that locate bifurcation points along
    the curve that is being continued

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

    Last modification: AMdR - Jan 12, 2018
***/
#ifndef BIFTEST
#define BIFTEST
#endif
#include "globals.h"


/*
 *====================================================================================================================================
 *      Routines to locate bifurcation points
 *====================================================================================================================================
 */

int LocateLP(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol)

/*
 * LocateLP - Routine locates a limit point in an equilibrium bifurcation curve.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y = (p,x)'.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables x.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              ytol    : Tolerance determining when change in y equals zero.
 *              rhstol  : Tolerance determining when RHS equals zero.
 */
{
  int     i, retval, outmax;
  double  *basemem, *lppoint;

  lppoint = basemem = calloc(pntdim, sizeof(double));
  if (!basemem) return ReportMemError("LocateLP");

  COPY(pntdim, y, 1, lppoint, 1);

  ReportMsg("\nStarting LP location from :\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", lppoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  LocalizeType = LP;
  retval       = FindPoint(pntdim, lppoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  LocalizeType = UNDEFINED;

  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New LP point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", lppoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      for (i = 0; i < pntdim; i++)
        {
#if (defined(R_PACKAGE))
          if (i) 
            STDOUT(",%15.8E", lppoint[i]*pnt_scale[i]);
          else
#endif
            STDOUT("%16.8E",  lppoint[i]*pnt_scale[i]);          
        }
      STDOUT("  ****  LP      ****");
      STDOUT("\n");
      outmax = DefineOutput(lppoint, Output);
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  LP      ****");
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
    }
  else
    {
      ReportMsg("\nFailed to locate LP bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate LP bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate LP bifurcation point ****\n");
#endif
    }
  fflush(NULL);
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  free(basemem);

  return retval;
}


/*==================================================================================================================================*/

int LocateBP(int *dimpntr, double *y, int (*fnc)(double *, double *), double dytol, double rhstol, int nr, int indx)

/*
 * LocateBP - Routine locates a transcritical bifurcation where the birth rate of the
 *            structured population with index nr crosses 0
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              ytol    : Tolerance determining when change in y equals zero.
 *              rhstol  : Tolerance determining when RHS equals zero.
 *              par2    : Pointer to the value of the 2nd bifurcation parameter
 *              nr      : Number of the population in which the bifurcation occurs
 *              indx    : Index in the solution vector y of birth rate of the
 *                        structured population in which the bifurcation occurs
 *                        (-1 in case of trivial equilibrium continuation)
 */
{
  int     i, j, retval, outmax;
  int     bppntdim;
  int     oldpntdim, oldPopBPIndex;
  double  *basemem, *bppoint;
  double  *tmpscales, *savedpntr;

  // Save the old pntdim value
  oldpntdim =*dimpntr;

  if (indx < 0)                                                                      // Detection from trivial equilibrium
    bppntdim = oldpntdim;
  else
    bppntdim = oldpntdim - 1;

  // Keep the vectors the same size as the original point to be able to revert them back
  // to the proper dimension for output generation
  bppoint = basemem = calloc(2*oldpntdim, sizeof(double));
  if (!basemem) return ReportMemError("LocateBP");

  tmpscales = bppoint + oldpntdim;

  // Copy the vectors, cutting out the BP variable
  for (i = 0, j = 0; i < oldpntdim; i++)
    {
      if (i == indx) continue;
      tmpscales[j] = pnt_scale[i];
      bppoint[j]   = y[i];
      j++;
    }

  savedpntr = pnt_scale;
  pnt_scale = tmpscales;

  *dimpntr = bppntdim;

  oldPopBPIndex = PopBPIndex;
  PopBPIndex    = nr;

  ReportMsg("\nStarting BP location from :\t");
  for (i = 0; i < bppntdim; i++) ReportMsg("%10.5E\t", bppoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  retval = FindPoint(bppntdim, bppoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New BP point :\t");
      for (i = 0; i < bppntdim; i++) ReportMsg("%16.8E  ", bppoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      // Now cast back to the original vector length
      for (i = oldpntdim - 1, j = bppntdim - 1; i >= 0; i--)
        {
          if (i == indx)
            bppoint[i] = 0.0;
          else
            {
              bppoint[i] = bppoint[j];
              j--;
            }
        }

      // And reset the values
      pnt_scale  = savedpntr;
      *dimpntr   = oldpntdim;
      PopBPIndex = oldPopBPIndex;

      for (i = 0; i < oldpntdim; i++) 
        {
#if (defined(R_PACKAGE))
          if (i) 
            STDOUT(",%15.8E", bppoint[i]*pnt_scale[i]);
          else
#endif
            STDOUT("%16.8E",  bppoint[i]*pnt_scale[i]);          
        }
      STDOUT("  ****  BP #%d   ****", nr);
      STDOUT("\n");

      outmax = DefineOutput(bppoint, Output);
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  BP #%d   ****", nr);
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }
  else
    {
      ReportMsg("\nFailed to locate BP bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate BP bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate BP bifurcation point ****\n");
#endif
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  pnt_scale  = savedpntr;
  *dimpntr   = oldpntdim;
  PopBPIndex = oldPopBPIndex;

  free(basemem);

  return retval;
}


/*==================================================================================================================================*/

int LocateBPE(int *dimpntr, double *y, int (*fnc)(double *, double *), double dytol, double rhstol, int nr, int indx)

/*
 * LocateBPE -  Routine locates a transcritical bifurcation where the environment
 *              variable with index nr crosses 0
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              ytol    : Tolerance determining when change in y equals zero.
 *              rhstol  : Tolerance determining when RHS equals zero.
 *              par2    : Pointer to the value of the 2nd bifurcation parameter
 *              nr      : Number of the environment variable in which the bifurcation occurs
 *              indx    : Index in the solution vector y of the environment variable
 *                        in which the bifurcation occurs (-1 in case of trivial equilibrium continuation)
 */
{
  int     i, j, retval, outmax;
  int     bppntdim, oldEnvBPIndex;
  int     oldpntdim;
  double  *basemem, *bppoint;
  double  *tmpscales, *savedpntr;

  // Save the old pntdim value
  oldpntdim =*dimpntr;

  if (indx < 0)                                                                      // Detection from trivial equilibrium
    bppntdim =*dimpntr;
  else
    bppntdim =*dimpntr - 1;

  // Keep the vectors the same size as the original point to be able to revert them back
  // to the proper dimension for output generation
  bppoint = basemem = calloc(2*oldpntdim, sizeof(double));
  if (!basemem) return ReportMemError("LocateBPE");

  tmpscales = bppoint + oldpntdim;

  // Copy the vectors, cutting out the BP variable
  for (i = 0, j = 0; i < oldpntdim; i++)
    {
      if (i == indx) continue;
      tmpscales[j] = pnt_scale[i];
      bppoint[j]   = y[i];
      j++;
    }

  savedpntr = pnt_scale;
  pnt_scale = tmpscales;

  *dimpntr = bppntdim;

  oldEnvBPIndex = EnvBPIndex;
  EnvBPIndex    = nr;

  ReportMsg("\nStarting BPE location from :\t");
  for (i = 0; i < bppntdim; i++) ReportMsg("%10.5E\t", bppoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  retval = FindPoint(bppntdim, bppoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New BPE point :\t");
      for (i = 0; i < bppntdim; i++) ReportMsg("%16.8E  ", bppoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      // Now cast back to the original vector length
      for (i = oldpntdim - 1, j = bppntdim - 1; i >= 0; i--)
        {
          if (i == indx)
            bppoint[i] = 0.0;
          else
            {
              bppoint[i] = bppoint[j];
              j--;
            }
        }

      // And reset the values
      pnt_scale  = savedpntr;
      *dimpntr   = oldpntdim;
      EnvBPIndex = oldEnvBPIndex;

      for (i = 0; i < oldpntdim; i++) 
        {
#if (defined(R_PACKAGE))
          if (i) 
            STDOUT(",%15.8E", bppoint[i]*pnt_scale[i]);
          else
#endif
            STDOUT("%16.8E",  bppoint[i]*pnt_scale[i]);          
        }
      STDOUT("  ****  BPE #%d  ****", nr);
      STDOUT("\n");

      outmax = DefineOutput(bppoint, Output);
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  BPE #%d  ****", nr);
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }
  else
    {
      ReportMsg("\nFailed to locate BPE bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate BPE bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate BPE bifurcation point ****\n");
#endif
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  pnt_scale  = savedpntr;
  *dimpntr   = oldpntdim;
  EnvBPIndex = oldEnvBPIndex;

  free(basemem);

  return retval;
}


/*==================================================================================================================================*/

int LocateESS(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol, const int R0index)

/*
 * LocateESS -  Routine locates an ESS point in an equilibrium bifurcation curve.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              ytol    : Tolerance determining when change in y equals zero.
 *              rhstol  : Tolerance determining when RHS equals zero.
 *              R0index : Index of the equilibrium condition (R0-1 = 0) in the result vector
 *                        returned by (*fnc)().
 */
{
  int     i, retval, outmax;
  double  *basemem, *esspoint;
  double  evJ, evH, zCz;

  esspoint = basemem = calloc(pntdim, sizeof(double));
  if (!basemem) return ReportMemError("LocateESS");

  // Initialize
  COPY(pntdim, y, 1, esspoint, 1);

  ReportMsg("\nStarting ESS location from :\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", esspoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  LocalizeType = ESS;
  retval       = FindPoint(pntdim, esspoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  LocalizeType = UNDEFINED;

  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New ESS point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", esspoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      for (i = 0; i < pntdim; i++)
        {
#if (defined(R_PACKAGE))
          if (i) 
            STDOUT(",%15.8E", esspoint[i]*pnt_scale[i]);
          else
#endif
            STDOUT("%16.8E",  esspoint[i]*pnt_scale[i]);          
        }
      // ESS located w.r.t. the bifurcation parameter
      retval = ESSclassify(pntdim, esspoint, fnc, dytol, rhstol, R0index, &evJ, &evH, &zCz, 1);
      if (retval == SUCCES)
        {
          if (evJ > 0)
            STDOUT("  ****   ERP #%d  ****", PopEVOIndex);
          else if ((evJ < 0) && (evH < 0))
            STDOUT("  ****   CSS #%d  ****", PopEVOIndex);
          else if ((evJ < 0) && (!evoParsDim || (zCz < 0)))
            STDOUT("  ****   EBP #%d  ****", PopEVOIndex);
          else if (evJ < 0)
            STDOUT("  ****  ENBP #%d  ****", PopEVOIndex);
          else
            STDOUT("  ****   ESS #%d  ****", PopEVOIndex);
        }
      else
        STDOUT("  ****  ESS #%d  ****", PopEVOIndex);
      STDOUT("\n");
      fflush(NULL);

      outmax = DefineOutput(esspoint, Output);
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          if (retval == SUCCES)
            {
              if (evJ > 0)
                fprintf(biffile, "  ****   ERP #%d  ****", PopEVOIndex);
              else if ((evJ < 0) && (evH < 0))
                fprintf(biffile, "  ****   CSS #%d  ****", PopEVOIndex);
              else if ((evJ < 0) && (!evoParsDim || (zCz < 0)))
                fprintf(biffile, "  ****   EBP #%d  ****", PopEVOIndex);
              else if (evJ < 0)
                fprintf(biffile, "  ****  ENBP #%d  ****", PopEVOIndex);
              else
                fprintf(biffile, "  ****   ESS #%d  ****", PopEVOIndex);
            }
          else
            fprintf(biffile, "  ****  ESS #%d  ****", PopEVOIndex);
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }
  else
    {
      ReportMsg("\nFailed to locate ESS bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate ESS bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate ESS bifurcation point ****\n");
#endif
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  free(basemem);

  return retval;
}


/*==================================================================================================================================*/
// Computing the Hessian instead of C01 seems numerically less stable
#ifndef COMPUTEHESSIAN
#define COMPUTEHESSIAN 1
#else
#if (COMPUTEHESSIAN == 0)
#warning COMPUTEHESSIAN overruled and set to 0
#elif (COMPUTEHESSIAN == 1)
#warning COMPUTEHESSIAN overruled and set to 1
#else
#error COMPUTEHESSIAN should be either set to 0 or to 1
#endif
#endif
#ifndef COMPUTEJACOBIAN
#define COMPUTEJACOBIAN       0
#else
#if (COMPUTEJACOBIAN == 0)
#warning COMPUTEJACOBIAN overruled and set to 0
#elif (COMPUTEJACOBIAN == 1)
#warning COMPUTEJACOBIAN overruled and set to 1
#else
#error COMPUTEJACOBIAN should be either set to 0 or to 1
#endif
#endif

#if (COMPUTEJACOBIAN == 1) && (COMPUTEHESSIAN == 1)
#define COMPUTECROSS          0
#else
#define COMPUTECROSS          1
#if (COMPUTEJACOBIAN == 0) && (COMPUTEHESSIAN == 0)
#warning COMPUTEJACOBIAN set to 1 as both COMPUTEJACOBIAN and COMPUTEHESSIAN are defined equal to 0
#undef COMPUTEJACOBIAN
#define COMPUTEJACOBIAN       1
#endif
#endif


int ESSclassify(const int pntdim, double *pnt, int (*fnc)(double *, double *), double dytol, double rhstol, const int R0index, double *EVmaxJ,
                double *EVmaxH, double *zC01z, const int detecting)

/*
 * ESSclassify -  Routine assess the stability properties of an ESS point by computing
 *                the second derivatives of R0(E(x),y) w.r.t. life history parameters, both
 *                the resident and mutant parameter value.
 *
 *                General equations for central difference approximation to the first derivative:
 *
 *                f'(x) = (f(x+h) - f(x-h))/(2h)
 *
 *                and to the second derivative:
 *
 *                f''(x) = (f(x+h) - 2f(x) + f(x-h))/h^2
 *
 *                Notice that this equation is derived from the first order difference approximation
 *                with step size h/2
 *
 *                f''(x) = (f'(x+h/2) - f'(x-h/2))/h
 *
 *                by substituting both f'(x+h/2) and f'(x-h/2) with their respective difference approximations
 *                with step size h/2:
 *
 *                f''(x) = ((f(x+h)-f(x))/h - (f(x)-f(x-h))/h)/h = (f(x+h) - 2f(x) + f(x-h))/h^2
 *
 *                In case of a 1-dimensional ESS point (a single ESS parameter):
 *
 *                R0_xx = c11 = (R0(E(x*+h), x*) - 2 R0(E(x*), x*) + R0(E(x*-h), x*))/h^2
 *                R0_yy = c00 = (R0(E(x*), x*+h) - 2 R0(E(x*), x*) + R0(E(x*), x*-h))/h^2
 *
 *                In the singular point moreover holds:
 *
 *                R0(E(x*),x*)=1
 *                R0(E(x*), x*-h) = R0(E(x*), x*+h)        (because dR0/dx = 0)
 *                R0(E(x*-h), x*) = R0(E(x*+h), x*)        (and hence also dR0/dE dE/dx = 0)
 *
 *                Hence, we only have to compute
 *
 *                R0_xx = c11 = 2*(R0(E(x*+h), x*)-1) / h^2
 *                R0_yy = c00 = 2*(R0(E(x*), x*+h)-1) / h^2
 *
 *                Moreover, since only the signs and the ratio of R0_xx and R0_yy are important
 *                the factor 2/h^2 could in principle be ignored
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the parameters) plus the dimension of the vector of
 *                        state variables.
 *              pnt     : Pointer to an array containing as first element the value
 *                        of the bifurcation parameter p and as subsequent elements
 *                        the values of the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              ytol    : Tolerance determining when change in y equals zero.
 *              rhstol  : Tolerance determining when RHS equals zero.
 *              R0index : Index of the equilibrium condition (R0-1 = 0) in the result vector
 *                        returned by (*fnc)().
 */
{
  int     i, j, ii, jj, indx1, indx2, evoparsdim;
  int     retval, OldCurveType;
  double  oldvalue1, oldvalue2;
  double  *basemem, *y, *tv, *rhs0, *rhs, *Jac, *pardif, *eigvecH, *eigvecJ, *CEhes, *CEjac, *CEC01;

/* Allocate the necessary memory:
 * local copy of solution point, tangent vector, right-hand side and new environment
 * Jacobian, Hessian and C_01 matrices of the canonical equation
 */
  evoparsdim = evoParsDim + detecting;                                              // Include Bifparone in the evolving parameter list if detecting
  y = basemem = calloc(4*pntdim + pntdim*pntdim + 3*evoparsdim + 3*evoparsdim*evoparsdim, sizeof(double));
  if (!basemem) return ReportMemError("ESSclassify");

  ReportMsg("\nClassifying ESS point:\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", pnt[i]*pnt_scale[i]);
  ReportMsg("\n");

  tv      = y + pntdim;
  rhs0    = tv + pntdim;
  rhs     = rhs0 + pntdim;
  Jac     = rhs + pntdim;
  pardif  = Jac + pntdim*pntdim;
  eigvecH = pardif + evoparsdim;
  eigvecJ = eigvecH + evoparsdim;
  CEhes   = eigvecJ + evoparsdim;
  CEjac   = CEhes + evoparsdim*evoparsdim;
  CEC01   = CEjac + evoparsdim*evoparsdim;

  // Initialize
  COPY(pntdim, pnt, 1, y, 1);
  memset((void *)CEhes, 0, evoparsdim*evoparsdim*sizeof(double));
  memset((void *)CEjac, 0, evoparsdim*evoparsdim*sizeof(double));
  memset((void *)CEC01, 0, evoparsdim*evoparsdim*sizeof(double));
  memset((void *)Jac, 0, pntdim*pntdim*sizeof(double));
  memset((void *)tv, 0, pntdim*sizeof(double));

  // Compute the Jacobian for general use in FindPoint() in the loop below
  Jacobian(pntdim, y, pntdim - 1, Jac, fnc, FORWARD);

  // FindPoint() below only needs Jacobian of system without evolutionary parameters
  // Move appropriate elements of Jacobian into right place, overwriting derivatives w.r.t. evolutionary parameters
  // memmove() used instead of memcpy() because src and dest overlap (undefined behavior in memcpy)
  for (i = 1; i < (pntdim - evoParsDim); i++)
    memmove(Jac + i*(pntdim - evoParsDim - 1), Jac + i*(pntdim - 1), (pntdim - evoParsDim - 1)*sizeof(double));

  // Function call to set all parameters to their appropriate values, preserve these values in rhs0
  if ((*fnc)(y, rhs0) == FAILURE)
    {
      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
      free(basemem);
      return FAILURE;
    }

  // Preparatory computations
  OldCurveType = CurveType;                                                         // Save the curve type
  CurveType    = EQ;                                                                // This prevents setting evolutionary parameters
  setBifParVal = 0;                                                                 // Freeze the value of parameter[Bifparone]

  // Compute all the step sizes
  for (i = 0; i < evoparsdim; i++)
    {
      if (detecting && !i)
        indx1 = Bifparone;                                                          // Classification after detecting an ESS in parameter[Bifparone]
      else
        indx1 = evoParsIndexPntr[i - detecting];                                    // Runs to i=0..evoParsDim in detecting and non-detecting case

      oldvalue1      = parPntr[indx1];
      pardif[i]      = max(fabs(Jacobian_Step*parPntr[indx1]), Jacobian_Min_Step);
      parPntr[indx1] = oldvalue1 + pardif[i];

      // Trick to reduce precision errors. See Num. Recipes 9.7, pg. 388
      pardif[i]      = parPntr[indx1] - oldvalue1;
      parPntr[indx1] = oldvalue1;
    }

/* Keep E(x*) as computed for the current point and compute the Hessian elements for the evolving parameters:
 *
 * ((R0(x*+dpi+dpj, E(x*))-1) - (R0(x*-dpi+dpj, E(x*))-1) - (R0(x*+dpi-dpj, E(x*))-1) + (R0(x*-dpi-dpj, E(x*))-1))/(4*eps_pi*eps_pj)        (*)
 *
 * where dpi = (0,.....eps_pi,....)^T and dpj = (0,.....eps_pj,....)^T are small perturbation vectors in the direction
 * of the parameter axis i and j, respectively.
 *
 * In case evoparsdim = 1 the loop below computes
 *
 *        R0_yy = c00 = (R0(E(x*), x*+h) - 2 R0(E(x*), x*) + R0(E(x*), x*-h))/h^2
 *
 * Given that in a singular point:
 *
 *        R0(E(x*),x*)=1
 *        R0(E(x*), x*-h) = R0(E(x*), x*+h)        (because dR0/dx = 0)
 */
#if (COMPUTEHESSIAN == 1)
  for (i = 0; i < evoparsdim; i++)
    {
      if (detecting && !i)
        indx1 = Bifparone;                                                          // Classification after detecting an ESS in parameter[Bifparone]
      else
        indx1 = evoParsIndexPntr[i - detecting];                                    // Runs to i=0..evoParsDim in detecting and non-detecting case

      oldvalue1 = parPntr[indx1];
      /*
       * For the diagonal elements we calculate
       * (R0(x*+dpi, E(x*)) - 2*R0(x*, E(x*)) + R0(x*-dpi, E(x*)))/(eps_pi*eps_pi)
       *
       * This should equal
       * ((R0(x*+dpi, E(x*))-1) + (R0(x*-dpi, E(x*))-1))/(eps_pi*eps_pi)
       */
      CEhes[i*evoparsdim + i] = -2*(rhs0[R0index] + 1.0);
      for (ii = 0; ii < 2; ii++)
        {
          if (ii == 0)
            parPntr[indx1] = oldvalue1 + pardif[i];
          else
            parPntr[indx1] = oldvalue1 - pardif[i];

          if ((*fnc)(y, rhs) == FAILURE)
            {
              ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
              setBifParVal = 1;
              CurveType    = OldCurveType;                                          // Reset the curve type
              free(basemem);
              return FAILURE;
            }
          CEhes[i*evoparsdim + i] += (rhs[R0index] + 1.0);

          // Compute the below diagonal elements
          for (j = 0; j < i; j++)
            {
              if (detecting && !j)
                indx2 = Bifparone;                                                  // Classification after detecting an ESS in parameter[Bifparone]
              else
                indx2 = evoParsIndexPntr[j - detecting];                            // Runs to i=0..evoParsDim in detecting and non-detecting case

              oldvalue2 = parPntr[indx2];
              for (jj = 0; jj < 2; jj++)
                {
                  if (jj == 0)
                    parPntr[indx2] = oldvalue2 + pardif[j];
                  else
                    parPntr[indx2] = oldvalue2 - pardif[j];

                  if ((*fnc)(y, rhs) == FAILURE)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                      setBifParVal = 1;
                      CurveType    = OldCurveType;                                  // Reset the curve type
                      free(basemem);
                      return FAILURE;
                    }
                  // The factor (1-2*ii)*(1-2*jj) gives the appropriate sign in the expression (*) above
                  CEhes[i*evoparsdim + j] += (1 - 2*ii)*(1 - 2*jj)*(rhs[R0index] + 1.0);
                }
              parPntr[indx2] = oldvalue2;                                           // Reset parameter
            }
        }
      parPntr[indx1] = oldvalue1;                                                   // Reset parameter

      // Scale all elements
      CEhes[i*evoparsdim + i] /= pardif[i]*pardif[i];                               // This like in the 1-dimensional case
      for (j = 0; j < i; j++)
        {
          CEhes[i*evoparsdim + j] /= 4*pardif[i]*pardif[j];
          CEhes[j*evoparsdim + i] = CEhes[i*evoparsdim + j];                        // Hessian is symmetric
        }
    }
#endif

  // Compute the Jacobian matrix and/or the matrix of cross-derivatives
  tv[0] = 1.0;                                                                      // Keep the first (curve) parameter constant, adjust only E and b
  for (i = 0; i < evoparsdim; i++)
    {
      if (detecting && !i)
        indx1 = Bifparone;                                                          // Classification after detecting an ESS in parameter[Bifparone]
      else
        indx1 = evoParsIndexPntr[i - detecting];                                    // Runs to i=0..evoParsDim in detecting and non-detecting case

      oldvalue1 = parPntr[indx1];                                       
      for (ii = 0; ii < 2; ii++)
        {
          if (ii == 0)
            parPntr[indx1] = oldvalue1 + pardif[i];
          else
            parPntr[indx1] = oldvalue1 - pardif[i];

          // Initialize anew
          COPY(pntdim, pnt, 1, y, 1);
          retval = FindPoint(pntdim - evoParsDim, y, Jac, tv, dytol, rhstol, MAXITER, fnc);
          if (retval != SUCCES)
            {
              ErrorMsg(__FILE__, __LINE__, "Location of fixed point for 2nd derivatives failed");
              setBifParVal = 1;
              CurveType    = OldCurveType;                                             // Reset the curve type
              free(basemem);
              return FAILURE;
            }
          /*
           *  y now contains the correct environment E(x*+dpi) or E(x*-dpi), which hence satisfies the condition
           *
           *  R0(x*+dpi, E(x*+dpi)) = 1         or      R0(x*-dpi, E(x*-dpi)) = 1
           *
           *  with all parameters fixed at their value of the ESS, except for the single perturbed parameter
           */

#if (COMPUTEJACOBIAN == 1)
          /*
           * With y containing the correct environment E(x*+dpi) or E(x*-dpi), compute the Jacobian matrix of the canonical equation
           * using the following expression for its elements:
           *
           * (R0(x*+dpi+dpj,E(x*+dpi)) - R0(x*+dpi-dpj,E(x*+dpi)) - R0(x*-dpi+dpj,E(x*-dpi)) + R0(x*-dpi-dpj,E(x*-dpi)))/(4*eps_pi*eps_pj)              (**)
           *
           * where dpi = (0,.....eps_pi,....)^T and dpj = (0,.....eps_pj,....)^T are small perturbation vectors in the direction
           * of the parameter axis i and j, respectively.
           */
          for (j = 0; j < evoparsdim; j++)
            {
              if (detecting && !j)
                indx2 = Bifparone;                                                  // Classification after detecting an ESS in parameter[Bifparone]
              else
                indx2 = evoParsIndexPntr[j - detecting];                            // Runs to i=0..evoParsDim in detecting and non-detecting case

              oldvalue2 = parPntr[indx2];
              for (jj = 0; jj < 2; jj++)
                {
                  if (jj == 0)
                    parPntr[indx2] = oldvalue2 + pardif[j];
                  else
                    parPntr[indx2] = oldvalue2 - pardif[j];

                  if ((*fnc)(y, rhs) == FAILURE)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                      setBifParVal = 1;
                      CurveType    = OldCurveType;                                  // Reset the curve type
                      free(basemem);
                      return FAILURE;
                    }
                  // The factor (1-2*ii)*(1-2*jj) gives the appropriate sign in the expression (**) above
                  CEjac[i*evoparsdim + j] += (1 - 2*ii)*(1 - 2*jj)*(rhs[R0index] + 1.0);
                }
              parPntr[indx2] = oldvalue2;                                           // Reset parameter
            }
#endif

          // Reset x*+dpi to its ESS value x*
          parPntr[indx1] = oldvalue1;

#if (COMPUTECROSS == 1)
          /*
           * With y containing the correct environment E(x*+dpi) or E(x*-dpi), compute the cross matrix of the canonical equation
           * using the following expression for its elements:
           *
           * (R0(x*+dpj, E(x*+dpi)) - R0(x*-dpj, E(x*+dpi)) - R0(x*+dpj, E(x*-dpi)) + R0(x*-dpj, E(x*-dpi)))/(4*eps_pi*eps_pj)        (**)
           *
           * where dpi = (0,.....eps_pi,....)^T and dpj = (0,.....eps_pj,....)^T are small perturbation vectors in the direction
           * of the parameter axis i and j, respectively.
           */
          for (j = 0; j < evoparsdim; j++)
            {
              if (detecting && !j)
                indx2 = Bifparone;                                                  // Classification after detecting an ESS in parameter[Bifparone]
              else
                indx2 = evoParsIndexPntr[j - detecting];                            // Runs to i=0..evoParsDim in detecting and non-detecting case

              oldvalue2 = parPntr[indx2];
              for (jj = 0; jj < 2; jj++)
                {
                  if (jj == 0)
                    parPntr[indx2] = oldvalue2 + pardif[j];
                  else
                    parPntr[indx2] = oldvalue2 - pardif[j];
                  if ((*fnc)(y, rhs) == FAILURE)
                    {
                      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
                      setBifParVal = 1;
                      CurveType    = OldCurveType;                                     // Reset the curve type
                      free(basemem);
                      return FAILURE;
                    }
                  // The factor (1-2*ii)*(1-2*jj) gives the appropriate sign in the expression (**) above
                  CEC01[i*evoparsdim + j] += (1 - 2*ii)*(1 - 2*jj)*(rhs[R0index] + 1.0);
                }
              parPntr[indx2] = oldvalue2;                                           // Reset parameter
            }
#endif
        }

      // Scale the elements
      for (j = 0; j < evoparsdim; j++)
        {
          if (detecting && !j)
            indx2 = Bifparone;                                                      // Classification after detecting an ESS in parameter[Bifparone]
          else
            indx2 = evoParsIndexPntr[j - detecting];                                // Runs to i=0..evoParsDim in detecting and non-detecting case

#if (COMPUTEJACOBIAN == 1)
          CEjac[i*evoparsdim + j] /= 4*pardif[i]*pardif[j];
#endif

#if (COMPUTECROSS == 1)
          CEC01[i*evoparsdim + j] /= 4*pardif[i]*pardif[j];
#endif
#if (COMPUTECROSS == 0)
          CEC01[i*evoparsdim + j] = CEjac[i*evoparsdim + j] - CEhes[i*evoparsdim + j];
#endif

#if (COMPUTEJACOBIAN == 0)
          CEjac[i*evoparsdim + j] = CEhes[i*evoparsdim + j] + CEC01[i*evoparsdim + j];
#endif
#if (COMPUTEHESSIAN == 0)
          CEhes[i*evoparsdim + j] = CEjac[i*evoparsdim + j] - CEC01[i*evoparsdim + j];
#endif
        }
    }

  if ((evoparsdim == 1) && !detecting)
    {
      *EVmaxH = CEhes[0];                                                          // Returns R0_yy
#if (COMPUTEHESSIAN == 1)                                                          // Returns R0_xx = C11 = -(C00 + C01 + C10) = -(J + C01) =
      *EVmaxJ = CEhes[0] - 2*CEjac[0];                                             //                       -(J + (J - H)) = -(2*J - H) = (H - 2*J)
#else
      *EVmaxJ = -(CEjac[0] + CEC01[0]);
#endif
      *zC01z = CEC01[0];
    }
  else if (evoparsdim == 1)
    {
      *EVmaxH = CEhes[0];                                                          // Returns R0_yy
      *EVmaxJ = CEjac[0];
#if (COMPUTEHESSIAN == 1)
      *zC01z = (CEhes[0] - 2*CEjac[0])/CEhes[0];                                   // Returns R0_xx/R0_yy, the direction of the curve in the PIP
#else
      *zC01z = -(CEjac[0] + CEC01[0])/(CEjac[0] - CEC01[0]);                       // Returns R0_xx/R0_yy, the direction of the curve in the PIP
#endif
    }
  else
    {
      Eigenval(evoparsdim, CEhes, 1, EVmaxH, DOMINANT, eigvecH, NULL, rhstol);
      Eigenval(evoparsdim, CEjac, 0, EVmaxJ, DOMINANT, eigvecJ, NULL, rhstol);
      for (i = 0; i < evoparsdim; i++)
        {
          /*
          * To check whether branches can coexist near an EBP we have to check Z^T C_01 Z,
          * C_01 the matrix of cross-derivatives to mutant and resident,
          * Z the dominant eigenvector of the Hessiaan. Only when the resulting expression
          * is negative coexistence is possible
          */

          // Compute C_00 Z
          eigvecJ[i] = 0;                                                           // We do not need the eigenvector of the Jacobiaan, use it store the result
          for (j = 0; j < evoparsdim; j++) eigvecJ[i] += CEC01[i*evoparsdim + j]*eigvecH[j];
        }
      // Compute Z^T C_01 Z
      *zC01z = 0.0;
      for (j = 0; j < evoparsdim; j++) *zC01z += eigvecJ[j]*eigvecH[j];
    }
  setBifParVal = 1;
  CurveType    = OldCurveType;                                                      // Reset the curve type

  // Additional call to reset global variables
  if ((*fnc)(y, rhs) == FAILURE)
    {
      ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
      free(basemem);
      return FAILURE;
    }
  ReportMsg("ESS point classification succeeded\n");

  free(basemem);

  return SUCCES;
}


/*==================================================================================================================================*/
