/***
  NAME
    io
  DESCRIPTION
    Implements all low-level I/O routines

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
#ifndef IO
#define IO
#endif
#include "globals.h"

/*==================================================================================================================================*/
/*
 * Start of function implementations.
 */

void ReportMsg(const char *fmt, ...)

{
  va_list argpnt;

  va_start(argpnt, fmt);
  if (errfile)
    {
      vfprintf(errfile, fmt, argpnt);
      (void)fflush(errfile);
    }
  va_end(argpnt);

  return;
}


/*==================================================================================================================================*/

void ErrorMsg(const char *name, const int line, const char *fmt, ...)

{
  va_list argpnt;

#if (defined(R_PACKAGE))
  REprintf("\n** %-12s (%3d): ", name, line);
#elif defined(MATLAB_MEX_FILE)
  mexPrintf("\n** %-12s (%3d): ", name, line);
#else
  (void)fprintf(stderr, "\n** %-12s (%3d): ", name, line);
#endif
  va_start(argpnt, fmt);
#if (defined(R_PACKAGE))
  REvprintf(fmt, argpnt);
#elif defined(MATLAB_MEX_FILE)
  mexPrintf(fmt, argpnt);
#else
  vfprintf(stderr, fmt, argpnt);
#endif
  va_end(argpnt);

#if (defined(R_PACKAGE))
  REprintf("\n");
#elif defined(MATLAB_MEX_FILE)
  mexPrintf("\n");
#else
  (void)fprintf(stderr, "\n");
#endif

  if (errfile)
    {
      (void)fprintf(errfile, "\n** %-12s (%3d): ", name, line);
      va_start(argpnt, fmt);
      vfprintf(errfile, fmt, argpnt);
      va_end(argpnt);
      (void)fprintf(errfile, "\n\n");
    }
  (void)fflush(NULL);
#ifdef MATLAB_MEX_FILE
  mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  return;
}

/*==================================================================================================================================*/

int ReportMemError(const char *name)

{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexWarnMsgIdAndTxt("MATLAB:MemoryError", "\nMemory allocation error in %s(). Exiting MEX module.....\n\n", name);
#else
  warning("\nMemory allocation error in %s(). Exiting module.....\n\n", name);
  R_FlushConsole();
  R_ProcessEvents();
#endif
  CtrlCPressed = TRUE;
#else
  ErrorMsg(__FILE__, __LINE__, "Memory allocation error in %s()!", name);
  exit(1);                                                                          // Only executed when in command-line mode
#endif

  return FAILURE;
}


/*==================================================================================================================================*/

void NumProcError(const char *name, const int line, const int NumProcErrorCode)

{
  switch (NumProcErrorCode)
    {
      case SINGULARITY:
        ErrorMsg(name, line, "%-45s", "Singular matrix encountered");
        break;
      case NORM_OVERFLOW:
        ErrorMsg(name, line, "%-45s", "Norm overflow in Newton iteration");
        break;
      case NO_CONVERGENCE:
        ErrorMsg(name, line, "%-45s", "No convergence in Newton iteration");
        break;
      case ILLEGAL_INPUT:
      case FAILED_EVALUATION:
        break;
      default:
        ErrorMsg(name, line, "%-45s%d", "Unknown numerical error code: ", NumProcErrorCode);
        break;
    }

  return;
}


/*==================================================================================================================================*/

void PrettyPrintArray(FILE *fp, const int dim, double *vec)

{
  register int i;
  double       tmp;

  for (i = 0; i < dim; i++)
    {
      tmp = vec[i];
      if (((fabs(tmp) <= 1.0E4) && (fabs(tmp) >= 1.0E-3)) || (tmp == 0.))
        {
          if (fp)
            (void)fprintf(fp, "%16.8f", tmp);
          else
            (void)STDOUT("%16.8f", tmp);
        }
      else
        {
          if (fp)
            (void)fprintf(fp, "%16.8E", tmp);
          else
            (void)STDOUT("%16.8E", tmp);
        }
    }
  if (fp)
    (void)fprintf(fp, "\n");
  else
    (void)STDOUT("\n");
  (void)fflush(fp);
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  return;
}


/*==================================================================================================================================*/
// Ctrl-C detection

#if defined(MATLAB_MEX_FILE)
int checkInterrupt()
{
  int         pressed;
  extern bool utIsInterruptPending();
  extern void utSetInterruptPending(bool);

  // check for a Ctrl-C event
  pressed = utIsInterruptPending();
  if (pressed)
    {
      utSetInterruptPending(false);
      mexPrintf("\n\nCtrl-C detected. Stopping computation\n\n");
      CtrlCPressed = TRUE;
    }

  return (pressed || CtrlCPressed);
}

#elif defined(OCTAVE_MEX_FILE)

int checkInterrupt()
{
  int pressed = 0;

  // checking of Ctrl-C event not implemented

  return pressed;
}

#elif defined(R_PACKAGE)

static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }

// this will call the above in a top-level context so it won't longjmp-out of your context
int checkInterrupt()
{
  Rboolean R_ToplevelExec(void (*fun)(void *), void *data);

  if (CtrlCPressed) return CtrlCPressed;

#if (defined(OPENMP))
  if (omp_get_thread_num() == 0)
#endif
    {
      CtrlCPressed = (R_ToplevelExec(chkIntFn, NULL) == FALSE);
      if (CtrlCPressed) REprintf("\n\nUser interrupt detected. Stopping computation\n\n");
    }

  return (CtrlCPressed);
}

#endif

/*==================================================================================================================================*/

char *ReadDouble(double *val, char *cpnt)

/*
   * ReadDouble - Routine reads a double value from the string pointed to
   *              by "cpnt". Invalid characters are skipped. It returns a
   *              pointer to the rest of the string or NULL on error.
   */

{
  register char *ch, *end = NULL;
  int           dot_start = 0;

  ch = cpnt;
  while ((!isdigit(*ch)) &&*ch) ch++;                                               // Skip non digits
  if (isdigit(*ch))
    {
      end = ch;
      if ((ch != cpnt) && (*(ch - 1) == '.'))                                       // Is previous a dot?
        {
          ch--;
          dot_start = 1;
        }
      if ((ch != cpnt) && (*(ch - 1) == '-')) ch--;                                 // Is previous a minus?

      while (isdigit(*end)) end++;                                                  // Skip digits
      if ((!dot_start) && (*end == '.'))                                            // Dot starts mantissa
        {
          end++;
          while (isdigit(*end)) end++;
        }

      if ((*end == 'e') || (*end == 'E'))                                           // Possible exponent
        {
          end++;
          if ((*end == '+') || (*end == '-')) end++;
          while (isdigit(*end)) end++;
        }
      *val = 0.0;
      (void)sscanf(ch, "%lg", val);
    }

  return end;
}


/*==================================================================================================================================*/

int ScanLine(FILE *infile, char *descr, double *value, int var_nr)

/*
   * ScanLine - Routine reads one line from the control variable file,
   *            splits it up in the description and the value part and
   *            returns the number of values read. "var_nr" is the number
   *            of variables to be read.
   */

{
  register char *ch = NULL, *dsp;
  char           input[1024], tmp_str[1024];
  int            read_no = 0;

  dsp = descr;                                                                      // Input line
  while ((!feof(infile)) && (!ferror(infile)) && (!ch))                             // Skip empty lines
    {
      if ((ch = fgets(input, 1024, infile)))
        {
          while (*ch == ' ') ch++;
          if (*ch == '\n') ch = NULL;
        }
    }
  if (feof(infile) || ferror(infile)) return 0;

  if (*ch == '"')                                                                   // Comment is in quotes
    {
      ch++;                                                                         // Skip the first quote
      while (*ch != '"')                                                            // Search for the closing
        {                                                                           // quotes, storing string
          (*dsp) = (*ch);
          dsp++;
          ch++;                                                                     // between them
        }
      ch++;                                                                         // Skip closing quotes
    }
  else                                                                              // Comment is not quoted
    {
      while ((!isdigit(*ch)) &&*ch)                                                 // Skip non-digits and put
        {                                                                           // them in description
          (*dsp) = (*ch);
          dsp++;
          ch++;
        }                                                                           // Remove trailing blanks
      // and tabs
      while ((*(dsp - 1) == ' ') || (*(dsp - 1) == '\t')) *(--dsp) = '\0';
    }

  ch = strcpy(tmp_str, ch);                                                         // Copy rest to string

  while (ch)
    if ((ch = ReadDouble(value + read_no, ch)) != NULL)
      if ((++read_no) == var_nr) return read_no;

  return read_no;
}


/*==================================================================================================================================*/

int ReadCvfFile(const char *fname)

// Routine reads parameter values from the .cvf file.

{
  double        tmp[1];
  register int  j;
  char          input[1024];
  int           read_no;
  FILE          *infile;

  // Add cvf extension to file name and open
  strcpy(input, fname);
  if (!strstr(input, ".cvf")) (void)strcat(input, ".cvf");
  infile = fopen(input, "r");
  if (!infile) return FAILURE;

  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD accuracy
  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD cohort_limit

  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD equal to 0

  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD max_time
  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD delt_out
  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD state_out

  read_no = ScanLine(infile, input, tmp, 1);                                        // DISCARD extinction treshold

  for (j = 0; j < IStateDim; j++) read_no = ScanLine(infile, input, tmp, 1);        // DISCARD rel_tols[][j]
  for (j = 0; j < IStateDim; j++) read_no = ScanLine(infile, input, tmp, 1);        // DISCARD abs_tols[][j]

  for (j = 0; j < ParameterNr; j++)
    {
      read_no = ScanLine(infile, input, parPntr + j, 1);
      if (read_no != 1)
        {
          fprintf(stderr, "Error while reading parameter %d from CVF file!", j);
          return FAILURE;
        }
    }
  (void)fclose(infile);

  return SUCCES;
}


/*==================================================================================================================================*/

/* In principle, it could be of interest to have the full density profile, average i-state & reproductive value written to the
 * state output file for every birth type separately. With large numbers of birth types this will increase the size of the output
 * a lot. Hence, the average i-state & reproductive value is calculated by weighing the contribution of the different birth types
 * with the right eigenvector
 */

static void CondenseStableDist(void)

{
  register int i, b, j, p, iloopmax;
  double       tmpval, cohortdens;

  iloopmax = IStateDim + 2*(CurveType == PGR);
  for (p = 0; p < CurPopulationNr; p++)
    {
      for (j = 0; j < CohortNr; j++)
        {
          // Compute the average cohort density
          for (b = 0, cohortdens = 0; b < birthStateNr[p]; b++)
            {
              if (RightEigenvec(p, b) < epsMach) continue;
              cohortdens += RightEigenvec(p, b)*PopDens(p, b, IStateDim, j);
            }

          for (i = 0; i < iloopmax; i++)
            {
              if (i == IStateDim) continue;                                         // Skip the density element
              if (cohortdens)
                {
                  for (b = 0, tmpval = 0; b < birthStateNr[p]; b++)
                    {
                      if (RightEigenvec(p, b) < epsMach) continue;
                      tmpval += RightEigenvec(p, b)*PopDens(p, b, IStateDim, j)*PopDens(p, b, i, j);
                    }
                  PopDens(p, 0, i, j) = tmpval/cohortdens;
                }
              else
                PopDens(p, 0, i, j) = 0.0;
            }

          // Now store the average cohort density
          PopDens(p, 0, IStateDim, j) = cohortdens;
        }
      for (b = 1; b < birthStateNr[p]; b++)
        {
          if (RightEigenvec(p, b) < epsMach) continue;
          Cohorts(p, 0) = max(Cohorts(p, 0), Cohorts(p, b));
        }
    }

  return;
}


/*==================================================================================================================================*/

#define FIELDNAMELN               40

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)                            // Matlab or Octave on Windows or Mac OS

#if defined(MATLAB_MEX_FILE)
#include "matrix.h"
#include "mat.h"

MATFile                           *pmat = NULL;
#endif

void WriteStateToFile(const int fullstateoutput)
{
  int         i, b, j, p, nalloc = 0, status;
  short int   bmag, pmag, num;
  static int  skip = 0;
  mxArray     *mystruct;
  mxArray     **mydata;
  const char  **fieldnames;
  char        matfname[MAXPATHLEN];
  char        tmpstr[FIELDNAMELN], tmpstr2[FIELDNAMELN];
  double      *dblptr;
#if defined(OCTAVE_MEX_FILE)
  mxArray     *prhs[4];
#endif

  if (skip || !fullstateoutput) return;

  if (fullstateoutput == 2)
    {
      mydata     = malloc((5 + CurPopulationNr*(MaxStatesAtBirth*CohortDim + 1))*sizeof(mxArray *));
      fieldnames = malloc((5 + CurPopulationNr*(MaxStatesAtBirth*CohortDim + 1))*sizeof(char *));

#if ((PSPMDEMO != 1) && (PSPMIND != 1))
      // Scale the cohort densities with the stable birth rate distribution over the birth states (right eigenvector)
      for (p = 0; p < CurPopulationNr; p++)
        {
          if (birthStateNr[p] == 1) continue;
          for (i = 0; i < CohortNr; i++)
            for (b = 0; b < birthStateNr[p]; b++)
              PopDens(p, b, IStateDim, i) *= RightEigenvec(p, b);
        }
#endif
    }
  else
    {
      if (MaxStatesAtBirth > 1) CondenseStableDist();

      mydata     = malloc((10 + CurPopulationNr*(CohortDim + 1))*sizeof(mxArray *));
      fieldnames = malloc((10 + CurPopulationNr*(CohortDim + 1))*sizeof(char *));
    }
  if (!mydata || !fieldnames)
    {
      ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
      if (mydata) free(mydata);
      if (fieldnames) free(fieldnames);
      skip = 1;
      return;
    }

  nalloc = 0;
  if (CurveType == EVODYN)
    {
      mydata[nalloc]     = mxCreateDoubleMatrix(1, 1, mxREAL);
      fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
      if (!mydata[nalloc] || !fieldnames[nalloc])
        {
          ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
          for (j = 0; j <= nalloc; j++)
            {
              if (mydata[j]) mxFree(mydata[j]);
              if (fieldnames[j]) mxFree((char *)fieldnames[j]);
            }
          skip = 1;
          free(mydata);
          free(fieldnames);
          return;
        }
      memcpy((void *)fieldnames[nalloc], "EvoTime", sizeof("EvoTime"));
      dblptr = mxGetPr(mydata[nalloc]);
      memcpy(dblptr, timePntr, sizeof(double));
      nalloc++;
    }

  if (Bifparone != -1)
    {
      if ((CurveType == PGR) || (CurveType == EQ) || (CurveType == ESS))
        mydata[nalloc] = mxCreateDoubleMatrix(1, 1, mxREAL);
      else
        mydata[nalloc]   = mxCreateDoubleMatrix(1, 2, mxREAL);
      fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
      if (!mydata[nalloc] || !fieldnames[nalloc])
        {
          ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
          for (j = 0; j <= nalloc; j++)
            {
              if (mydata[j]) mxFree(mydata[j]);
              if (fieldnames[j]) mxFree((char *)fieldnames[j]);
            }
          skip = 1;
          free(mydata);
          free(fieldnames);
          return;
        }
      memcpy((void *)fieldnames[nalloc], "BifPars", sizeof("BifPars"));
      dblptr    = mxGetPr(mydata[nalloc]);
      dblptr[0] = parPntr[Bifparone];

      if (CurveType == PIP)
        dblptr[1] = MutantParVal;
      else if ((CurveType == BP) || (CurveType == LP) || (CurveType == BPE))
        dblptr[1] = parPntr[Bifpartwo];
      nalloc++;
    }

  if ((CurveType == EVODYN) || (CurveType == ESS))
    {
      mydata[nalloc]     = mxCreateDoubleMatrix(1, evoParsDim, mxREAL);
      fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
      if (!mydata[nalloc] || !fieldnames[nalloc])
        {
          ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
          for (j = 0; j <= nalloc; j++)
            {
              if (mydata[j]) mxFree(mydata[j]);
              if (fieldnames[j]) mxFree((char *)fieldnames[j]);
            }
          skip = 1;
          free(mydata);
          free(fieldnames);
          return;
        }
      memcpy((void *)fieldnames[nalloc], "EvoPars", sizeof("EvoPars"));
      dblptr = mxGetPr(mydata[nalloc]);
      for (i = 0; i < evoParsDim; i++) dblptr[i] = parPntr[evoParsIndexPntr[i]];
      nalloc++;
    }

  mydata[nalloc]     = mxCreateDoubleMatrix(1, ParameterNr, mxREAL);
  fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
  if (!mydata[nalloc] || !fieldnames[nalloc])
    {
      ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
      for (j = 0; j <= nalloc; j++)
        {
          if (mydata[j]) mxFree(mydata[j]);
          if (fieldnames[j]) mxFree((char *)fieldnames[j]);
        }
      skip = 1;
      free(mydata);
      free(fieldnames);
      return;
    }
  memcpy((void *)fieldnames[nalloc], "Parameters", sizeof("Parameters"));
  dblptr = mxGetPr(mydata[nalloc]);
  memcpy(dblptr, parPntr, ParameterNr*sizeof(double));
  nalloc++;

  mydata[nalloc]     = mxCreateDoubleMatrix(1, EnvironDim, mxREAL);
  fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
  if (!mydata[nalloc] || !fieldnames[nalloc])
    {
      ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
      for (j = 0; j <= nalloc; j++)
        {
          if (mydata[j]) mxFree(mydata[j]);
          if (fieldnames[j]) mxFree((char *)fieldnames[j]);
        }
      skip = 1;
      free(mydata);
      free(fieldnames);
      return;
    }
  if (CurveType == PGR)
    memcpy((void *)fieldnames[nalloc], "PGR", sizeof("PGR"));
  else
    memcpy((void *)fieldnames[nalloc], "Environment", sizeof("Environment"));
  dblptr = mxGetPr(mydata[nalloc]);
  memcpy(dblptr, eVarPntr, EnvironDim*sizeof(double));
  nalloc++;

  pmag = 1;
  num  = CurPopulationNr;
  while (num > 0)
    {
      pmag++;
      num = num/10;
    }
  for (p = 0, b = 0; p < CurPopulationNr; p++, b = 0)
    {
      if (CurveType != IND)
        {
          // Store the right eigenvector as a separate population.
          if (birthStateNr[p] > 1)
            {
              mydata[nalloc]     = mxCreateDoubleMatrix(1, birthStateNr[p], mxREAL);
              fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
              if (!mydata[nalloc] || !fieldnames[nalloc])
                {
                  ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
                  for (j = 0; j <= nalloc; j++)
                    {
                      if (mydata[j]) mxFree(mydata[j]);
                      if (fieldnames[j]) mxFree((char *)fieldnames[j]);
                    }
                  free(mydata);
                  free(fieldnames);
                  skip = 1;
                  return;
                }
              dblptr = mxGetPr(mydata[nalloc]);
              memcpy(dblptr, rightEigenvecMem + p*MaxStatesAtBirth, birthStateNr[p]*sizeof(double));
              sprintf(tmpstr, "Pop%0*d_StableBirthDist", pmag, p);
              memcpy((void *)fieldnames[nalloc], tmpstr, strlen(tmpstr)*sizeof(char));
              nalloc++;
            }
        }
      // Store all birth states of the population
      mydata[nalloc]     = mxCreateDoubleMatrix(birthStateNr[p], IStateDim, mxREAL);
      fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
      if (!mydata[nalloc] || !fieldnames[nalloc])
        {
          ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
          for (j = 0; j <= nalloc; j++)
            {
              if (mydata[j]) mxFree(mydata[j]);
              if (fieldnames[j]) mxFree((char *)fieldnames[j]);
            }
          free(mydata);
          free(fieldnames);
          skip = 1;
          return;
        }
      dblptr = mxGetPr(mydata[nalloc]);
      for (i = 0; i < IStateDim; i++)
        for (b = 0; b < birthStateNr[p]; b++) memcpy(dblptr++, BirthStatePnt(p, b, i), sizeof(double));
      sprintf(tmpstr, "Pop%0*d_BirthStates", pmag, p);
      memcpy((void *)fieldnames[nalloc], tmpstr, strlen(tmpstr)*sizeof(char));
      nalloc++;

      // Store all cohorts originating from all birth states as separate populations
      bmag = 1;
      num  = birthStateNr[p];
      while (num > 0)
        {
          bmag++;
          num = num/10;
        }
      for (b = 0; b < birthStateNr[p]; b++)
        {
          if (!Cohorts(p, b)) continue;
          if (CurveType == IND)
            mydata[nalloc] = mxCreateDoubleMatrix(Cohorts(p, b), CohortDim + birthStateNr[p], mxREAL);
          else
            mydata[nalloc]   = mxCreateDoubleMatrix(Cohorts(p, b), CohortDim, mxREAL);
          fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
          if (!mydata[nalloc] || !fieldnames[nalloc])
            {
              ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteStateToFile(). Further state output suppressed!");
              for (j = 0; j <= nalloc; j++)
                {
                  if (mydata[j]) mxFree(mydata[j]);
                  if (fieldnames[j]) mxFree((char *)fieldnames[j]);
                }
              skip = 1;
              free(mydata);
              free(fieldnames);
              return;
            }
          dblptr = mxGetPr(mydata[nalloc]);
          memcpy(dblptr, &(PopDens(p, b, IStateDim, 0)), Cohorts(p, b)*sizeof(double));
          dblptr += Cohorts(p, b);
          for (i = 0; i < IStateDim; i++, dblptr += Cohorts(p, b)) memcpy(dblptr, &(PopDens(p, b, i, 0)), Cohorts(p, b)*sizeof(double));
          if (CurveType == PGR)
            memcpy(dblptr, &(PopDens(p, b, IStateDim + 1, 0)), Cohorts(p, b)*sizeof(double));
          else if (CurveType == IND)
            for (i = IStateDim + 1; i < (CohortDim + birthStateNr[p]); i++, dblptr += Cohorts(p, b))
              memcpy(dblptr, &(PopDens(p, b, i, 0)), Cohorts(p, b)*sizeof(double));
          if ((fullstateoutput == 2) && (MaxStatesAtBirth > 1))
            sprintf(tmpstr, "Pop%0*d_Bstate%0*d", pmag, p, bmag, b);
          else
            sprintf(tmpstr, "Pop%0*d", pmag, p);
          memcpy((void *)fieldnames[nalloc], tmpstr, strlen(tmpstr)*sizeof(char));
          nalloc++;
          if (fullstateoutput == 1) break;
        }
    }

  mystruct = mxCreateStructMatrix(1, 1, nalloc, (const char **)fieldnames);

  if (mystruct)
    {
      for (j = 0; j < nalloc; j++) mxSetFieldByNumber(mystruct, 0, j, mydata[j]);
      if ((Bifparone != -1) || (CurveType == EVODYN))
        {
          if (CurveType == EVODYN)
            sprintf(tmpstr, "State_%.6E", *timePntr);
          else
            sprintf(tmpstr, "State_%.6E", parPntr[Bifparone]);
          for (i = 0, j = 0; i < strlen(tmpstr); i++)
            {
              if (tmpstr[i] == '.')
                tmpstr2[j++] = '_';
              else if (tmpstr[i] == '+')
                ;
              else if (tmpstr[i] == '-')
                tmpstr2[j++] = '_';
              else
                tmpstr2[j++] = tmpstr[i];
            }
          tmpstr2[j] = '\0';
        }
      else if (CurveType == IND)
        sprintf(tmpstr2, "LifeHistory");
      else
        sprintf(tmpstr2, "State");
      sprintf(matfname, "%s.mat", runname);

#if defined(MATLAB_MEX_FILE)
      if (pmat == NULL) pmat = matOpen(matfname, "wz");
      if (pmat == NULL)
        {
          ErrorMsg(__FILE__, __LINE__, "Can not open file in WriteStateToFile(). Further state output suppressed!");
          for (j = 0; j < nalloc; j++)
            {
              if (mydata[j]) mxFree(mydata[j]);
              if (fieldnames[j]) mxFree((char *)fieldnames[j]);
            }
          skip = 1;
          free(mydata);
          free(fieldnames);
          return;
        }

      status = matPutVariableAsGlobal(pmat, tmpstr2, mystruct);
      if (status != 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Failed to write to file in WriteStateToFile(). Further state output suppressed!");
          skip = 1;
        }
#else
      status = mexPutVariable("caller", tmpstr2, mystruct);                         // Put the structure in the workspace
      if (status != 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Failed to put state into global workspace WriteStateToFile(). Further state output suppressed!");
          skip = 1;
        }
      else
        {
          prhs[0] = mxCreateString("-append");
          prhs[1] = mxCreateString("-V7");
          prhs[2] = mxCreateString(matfname);
          prhs[3] = mxCreateString(tmpstr2);
          mexCallMATLAB(0, NULL, 4, prhs, "save");                                  // Save the structure to file
          mexCallMATLAB(0, NULL, 1, prhs + 3, "clear");                             // Clear the structure from the workspace
          mxFree(prhs[0]);
          mxFree(prhs[1]);
          mxFree(prhs[2]);
          mxFree(prhs[3]);
        }
#endif
    }
  else
    {
      ErrorMsg(__FILE__, __LINE__, "Failed to create output structure in WriteStateToFile(). Further state output suppressed!");
      skip = 1;
    }

  for (j = 0; j < nalloc; j++)
    {
      // Matlab deallocates the data arrays when calling mxDestroyArray(). Deallocating it on beforehand makes Matlab crash (but not Octave)
      // if (mydata[j]) mxFree(mydata[j]);
      if (fieldnames[j]) mxFree((char *)fieldnames[j]);
    }
  mxDestroyArray(mystruct);
  free(mydata);
  free(fieldnames);

  return;
}


/*==================================================================================================================================*/

#else

// Magic key of the type of CSB file written by this module
#define CSB_MAGIC_KEY             20030509
#define statelabels(p)            (statelblmem + p*MAX_LBL_LEN*sizeof(double))

typedef struct envdim
{
  double   timeval;
  int      columns;
  int      data_offset;
  uint32_t memory_used;
} Envdim;

typedef struct popdim
{
  double timeval;
  int    population;
  int    cohorts;
  int    columns;
  int    data_offset;
  int    lastpopdim;
} Popdim;

void WriteStateToFile(const int fullstateoutput)

/*
 * WriteBinStateToFile  - Routine writes the entire state of the populations
 *                        and the environment in binary format to the named
 *                        file (with csb extension) after updating the
 *                        population state using the procedure (*updatepop)().
 */

{
  register int  i, b, p, totpop;
  size_t        hdrdbls;
  Envdim        cenv[2];
  Popdim        cpop[2];
  double        zero = 0.0, tmpval;
  uint32_t      tmpuint32;
  int           tmpint, num;
  char          input[MAX_STR_LEN], tmpstr[MAX_STR_LEN];
  char          *statelblmem;
  char          envlabel[MAX_LBL_LEN*sizeof(double)];
  FILE          *fp;
  struct stat   st;
  static int    skip = 0;
  short int     pmag, bmag;

  if (skip || !fullstateoutput) return;

  statelblmem = (char *)malloc(CurPopulationNr*MAX_LBL_LEN*sizeof(double));
  if (!statelblmem)
    {
      ErrorMsg(__FILE__, __LINE__, "Memory allocation error in WriteBinStateToFile(). Further state output suppressed!");
      skip = 1;
      return;
    }

  if ((fullstateoutput == 1) && (MaxStatesAtBirth > 1)) CondenseStableDist();
#if ((PSPMDEMO != 1) && (PSPMIND != 1))
  else if (fullstateoutput == 2)
    {
      // Scale the cohort densities with the stable birth rate distribution over the birth states (right eigenvector)
      for (p = 0; p < CurPopulationNr; p++)
        {
          if (birthStateNr[p] == 1) continue;
          for (i = 0; i < CohortNr; i++)
            for (b = 0; b < birthStateNr[p]; b++)
              PopDens(p, b, IStateDim, i) *= RightEigenvec(p, b);
        }
    }
#endif

  // Add csb extension to file name if required and open the file
  strcpy(input, runname);
  if (strcmp(input + strlen(input) - 4, ".csb")) (void)strcat(input, ".csb");

  if ((stat(input, &st) != 0) || (st.st_size == 0))
    {
      // File does not exist yet. Write magic key and parameters
      fp = fopen(input, "wb");
      if (!fp)
        {
          ErrorMsg(__FILE__, __LINE__, "Unable to open CSB file!");
          skip = 1;
          free(statelblmem);
          return;
        }

      tmpuint32 = CSB_MAGIC_KEY;
      fwrite((void *)(&tmpuint32), 1, sizeof(uint32_t), fp);

      if (((CurveType == PGR) && (Bifparone == -1)) || (CurveType == IND))
        {
          tmpint = ParameterNr + 1;
          fwrite((void *)(&tmpint), 1, sizeof(int), fp);
          fwrite((void *)(parPntr), ParameterNr, sizeof(double), fp);

          tmpval = 0;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
        }
      else if ((CurveType == PGR) || (CurveType == EQ))
        {
          tmpint = ParameterNr + 2;
          fwrite((void *)(&tmpint), 1, sizeof(int), fp);
          fwrite((void *)(parPntr), ParameterNr, sizeof(double), fp);

          tmpval = Bifparone;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
          tmpval = 1;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
        }
      else if (CurveType == ESS)
        {
          tmpint = ParameterNr + evoParsDim + 2;
          fwrite((void *)(&tmpint), 1, sizeof(int), fp);
          fwrite((void *)(parPntr), ParameterNr, sizeof(double), fp);

          tmpval = Bifparone;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
          for (i = 0; i < evoParsDim; i++)
            {
              tmpval = evoParsIndexPntr[i];
              fwrite((void *)(&tmpval), 1, sizeof(double), fp);
            }
          tmpval = evoParsDim + 1;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
        }
      else if (CurveType == EVODYN)
        {
          tmpint = ParameterNr + evoParsDim + 1;
          fwrite((void *)(&tmpint), 1, sizeof(int), fp);
          fwrite((void *)(parPntr), ParameterNr, sizeof(double), fp);

          for (i = 0; i < evoParsDim; i++)
            {
              tmpval = evoParsIndexPntr[i];
              fwrite((void *)(&tmpval), 1, sizeof(double), fp);
            }
          tmpval = evoParsDim;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
        }
      else
        {
          tmpint = ParameterNr + 3;
          fwrite((void *)(&tmpint), 1, sizeof(int), fp);
          fwrite((void *)(parPntr), ParameterNr, sizeof(double), fp);

          tmpval = Bifparone;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
          tmpval = Bifpartwo;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
          tmpval = 2;
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
        }
    }
  else
    {
      fp = fopen(input, "ab");
      if (!fp)
        {
          ErrorMsg(__FILE__, __LINE__, "Unable to open CSB file!");
          skip = 1;
          free(statelblmem);
          return;
        }
    }

  hdrdbls = (sizeof(Envdim)/sizeof(double)) + 1;
  (void)memset((void *)cenv, 0, hdrdbls*sizeof(double));

  if ((CurveType == PGR) && (Bifparone == -1))
    {
      cenv->timeval = eVarPntr[0];
      cenv->columns = (CurPopulationNr);
      sprintf(envlabel, "PGR: Population growth rates");
    }
  else if (CurveType == PGR)
    {
      cenv->timeval = parPntr[Bifparone];
      cenv->columns = (CurPopulationNr + 1);
      sprintf(envlabel, "PGR: Population growth rates for parameter value %G", parPntr[Bifparone]);
    }
  else if (CurveType == IND)
    {
      cenv->timeval = eVarPntr[0];
      cenv->columns = (EnvironDim);
      sprintf(envlabel, "Individual life history dynamics under given environmental conditions");
    }
  else if (CurveType == EQ)
    {
      cenv->timeval = parPntr[Bifparone];
      cenv->columns = (EnvironDim + 1);
      sprintf(envlabel, "Environment: Environment variables for parameter value %G", parPntr[Bifparone]);
    }
  else if (CurveType == ESS)
    {
      cenv->timeval = parPntr[Bifparone];
      cenv->columns = (EnvironDim + 1 + evoParsDim);
      sprintf(envlabel, "Environment: Environment variables and evolutionary parameters for parameter value %G", parPntr[Bifparone]);
    }
  else if (CurveType == EVODYN)
    {
      cenv->timeval =*timePntr;
      cenv->columns = (EnvironDim + evoParsDim);
      sprintf(envlabel, "Environment: Environment variables for evolutionary time value %G", *timePntr);
    }
  else if (CurveType == PIP)
    {
      cenv->timeval = parPntr[Bifparone];
      cenv->columns = (EnvironDim + 2);
      sprintf(envlabel, "Environment: Environment variables for parameter values %G and %G", parPntr[Bifparone], MutantParVal);
    }
  else
    {
      cenv->timeval = parPntr[Bifparone];
      cenv->columns = (EnvironDim + 2);
      sprintf(envlabel, "Environment: Environment variables for parameter values %G and %G", parPntr[Bifparone], parPntr[Bifpartwo]);
    }
  cenv->data_offset = hdrdbls + MAX_LBL_LEN;

  // Compute the memory to be written to file
  cenv->memory_used = (hdrdbls + MAX_LBL_LEN)*sizeof(double);
  cenv->memory_used += (cenv->columns)*sizeof(double);

  hdrdbls = (sizeof(Popdim)/sizeof(double)) + 1;
  for (p = 0; p < CurPopulationNr; p++)
    {
      if (CurveType != IND)                                                         // In case of life history dynamics stable birth distribution not computed
        {
          if (birthStateNr[p] > 1)
            {
              cenv->memory_used += (hdrdbls + MAX_LBL_LEN)*sizeof(double);
              cenv->memory_used += birthStateNr[p]*sizeof(double);                  // Store the stable birth distribution of every population as a separate population
            }
        }

      cenv->memory_used += (hdrdbls + MAX_LBL_LEN)*sizeof(double);
      cenv->memory_used += birthStateNr[p]*IStateDim*sizeof(double);                // Store the birth states of every population as a separate population
      for (b = 0; b < birthStateNr[p]; b++)
        {
          cenv->memory_used += (hdrdbls + MAX_LBL_LEN)*sizeof(double);
          if (CurveType == IND)
            cenv->memory_used += max(Cohorts(p, b), 1)*(CohortDim + birthStateNr[p])*sizeof(double);
          else
            cenv->memory_used += max(Cohorts(p, b), 1)*CohortDim*sizeof(double);
          if (fullstateoutput == 1) break;                                          // Store only b=0 in case of condensed output
        }
    }

  hdrdbls = (sizeof(Envdim)/sizeof(double)) + 1;
  fwrite((void *)cenv, 1, hdrdbls*sizeof(double), fp);
  fwrite((void *)envlabel, 1, MAX_LBL_LEN*sizeof(double), fp);

  for (i = 0; i < EnvironDim; i++)
    {
      tmpval = (double)eVarPntr[i];
      fwrite((void *)(&tmpval), 1, sizeof(double), fp);
    }

  if (Bifparone != -1)
    {
      tmpval = (double)parPntr[Bifparone];
      fwrite((void *)(&tmpval), 1, sizeof(double), fp);
    }
  if (CurveType == PIP)
    {
      tmpval = MutantParVal;
      fwrite((void *)(&tmpval), 1, sizeof(double), fp);
    }
  else if ((CurveType == ESS) || (CurveType == EVODYN))
    {
      for (i = 0; i < evoParsDim; i++)
        {
          tmpval = (double)parPntr[evoParsIndexPntr[i]];
          fwrite((void *)(&tmpval), 1, sizeof(double), fp);
        }
    }
  else if ((CurveType == BP) || (CurveType == LP) || (CurveType == BPE))
    {
      tmpval = (double)parPntr[Bifpartwo];
      fwrite((void *)(&tmpval), 1, sizeof(double), fp);
    }
  pmag = 1;
  num  = CurPopulationNr;
  while (num > 0)
    {
      pmag++;
      num = num/10;
    }

  for (p = 0, totpop = 0; p < CurPopulationNr; p++)
    {
      if (CurveType != IND)                                                         // In case of life history dynamics only store time course of all variables
        {
          // Store the right eigenvector as a separate population.
          if (birthStateNr[p] > 1)
            {
              hdrdbls = (sizeof(Popdim)/sizeof(double)) + 1;
              (void)memset((void *)cpop, 0, hdrdbls*sizeof(double));
              sprintf(tmpstr, "Pop%0*d_StableBirthDist", pmag, p);
              if ((CurveType == PGR) && (Bifparone == -1))
                sprintf(statelabels(p), "%s: Stable birth state distribution of population #%d", tmpstr, p);
              else if ((CurveType == PGR) || (CurveType == EQ) || (CurveType == ESS))
                sprintf(statelabels(p), "%s: Stable birth state distribution of population #%d for parameter value %G", tmpstr, p,
                        parPntr[Bifparone]);
              else if (CurveType == EVODYN)
                sprintf(statelabels(p), "%s: Stable birth state distribution of population #%d for evolutionary time value %G", tmpstr, p, *timePntr);
              else if (CurveType == PIP)
                sprintf(statelabels(p), "%s: Stable birth state distribution of population #%d for parameter values %G and %G", tmpstr, p,
                        parPntr[Bifparone], MutantParVal);
              else
                sprintf(statelabels(p), "%s: Stable birth state distribution of population #%d for parameter values %G and %G", tmpstr, p,
                        parPntr[Bifparone], parPntr[Bifpartwo]);

              if (Bifparone != -1)
                cpop->timeval = parPntr[Bifparone];
              else if (CurveType == EVODYN)
                cpop->timeval =*timePntr;
              else
                cpop->timeval = 0.0;

              cpop->population  = totpop++;
              cpop->columns     = birthStateNr[p];
              cpop->cohorts     = 1;
              cpop->data_offset = hdrdbls + MAX_LBL_LEN;
              cpop->lastpopdim  = 0;
              fwrite((void *)cpop, 1, hdrdbls*sizeof(double), fp);
              fwrite((void *)statelabels(p), 1, MAX_LBL_LEN*sizeof(double), fp);
              fwrite((void *)(&(RightEigenvec(p, 0))), birthStateNr[p], sizeof(double), fp);
            }
        }

      // Store all birth states as a separate population
      hdrdbls = (sizeof(Popdim)/sizeof(double)) + 1;
      (void)memset((void *)cpop, 0, hdrdbls*sizeof(double));
      sprintf(tmpstr, "Pop%0*d_BirthStates", pmag, p);
      if (((CurveType == PGR) && (Bifparone == -1)) || (CurveType == IND))
        sprintf(statelabels(p), "%s: Birth states of population #%d", tmpstr, p);
      else if ((CurveType == PGR) || (CurveType == EQ) || (CurveType == ESS))
        sprintf(statelabels(p), "%s: Birth states of population #%d for parameter value %G", tmpstr, p, parPntr[Bifparone]);
      else if (CurveType == EVODYN)
        sprintf(statelabels(p), "%s: Birth states of population #%d for evolutionary time value %G", tmpstr, p, *timePntr);
      else if (CurveType == PIP)
        sprintf(statelabels(p), "%s: Birth states of population #%d for parameter values %G and %G", tmpstr, p, parPntr[Bifparone], MutantParVal);
      else
        sprintf(statelabels(p), "%s: Birth states of population #%d for parameter values %G and %G", tmpstr, p, parPntr[Bifparone],
                parPntr[Bifpartwo]);

      if (Bifparone != -1)
        cpop->timeval = parPntr[Bifparone];
      else if (CurveType == EVODYN)
        cpop->timeval =*timePntr;
      else
        cpop->timeval = 0.0;

      cpop->population  = totpop++;
      cpop->columns     = IStateDim;
      cpop->cohorts     = birthStateNr[p];
      cpop->data_offset = hdrdbls + MAX_LBL_LEN;
      cpop->lastpopdim  = 0;
      fwrite((void *)cpop, 1, hdrdbls*sizeof(double), fp);
      fwrite((void *)statelabels(p), 1, MAX_LBL_LEN*sizeof(double), fp);
      for (i = 0; i < IStateDim; i++)
        for (b = 0; b < birthStateNr[p]; b++) fwrite((void *)BirthStatePnt(p, b, i), 1, sizeof(double), fp);

      // Store all cohorts originating from all birth states as separate populations
      bmag = 1;
      num  = birthStateNr[p];
      while (num > 0)
        {
          bmag++;
          num = num/10;
        }
      for (b = 0; b < birthStateNr[p]; b++)
        {
          hdrdbls = (sizeof(Popdim)/sizeof(double)) + 1;
          (void)memset((void *)cpop, 0, hdrdbls*sizeof(double));
          if ((fullstateoutput == 2) && (MaxStatesAtBirth > 1))
            {
              sprintf(tmpstr, "Pop%0*d_Bstate%0*d", pmag, p, bmag, b);
              if (CurveType == IND)
                sprintf(statelabels(p), "%s: Life history dynamics in population #%d with birth state #%0*d", tmpstr, p, bmag, b);
              else if ((CurveType == PGR) && (Bifparone == -1))
                sprintf(statelabels(p), "%s: State of population #%d with birth state #%0*d", tmpstr, p, bmag, b);
              else if ((CurveType == PGR) || (CurveType == EQ) || (CurveType == ESS))
                sprintf(statelabels(p), "%s: State of population #%d with birth state #%0*d for parameter value %G", tmpstr, p, bmag, b,
                        parPntr[Bifparone]);
              else if (CurveType == EVODYN)
                sprintf(statelabels(p), "%s: State of population #%d with birth state #%0*d for evolutionary time value %G", tmpstr, p, bmag, b,
                        *timePntr);
              else if (CurveType == PIP)
                sprintf(statelabels(p), "%s: State of population #%d with birth state #%0*d for parameter values %G and %G", tmpstr, p, bmag, b,
                        parPntr[Bifparone], MutantParVal);
              else
                sprintf(statelabels(p), "%s: State of population #%d with birth state #%0*d for parameter values %G and %G", tmpstr, p, bmag, b,
                        parPntr[Bifparone], parPntr[Bifpartwo]);
            }
          else
            {
              sprintf(tmpstr, "Pop%0*d", pmag, p);
              if ((CurveType == PGR) && (Bifparone == -1))
                sprintf(statelabels(p), "%s: State of population #%d", tmpstr, p);
              else if ((CurveType == PGR) || (CurveType == EQ) || (CurveType == ESS))
                sprintf(statelabels(p), "%s: State of population #%d for parameter value %G", tmpstr, p, parPntr[Bifparone]);
              else if (CurveType == EVODYN)
                sprintf(statelabels(p), "%s: State of population #%d for evolutionary time value %G", tmpstr, p, *timePntr);
              else if (CurveType == PIP)
                sprintf(statelabels(p), "%s: State of population #%d for parameter values %G and %G", tmpstr, p, parPntr[Bifparone], MutantParVal);
              else
                sprintf(statelabels(p), "%s: State of population #%d for parameter values %G and %G", tmpstr, p, parPntr[Bifparone],
                        parPntr[Bifpartwo]);
            }

          if (Bifparone != -1)
            cpop->timeval = parPntr[Bifparone];
          else if (CurveType == EVODYN)
            cpop->timeval =*timePntr;
          else if (CurveType == IND)
            cpop->timeval = IStateDim;
          else
            cpop->timeval = 0.0;

          cpop->population = totpop++;
          if (CurveType == IND)
            cpop->columns = CohortDim + birthStateNr[p];
          else
            cpop->columns   = CohortDim;
          cpop->cohorts     = max(Cohorts(p, b), 1);
          cpop->data_offset = hdrdbls + MAX_LBL_LEN;
          if (fullstateoutput == 2)
            cpop->lastpopdim = (p == (CurPopulationNr - 1)) && (b == (birthStateNr[CurPopulationNr - 1] - 1));
          else
            cpop->lastpopdim = (p == (CurPopulationNr - 1));
          fwrite((void *)cpop, 1, hdrdbls*sizeof(double), fp);
          fwrite((void *)statelabels(p), 1, MAX_LBL_LEN*sizeof(double), fp);

          if (Cohorts(p, b))
            {
              fwrite((void *)(&(PopDens(p, b, IStateDim, 0))), Cohorts(p, b), sizeof(double), fp);
              for (i = 0; i < IStateDim; i++) fwrite((void *)(&(PopDens(p, b, i, 0))), Cohorts(p, b), sizeof(double), fp);
              if (CurveType == PGR) fwrite((void *)(&(PopDens(p, b, IStateDim + 1, 0))), Cohorts(p, b), sizeof(double), fp);
              if (CurveType == IND)
                for (i = IStateDim + 1; i < (CohortDim + birthStateNr[p]); i++)
                  fwrite((void *)(&(PopDens(p, b, i, 0))), Cohorts(p, b), sizeof(double), fp);
            }
          else
            {
              for (i = 0; i < CohortDim; i++) fwrite((void *)(&zero), 1, sizeof(double), fp);
            }
          if (fullstateoutput == 1) break;                                          // Store only b=0 in case of condensed output
        }
    }

  (void)fflush(fp);                                                                 // Flush the file buffer
  (void)fclose(fp);
  free(statelblmem);

  return;
}

#endif


/*==================================================================================================================================*/
