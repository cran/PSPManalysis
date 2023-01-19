/*
  NAME
     ebtutils.c

     This file contains all routines that are part of the original Escalator Boxcar
     Train software package.

  Last modification: AMdR - Jan 19, 2023
*/
#define EBTUTILS_C
#define EBTLIB

#include "escbox.h"
#include "ebtcohrt.h"
#include "ebtdopri5.h"
#include "ebtmain.h"
#include "ebtutils.h"


/*==================================================================================================================================*/
/*
 * The error messages that occur in the routines in the present file.
 */

#define ECSB "Error writing to CSB file. Further state output will be disabled!"
#define ESF  "Unable to open ESF file for storing end state of populations!"
#define MAFC "Memory allocation failure for cohort variables!"
#define MAFI "Memory allocation failure for cohort constants!"


/*==================================================================================================================================*/
/*
 * Start of function implementations.
 */
/*==================================================================================================================================*/

int imin(int a, int b) { return (a < b) ? a : b; }


int imax(int a, int b) { return (a > b) ? a : b; }


double min(double a, double b) { return (a < b) ? a : b; }


double max(double a, double b) { return (a > b) ? a : b; }


int	  ismissing(double a) { return (fabs(a) > 0.95*MISSING_VALUE); }

int isequal2zero(double a){ return (fabs(a) < identical_zero); }

int isequal(double a, double b)
{
  double diff = fabs(a - b);
  return ((diff < identical_zero) || (diff < 0.5*identical_zero*(fabs(a) + fabs(b))));
}


/*==================================================================================================================================*/

void ShutDown(int exitcode)

/* 
   * ShutDown - Routine closes the output file and saves the entire end
   *		state of the population and the environment in the end
   *		state file. 
   */

{
  static int                      first = 1;

  if (!first) return;

  if (outfile) (void)fclose(outfile);                                               // Close result file
  outfile = NULL;
  if (csbfile) (void)fclose(csbfile);                                               // Close binary state file
  csbfile = NULL;
  if (dbgfile) (void)fclose(dbgfile);                                               // Close debugging file
  dbgfile = NULL; 

  if (BifurcationRun)
    {
      if (averages) (void)fclose(averages);                                         // Close bifurcation output files
      averages  = NULL;
      if (gaverages) (void)fclose(gaverages);
      gaverages = NULL;
      if (variances) (void)fclose(variances);
      variances = NULL;
      if (extrema) (void)fclose(extrema);
      extrema   = NULL;
    }

  (void)fprintf(stderr, "\n\nRUN %-s COMPLETED at T = %.2f:\n", runname, env[0]);
  (void)fprintf(stderr, "** %-70s **\n\n", "Program terminated. Normal closure of output files succeeded.");

  ResetDopri5();
  ResetCohorts();
  ResetBifStats();

  (void)fflush(stdout);
  (void)fflush(stderr);

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
  if (exitcode >= 0)
    {
#if (DEBUG > 0)
      if (exitcode)
        abort();
      else
#endif // DEBUG
        exit(exitcode);
    }
#endif

  first = 0;

  return;
}


/*==================================================================================================================================*/

void ErrorAbort(CONST char *mes)

/*
   * ErrorAbort - Routine issues an error message and exits immediately
   *		  without trying to save the state of the system.
   */

{
  char rn[MAXFILENAMELEN];

  (void)strcpy(rn, runname);
  rn[strlen(rn) - 1] = '\0';                                                        // Remove trainling dot

  (void)fprintf(stderr, "\nRUN %-s: ERROR at T = %.2f:\n", rn, env[0]);
  (void)fprintf(stderr, "** %-70s **\n", mes);
  (void)fprintf(stderr, "** %-70s **\n\n", "Normal closure of output files failed.");
  exit(1);

  return;
}


/*==================================================================================================================================*/

void ErrorExit(CONST int exitcode, CONST char *mes)

/*
   * ErrorExit - Routine issues an error message and tries to call the
   *		 ShutDown() routine to save the state of the system. The
   *		 program subsequently terminates.
   */

{
  char rn[MAXFILENAMELEN];

  (void)strcpy(rn, runname);
  rn[strlen(rn) - 1] = '\0';                                                        // Remove trainling dot

  (void)fprintf(stderr, "\nRUN %-s: ERROR at T = %.2f:\n", rn, env[0]);
  (void)fprintf(stderr, "** %-70s **\n", mes);
  (void)fprintf(stderr, "** %-70s **\n\n", "Normal closure of output files attempted.");
  ShutDown(exitcode);

  return;
}


/*==================================================================================================================================*/

void Warning(CONST char *mes)

/*
   * Warning - Routine issues a warning message pertaining to the current
   *	       state of the program. The program continues normally.
   */

{
  char rn[MAXFILENAMELEN];

  (void)strcpy(rn, runname);
  rn[strlen(rn) - 1] = '\0';                                                        // Remove trainling dot

  (void)fprintf(stderr, "\nRUN %-s: WARNING at T = %.2f:\n", rn, env[0]);
  (void)fprintf(stderr, "** %-70s **\n", mes);
  (void)fprintf(stderr, "** %-70s **\n\n", "Program continues normally.");

  return;
}


/*==================================================================================================================================*/

void *Myalloc(void *pnt, size_t count, size_t eltsize)

/*
   * Myalloc - Replacement routine for the 'calloc()' and 'realloc()'
   *	       functions.  Allocates or reallocates memory and sets it to 0.
   *	       Implemented to watch memory use.
   */

{
  size_t size;
  void * value;

  size = count*eltsize;
  if (pnt)                                                                          // Reallocation call
    value = (void *)realloc((DEF_TYPE *)pnt, (SIZE_TYPE)size);
  else                                                                              // Allocation call
    {
      value = (void *)malloc((SIZE_TYPE)size);
      if (value != NULL) (void)memset((DEF_TYPE *)value, 0, size);
    }

  return value;
}


/*==================================================================================================================================*/

void PrettyPrint(FILE *fp, double value)

/*
 * PrettyPrint - Formatted print of output to the file pointed to by fp.
 */

{
#if !defined(_MSC_VER)
  int fpclass = fpclassify(value);

  if (fpclass == FP_NORMAL)
#endif
    {
      if (((fabs(value) <= 1.0E4) && (fabs(value) >= 1.0E-4)) || (value == 0))
        (void)fprintf(fp, "%.10f", value);
      else
        (void)fprintf(fp, "%.6E", value);
    }
#if !defined(_MSC_VER)
  else
    (void)fprintf(fp, "%E", value);
#endif

  return;
}
/*==================================================================================================================================*/

void WriteStateToFile(FILE *fp, double *data)

/*
 * WriteStateToFile - Routine writes the entire state of the populations
 *		      and the environment to the file pointed to by fp.
 */

{
  register int i;

  // Write environment state
  for (i = 0; i < ENVIRON_DIM_EBT; i++)
    {
      if (i) (void)fprintf(fp, "\t");
      if (data)
        PrettyPrint(fp, data[i]);
      else
        PrettyPrint(fp, env[i]);
    }
  (void)fprintf(fp, "\n\n");

  register int j, k;
  int          len  = ENVIRON_DIM_EBT;
  population   lpop = NULL;

  for (i = 0; i < POPULATION_NR; i++)
    {
      // Write boundary cohorts
      if (data) lpop = (population)(data + len + CohortNo[i]*COHORT_SIZE);
      for (j = BpointNo[i] - 1; (j >= 0); j--)
        {
          for (k = 0; k < COHORT_SIZE; k++)
            {
              if (k) (void)fprintf(fp, "\t");
              if (data)
                PrettyPrint(fp, lpop[j][k]);
              else
                PrettyPrint(fp, ofs[i][j][k]);
            }
          for (k = 0; k < I_CONST_DIM; k++)
            {
              (void)fprintf(fp, "\t");
              PrettyPrint(fp, ofsIDcard[i][j][k]);
            }
          (void)fprintf(fp, "\n");
        }

      // Write internal cohorts
      if (data) lpop = (population)(data + len);
      for (j = CohortNo[i] - 1; (j >= 0); j--)
        {
          for (k = 0; k < COHORT_SIZE; k++)
            {
              if (k) (void)fprintf(fp, "\t");
              if (data)
                PrettyPrint(fp, lpop[j][k]);
              else
                PrettyPrint(fp, pop[i][j][k]);
            }
          for (k = 0; k < I_CONST_DIM; k++)
            {
              (void)fprintf(fp, "\t");
              PrettyPrint(fp, popIDcard[i][j][k]);
            }
          (void)fprintf(fp, "\n");
        }

      if (!CohortNo[i])
        {
          for (j = 0; j < (COHORT_SIZE + I_CONST_DIM); j++)
            {
              if (j) (void)fprintf(fp, "\t");
              (void)fprintf(fp, "%f", 0.0);
            }
          (void)fprintf(fp, "\n");
        }
      else
        len += (CohortNo[i] + BpointNo[i])*COHORT_SIZE;

      (void)fprintf(fp, "\n");
    }

  return;
}


/*==================================================================================================================================*/

void FileOut(void)

/*
   * FileOut - Routine sends output to file. All output statistics are
   *	       defined by the user, except for the first one which is the
   *	       current time value.
   */

{
  register int i;

  for (i = 0; i < OUTPUT_VAR_NR; i++) output[i] = 0.0;
  for (i = 0; i < POPULATION_NR; i++) cohort_no[i] = CohortNo[i];

  DefineOutput(env, pop, output); /* User output values       */

  // Add bifurcation parameter as last column
  if (BifurcationRun) output[OUTPUT_VAR_NR] = parameter[BifParIndex];

  if (env[0] < 1.0E6)
    (void)fprintf(outfile, "%10.1f", env[0]);
  else
    (void)fprintf(outfile, "%10.4E", env[0]);
  
  for (i = 0; i < output_var_nr; i++)
    {
      (void)fprintf(outfile, "\t");
      PrettyPrint(outfile, output[i]);
    }
  (void)fprintf(outfile, "\n");
  (void)fflush(outfile);

  for (i = output_var_nr; i > 0; i--) output[i] = output[i - 1];
  output[0]                                           = env[0];

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

void FileState(void)

/*
 * FileState - Routine writes the entire state of the environment and all populations in binary format to the '.csb' file.
 */

{
  int         i, j, p, nalloc = 0, status;
  int         pmag, num;
  static int  skip = 0;
  mxArray     *mystruct;
  mxArray     **mydata;
  const char  **fieldnames;
  char        matfname[MAXFILENAMELEN];
  char        tmpstr[FIELDNAMELN], tmpstr2[FIELDNAMELN];
  double      *dblptr;
#if defined(OCTAVE_MEX_FILE)
  mxArray     *prhs[4];
#endif

  if (skip) return;

  mydata     = malloc((10 + POPULATION_NR*(COHORT_SIZE + I_CONST_DIM))*sizeof(mxArray *));
  fieldnames = malloc((10 + POPULATION_NR*(COHORT_SIZE + I_CONST_DIM))*sizeof(char *));

  if (!mydata || !fieldnames)
    {
      Warning("Memory allocation error in FileState(). Further state output suppressed!");
      if (mydata) free(mydata);
      if (fieldnames) free(fieldnames);
      skip = 1;
      return;
    }

  nalloc = 0;
  mydata[nalloc]     = mxCreateDoubleMatrix(1, 1, mxREAL);
  fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
  if (!mydata[nalloc] || !fieldnames[nalloc])
    {
      Warning("Memory allocation error in FileState(). Further state output suppressed!");
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
  memcpy((void *)fieldnames[nalloc], "Time", sizeof("Time"));
  dblptr = mxGetPr(mydata[nalloc]);
  memcpy(dblptr, env, sizeof(double));
  nalloc++;

  mydata[nalloc]     = mxCreateDoubleMatrix(1, PARAMETER_NR, mxREAL);
  fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
  if (!mydata[nalloc] || !fieldnames[nalloc])
    {
      Warning("Memory allocation error in FileState(). Further state output suppressed!");
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
  memcpy(dblptr, parameter, PARAMETER_NR*sizeof(double));
  nalloc++;

  mydata[nalloc]     = mxCreateDoubleMatrix(1, ENVIRON_DIM_EBT, mxREAL);
  fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
  if (!mydata[nalloc] || !fieldnames[nalloc])
    {
      Warning("Memory allocation error in FileState(). Further state output suppressed!");
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
  memcpy((void *)fieldnames[nalloc], "Environment", sizeof("Environment"));
  dblptr = (double *)mxGetPr(mydata[nalloc]);
  memcpy(dblptr, env, ENVIRON_DIM_EBT*sizeof(double));
  nalloc++;

  pmag = 1; num  = POPULATION_NR; while (num > 0) { pmag++; num = num/10; }

  for (p = 0; p < POPULATION_NR; p++)
    {
      if (!CohortNo[p]) continue;

      mydata[nalloc]   = mxCreateDoubleMatrix(CohortNo[p], (COHORT_SIZE + I_CONST_DIM), mxREAL);
      fieldnames[nalloc] = (char *)mxCalloc(FIELDNAMELN, sizeof(char));
      if (!mydata[nalloc] || !fieldnames[nalloc])
        {
          Warning("Memory allocation error in FileState(). Further state output suppressed!");
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
      dblptr = (double *)mxGetPr(mydata[nalloc]);
      for (j = CohortNo[p] - 1; j >= 0; j--, dblptr++) *dblptr = popIDcard[p][j][number]*exp(-(pop[p][j][number] - 1.0));
      for (i = 1; i < COHORT_SIZE; i++)
        for (j = CohortNo[p] - 1; j >= 0; j--, dblptr++) *dblptr = pop[p][j][i];

      for (i = 0; i < I_CONST_DIM; i++)
        for (j = CohortNo[p] - 1; j >= 0; j--, dblptr++) *dblptr = popIDcard[p][j][i];
        
      snprintf(tmpstr, sizeof(tmpstr), "Pop%0*d", pmag, p);
      memcpy((void *)fieldnames[nalloc], tmpstr, strlen(tmpstr)*sizeof(char));
      nalloc++;
    }

  mystruct = mxCreateStructMatrix(1, 1, nalloc, (const char **)fieldnames);

  if (mystruct)
    {
      for (j = 0; j < nalloc; j++) mxSetFieldByNumber(mystruct, 0, j, mydata[j]);

      snprintf(tmpstr, sizeof(tmpstr), "State_%.6E", env[0]);
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
      snprintf(matfname, sizeof(matfname), "%s.mat", runname);

#if defined(MATLAB_MEX_FILE)
      if (pmat == NULL) pmat = matOpen(matfname, "wz");
      if (pmat == NULL)
        {
          Warning("Memory allocation error in FileState(). Further state output suppressed!");
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
          Warning("Failed to write to file in FileState(). Further state output suppressed!");
          skip = 1;
        }
#else
      status = mexPutVariable("caller", tmpstr2, mystruct);                         // Put the structure in the workspace
      if (status != 0)
        {
          Warning("Failed to put population state into global workspace in FileState(). Further state output suppressed!");
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
      Warning("Failed to create output structure in FileState(). Further state output suppressed!");
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

// Definition of Envdim and Popdim structure types

#if defined(_MSC_VER) && (_MSC_VER <= 1600)
typedef __int32 uint32_t;
#else
#include <stdint.h>
#endif

typedef struct envdim {
                        double   timeval;
                        int      columns;
                        int      data_offset;
                        uint32_t memory_used;
                      } Envdim;

typedef struct popdim {
                        double timeval;
                        int    population;
                        int    cohorts;
                        int    columns;
                        int    data_offset;
                        int    lastpopdim;
                      } Popdim;


// Magic key of the type of CSB file written
const uint32_t                    CSB_MAGIC_KEY = 20030509;

#define MAX_LBL_LEN               32                                                // Measured in sizeof(double)

void FileState(void)

/*
 * FileState - Routine writes the entire state of the environment and all populations in binary format to the '.csb' file.
 */

{
  uint32_t      tmpint32;
  int           writeOK = 1;
  register int  i, j, k;
  int           tmpint, pmag, num;
  size_t        hdrdbls;
  Envdim        cenv[2];
  Popdim        cpop[2];
  double        zero = 0.0, cohortDensity, tmpval;
  char          statelabels[POPULATION_NR][MAX_LBL_LEN*sizeof(double)];
  char          envlabel[MAX_LBL_LEN*sizeof(double)];

  if (csbfile && csbnew)
    {                                                                               // New CSB file: Write magic key and parameters
      tmpint32             = CSB_MAGIC_KEY;
      writeOK              = (fwrite((void *)(&tmpint32), 1, sizeof(uint32_t), csbfile) == sizeof(uint32_t));

      // Add an element (value: 0) to the parameter array to indicate that no bifurcation parameters are changed 
      tmpint               = PARAMETER_NR + 1;
      if (writeOK) writeOK = (fwrite((void *)(&tmpint), 1, sizeof(int), csbfile) == sizeof(int));
      if (writeOK) writeOK = (fwrite((void *)(&parameter), PARAMETER_NR, sizeof(double), csbfile) == sizeof(double));

      tmpval = 0;
      fwrite((void *)(&tmpval), 1, sizeof(double), csbfile);
      csbnew               = 0;
    }

  if (!writeOK)
    {
      Warning(ECSB);
      (void)fclose(csbfile);
      csbfile = NULL;
    }

  hdrdbls = (sizeof(Envdim)/sizeof(double)) + 1;
  (void)memset((DEF_TYPE *)cenv, 0, hdrdbls*sizeof(double));

  cenv->timeval       = env[0];
  cenv->columns       = ENVIRON_DIM;
  snprintf(envlabel, sizeof(envlabel), "Environment: Environment variables at time %G", env[0]);
  cenv->data_offset   = hdrdbls + MAX_LBL_LEN;
  cenv->data_offset   = hdrdbls + MAX_LBL_LEN;
  
  // Compute the memory to be written to file
  cenv->memory_used   = (hdrdbls + MAX_LBL_LEN)*sizeof(double);
  cenv->memory_used  += (cenv->columns)*sizeof(double);

  for (i = 0; i < POPULATION_NR; i++)
    {
      hdrdbls = (sizeof(Popdim)/sizeof(double)) + 1;
      cenv->memory_used += (hdrdbls + MAX_LBL_LEN)*sizeof(double);
      cenv->memory_used += (imax(CohortNo[i], 1)*(COHORT_SIZE + I_CONST_DIM)*sizeof(double));
    }
  
  hdrdbls = (sizeof(Envdim)/sizeof(double)) + 1;
  writeOK = (fwrite((void *)cenv, 1, hdrdbls*sizeof(double), csbfile) == (hdrdbls*sizeof(double)));
  if (writeOK) writeOK = (fwrite((void *)envlabel, 1, MAX_LBL_LEN*sizeof(double), csbfile) == MAX_LBL_LEN*sizeof(double));
  if (writeOK) writeOK = (fwrite((void *)(env + 1), 1, ENVIRON_DIM*sizeof(double), csbfile) == (ENVIRON_DIM*sizeof(double)));

  pmag = 1; num  = POPULATION_NR; while (num > 0) { pmag++; num = num/10;}

  for (i = 0; (i < POPULATION_NR) && (writeOK); i++)
    {
      hdrdbls = (sizeof(Popdim)/sizeof(double)) + 1;
      (void)memset((void *)cpop, 0, hdrdbls*sizeof(double));
      snprintf(statelabels[i], MAX_LBL_LEN * sizeof(double), "Pop%0*d: State of population #%d at time %G", pmag, i, i, env[0]);
      cpop->timeval     = env[0];
      cpop->population  = i;
      cpop->columns     = (COHORT_SIZE + I_CONST_DIM);
      cpop->cohorts     = imax(CohortNo[i], 1);
      cpop->data_offset = hdrdbls + MAX_LBL_LEN;
      cpop->lastpopdim  = (i == (POPULATION_NR - 1));

      writeOK = (fwrite((void *)cpop, 1, hdrdbls*sizeof(double), csbfile) == (hdrdbls*sizeof(double)));
      if (!writeOK) break;

      writeOK = (fwrite((void *)statelabels[i], 1, MAX_LBL_LEN*sizeof(double), csbfile) == (MAX_LBL_LEN*sizeof(double)));
      if (!writeOK) break;

      if (CohortNo[i])
        {
          for (k = 0; (k < COHORT_SIZE) && writeOK; k++)
            for (j = CohortNo[i] - 1; (j >= 0) && writeOK; j--)
              {
                if (!k)
                  {
                    cohortDensity = popIDcard[i][j][number]*exp(-(pop[i][j][number] - 1.0));
                    writeOK = (fwrite((void *)(&cohortDensity), sizeof(double), 1, csbfile) == 1);
                  }
                else
                  writeOK = (fwrite((void *)(pop[i][j] + k), sizeof(double), 1, csbfile) == 1);
              }
          for (k = 0; (k < I_CONST_DIM) && writeOK; k++)
            for (j = CohortNo[i] - 1; (j >= 0) && writeOK; j--)
              writeOK = (fwrite((void *)(popIDcard[i][j] + k), sizeof(double), 1, csbfile) == 1);
        }
      else
        for (k = 0; (k < (COHORT_SIZE + I_CONST_DIM)) && writeOK; k++) writeOK = (fwrite((void *)&zero, sizeof(double), 1, csbfile) == 1);
    }

  (void)fflush(csbfile);

  if (!writeOK)
    {
      Warning(ECSB);
      (void)fclose(csbfile);
      csbfile = NULL;
    }

  return;
}
#endif


/*==================================================================================================================================*/

static double next_bif_output;
static int    DoBifOutput = 0;

#include "ebtbifstats.c"

void SetBifOutputTimes(double *env)

{
  static int first = 1;

  // At start no (state) output
  if (first)
    {
      BifPeriod = max_time;

      max_time = floor(((BifParLastVal - parameter[BifParIndex])/BifParStep) + 1.0 + BIFTINY);

      max_time *= BifPeriod;

      next_state_output = ((floor((BIFTINY + env[0])/BifPeriod) + 1)*BifPeriod - BifStateOutput);
      next_output       = ((floor((BIFTINY + env[0])/BifPeriod) + 1)*BifPeriod - BifOutput);
      next_bif_output   = (floor((BIFTINY + env[0])/BifPeriod) + 1)*BifPeriod;

      first = 0;

      return;
    }

  if (env[0] >= (next_bif_output - identical_zero))
    {
      // If env[0] equal to integer multiple of BifPeriod both state output
      // and regular output, unless this already last time we were here

      next_state_output = env[0];
      next_output       = env[0];
      next_bif_output += BifPeriod;

      DoBifOutput = 1;
    }
  else
    {
      // Otherwise schedule next state output for the end of the
      // current period with constant bifurcation parameter and
      // regular output for the last BifOutput timesteps of this
      // period
      //
      next_state_output = max(next_state_output, next_bif_output - BifStateOutput);
      next_output       = max(next_output, next_bif_output - BifOutput);
      DoBifOutput       = 0;
    }

  return;
}


/*==================================================================================================================================*/
// Ctrl-C detection

#if defined(MATLAB_MEX_FILE)
int checkInterrupt(void)
{
  int         pressed;
  extern bool utIsInterruptPending(void);
  extern void utSetInterruptPending(bool);

  // check for a Ctrl-C event
  pressed = utIsInterruptPending();
  if (pressed)
    {
      utSetInterruptPending(false);
      mexPrintf("\n\nCtrl-C detected. Stopping computation\n\n");
      CtrlCPressed = true;
    }

  return (pressed || CtrlCPressed);
}

#elif defined(OCTAVE_MEX_FILE)

int checkInterrupt(void)
{
  int pressed = 0;

  // checking of Ctrl-C event not implemented

  return pressed;
}

#elif defined(R_PACKAGE)

static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }

// this will call the above in a top-level context so it won't longjmp-out of your context
int checkInterrupt(void)
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
