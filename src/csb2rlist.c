/***
   NAME
     csb2rlist
   DESCRIPTION
     Program converts a binary CSB file to text. Either the entire file is
     converted, or just the single state indicated by the time value that is
     passed to the program as command-line

   Last modification: AMdR - Jan 19, 2023
***/
#ifndef CSB2RLIST
#define CSB2RLIST
#endif

#include "stdio.h"
#include "stdint.h"
#include "stdlib.h"
#include "string.h"
#include <sys/param.h>

#include <R.h>
#include <Rdefines.h>

// Type definition for environment and population
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

/*==================================================================================================================================*/

SEXP csb2rlist(SEXP filename, SEXP command, SEXP stateindex, SEXP statetime)

{
  char      tmpstr[MAXPATHLEN], fname[MAXPATHLEN], cmd[128], *cp;
  int       parDim, file_ok = 1, popinfile = 0, popindex;
  double    timeval = -1.0, closest_time = HUGE_VAL;
  FILE      *   fp      = NULL;
  uint32_t  filepos = 0, closest_fpos = 0, magic;
  Envdim    cur_env;
  void      *mem_base     = NULL;
  long      MemAllocated = 0L;
  double    *parVals;

  strcpy(fname, CHAR(STRING_ELT(filename, 0)));
  strcpy(cmd, CHAR(STRING_ELT(command, 0)));
  popindex = INTEGER(stateindex)[0];
  timeval  = REAL(statetime)[0];

  fp = fopen(fname, "rb");                                                          // Open CSB file, error checking occurs in R

  if (fread(&magic, sizeof(uint32_t), 1, fp) != 1) error("\nCSB file does not contain any data!\n\n");

  file_ok = (fread(&parDim, sizeof(int), 1, fp) == 1);
  if (file_ok)
    {
      parVals = (double *)R_alloc(parDim, sizeof(double));
      file_ok = (fread(parVals, sizeof(double), (size_t)parDim, fp) == (size_t)parDim);
    }
  if (!file_ok)
    {
#if (DEBUG == 1)
      for (int i = 0; i < parDim; i++) printf("parameter[%d] = %G\n", i, parVals[i]);
#endif
      fclose(fp);
      error("\nCSB file seems to be invalid or corrupt!\n\n");
    }
  filepos = ftell(fp);

  if (!strcmp(cmd, "list")) Rprintf("\nStates in file %s:\n\n", fname);
  while (fread(&cur_env, sizeof(Envdim), 1, fp) == 1)
    {
      fseek(fp, filepos, SEEK_SET);

      if (cur_env.memory_used > MemAllocated)
        {
          mem_base = (void *)realloc(mem_base, (size_t)cur_env.memory_used);
          if (!mem_base) error("\nCould not allocate memory to read in data!\n\n");
          MemAllocated = cur_env.memory_used;
        }

      if (fread(mem_base, 1, (size_t)cur_env.memory_used, fp) != (size_t)cur_env.memory_used)
        error("\nFailed to read values from file: %s!\n\n", fname);

      popinfile++;
      if (!strcmp(cmd, "read"))
        {
          if (popindex >= 0)
            {
              if (popindex == popinfile)
                {
                  closest_fpos = filepos;
                  break;
                }
            }
          else if (fabs(cur_env.timeval - timeval) < fabs(closest_time - timeval))
            {
              closest_time = cur_env.timeval;
              closest_fpos = filepos;
            }
        }
      else if (!strcmp(cmd, "list"))
        {
          snprintf(tmpstr, sizeof(tmpstr), "%5d: State-%.6E", popinfile, cur_env.timeval);
          Rprintf("%s\n", tmpstr);
        }

      fseek(fp, filepos, SEEK_SET);
      fseek(fp, cur_env.memory_used, SEEK_CUR);
      filepos = ftell(fp);
    }

  if (!strcmp(cmd, "list"))
    {
      SEXP         csbdim;
    
      PROTECT(csbdim = allocVector(INTSXP, 1));
      INTEGER(csbdim)[0] = popinfile;
      UNPROTECT(1);
      
      if (mem_base) free(mem_base);
      fclose(fp);
      Rprintf("\n");

      return (csbdim);
    }
  else if ((!strcmp(cmd, "read")) && (popindex > popinfile))
    error("No such population state in file %s: index out of bounds. Index should be in the range 1-%d\n\n", fname, popinfile);
  else if ((!strcmp(cmd, "read")) && (closest_fpos > (uint32_t)0))
    {
      register int  n, i, j, k;
      double *      cdbl, *base_pnt, *poppnt;
      Envdim *      cenv;
      Popdim *      cpop;
      int           npops = 0, nprotect = 0, nsexpel = 0, num, indx, IStateDim;
      int           changedParsDim;
      int           demo = 0, ecodyn = 0, evodyn = 0, lifedyn = 0, evopars = 0, bifpar1 = 0, bifpar2 = 0;
      unsigned char mag;
      SEXP          tvalvec, parvec, bifparvec, envvec, *popmat;
      SEXP *        popcolnames, *dimnames;
      SEXP          result, ret_names;
      char          fmtstr[MAXPATHLEN];

      fseek(fp, closest_fpos, SEEK_SET);
      if (fread(&cur_env, sizeof(Envdim), 1, fp) != 1) error("\nFailed to read values from file: %s!\n\n", fname);

      fseek(fp, closest_fpos, SEEK_SET);
      if (fread(mem_base, 1, (size_t)cur_env.memory_used, fp) != (size_t)cur_env.memory_used)
        error("\nFailed to read values from file: %s!\n\n", fname);

      cenv     = (Envdim *)mem_base;
      base_pnt = ((double *)mem_base + cenv->data_offset);

      // Count the populations
      cdbl = ((double *)mem_base + cenv->data_offset + cenv->columns);
      cpop = (Popdim *)cdbl;
      while (1)
        {
          npops++;
          if (cpop->lastpopdim) break;
          cdbl = (((double *)cpop) + cpop->data_offset + (cpop->cohorts*cpop->columns));
          cpop = (Popdim *)cdbl;
        }
      cdbl = ((double *)mem_base + cenv->data_offset + cenv->columns);
      cpop = (Popdim *)cdbl;

      demo    = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "PGR") != NULL);
      ecodyn  = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "Environment: Environment variables at time") != NULL);
      evodyn  = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "evolutionary time") != NULL);
      lifedyn = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "Individual life history dynamics") != NULL);
      evopars = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "evolutionary parameters for parameter value") != NULL);
      bifpar1 = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "for parameter value ") != NULL);
      bifpar2 = (strstr((char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)), "for parameter values") != NULL);

      // Allocate the output list with the appropriate dimension, a list to hold the names of its elements
      PROTECT(result = allocVector(VECSXP, npops + 2 + bifpar1 + bifpar2 + evopars + 2*evodyn + ecodyn));
      nprotect++;
      PROTECT(ret_names = allocVector(STRSXP, npops + 2 + bifpar1 + bifpar2 + evopars + 2*evodyn + ecodyn));
      nprotect++;

      // Adjust the values of the bifurcation parameters
      changedParsDim = (int)floor(parVals[parDim - 1] + 0.5);
      for (i = 0; i < changedParsDim; i++)
        {
          indx = (int)floor(parVals[parDim - 1 - changedParsDim + i] + 0.5);
          if ((changedParsDim == 2) && (i == 1) && (indx == ((int)floor(parVals[parDim - 1 - changedParsDim] + 0.5))))
            continue;                                                                // PIP results
          parVals[indx] = base_pnt[cenv->columns - changedParsDim + i];
        }
      parDim -= changedParsDim + 1;

      if (ecodyn)
        {
          PROTECT(tvalvec = allocVector(REALSXP, 1));
          nprotect++;
          REAL(tvalvec)[0] = cenv->timeval;
          SET_VECTOR_ELT(result, nsexpel, tvalvec);
          SET_STRING_ELT(ret_names, nsexpel++, mkChar("Time"));
        }
      else if (evodyn)
        {
          PROTECT(tvalvec = allocVector(REALSXP, 1));
          nprotect++;
          REAL(tvalvec)[0] = cenv->timeval;
          SET_VECTOR_ELT(result, nsexpel, tvalvec);
          SET_STRING_ELT(ret_names, nsexpel++, mkChar("EvoTime"));
        }
      else if (bifpar1)
        {
          PROTECT(bifparvec = allocVector(REALSXP, 1));
          nprotect++;
          REAL(bifparvec)[0] = base_pnt[cenv->columns - changedParsDim];
          SET_VECTOR_ELT(result, nsexpel, bifparvec);
          SET_STRING_ELT(ret_names, nsexpel++, mkChar("BifPars"));
        }
      else if (bifpar2)
        {
          PROTECT(bifparvec = allocVector(REALSXP, changedParsDim));
          nprotect++;
          for (i = 0; i < changedParsDim; i++) REAL(bifparvec)[i] = base_pnt[cenv->columns - changedParsDim + i];
          SET_VECTOR_ELT(result, nsexpel, bifparvec);
          SET_STRING_ELT(ret_names, nsexpel++, mkChar("BifPars"));
        }
      if (evodyn || evopars)
        {
          PROTECT(bifparvec = allocVector(REALSXP, changedParsDim - evopars));
          nprotect++;
          for (i = evopars; i < changedParsDim; i++)                                // In case of EVODYN start at i=0, otherwise skip the first (bifurcation) parameter
            REAL(bifparvec)[i - evopars] = base_pnt[cenv->columns - changedParsDim + i];
          SET_VECTOR_ELT(result, nsexpel, bifparvec);
          SET_STRING_ELT(ret_names, nsexpel++, mkChar("EvoPars"));
        }

      PROTECT(parvec = allocVector(REALSXP, parDim));
      nprotect++;
      for (i = 0; i < parDim; i++) REAL(parvec)[i] = parVals[i];
      SET_VECTOR_ELT(result, nsexpel, parvec);
      SET_STRING_ELT(ret_names, nsexpel++, mkChar("Parameters"));

      PROTECT(envvec = allocVector(REALSXP, cenv->columns - changedParsDim));
      nprotect++;
      for (i = 0; i < cenv->columns - changedParsDim; i++) REAL(envvec)[i] = base_pnt[i];
      SET_VECTOR_ELT(result, nsexpel, envvec);

      if (lifedyn)
        SET_STRING_ELT(ret_names, nsexpel++, mkChar("Environment"));
      else
        {
          strcpy(tmpstr, (char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)));
          cp          = strstr(tmpstr, ": ");
          if (cp) *cp = '\0';
          SET_STRING_ELT(ret_names, nsexpel++, mkChar(tmpstr));
        }

      // Allocate arrays (persisting as long as this C code is running) to hold the SEXP object
      // representing the population matrices, the names of their column headers and lists of both column and row headers
      popmat      = (SEXP *)R_alloc(npops, sizeof(SEXP));
      popcolnames = (SEXP *)R_alloc(npops, sizeof(SEXP));
      dimnames    = (SEXP *)R_alloc(npops, sizeof(SEXP));

      for (n = 0; n < npops; n++)
        {
          strcpy(tmpstr, (char *)(((double *)cpop) + (sizeof(Popdim)/sizeof(double) + 1)));
          cp          = strstr(tmpstr, ": ");
          if (cp) *cp = '\0';
          SET_STRING_ELT(ret_names, nsexpel, mkChar(tmpstr));

          base_pnt  = ((double *)cpop + cpop->data_offset);
          popmat[n] = PROTECT(allocMatrix(REALSXP, cpop->cohorts, cpop->columns));
          nprotect++;
          poppnt = REAL(popmat[n]);
          for (j = 0; j < cpop->cohorts; j++)
            for (k = 0; k < cpop->columns; k++) poppnt[k*cpop->cohorts + j] =*(base_pnt + k*cpop->cohorts + j);

          // Define the column labels of the different populations here
          mag = 1;
          num = cpop->columns;
          while (num > 0)
            {
              mag++;
              num = num/10;
            }
          PROTECT(popcolnames[n] = allocVector(STRSXP, cpop->columns));
          nprotect++;
          if (strstr(tmpstr, "_StableBirthDist") != NULL)
            {
              for (k = 0; k < cpop->columns; k++)
                {
                  snprintf(fmtstr, sizeof(fmtstr), "Bstate%%0%dd", mag);
                  snprintf(tmpstr, sizeof(tmpstr), fmtstr, k);
                  SET_STRING_ELT(popcolnames[n], k, mkChar(tmpstr));
                }
            }
          else if (strstr(tmpstr, "_BirthStates") != NULL)
            {
              for (k = 0; k < cpop->columns; k++)
                {
                  snprintf(fmtstr, sizeof(fmtstr), "Istate%%0%dd", mag);
                  snprintf(tmpstr, sizeof(tmpstr), fmtstr, k);
                  SET_STRING_ELT(popcolnames[n], k, mkChar(tmpstr));
                }
            }
          else
            {
              if (lifedyn) IStateDim = (int)(cpop->timeval + 0.01);
              for (k = 0; k < (cpop->columns - demo); k++)
                {
                  if (demo && (k == 0))
                    SET_STRING_ELT(popcolnames[n], k, mkChar("StableDist"));
                  else if (lifedyn)
                    {
                      if (k == 0)
                        SET_STRING_ELT(popcolnames[n], k, mkChar("Survival"));
                      else if (k == cpop->columns - 1)
                        SET_STRING_ELT(popcolnames[n], k, mkChar("R0"));
                      else
                        {
                          if ((k - 1) < IStateDim)
                            {
                              snprintf(fmtstr, sizeof(fmtstr), "Istate%%0%dd", mag);
                              snprintf(tmpstr, sizeof(tmpstr), fmtstr, k - 1);
                            }
                          else
                            {
                              snprintf(fmtstr, sizeof(fmtstr), "Impact%%0%dd", mag);
                              snprintf(tmpstr, sizeof(tmpstr), fmtstr, k - 1  - IStateDim);
                            }
                          SET_STRING_ELT(popcolnames[n], k, mkChar(tmpstr));
                        }
                    }
                  else if (k == 0)
                    SET_STRING_ELT(popcolnames[n], k, mkChar("Density"));
                  else
                    {
                      snprintf(fmtstr, sizeof(fmtstr), "Istate%%0%dd", mag);
                      snprintf(tmpstr, sizeof(tmpstr), fmtstr, k - 1);
                      SET_STRING_ELT(popcolnames[n], k, mkChar(tmpstr));
                    }
                }
              if (demo) SET_STRING_ELT(popcolnames[n], k, mkChar("ReproVal"));
            }

          // Now assign them
          PROTECT(dimnames[n] = allocVector(VECSXP, 2));
          nprotect++;
          SET_VECTOR_ELT(dimnames[n], 1, popcolnames[n]);
          setAttrib(popmat[n], R_DimNamesSymbol, dimnames[n]);

          // Add the matrix to the output list
          SET_VECTOR_ELT(result, nsexpel++, popmat[n]);
          if (cpop->lastpopdim) break;

          cdbl = (((double *)cpop) + cpop->data_offset + (cpop->cohorts*cpop->columns));
          cpop = (Popdim *)cdbl;
        }

#if (DEBUG == 1)
      double output;

      cdbl = ((double *)mem_base + cenv->data_offset + cenv->columns);
      cpop = (Popdim *)cdbl;

      // The description
      (void)printf("# %s\n", (char *)(((double *)cenv) + (sizeof(Envdim)/sizeof(double) + 1)));
      for (i = 0; i < cenv->columns; i++) /* Write environment state  */
        {
          if (i == 0)
            output = base_pnt[i];
          else
            {
              (void)printf("\t");
              output = base_pnt[i];
            }
          if (((fabs(output) <= 1.0E4) && (fabs(output) >= 1.0E-4)) || (output == 0))
            (void)printf("%.10f", output);
          else
            (void)printf("%.6E", output);
        }
      (void)printf("\n\n");

      // Write population state
      while (1)
        {
          // Description
          (void)printf("# %s\n", (char *)(((double *)cpop) + (sizeof(Popdim)/sizeof(double) + 1)));
          base_pnt = ((double *)cpop + cpop->data_offset);
          for (j = 0; j < cpop->cohorts; j++, (void)printf("\n"))
            for (k = 0; k < cpop->columns; k++)
              {
                if (k > 0) (void)printf("\t");
                output =*(base_pnt + k*cpop->cohorts + j);
                if (((fabs(output) <= 1.0E4) && (fabs(output) >= 1.0E-4)) || (output == 0))
                  (void)printf("%.10f", output);
                else
                  (void)printf("%.6E", output);
              }
          (void)printf("\n");
          if (cpop->lastpopdim) break;

          cdbl = (((double *)cpop) + cpop->data_offset + (cpop->cohorts*cpop->columns));
          cpop = (Popdim *)cdbl;
        }
#endif
      if (mem_base) free(mem_base);
      setAttrib(result, R_NamesSymbol, ret_names);
      UNPROTECT(nprotect);
      fclose(fp);
      return (result);
    }

  if (mem_base) free(mem_base);
  fclose(fp);
  return R_NilValue;
}


/*==================================================================================================================================*/
