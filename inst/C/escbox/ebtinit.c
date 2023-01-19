/*
  NAME
     ebtinit.c

     This file contains the Initialize() routine and all functions that are
     used by it to initialize global variables, to read the input file with
     the initial state values of the environment and the populations and to
     read the control variable file. 
     
   Last modification: AMdR - Jan 19, 2023
*/

#define EBTINIT_C 				                                                          // Identification of file
#define EBTLIB					                                                            // and file grouping

#include "escbox.h"
#include "ebtinit.h"
#include "ebtmain.h"
#include "ebtcohrt.h"
#include "ebtutils.h"


/*==================================================================================================================================*/
/*
 * Defining all constants that are local to this specific file.
 */
#define MAX_INPUT_LINE 2048			                                                    // The maximum length of a line of input




/*==================================================================================================================================*/
/*
 * The error messages that occur in the routines in the present file.
 */

#define EENV "Unexpected end/error while reading environment from ISF file!"
#define EISF "Unexpected end/error while reading populations from ISF file!"
#define ICS  "Incomplete cohort specification(s) encountered in ISF file!"
#define ISF  "Unable to open ISF file! Expecting initialization in UserInit()!"
#define MAFC "Memory allocation failure for cohort variables!"
#define MAFI "Memory allocation failure for cohort constants!"
#define WNEV "Incomplete environment specification encountered in ISF file! Expecting initialization in UserInit()!"


/*==================================================================================================================================*/
/*
 * Start of function implementations.
 */
/*==================================================================================================================================*/

static char	*ReadDouble(double *val, char *cpnt)

  /* 
   * ReadDouble - Routine reads a double value from the string pointed to
   *		  by "cpnt". Invalid characters are skipped. It returns a 
   *		  pointer to the rest of the string or NULL on error.
   */

{
  register char *ch, *end = NULL;
  int            dot_start = 0;

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
      if ((ch != cpnt) && (*(ch - 1) == '-'))
        ch--;                                                                       // Is previous a minus?

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

static int ReadInputEnv(FILE *infile)

/* 
   * ReadInputEnv - Read the initial values of the environment variables 
   *		    from the already opened .isf file.
   */

{
  char *ch, input[MAX_INPUT_LINE];
  int   read_no;
  // Number of variables already read into array
  read_no = 0;
  while ((!feof(infile)) && (!ferror(infile)))
    {
      // Input line
      ch = fgets(input, MAX_INPUT_LINE, infile);
      while (isspace(*ch)) ch++;
      if (*ch == '#') continue;

      // Read double from string, stop on line end
      while (ch)
        if ((ch = ReadDouble(env + read_no, ch)) != NULL)
          if ((++read_no) == ENVIRON_DIM_EBT) return read_no;
    }

  // Stop loop if array contains enough  data
  if (feof(infile) || ferror(infile)) Warning(EENV);

  return read_no;
}


/*==================================================================================================================================*/

static void ReadInputPop(FILE *infile)

/* 
   * ReadInputPop - Routine reads the initial values for the state of all 
   *		    populations from the already opened .isf file.
   */

{
  register int i, j;
  char *       ch, input[MAX_INPUT_LINE];
  int          done, read_no, warnics = 1;
  double       val_tmp[COHORT_SIZE + I_CONST_DIM];
  long         mem_req;

  for (i = 0; i < POPULATION_NR; i++)
    {
      done = 0;                                                                     // Flag indicating end of population data
      while ((!feof(infile)) && (!ferror(infile)) && (!done))
        {
          // Input line
          ch = fgets(input, MAX_INPUT_LINE, infile);
          while (isspace(*ch)) ch++;
          if (*ch == '#') continue;

          // Initialize all data to default: MISSING_VALUE
          for (j = 0; j < (COHORT_SIZE + I_CONST_DIM); j++) val_tmp[j] = MISSING_VALUE;

          // Read complete cohort	from one line. Stop on line end or full cohort
          for (j = 0, read_no = 0; (j < (COHORT_SIZE + I_CONST_DIM)) && (ch); j++)
            if ((ch = ReadDouble(val_tmp + j, ch)) != NULL) read_no++;

          // End of population reached
          if ((read_no == 0) && CohortNo[i]) done = 1;
          else if (read_no > 0)                                                     // Store cohort read
            {
              // Warn of incomplete cohort
              if ((read_no != (COHORT_SIZE + I_CONST_DIM)) && warnics)
                {
                  Warning(ICS);
                  warnics = 0;
                }
              mem_req = (CohortNo[i] + 1)*COHORT_SIZE;
              if (!(mem_req < DataMemAllocated[i]))
                {
                  DataMemAllocated[i] = MemBlocks(mem_req);
                  pop[i]              = (population)Myalloc((void *)pop[i], (size_t)DataMemAllocated[i], sizeof(double));
                  if (!(pop[i])) ErrorAbort(MAFC);
                }
              mem_req = (CohortNo[i] + 1)*I_CONST_DIM;
              if (!(mem_req < IDMemAllocated[i]))
                {
                  IDMemAllocated[i] = MemBlocks(mem_req);
                  popIDcard[i]      = (popID)Myalloc((void *)popIDcard[i], (size_t)IDMemAllocated[i], sizeof(double));
                  if (!(popIDcard[i])) ErrorAbort(MAFI);
                }
              if (CohortNo[i])
                CohortNo[i] = InsCohort(val_tmp, val_tmp + COHORT_SIZE, i);
              else
                {
                  for (j = 0; j < COHORT_SIZE; j++) pop[i][CohortNo[i]][j] = val_tmp[j];
                  for (j = 0; j < I_CONST_DIM; j++) popIDcard[i][CohortNo[i]][j] = val_tmp[j + COHORT_SIZE];
                  CohortNo[i]++;
                }
            }
        }
      if (ferror(infile) || (feof(infile) && !CohortNo[POPULATION_NR - 1]))
        Warning(EISF);                                                              // On read error exit
    }

  for (i = 0; i < POPULATION_NR; i++)
    {
      (void)snprintf(input, sizeof(input), "Initial population %d is empty! %s", i, "Expecting initialization in UserInit()!");
      if (!CohortNo[i]) Warning(input);
      cohort_no[i] = CohortNo[i];
    }

  return;
}


/*==================================================================================================================================*/

static void Usage(char *progname)

{
  fprintf(stderr, "Usage:\t%s %s %s %s %s", progname, "[-info <1|2|3|4>]", "<ISF filename>",
          "<Cohort cycle time> <Output time interval> <State output time interval> <Maximum integration time>",
          "[<Index> <First> <Step> <Last> <Output> <State output>]");
  fprintf(stderr, "\n\n%s\n\n", "Aim:\tSimulating ecological dynamics of a structured population using the Escalator Boxcar Train");
  fprintf(stderr, "Command-line arguments:\n\n");
  fprintf(stderr, "\t<Cohort cycle time>        : Time interval between starts of new boundary cohorts\n");
  fprintf(stderr, "\t<Output time interval>     : Time interval between data output to .out file\n");
  fprintf(stderr, "\t<State output interval>    : Time interval between complete state output to .csb file\n");
  fprintf(stderr, "\t<Maximum integration time> : Maximum time value until which to continue the integration\n");

  fprintf(stderr, "\nWhen the following command-line arguments are defined a bifurcation simulation is carried out:\n\n");
  fprintf(stderr, "\t<Index>          : Index of the bifurcation parameter\n");
  fprintf(stderr, "\t<First>          : Starting value of the bifurcation parameter\n");
  fprintf(stderr, "\t<Step>           : Step size in the bifurcation parameter\n");
  fprintf(stderr, "\t<Last>           : Final value of the bifurcation parameter\n");
  fprintf(stderr, "\t<Output>         : Period of producing data output during each bifurcation interval\n");
  fprintf(stderr, "\t<State output>   : Period of producing state output during each bifurcation interval\n");
  
  fprintf(stderr, "\nOptional command-line arguments:\n\n");
  fprintf(stderr, "\t-info   <1|2|3|4>       : Level of performance information on the DOPRI5 integrator written to .err file\n");
  fprintf(stderr, "\t-report <1, 2, 3, ...>  : Interval between reporting of data output to console ( > 0)\n");
  
  fprintf(stderr, "\n%s, Copyright (C) 2015, Andre M. de Roos, University of Amsterdam\n\n", progname);
  fprintf(stderr, "This program comes with ABSOLUTELY NO WARRANTY; without even the implied warranty of\n");
  fprintf(stderr, "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public\n");
  fprintf(stderr, "License (<http://www.gnu.org/philosophy/why-not-lgpl.html>) for more details\n\n");

  exit(1);                                                                           // Only executed when in command-line mode

  return;
}



/*==================================================================================================================================*/

void	  Initialize(int argc, char **argv)

  /* 
   * Initialize - Routine initializes the global variables, reads the 
   *		  constants from the command-line and the initial state from
   *		  the .isf file and takes care of the output at start up.
   */

{
  int      argnr = 1;
  char     filename[MAXFILENAMELEN];
  FILE    *isf;

  if (argc < 2)
    {
      fprintf(stderr, "\nNot enough command-line arguments\n\n");
      Usage(argv[0]);
    }
  if (!strcmp(argv[argnr], "-info"))
    {
      argnr++;
      if (!argv[argnr])
        {
          fprintf(stderr, "\nOutput level of the DOPRI5 time integrator not specified!\n\n");
          Usage(argv[0]);
        }
      debug_level = atoi(argv[argnr]);
      if ((debug_level < 1) || (debug_level > 4))
        {
          fprintf(stderr, "\nOutput level of the DOPRI5 time integrator (%d) not in the appropriate range (0 < i < 4)!\n\n", debug_level);
          Usage(argv[0]);
        }
      argnr++;
    }
  else if (!strcmp(argv[argnr], "-report"))
    {
      argnr++;
      if (!argv[argnr])
        {
          fprintf(stderr, "\nInterval for reporting of data output to console not specified!\n\n");
          Usage(argv[0]);
        }
      report_level = atoi(argv[argnr]);
      if (report_level < 1)
        {
          fprintf(stderr, "\nInterval for reporting of data output to console (%d) should be at least 1!\n\n", report_level);
          Usage(argv[0]);
        }
      argnr++;
    }
  else if (!strncmp(argv[argnr], "-", 1))
    {
      fprintf(stderr, "\nUnknown command line option: %s\n", argv[argnr]);
      Usage(argv[0]);
    }

  strcpy(filename, argv[argnr++]);                                                  // Store name of the run
  if (strcmp(filename + strlen(filename) - 4, ".isf"))
    strcat(filename, ".isf");                                                       // Add .isf if not there

  if (((argc - argnr) != 4) && ((argc - argnr) != 10))
    {
      fprintf(stderr, "\nWrong number of command-line arguments!\n\n");
      Usage(argv[0]);
    }

  BifurcationRun = ((argc - argnr) == 10);
  cohort_limit = atof(argv[argnr]); argnr++;
  delt_out     = atof(argv[argnr]); argnr++;
  state_out    = atof(argv[argnr]); argnr++;
  max_time     = atof(argv[argnr]); argnr++;
  if (cohort_limit <= 0.0)
    {
      fprintf(stderr, "\nCohort cylce time should be positive!\n\n");
      Usage(argv[0]);
    }
  if (delt_out < cohort_limit)
    {
      fprintf(stderr, "\nInterval for data output to .out file should be larger than or equal to cohort cylce time!\n\n");
      Usage(argv[0]);
    }
  if ((state_out != 0) && (state_out < cohort_limit))
    {
      fprintf(stderr, "\nInterval for complete state output to .csb file should either be 0 or larger than or equal to cohort cylce time!\n\n");
      Usage(argv[0]);
    }
  if (max_time < cohort_limit)
    {
      fprintf(stderr, "\nMaximum integration time should be larger than cohort or equal to cylce time!\n\n");
      Usage(argv[0]);
    }

  if (BifurcationRun)
    {
      BifParIndex             = atoi(argv[argnr]); argnr++;
      if ((BifParIndex < 0) || (BifParIndex > PARAMETER_NR))
        {
          fprintf(stderr, "\nIndex of bifurcation parameter (%d) not in the appropriate range (0 <= i < %d)!\n\n", BifParIndex, PARAMETER_NR);
          Usage(argv[0]);
        }

      parameter[BifParIndex]  = atof(argv[argnr]); argnr++;

      BifParStep              = atof(argv[argnr]); argnr++;
      // Sanitize the bifurcation control variable: Make it absolute
      BifParStep = fabs(BifParStep);
      if (BifParStep < Odesolve_Min_Step)
        {
          fprintf(stderr, "\nStep size in bifurcation parameter (%G) should be larger than %G!\n\n", BifParStep, Odesolve_Min_Step);
          Usage(argv[0]);
        }

      BifParLastVal           = atof(argv[argnr]); argnr++;
      // Set the sign of bifurcation step size
      if (parameter[BifParIndex] > BifParLastVal) BifParStep *= -1;

      BifOutput               = atof(argv[argnr]); argnr++;
      // Sanitize the bifurcation control variable
      BifOutput = max(BifOutput, 0.0);
      BifOutput = min(BifOutput, max_time);

      BifStateOutput          = atof(argv[argnr]); argnr++;
      // Sanitize the bifurcation control variable
      BifStateOutput = max(BifStateOutput, 0.0);
      BifStateOutput = min(BifStateOutput, max_time);
    }

  // Open ISF file with lower case extension
  isf = fopen(filename, "r");
  if (isf)                                                                          // Read initial state of environment
    {
      if (ReadInputEnv(isf) != ENVIRON_DIM_EBT) Warning(WNEV);
      ReadInputPop(isf);                                                            // Read initial state
      (void)fclose(isf);
    }
  else Warning(ISF);                                                                // UserInit()

  return;
}


/*==================================================================================================================================*/
