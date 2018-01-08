#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP csb2rlist(SEXP, SEXP, SEXP, SEXP);
/*
extern SEXP PSPMdemo(SEXP, SEXP, SEXP, SEXP);
extern SEXP PSPMecodyn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PSPMequi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PSPMevodyn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PSPMind(SEXP, SEXP, SEXP, SEXP);
*/

static const R_CallMethodDef CallEntries[] = {
    {"csb2rlist",  (DL_FUNC) &csb2rlist,  4},
/*
    {"PSPMdemo",   (DL_FUNC) &PSPMdemo,   4},
    {"PSPMecodyn", (DL_FUNC) &PSPMecodyn, 6},
    {"PSPMequi",   (DL_FUNC) &PSPMequi,   7},
    {"PSPMevodyn", (DL_FUNC) &PSPMevodyn, 7},
    {"PSPMind",    (DL_FUNC) &PSPMind,    4},
*/
    {NULL, NULL, 0}
};

void R_init_PSPManalysis(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

