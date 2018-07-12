#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP poisbinom_dpoisbinom(SEXP, SEXP, SEXP);
extern SEXP poisbinom_ppoisbinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP poisbinom_qpoisbinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP poisbinom_rpoisbinom(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"poisbinom_dpoisbinom", (DL_FUNC) &poisbinom_dpoisbinom, 3},
    {"poisbinom_ppoisbinom", (DL_FUNC) &poisbinom_ppoisbinom, 4},
    {"poisbinom_qpoisbinom", (DL_FUNC) &poisbinom_qpoisbinom, 4},
    {"poisbinom_rpoisbinom", (DL_FUNC) &poisbinom_rpoisbinom, 2},
    {NULL, NULL, 0}
};

void R_init_poisbinom(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
