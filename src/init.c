#include "pedigree.h"
#include <R_ext/Rdynload.h>

extern SEXP _pedigreeTools_inbreeding_meta(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"pedigree_chol", (DL_FUNC) &pedigree_chol, 2},
    {"pedigree_inbreeding", (DL_FUNC) &pedigree_inbreeding, 1},
    {"_pedigreeTools_inbreeding_meta", (DL_FUNC) &_pedigreeTools_inbreeding_meta, 2},
    {NULL, NULL, 0}
};

void R_init_pedigreeTools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
