#include <R.h>
#include <R_ext/Rdynload.h>
#include "crlmm.h"
#include <Rinternals.h>

static const R_CallMethodDef CallEntries[] = {
    {"gtypeCallerPart1", (DL_FUNC)&gtypeCallerPart1, 17},
    {"gtypeCallerPart2", (DL_FUNC)&gtypeCallerPart2, 19},
    {"normalizeBAF", (DL_FUNC)&normalizeBAF, 2},
    {"krlmmComputeM", (DL_FUNC)&krlmmComputeM, 2},
    {"krlmmComputeS", (DL_FUNC)&krlmmComputeS, 2},
    {"readGenCallOutputCFunc", (DL_FUNC)&readGenCallOutputCFunc, 8},
    {"krlmmConfidenceScore", (DL_FUNC)&krlmmConfidenceScore, 2},
    {"krlmmHardyweinberg", (DL_FUNC)&krlmmHardyweinberg, 1},
    {"countFileLines", (DL_FUNC)&countFileLines, 1},
    {NULL, NULL, 0}
};

void R_init_crlmm(DllInfo *dll){
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}

SEXP subColSummarizeMedianPP(SEXP RMatrix, SEXP R_rowIndexList){
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL)
    fun =  (SEXP(*)(SEXP, SEXP))R_GetCCallable("preprocessCore","R_subColSummarize_median");
  return fun(RMatrix, R_rowIndexList);
}
