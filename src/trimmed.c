#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "utils.h"

SEXP R_trimmed_stats(SEXP X, SEXP Y, SEXP trim){
  SEXP dim1;
  SEXP estimates1, estimates2, estimates3, output;
  double *Xptr, *Mptr1, *Mptr2, *Mptr3, *Tptr;
  int *Yptr;
  int rows, cols;

  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];
  
  Xptr = NUMERIC_POINTER(AS_NUMERIC(X));
  Yptr = INTEGER_POINTER(AS_INTEGER(Y));
  Tptr = NUMERIC_POINTER(AS_NUMERIC(trim));

  PROTECT(estimates1 = allocMatrix(REALSXP, rows, 3));
  PROTECT(estimates2 = allocMatrix(REALSXP, rows, 3));
  PROTECT(estimates3 = allocMatrix(REALSXP, rows, 3));
  
  Mptr1 = NUMERIC_POINTER(estimates1);
  Mptr2 = NUMERIC_POINTER(estimates2);
  Mptr3 = NUMERIC_POINTER(estimates3);
  
  trimmed_stats(Xptr, Mptr1, Mptr2, Mptr3, Yptr, rows, cols, Tptr);

  PROTECT(output = allocVector(VECSXP,3));
  SET_VECTOR_ELT(output, 0, estimates1);
  SET_VECTOR_ELT(output, 1, estimates2);
  SET_VECTOR_ELT(output, 2, estimates3);

  UNPROTECT(5);
  
  return output;

}
