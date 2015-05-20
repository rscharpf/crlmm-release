#include <R.h>
#include <Rinternals.h>

SEXP gtypeCallerPart1(SEXP *, SEXP *, SEXP *, SEXP *,SEXP *, SEXP *,
		      SEXP *, SEXP *, SEXP *, SEXP *,SEXP *, SEXP *,
		      SEXP *, SEXP *, SEXP *, SEXP *,SEXP *);

SEXP gtypeCallerPart2(SEXP *, SEXP *, SEXP *, SEXP *,
		      SEXP *, SEXP *, SEXP *, SEXP *,
                      SEXP *, SEXP *, SEXP *, SEXP *,
		      SEXP *, SEXP *, SEXP *, SEXP *,
                      SEXP *, SEXP *, SEXP *);

SEXP normalizeBAF(SEXP *, SEXP *);

SEXP krlmmComputeM(SEXP *, SEXP *);

SEXP krlmmComputeS(SEXP *, SEXP *);

SEXP readGenCallOutputCFunc(SEXP *, SEXP *, SEXP *, SEXP *, SEXP *, SEXP *, SEXP *, SEXP *);

SEXP krlmmConfidenceScore(SEXP *, SEXP *);

SEXP krlmmHardyweinberg(SEXP *);

SEXP countFileLines(SEXP *);
