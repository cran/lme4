#ifndef NLME_UTILITIES_H
#define NLME_UTILITIES_H
#include <Rdefines.h>
#include <R_ext/Lapack.h>

typedef struct SEXPREC SEXPREC;

double *nlme_symmetrize(double *a, const int nc);

double *LMEgradient(const double* factor, const double* A,
		    const int nlev, const int nc,
		    const double* pdgradient, const int plen,
		    double* value);

SEXP nlme_replaceSlot(SEXP obj, const SEXPREC* names, const SEXPREC* value);

void nlme_check_Lapack_error(int info, const char *laName);

#endif
