#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "bitfield.h"
#include "utilities.h"

double
nlme_decomposeChunk(const double* srcmat, int nrowSrc, int ncolSrc,
                    double* destmat, int nrowDest, int ncolDest, SEXP storemat,
                    const SEXPREC* Delta, const SEXPREC* srcRowIndx,
                    const SEXPREC* storeRowIndx, int startCol, int ncol,
                    int q, int ncolRest, nlme_bitfield* srcValidRows,
                    nlme_bitfield* destValidRows, double* scratch,
                    double* qraux, int* pivot, double* work, int lwork);

int nlme_setParameters(SEXP* reStruct, const SEXPREC* pars);

double nlme_logLikelihood_internal(SEXP reStruct, int store,
                                   int GLMMLaplace2);

void nlme_estimate_level(SEXP stored, const SEXPREC* lmeLevel, SEXP bbetas);

void nlme_invert_level(SEXP stored, const SEXPREC* lmeLevel, int ncol_levels);

void nlme_factor_level(SEXP stored, SEXP lmeLevel, SEXP bbetas,
                       double sigmainv, int fixedStartCol);

SEXP nlme_solveOnly(SEXP reStruct);

SEXP nlme_reStruct_fitted(const SEXPREC* reStruct, const SEXPREC* level);

SEXP nlme_reStruct_fitted_internal(const SEXPREC* reStruct, SEXP ans,
                                   const SEXPREC* level);

int nlme_scratchRowSize(const SEXPREC* vecList);

void nlme_copyRows(const double* src, int nrowSrc, int ncolSrc,
                   double* dest, int nrowDest, int ncolDest,
                   int srcStartCol, int destStartCol, int ncol,
                   int srcStartRow, int m,
                   int destStartRow, int mdest,
                   int upper,
                   nlme_bitfield *srcNonZero,
                   nlme_bitfield *destNonZero);

SEXP nlme_predecompose(SEXP reStruct);

SEXP nlme_commonDecompose(SEXP reStruct, const SEXPREC* pars);

SEXP nlme_getFixDF(const SEXPREC* reStruct);
