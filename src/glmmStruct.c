/**
 * @file   glmmStruct.c
 * @author Saikat DebRoy <saikat@stat.wisc.edu>
 * @author Douglas Bates <bates@stat.wisc.edu>
 * @date   $Date: 2004/01/14 20:38:44 $
 * 
 * @brief  functions for handling glmmStruct objects.
 * 
 */

#include <limits.h>
/* #include "bitfield.h"*/
#include "utilities.h"
#include "reStruct.h"

#ifndef max
#define max(m, n) (((m) > (n))?(m):(n))
#define min(m, n) (((m) < (n))?(m):(n))
#endif

/** 
 * Calculate the penalty term in the glmm log-likelihood.
 * 
 * Argument must be protected.
 *
 * @param glmmStruct A glmmStruct object
 * 
 * @return The logarithm of the penalty term in the log-likelihood
 */
static double
nlme_glmmLa2_penalty(SEXP glmmStruct)
{
    double ans = 0.;
    const SEXPREC* random = GET_SLOT((SEXP) glmmStruct, install("random"));
    double* bbetaVals = REAL(GET_SLOT((SEXP) glmmStruct, install("bbetas")));
    SEXP precision_sym = install("precision");
    SEXP storedRows_sym = install("storedRows");
    SEXP factor_sym = install("factor");
    SEXP columns_sym = install("columns");
    int nlevel = LENGTH((SEXP) random) - 2;
    int lev;
    
    for (lev = 0; lev < nlevel; lev++) {
        const SEXPREC* lmeLevel = VECTOR_ELT((SEXP)random, lev);
        const SEXPREC* storedRows = GET_SLOT((SEXP)lmeLevel,
                                             storedRows_sym);
        SEXP precision = GET_SLOT((SEXP)lmeLevel, precision_sym);
        double* Delta = REAL(GET_SLOT(precision, factor_sym));
        int q = LENGTH(GET_SLOT((SEXP)lmeLevel, columns_sym));
        int M = LENGTH((SEXP)storedRows);
        double *scratch = Calloc(q, double);
        double one = 1., zero = 0.;
        int i_one = 1;
        double tmp = 0.0;
        int j, k;
        
        for (j = 0; j < M; j++) {
            int startRow = INTEGER(VECTOR_ELT((SEXP)storedRows, j))[0];
            F77_CALL(dgemv)("N", &q, &q, &one, Delta, &q,
                            bbetaVals + startRow - 1, &i_one,
                            &zero, scratch, &i_one);

            for (k = 0; k < q; k++) tmp += scratch[k] * scratch[k];
        }
        Free(scratch);
        ans += tmp;
    }
    return -0.5*ans;
}

/* static void */
/* nlme_glmmLa2_calculateCorrection(SEXP original, SEXP weighted, SEXP stored, */
/*                                  SEXP random, const SEXPREC* eta, SEXP fam, */
/*                                  SEXP w) */
/* { */
/*     SEXP origRows_sym = install("originalRows"); */
/*     SEXP storedRows_sym = install("storedRows"); */
/*     SEXP columns_sym = install("columns"); */
/*     int nlevel = LENGTH(random)-2; */
/*     double* dorig = REAL(original); */
/*     double* dstored = REAL(stored); */
/*     int* dim = INTEGER(GET_DIM(original)); */
/*     int nrow = dim[0]; */
/*     int ncol = dim[1]; */
/*     int ncol_ranef = */
/*         asInteger(GET_SLOT(VECTOR_ELT(random, nlevel), columns_sym))-1; */
/*     double* dorig_last = dorig+nrow*(ncol-1); */
/*     int nrowStored = INTEGER(GET_DIM(stored))[0]; */
/*     double* dstored_last = dstored+nrowStored*(ncol-1); */
/*     const double zero = 0.0; */
/*     const double one = 1.0; */
/*     const double neg_one = -1.0; */
/*     const int i_one = 1; */
/*     double* wt2; */
/*     double* dorig_copy = Calloc(nrow*ncol_ranef, double); */
/*     int lev, i; */

/*     memset(dorig_last, 0, nrow*sizeof(double)); */
/*     memset(dstored_last, 0, nrowStored*sizeof(double)); */
/*     memcpy(dorig_copy, dorig, nrow*ncol_ranef*sizeof(double)); */
/*     for (lev = 0; lev < nlevel; lev++) { */
/*         const SEXPREC* lmeLevel = VECTOR_ELT((SEXP)random, lev); */
/*         const SEXPREC* origLevelRows = GET_SLOT((SEXP)lmeLevel, */
/*                                                 origRows_sym); */
/*         const SEXPREC* storedLevelRows = GET_SLOT((SEXP)lmeLevel, */
/*                                                   storedRows_sym); */
/*         SEXP columns = GET_SLOT((SEXP)lmeLevel, columns_sym); */
/*         int startCol = INTEGER(columns)[0]-1; */
/*         int q = LENGTH(columns); */
/*         int M = LENGTH((SEXP) origLevelRows); */
/*         int maxRow = nlme_scratchRowSize(origLevelRows); */
/*         int ncolRest = ncol_ranef - startCol - q; */
/*         for (i = 0; i < M; i++) { */
/*             SEXP origChunkRows = VECTOR_ELT((SEXP)origLevelRows, i); */
/*             int n = LENGTH(origChunkRows); */
/*             int origStartRow = INTEGER(origChunkRows)[0]-1; */
/*             int storedStartRow = INTEGER(VECTOR_ELT((SEXP)storedLevelRows, i))[0]-1; */
/*             double* a = dorig_last+origStartRow; */
/*             double* dZ = dorig_copy + nrow*startCol+origStartRow; */
/*             double* dR = dstored+nrowStored*startCol+storedStartRow; */
/*             int j; */

/*             F77_CALL(dtrsm)("R", "U", "N", "N", &n, &q, &one, dR, */
/*                             &nrowStored, dZ, &nrow); */
/*             for (j = 0; j < n; j++) { */
/*                 double* tmp2 = dZ + j; */
/*                 double tmp1 = 0.0; */
/*                 int k; */

/*                 for (k = 0; k < q; k++, tmp2 += n) { */
/*                     tmp1 += *tmp2 * *tmp2; */
/*                 } */
/*                 a[j] += tmp1; */
/*             } */
/*             if (ncolRest > 0) */
/*                 F77_CALL(dgemm)("N", "N", &n, &ncolRest, &q, &neg_one, */
/*                                 dZ, &nrow, dR+nrowStored*q, &nrowStored, */
/*                                 &one, dZ+nrow*q, &nrow); */
/*         } */
/*     } */
/*     Free(dorig_copy); */

/*     wt2 = REAL(eval(lang4(install("glmmLa2Wt2"), fam, (SEXP)eta, w), R_GlobalEnv)); */
/*     for (i = 0; i < nrow; i++) { */
/*         dorig_last[i] *= wt2[i]; */
/*     } */
/*     for (lev = 0; lev <= nlevel; lev++) { */
/*         const SEXPREC* lmeLevel = VECTOR_ELT((SEXP)random, lev); */
/*         const SEXPREC* origLevelRows = GET_SLOT((SEXP)lmeLevel, */
/*                                                 origRows_sym); */
/*         const SEXPREC* storedLevelRows = GET_SLOT((SEXP)lmeLevel, */
/*                                                   storedRows_sym); */
/*         SEXP columns = GET_SLOT((SEXP)lmeLevel, columns_sym); */
/*         int startCol = INTEGER(columns)[0]-1; */
/*         int q = LENGTH(columns); */
/*         int M = LENGTH((SEXP) origLevelRows); */
/*         double* Z = dorig + startCol*nrow; */
/*         double* Rs = dstored+nrowStored*startCol; */

/*         for (i = 0; i < M; i++) { */
/*             SEXP origChunkRows = VECTOR_ELT((SEXP)origLevelRows, i); */
/*             int n = LENGTH(origChunkRows); */
/*             int origStartRow = INTEGER(origChunkRows)[0]-1; */
/*             int storedStartRow = INTEGER(VECTOR_ELT((SEXP)storedLevelRows, i))[0]-1; */
/*             double* R = Rs+storedStartRow; */
/*             double* u = dstored_last+storedStartRow; */

/*             F77_CALL(dgemv)("T", &n, &q, &one, Z+origStartRow, &nrow, */
/*                             dorig_last+origStartRow, &i_one, */
/*                             &zero, u, &i_one); */
/*             F77_CALL(dtrsv)("U", "T", "N", &q, R, &nrowStored, u, &i_one); */
/*             if (lev < nlevel) */
/*                 F77_CALL(dtrsv)("U", "N", "N", &q, R, &nrowStored, u, &i_one); */
/*         } */
/*     } */
/* } */

/** 
 * Decompose the decomposed matrix and populate the stored matrix
 * 
 * Both SEXP arguments must be protected.
 *
 * @param reStruct An reStruct object
 * 
 * @return A scalar real with the log-likelihood
 */
static void
nlme_glmm_conditionalDecompose(SEXP reStruct, SEXP random, SEXP stored,
                               int origRowNum, int nlevel)
{
    const SEXPREC* decomposed = GET_SLOT(nlme_predecompose(reStruct),
                                         install("decomposed"));
    SEXP precision_sym = install("precision");
    SEXP storedRows_sym = install("storedRows");
    SEXP decomposedRows_sym = install("decomposedRows");
    SEXP columns_sym = install("columns");
    SEXP factor_sym = install("factor");
    int ncol_levels = asInteger(GET_SLOT(VECTOR_ELT(random, nlevel-1),
                                         columns_sym));
    int lev;
    int* dim = INTEGER(GET_DIM((SEXP)decomposed));
    int nrow = dim[0];
    int ncol = dim[1];
    int startCol = 0;
    double* mat;
    SEXP tmp;
    int nprotect = 0;
    
    if (ncol_levels > ncol) {
        error("Incorreect number of columns in decomposed matrix: %d instead of >= %d",
              ncol, ncol_levels);
    }

    /* checking input types - we may skip some of this in the future */
    if (stored != NULL) {
        if (TYPEOF(stored) != REALSXP)
            error("stored must be of storage mode double");
        dim = INTEGER(GET_DIM(stored));
        if (dim[0] != nrow || dim[1] != ncol)
            error("Dimensions of decomposed matrix do not match stored");
        memset(REAL(stored), 0, ncol*nrow*sizeof(double));
    }
    
    /*
     * Make sure decomposed is of mode double.
     * If so we make a duplicate of it for internal use.
     * Otherwise, we use the coerced version of decomposed.
     */
    tmp = coerceVector((SEXP)decomposed, REALSXP);
    if (tmp == decomposed) {
        mat = Calloc(nrow*ncol, double);
        memcpy(mat, REAL((SEXP)decomposed), nrow*ncol*sizeof(double));
    } else {
        nprotect = 1;
        mat = REAL(PROTECT(tmp));
    }
    
    for (lev = 0; lev < nlevel; lev++) {
        const SEXPREC* lmeLevel = VECTOR_ELT(random, lev);
        const SEXPREC* storedLevelRows = GET_SLOT((SEXP)lmeLevel,
                                                  storedRows_sym);
        const SEXPREC* decLevelRows = GET_SLOT((SEXP)lmeLevel,
                                               decomposedRows_sym);
        SEXP precision = GET_SLOT((SEXP)lmeLevel, precision_sym);
        SEXP Delta = ((lev < nlevel-2)?
                      GET_SLOT(precision, factor_sym):R_NilValue);
        int q = LENGTH(GET_SLOT((SEXP)lmeLevel, columns_sym));
        int M = LENGTH((SEXP)storedLevelRows);
        int ncolRest = ncol_levels-startCol-q;
        int maxRow = nlme_scratchRowSize(decLevelRows);
        double* qraux = Calloc(q, double);
        int* pivot = Calloc(q, int);
        int scratchCols = q+ncolRest;
        double* scratch = Calloc(((Delta==R_NilValue)?maxRow:(maxRow+q))*
                                 scratchCols, double);
        double* work;
        int j, info;
        double tmp;
        int lwork = -1;
        
        for (j = 0; j < q; j++)
            pivot[j] = 1;
        
        F77_CALL(dgeqp3)(&maxRow, &q, scratch, &maxRow,
                         pivot, qraux, &tmp, &lwork, &info);
        nlme_check_Lapack_error(info, "dgeqp3");
        j = (int) tmp;
        F77_CALL(dormqr)("Left", "Trans", &maxRow, &ncolRest,
                         &q, scratch, &maxRow, qraux,
                         scratch+maxRow*q, &maxRow, &tmp, &lwork,
                         &info);
        nlme_check_Lapack_error(info, "dormqr");
        if (j < (int) tmp)
            lwork = (int) tmp;
        else lwork = j;
        work = Calloc(lwork, double);
        
        /*
         *  Decompose and store
         */
        for (j = 0; j < M; j++) {
            nlme_decomposeChunk(mat, nrow, ncol_levels, mat, nrow,
                                ncol_levels, stored, Delta,
                                VECTOR_ELT((SEXP)decLevelRows, j),
                                VECTOR_ELT((SEXP)storedLevelRows, j),
                                startCol, ncol_levels, q, ncolRest,
                                NULL, NULL, scratch, qraux,
                                pivot, work, lwork);
        }
        
        startCol += q;
        Free(qraux);
        Free(pivot);
        Free(work);
        Free(scratch);
    }
    if (nprotect == 0) {
        Free(mat);
    } else {
        UNPROTECT(nprotect);
    }
}

SEXP
nlme_glmm_ranefIRLS(SEXP glmm,
                    const SEXPREC* nIRLS,
                    const SEXPREC* fixedLevels)
{
  int niter = asInteger((SEXP)nIRLS);
  int nFixedLevs = asInteger((SEXP)fixedLevels);
  SEXP eta;
  if (NAMED(glmm) && !asLogical(GET_SLOT(GET_SLOT(glmm, install("reStruct")),
                                         install("dontCopy"))))
      glmm = duplicate(glmm);
  PROTECT(glmm);
  {
    double eps = sqrt(DOUBLE_EPS);
    SEXP reStruct = GET_SLOT(glmm, install("reStruct"));
    SEXP random = GET_SLOT(reStruct, install("random"));
    SEXP stored = GET_SLOT(reStruct, install("stored"));
    SEXP bbetas = GET_SLOT(reStruct, install("bbetas"));
    SEXP original = GET_SLOT(reStruct, install("original"));
    SEXP weighted = GET_SLOT(reStruct, install("weighted"));
    SEXP fam = GET_SLOT(glmm, install("family"));
    SEXP origy = GET_SLOT(glmm, install("origy"));
    SEXP w = GET_SLOT(glmm, install("prior.weights"));
    int* dim = INTEGER(GET_DIM(original));
    int nrow = dim[0];
    int ncol = dim[1];
    int nrowStored = INTEGER(GET_DIM(stored))[0];
    int nlevel = LENGTH(random);
    int respcol = INTEGER(GET_SLOT(VECTOR_ELT(random, nlevel-1),
                                   install("columns")))[0] - 1;
    int ranefLen =
        asInteger(VECTOR_ELT(GET_SLOT(VECTOR_ELT(random, nlevel-nFixedLevs-1),
                                      install("storedRows")), 0))-1;
    SEXP respWtCall;
    SEXP respWt;
    double* dorig = REAL(original);
    double* dwt;
    double* deta;
    double* detaold;
    double* dweighted;
    int newWeighted = 1;
    int i, j, lev;

    if (LENGTH(weighted) != nrow*ncol) {
        SET_SLOT(reStruct, install("weighted"),
                 allocMatrix(REALSXP, nrow, ncol));
        weighted = GET_SLOT(reStruct, install("weighted"));
    } else if (asLogical(GET_SLOT(reStruct, install("useWeighted")))) {
        newWeighted = 0;
    }
    dweighted = REAL(weighted);

    eta = PROTECT(nlme_reStruct_fitted(reStruct, R_NilValue));
    deta = REAL(eta);
    respWtCall =
        PROTECT(lcons(install("glmmLa2RespWt"),
                      cons(fam,
                           list4(eta, origy, w,
                                 GET_SLOT(reStruct,
                                          install("offset"))))));
    detaold = Calloc(nrow, double);

    for (i = 0; i < niter; i++) {
        double dist = 0.0;
        double etanorm = 0.0;
        
        memcpy(detaold, deta,
               nrow*sizeof(double));
        respWt = PROTECT(eval(respWtCall, R_GlobalEnv));
        /* Set the response */
        memcpy(dorig+nrow*respcol,
               REAL(VECTOR_ELT(respWt, 0)), nrow*sizeof(double));
        /* multiply by the weights */
        dwt = REAL(VECTOR_ELT(respWt, 1));
        for (j = 0; j < ncol-1; j++) {
            int k;
            for (k = 0; k < nrow; k++) {
                dweighted[k+j*nrow] = dwt[k] * dorig[k+j*nrow];
            }
        }
        UNPROTECT(1);
        SET_SLOT(reStruct, install("dirtyDecomposed"), ScalarLogical(1));
        
        nlme_glmm_conditionalDecompose(reStruct, random, stored, nrow,
                                       nlevel);
        memcpy(REAL(bbetas),
               REAL(stored) + nrowStored*respcol, ranefLen*sizeof(double));
        if (nFixedLevs > 0) {
            SEXP lmeLevel = VECTOR_ELT(random, nlevel-nFixedLevs-1);
            int nlev = asInteger(GET_SLOT((SEXP) lmeLevel, install("nlev")));
            const SEXPREC* decomposedRows = GET_SLOT((SEXP) lmeLevel,
                                                     install("decomposedRows"));
            const SEXPREC* storedRows = GET_SLOT((SEXP) lmeLevel,
                                                 install("storedRows"));
            const SEXPREC* columns = GET_SLOT((SEXP) lmeLevel,
                                              install("columns"));
            int ldstored = INTEGER(GET_DIM(stored))[0];
            double* dBbetas = REAL(bbetas);
            int ncolR = LENGTH((SEXP)columns);
            int startCol = INTEGER((SEXP)columns)[0]-1;
            double* dStored = REAL(stored) + ldstored*startCol;
            const int one = 1;
            const double one_d = 1.0;
            const double neg_one_d = -1.0;

            for (j = 0; j < nlev; j++) {
                int startRow = INTEGER(VECTOR_ELT((SEXP)storedRows, j))[0];
                int nrowUpdate =
                    startRow-INTEGER(VECTOR_ELT((SEXP)decomposedRows, j))[0];

                /*
                 * Update the bbetas values for upper levels
                 */
                if (nrowUpdate > 0) {
                    double* curr_coef = dBbetas + --startRow;
                    double* R = dStored + startRow;
                    F77_CALL(dgemm)("N", "N", &nrowUpdate, &one, &ncolR, &neg_one_d,
                                    R-nrowUpdate, &ldstored, curr_coef, &ncolR, &one_d,
                                    curr_coef-nrowUpdate, &nrowUpdate);
                }
            }
        }
        for (lev = nlevel-nFixedLevs-2; lev >= 0; lev--) {
            nlme_estimate_level(stored, VECTOR_ELT(random, lev), bbetas);
        }

        nlme_reStruct_fitted_internal(reStruct, eta, R_NilValue);
        for (j = 0; j < nrow; j++) {
            double tmp = deta[j] - detaold[j];
            dist += tmp*tmp;
            etanorm += deta[j]*deta[j];
        }

        if (dist < etanorm*eps) {
            break;
        }
    }
    Free(detaold);
    SET_SLOT(reStruct, install("logLik"),
             eval(PROTECT(lang3(install("glmmLa2LogLikComp"),
                                glmm, eta)), R_GlobalEnv));

    UNPROTECT(3);
  }
  UNPROTECT(1);
  return glmm;
}

static double
nlme_glmmLaplace_logDets(SEXP reStruct)
{
    double logDets = 0.0;
    SEXP random = GET_SLOT(reStruct, install("random"));
    SEXP stored = GET_SLOT(reStruct, install("stored"));
    int* dim = INTEGER(GET_DIM(stored));
    int nrow = dim[0];
    double* dstored = REAL(stored);
    int nlevel = LENGTH(random)-2;
    SEXP storedRows_sym = install("storedRows");
    SEXP columns_sym = install("columns");
    SEXP precision_sym = install("precision");
    SEXP logDet_sym = install("logDet");
    int lev;

    for (lev = 0; lev < nlevel; lev++) {
        const SEXPREC* lmeLevel = VECTOR_ELT(random, lev);
        const SEXPREC* storedLevelRows = GET_SLOT((SEXP)lmeLevel,
                                                  storedRows_sym);
        SEXP columns = GET_SLOT((SEXP)lmeLevel, columns_sym);
        int q = LENGTH(columns);
        int qp1 = q+1;
        int M = LENGTH((SEXP)storedLevelRows);
        double* storedMat = dstored+(asInteger(columns)-1)*nrow;
        double logDets_lev = 0.0;
        int j;
        
        for (j = 0; j < M; j++) {
            double* R =
                storedMat+asInteger(VECTOR_ELT((SEXP)storedLevelRows, j))-1;
            int k;
            double logDet = 0.0;
            for (k = 0; k < q; k++)
                logDet += log(fabs(R[k*qp1]));
            logDets_lev += logDet;
        }
        logDets += M*asReal(GET_SLOT(GET_SLOT((SEXP)lmeLevel,
                                              precision_sym),
                                     logDet_sym)) - logDets_lev;
    }
    return logDets;
}

SEXP
nlme_glmmLaplace_solveOnly(SEXP glmm,
                           const SEXPREC* nIRLS,
                           const SEXPREC* fixedLevels)
{
    double logLikVal;
    SEXP logLik;
    SEXP reStruct= GET_SLOT(glmm, install("reStruct"));
    int nFixedLevs = asInteger((SEXP)fixedLevels);
    SEXP random, stored;
    int lev, ncol_levels, nlevel;

    if (NAMED(glmm) && !asLogical(GET_SLOT(reStruct,
                                           install("dontCopy")))) {
        glmm = duplicate(glmm);
    }
    PROTECT(glmm);
    glmm =
        nlme_glmm_ranefIRLS(glmm, nIRLS,
                            fixedLevels);
    UNPROTECT(1);
    PROTECT(glmm);
    reStruct = GET_SLOT(glmm, install("reStruct"));
    logLik = GET_SLOT(reStruct, install("logLik"));
    logLikVal = asReal(logLik);
    logLikVal += nlme_glmmLaplace_logDets(reStruct);
    logLikVal += nlme_glmmLa2_penalty(reStruct);
    REAL(logLik)[0] = logLikVal;
    SET_SLOT(reStruct, install("logLik"), logLik);

    random = GET_SLOT(reStruct, install("random"));
    stored = GET_SLOT(reStruct, install("stored"));
    nlevel = LENGTH(random);
    ncol_levels = asInteger(GET_SLOT(VECTOR_ELT(random, nlevel-1),
                                     install("columns")));

    for (lev = nlevel-nFixedLevs-2; lev >= 0; lev--) {
        nlme_invert_level(stored, VECTOR_ELT(random, lev), ncol_levels);
    }

    LOGICAL(GET_SLOT(reStruct, install("dirtyBbetas")))[0] = 0;
				/* FIXME - part of stored is dirty */
    LOGICAL(GET_SLOT(reStruct, install("dirtyStored")))[0] = 0; 
    UNPROTECT(1);
    return glmm;
}

SEXP
nlme_glmmLaplace_logLikelihood(SEXP glmm,
                               const SEXPREC* pars,
                               const SEXPREC* nIRLS,
                               const SEXPREC* fixedLevels)
{
    SEXP reStruct= GET_SLOT(glmm, install("reStruct"));
    SEXP logLik = GET_SLOT(reStruct,
                           install("logLik"));
    int newPars = length((SEXP)pars) > 0 && nlme_setParameters(&reStruct,
                                                               pars);
    if (newPars) {
        PROTECT(reStruct);
        if (NAMED(glmm) && !asLogical(GET_SLOT(reStruct,
                                               install("dontCopy")))) {
            glmm = duplicate(glmm);
        }
        PROTECT(glmm);
        SET_SLOT(glmm, install("reStruct"), reStruct);
        UNPROTECT(2);
    }
    PROTECT(glmm);
    if (newPars || (ISNA(asReal(logLik)))) {
        logLik = GET_SLOT(GET_SLOT(nlme_glmmLaplace_solveOnly(glmm, nIRLS,
                                                              fixedLevels),
                                   install("reStruct")),
                          install("logLik"));
    }
    UNPROTECT(1);
    return duplicate(logLik);
}

/** 
 * Calculate the updateFactor slot for a level
 * 
 * Arguments must be protected.
 *
 * @param stored The stored slot of the reStruct object
 * @param lmeLevel lmeLevel object for the level
 * @param bbetas The bbetas slot of the reStruct object
 * @param sigmainv Estimate of \sigma^{-1}
 * @param fixedStartCol First column in the fixed level.
 */
void
nlme_glmmLa2_factor_level(SEXP stored, SEXP lmeLevel, SEXP bbetas,
                          double sigmainv, int fixedStartCol)
{
    int nlev = asInteger(GET_SLOT((SEXP) lmeLevel, install("nlev")));
    const SEXPREC* columns = GET_SLOT(lmeLevel, install("columns"));
    const SEXPREC* storedRows = GET_SLOT((SEXP) lmeLevel,
                                         install("storedRows"));
    int* dim = INTEGER(GET_DIM(stored));
    int ldstored = dim[0];
    int ncol = LENGTH((SEXP)columns);
    int srcStartCol = INTEGER((SEXP)columns)[0] - 1;
    double* dStored = REAL(stored) + ldstored*srcStartCol - 1;
    double* stored_last = REAL(stored) + ldstored*(dim[1]-1);
    double* dBbetas = REAL(bbetas) - 1;
    int nxcol = fixedStartCol - srcStartCol;
    int nrow = (nxcol + 1) * nlev;
    double* factor = REAL(GET_SLOT(lmeLevel, install("updateFactor")));
    double* scratch = Calloc(max(nrow, 2*ncol)*ncol, double);
    double* qraux = Calloc(ncol, double);
    int* pivot = Calloc(ncol, int);
    double* work;
    double tmp;
    const double one = 1.0;
    int lwork = -1;
    int info, i;

    for (i = 0; i < nlev; i++) {
        int srcStartRow = INTEGER(VECTOR_ELT((SEXP) storedRows, i))[0];
        double* from = dStored + srcStartRow;
        double* coefs = dBbetas + srcStartRow;
        double* u = stored_last + srcStartRow;
        double* to = scratch + (nxcol+1)*i;
        int j;

        for (j = 0; j < ncol; j++, from++, to += nrow) {
            int k;

            for (k = 0; k < nxcol; k++) {
                to[k] = from[k*ldstored];
            }
            to[nxcol] = (coefs[j] - u[j])*sigmainv;
        }
    }

    F77_CALL(dgeqp3)(&nrow, &ncol, scratch, &nrow,
                     pivot, qraux, &tmp, &lwork, &info);
    nlme_check_Lapack_error(info, "dgeqp3");
    lwork = (int) tmp;
    work = Calloc(lwork, double);
    F77_CALL(dgeqp3)(&nrow, &ncol, scratch, &nrow, pivot,
                     qraux, work, &lwork, &info);
    nlme_check_Lapack_error(info, "dgeqp3");
    nlme_copyRows(scratch, nrow, ncol, factor,
                  ncol, ncol,
                  0, 0, ncol, 0, ncol, 0, ncol, 1, NULL, NULL);
    for (i = 0; i < nlev; i++) {
        double* u = stored_last + INTEGER(VECTOR_ELT((SEXP) storedRows, i))[0];
        double* to = scratch +i;
        int j;

        for (j = 0; j < ncol; j++, to += nlev)
            *to = u[j]*sigmainv;
    }
    lwork = -1;
    info = 0;
    F77_CALL(dgeqp3)(&nlev, &ncol, scratch, &nlev,
                     pivot, qraux, &tmp, &lwork, &info);
    nlme_check_Lapack_error(info, "dgeqp3");
    if ((int) tmp > lwork) {
        lwork = (int) tmp;
        Free(work);
        work = Calloc(lwork, double);
    }
    F77_CALL(dgeqp3)(&nlev, &ncol, scratch, &nlev,
                     pivot, qraux, work, &lwork, &info);
    nlme_check_Lapack_error(info, "dgeqp3");
    nlme_copyRows(scratch, nrow, ncol, scratch+ncol*ncol,
                  ncol, ncol,
                  0, 0, ncol, 0, ncol, 0, ncol, 1, NULL, NULL);
    F77_CALL(dtrmm)("L", "U", "T", "N", &ncol, &ncol, &one, scratch,
                    &ncol, scratch+ncol*ncol, &nrow);

    nlme_copyRows(factor, ncol, ncol, scratch,
                  ncol, ncol,
                  0, 0, ncol, 0, ncol, 0, ncol, 1, NULL, NULL);
    F77_CALL(dtrmm)("L", "U", "T", "N", &ncol, &ncol, &one, factor, &ncol,
                    scratch, &ncol);
    for (i = 0; i < ncol; i++) {
        int j, k;
        for (j = 0, k = ncol; j < ncol; j++, k++) {
            scratch[i+j*ncol] -= scratch[i+k*ncol];
        }
    }
    info = 0;
    F77_CALL(dpotrf)("U", &ncol, scratch, &ncol, &info);
    if (info == 0) {
        nlme_copyRows(scratch, nrow, ncol, factor,
                      ncol, ncol,
                      0, 0, ncol, 0, ncol, 0, ncol, 1, NULL, NULL);
    }

    Free(scratch);
    Free(work);
    Free(qraux);
    Free(pivot);
}

/**
 * Decomposition routine for computing gradient of a glmmStruct
 *
 * Calculates the \hat{b_i}'s and R_{11(i)}^{-1}'s and A's from that.
 * 
 * Arguments must be protected.
 *
 * @param glmmStruct An glmmStruct object
 * @param pars If of positive length, use to set the parameter values
 *             of glmmStruct
 * 
 * @return glmmStruct
 */
/* SEXP */
/* nlme_glmmLa2_commonDecompose(SEXP glmmStruct, const SEXPREC* pars) */
/* { */
/*     SEXP logLik, random, stored, bbetas; */
/*     int newPars = length((SEXP)pars) > 0 && nlme_setParameters(&glmmStruct, pars); */
/*     int nlevel, ncol_levels, fixedStartCol, lev; */
/*     double sigmainv = 1.0; */
/* /\*         double sqrtDF; *\/ */

/*     if (!(newPars || ISNA(asReal(GET_SLOT(glmmStruct, install("logLik")))) || */
/*           asLogical(GET_SLOT(glmmStruct, install("dirtyBbetas"))) || */
/*           asLogical(GET_SLOT(glmmStruct, install("dirtyStored"))))) */
/*         return glmmStruct; */

/*     if (NAMED(glmmStruct)) */
/*         glmmStruct = duplicate(glmmStruct); */
/*     PROTECT(glmmStruct); */
/*     logLik = GET_SLOT(glmmStruct, install("logLik")); */
/*     if (newPars || ISNA(REAL(logLik)[0])) */
/*         REAL(logLik)[0] = nlme_glmmLa2_logLikelihood_internal(glmmStruct); */


/*     random = GET_SLOT(glmmStruct, install("random")); */
/*     stored = GET_SLOT(glmmStruct, install("stored")); */
/*     bbetas = GET_SLOT(glmmStruct, install("bbetas")); */
/*     nlevel = LENGTH(random) - 2; */
/*     ncol_levels = asInteger(GET_SLOT(VECTOR_ELT(random, nlevel+1), */
/*                                      install("columns"))); */
/*     fixedStartCol =  */
/*         asInteger(GET_SLOT(VECTOR_ELT(random, nlevel), install("columns")))-1; */

/* /\*         sqrtDF = sqrt((double) nrows); *\/ */

/*     for (lev = nlevel; lev >= 0; lev--) { */
/*         nlme_invert_level(stored, VECTOR_ELT(random, lev), ncol_levels); */
/*     } */

/* /\*         sigmainv = REAL(stored)[dim[0]*ncol_levels-1]; *\/ */
/* /\*         sigmainv = sqrtDF/((sigmainv < 0.0)?-sigmainv:sigmainv); *\/ */

/*     for (lev = 0; lev < nlevel; lev++) { */
/*         nlme_glmmLa2_factor_level(stored, VECTOR_ELT(random, lev), bbetas, */
/*                                   sigmainv, fixedStartCol); */
/*     } */
/*     SET_SLOT(glmmStruct, install("dirtyStored"), ScalarLogical(0)); */
/*     UNPROTECT(1); */
/*     return glmmStruct; */
/* } */
