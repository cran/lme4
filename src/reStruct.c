/**
 * @file   reStruct.c
 * @author Saikat DebRoy <saikat@stat.wisc.edu>
 * @author Douglas Bates <bates@stat.wisc.edu>
 * @date   $Date: 2003/09/30 16:27:39 $
 * 
 * @brief  functions for handling reStruct objects.
 * 
 */

#include <limits.h>
#include "utilities.h"
#include "reStruct.h"

#ifndef max
#define max(m, n) (((m) > (n))?(m):(n))
#define min(m, n) (((m) < (n))?(m):(n))
#endif

/** 
 * Calculate the maximum length of elements in vecList
 * 
 * Arguments need not be protected.
 *
 * @param vecList A list of objects on which LENGTH can be applied
 * 
 * @return the maximum length
 */
int
nlme_scratchRowSize(const SEXPREC* vecList)
{
    int maxlen = 0;
    int n = LENGTH((SEXP)vecList);
    int i;
    
    for (i = 0; i < n; i++) {
        int len = LENGTH(VECTOR_ELT((SEXP)vecList, i));
        if (maxlen < len)
            maxlen = len;
    }
    return maxlen;
}


/** 
 * Zero out some rows in mat.
 *
 * If validRows is Non-NULL, the appropriate bits in it are
 * unset. Otherwise, actual rows in mat are set to zero.
 * 
 * Arguments need not be protected.
 *
 * @param mat A matrix
 * @param ldmat Leading dimension of mat
 * @param ncol Number of columns in mat
 * @param startRow First row to set to zero
 * @param endRow 1+last row to set to zero
 * @param validRows (optional) an nlme_bitfield object of length same
 *                  as nrows(mat) 
 */
static void
nlme_zeroRows(double* mat, int ldmat, int ncol, int startRow,
              int endRow, nlme_bitfield* validRows)
{
    if (validRows == NULL) {
        int nzeroBytes = (endRow-startRow)*sizeof(double);
        int i;

        mat += startRow;
        for (i = ncol-1; i >= 0; i--)
            memset(mat+i*ldmat, 0, nzeroBytes);
    } else {
        nlme_bitfield_unset(validRows, startRow, endRow);
    }
}

/** 
 * Copy rows from src to dest
 * 
 * Arguments need not be protected.
 *
 * @param src A matrix
 * @param nrowSrc Number of rows in src
 * @param ncolSrc Number of columns in src
 * @param dest A matrix of type numeric
 * @param nrowDest Number of rows in dest
 * @param ncolDest Number of columns in dest
 * @param srcStartCol First column of src to copy
 * @param destStartCol First column of dest to copy to
 * @param ncol Number of columns to copy
 * @param srcStartRow First row of src to copy
 * @param m Number of rows to copy
 * @param destStartRow First row of dest to copy to
 * @param mdest Maximum number of rows in dest that can be used for copying
 * @param upper Should we only copy the upper triangular part of src?
 * @param srcNonZero If non-NULL, only rows corresponding to set bits
 *                   in this are copied. Rows corresponding to unset
 *                   bits are not copied (not even as zero rows).
 * @param destNonZero If non-NULL rows that were copied to are set in this.
 */
void
nlme_copyRows(const double* src, int nrowSrc, int ncolSrc,
              double* dest, int nrowDest, int ncolDest,
              int srcStartCol, int destStartCol, int ncol,
              int srcStartRow, int m,
              int destStartRow, int mdest,
              int upper,
              nlme_bitfield *srcNonZero,
              nlme_bitfield *destNonZero)
{
    int srcEndRow, destEndRow;
    nlme_range range;
    const char *uplo = upper?"U":"N";

    if (srcStartCol < 0)
        error("can not copy from negative column indices");
    if (destStartCol < 0)
        error("can not copy to negative column indices");
    if (srcStartCol + ncol > ncolSrc)
        error("Number of columns to copy is too large: %d, %d, %d",
              srcStartCol, ncol, ncolSrc);
    if (destStartCol + ncol > ncolDest)
        error("Number of columns to copy to is too large");
    if (m <= 0 || srcStartRow < 0 || srcStartRow+m > nrowSrc) {
        error("invalid source matrix: %d, %d, %d",
              srcStartRow, m, nrowSrc);
    }
    if (mdest <= 0 || destStartRow < 0 || destStartRow+mdest > nrowDest) {
        error("invalid destinition matrix: %d, %d, %d",
              destStartRow, mdest, nrowDest);
    }
    if (upper && m < ncol) {
        error("invalid source matrix");
    }
    if ((upper && mdest < ncol)||
        (!upper && mdest < m))
        error("insufficient space in destination matrix, %d, %d",
              mdest, ncol);

    destEndRow = destStartRow;
    range.end = srcStartRow;
    srcEndRow = srcStartRow+m;
    /* copy src in the upper part of dest */
    if (upper) {
        nlme_bitfield_next_range(srcNonZero, srcEndRow, &range);
        F77_CALL(dlacpy)(uplo, &ncol, &ncol,
                         src+srcStartCol*nrowSrc+range.beg,
                         &nrowSrc,
                         dest+destStartCol*nrowDest+destEndRow,
                         &nrowDest);
        destEndRow += ncol;
    } else {
        while (nlme_bitfield_next_range(srcNonZero, srcEndRow, &range)) {
            int nrow = range.end - range.beg;
            F77_CALL(dlacpy)(uplo, &nrow, &ncol,
                             src+srcStartCol*nrowSrc+range.beg,
                             &nrowSrc,
                             dest+destStartCol*nrowDest+destEndRow,
                             &nrowDest);
            destEndRow += nrow;
        }
    }
    if (destNonZero != NULL) {
        nlme_bitfield_set(destNonZero, destStartRow, destEndRow);
        nlme_bitfield_unset(destNonZero, destEndRow, destStartRow+mdest);
    } else if (destEndRow < destStartRow+mdest) {
        nlme_zeroRows(dest, nrowDest, ncolDest, destEndRow,
                      destStartRow+mdest, NULL);
    }
}

/** 
 * Decompose a single chunk in a particular level.
 * 
 * Arguments must be protected.
 *
 * @param srcmat Contains the matrix ZXy to be decomposed
 * @param nrowSrc Number of rows in source
 * @param ncolSrc Number of columns in source
 * @param destmat Store all but first q rows of Q'Xy here : often same as srcmat
 * @param nrowDest Number of rows in dest
 * @param ncolDest Number of columns in dest
 * @param storemat If non-NULL store first q rows of Q'ZXy here
 * @param Delta If non-NULL square matrix of dimension (q, q) to be appended below Z before decomposition
 * @param srcRowIndx Row indices for srcmat and destmat
 * @param storeRowIndx Row indices for storemat matrix
 * @param startCol First column of ZXy
 * @param ncol Number of columns in ZXy
 * @param q Number of columns in Z
 * @param ncolRest Number of columns in Xy
 * @param srcValidRows If non-NULL a bitfield vector indicating non-zero rows in srcmat
 * @param destValidRows If non-NULL a bitfield vector indicating non-zero rows in destmat
 * @param scratch A matrix of at least length(srcRowIndx)+q rows and q+ncolRest columns
 * @param qraux A double vector of length at least q
 * @param pivot An int vector of length at least q
 * @param work A double vector of length lwork
 * @param lwork length of work
 * 
 * @return log of the absolute value of determinant of R, which can be
 *         -Inf if the total number of rows is less than q.
 */
double
nlme_decomposeChunk(const double* srcmat, int nrowSrc, int ncolSrc,
                    double* destmat, int nrowDest, int ncolDest, SEXP storemat,
                    const SEXPREC* Delta, const SEXPREC* srcRowIndx,
                    const SEXPREC* storeRowIndx, int startCol, int ncol,
                    int q, int ncolRest, nlme_bitfield* srcValidRows,
                    nlme_bitfield* destValidRows, double* scratch,
                    double* qraux, int* pivot, double* work, int lwork)
{
    int m = LENGTH((SEXP)srcRowIndx);
    int srcStartRow = INTEGER((SEXP)srcRowIndx)[0]-1;
    int storeStartRow = INTEGER((SEXP)storeRowIndx)[0]-1;
    int nrow = (Delta==R_NilValue)?m:(m+q);
    int scratchCols = q+ncolRest;
    double ans = 0.0;
    int info, i;

    memset(scratch, 0, nrow*(q+ncolRest)*sizeof(double));

    nlme_copyRows(srcmat, nrowSrc, ncolSrc, scratch,
                  nrow, scratchCols,
                  startCol, 0, scratchCols,
                  srcStartRow, m, 0, m,
                  FALSE, srcValidRows, NULL);
    if (Delta != R_NilValue) {
        nlme_copyRows(REAL((SEXP)Delta), q, q, scratch, nrow,
                      scratchCols, 0, 0, q, 0, q,
                      m, q, TRUE,
                      NULL, NULL);
    }
    F77_CALL(dgeqp3)(&nrow, &q, scratch, &nrow, pivot,
                     qraux, work, &lwork, &info);
    nlme_check_Lapack_error(info, "dgeqp3");
    if (storemat != NULL) {
        nlme_copyRows(scratch, nrow, scratchCols, REAL(storemat),
                      INTEGER(GET_DIM(storemat))[0],
                      INTEGER(GET_DIM(storemat))[1], 0, startCol, q,
                      0, min(q, nrow), storeStartRow, q, TRUE, NULL, NULL);
    }
    if (ncolRest > 0) {
        startCol += q;
        F77_CALL(dormqr)("Left", "Trans", &nrow, &ncolRest, &q,
                         scratch, &nrow, qraux,
                         scratch+q*nrow, &nrow, work, &lwork,
                         &info);
        if (storemat != NULL) {
            nlme_copyRows(scratch, nrow, scratchCols, REAL(storemat),
                          INTEGER(GET_DIM(storemat))[0],
                          INTEGER(GET_DIM(storemat))[1], q, startCol,
                          ncolRest, 0, min(q, nrow), storeStartRow, q,
                          FALSE, NULL, NULL);
        }
        if (nrow > q) {
            nlme_copyRows(scratch, nrow, scratchCols, destmat,
                          nrowDest, ncolDest, q, startCol, ncolRest, q,
                          nrow-q, srcStartRow, m, FALSE, NULL,
                          destValidRows);
        } else {
            nlme_zeroRows(destmat, nrowDest, ncolDest, srcStartRow,
                          m+srcStartRow, destValidRows);
        }
    }
    if (nrow < q)
        return R_NegInf;
    for (i = (q-1)*(nrow+1); i >= 0; i-= nrow+1) {
        ans += log(fabs(scratch[i]));
    }
    return ans;
}

/** 
 * Pre-decompose the model matrix in the weighted slot, if present, otherwise
 * the model matrix in the original slot, and populate the decomposed slot
 * 
 * Argument must be protected.
 *
 * @param reStruct An reStruct object
 */
SEXP
nlme_predecompose(SEXP reStruct)
{
    int useWeighted = asLogical(GET_SLOT(reStruct,
                                         install("useWeighted")));
    const SEXPREC* original = GET_SLOT(reStruct,
                                       (useWeighted?install("weighted"):
                                        install("original")));
    SEXP decomposed = GET_SLOT(reStruct, install("decomposed"));
    const SEXPREC* random = GET_SLOT(reStruct, install("random"));
    int nlevel = LENGTH((SEXP)random);
    int* dim = INTEGER(GET_DIM((SEXP)original));
    int nrow = dim[0];
    int ncol = dim[1];
    if (LENGTH(decomposed) == 0) {
        int nrowDecomposed =
            LENGTH(VECTOR_ELT(GET_SLOT(VECTOR_ELT((SEXP)random, nlevel-1),
                                       install("decomposedRows")), 0));
        SET_SLOT(reStruct, install("decomposed"),
                 allocMatrix(REALSXP, nrowDecomposed, ncol));
        decomposed = GET_SLOT(reStruct, install("decomposed"));
        memset(REAL(decomposed), 0, ncol*nrowDecomposed*sizeof(double));
        SET_SLOT(reStruct, install("dirtyDecomposed"), ScalarLogical(1));
    }
    if (asLogical(GET_SLOT(reStruct, install("dirtyDecomposed")))) {
        SEXP originalRows_sym = install("originalRows");
        SEXP storedRows_sym = install("storedRows");
        SEXP columns_sym = install("columns");
        int ncol_levels = asInteger(GET_SLOT(VECTOR_ELT((SEXP)random, nlevel-1),
                                             columns_sym));
        int startCol = 0;
        nlme_bitfield* validRows = NULL;
        const double* srcmat;
        double* destmat;
        SEXP tmp;
        int lev;
        int nprotect = 0;

        /* various error checking - we may skip some of them in future if
         * we are confident enough */

        if (ncol_levels > ncol) {
            error("Incoreect number of columns in original matrix: %d instead of >= %d",
                  ncol, ncol_levels);
        }

        if (TYPEOF(decomposed) != REALSXP)
            error("decomposed must be of storage mode double");
        dim = INTEGER(GET_DIM(decomposed));
        if (dim[1] != ncol)
            error("Column dimension of decomposed matrix do not match original");
    
        /*
         * Make sure original is of mode double.
         * If so we allocate new memory for destmat. Otherwise, we use
         * the coerced version of original as destmat.
         */
        tmp = coerceVector((SEXP)original, REALSXP);
        if (tmp == original) {
            destmat = Calloc(nrow*ncol, double);
            srcmat = REAL((SEXP)original);
        } else {
            nprotect = 1;
            srcmat = destmat = REAL(PROTECT(tmp));
        }
        
        validRows = nlme_bitfield_alloc(nrow);
        nlme_bitfield_set(validRows, 0, nrow);

        for (lev = 0; lev < nlevel; lev++) {
            const SEXPREC* lmeLevel = VECTOR_ELT((SEXP)random, lev);
            const SEXPREC* origLevelRows = GET_SLOT((SEXP)lmeLevel,
                                                    originalRows_sym);
            const SEXPREC* decLevelRows = GET_SLOT((SEXP)lmeLevel,
                                                   storedRows_sym);
            int q = LENGTH(GET_SLOT((SEXP)lmeLevel, columns_sym));
            int M = LENGTH((SEXP)origLevelRows);
            int ncolRest = ncol_levels-startCol-q;
            int maxRow = nlme_scratchRowSize(origLevelRows);
            int scratchCols = q+ncolRest;
            double* work;
            int info, j;
            double tmp;
            int lwork = -1;
            double* scratch;
            double* qraux;
            int* pivot;


            scratch = Calloc(maxRow*scratchCols, double);
            qraux = Calloc(q, double);
            pivot = Calloc(q, int);
            for (j = 0; j < q; j++)
                pivot[j] = 1;

            /* determine the size of the work array and allocate it only once */
            F77_CALL(dgeqp3)(&maxRow, &q, scratch, &maxRow,
                             pivot, qraux, &tmp, &lwork, &info);
            nlme_check_Lapack_error(info, "dgeqp3");
            j = (int) tmp;
            F77_CALL(dormqr)("Left", "Trans", &maxRow, &ncolRest,
                             &q, scratch, &maxRow, qraux,
                             scratch+maxRow*q, &maxRow, &tmp, &lwork,
                             &info);
            nlme_check_Lapack_error(info, "dormqr");
            lwork = max((int) tmp, j);
            work = Calloc(lwork, double);

            /*
             *  Decompose srcmat and store in destmat
             */
            for (j = 0; j < M; j++)
                nlme_decomposeChunk(srcmat, nrow, ncol_levels,
                                    destmat, nrow, ncol_levels,
                                    decomposed, R_NilValue,
                                    VECTOR_ELT((SEXP)origLevelRows, j),
                                    VECTOR_ELT((SEXP)decLevelRows, j),
                                    startCol, ncol_levels, q, ncolRest,
                                    validRows, validRows, scratch,
                                    qraux, pivot, work, lwork);

            startCol += q;
            if (destmat != srcmat) {
                srcmat = destmat;
            }
            Free(qraux);
            Free(pivot);
            Free(work);
            Free(scratch);
        }
        nlme_bitfield_free(validRows);
        if (nprotect == 0) {
            Free(destmat);
        } else {
            UNPROTECT(nprotect);
        }
        SET_SLOT(reStruct, install("dirtyDecomposed"), ScalarLogical(0));
        REAL(GET_SLOT(reStruct, install("logLik")))[0] = NA_REAL;
    }
    return reStruct;
}

/** 
 * Set the parameters in an reStruct object if parameter value is new
 * 
 * The first argument must be protected
 *
 * @param reStruct An reStruct object
 * @param pars The vector used to set the parameters. The length
 *             of the vector is assumed to correct and is not
 *             checked. The vector is coerced to type REALSXP.
 * 
 * @return 1 if parameter was set and 0 otherwise.
 */
int
nlme_setParameters(SEXP* reStructPtr, const SEXPREC* pars)
{
    SEXP reStruct = *reStructPtr;
    SEXP tmp =
        PROTECT((length((SEXP)pars) > 0)?coerceVector((SEXP) pars,
                                                      REALSXP):R_NilValue);
    SEXP random = GET_SLOT(reStruct, install("random"));
    int nlevel = LENGTH(random)-2;
    SEXP parsInd_sym = install("parsInd");
    SEXP coefcall = PROTECT(lang2(install("coef"), reStruct));
    double* dpars = REAL(tmp);
    int lev;

    pars = (const SEXPREC*) tmp;
    tmp = eval(coefcall, R_GlobalEnv);
    UNPROTECT(1);
    for (lev = LENGTH(tmp)-1;
         lev >= 0 && REAL(tmp)[lev] == dpars[lev];
         lev--) {
    }
    if (lev == -1) {
        UNPROTECT(1);
        return 0;
    }
    dpars--;

    coefcall = PROTECT(lang3(install("coef<-"), R_NilValue,
                             R_NilValue));
    reStruct = PROTECT((NAMED(*reStructPtr) && !asLogical(GET_SLOT(*reStructPtr, install("dontCopy"))))?duplicate(*reStructPtr):*reStructPtr);
    random = GET_SLOT(reStruct, install("random"));
    for (lev = 0; lev < nlevel; lev++) {
        SEXP lmeLevel = VECTOR_ELT(random, lev);
        SEXP parsInd = GET_SLOT(lmeLevel, parsInd_sym);
        int parlen = LENGTH(parsInd);
        SEXP newpar = PROTECT(allocVector(REALSXP, parlen));
        SEXP tmp;
        
        SETCADR(coefcall, lmeLevel);
        SETCADDR(coefcall, newpar);
        memcpy(REAL(newpar), dpars+INTEGER(parsInd)[0],
               parlen*sizeof(double));
        tmp = eval(coefcall, R_GlobalEnv);
        if (tmp != lmeLevel) {
            SET_VECTOR_ELT(random, lev, tmp);
        }
        UNPROTECT(1);
    }
    REAL(GET_SLOT(reStruct, install("logLik")))[0] = NA_REAL;
    UNPROTECT(3);
    *reStructPtr = reStruct;
    return 1;
}

/** 
 * Decompose the decomposed matrix and optionally populate the stored matrix
 * returning the log-likelihood
 * 
 * Both SEXP arguments must be protected.
 *
 * @param reStruct An reStruct object
 * @param store Logical indicator of whether to populate the stored slot
 * @param GLMMLaplace2 If TRUE, calculate a component of the 2nd order Laplace
 *             approximation to the log-likelihood for a GLMM models
 * 
 * @return A scalar real with the log-likelihood
 */
double
nlme_logLikelihood_internal(SEXP reStruct, int store,
                            int GLMMLaplace2)
{
    double logLik = 0.0;
    SEXP stored = store?GET_SLOT(reStruct, install("stored")):NULL;
    const SEXPREC* decomposed = GET_SLOT(nlme_predecompose(reStruct),
                                         install("decomposed"));
    SEXP random = GET_SLOT(reStruct, install("random"));
    SEXP precision_sym = install("precision");
    SEXP storedRows_sym = install("storedRows");
    SEXP decomposedRows_sym = install("decomposedRows");
    SEXP columns_sym = install("columns");
    SEXP factor_sym = install("factor");
    SEXP logDet_sym = install("logDet");
    int REML = asLogical(GET_SLOT(reStruct, install("REML")));
    int origRowNum = INTEGER(GET_DIM(GET_SLOT(reStruct,
                                              install("original"))))[0];
    int nlevel = LENGTH(random);
    int ncol_levels = asInteger(GET_SLOT(VECTOR_ELT(random, nlevel-1),
                                         columns_sym));
    int lev, p = 0;
    int* dim = INTEGER(GET_DIM((SEXP)decomposed));
    int nrow = dim[0];
    int ncol = dim[1];
    int startCol = 0;
    double* mat;
    SEXP tmp;
    double logLikComp = 0.0;
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
        
        logLikComp = 0.0;
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
            logLikComp +=
                nlme_decomposeChunk(mat, nrow, ncol_levels, mat, nrow,
                                    ncol_levels, stored, Delta,
                                    VECTOR_ELT((SEXP)decLevelRows, j),
                                    VECTOR_ELT((SEXP)storedLevelRows, j),
                                    startCol, ncol_levels, q, ncolRest,
                                    NULL, NULL, scratch, qraux,
                                    pivot, work, lwork);
        }
        
        if (Delta != R_NilValue) {
            logLik += M*asReal(GET_SLOT(precision, logDet_sym))-logLikComp;
        } else if (lev == nlevel-2) { /* fixed-effects */
            if (!GLMMLaplace2 && REML)
                logLik -= logLikComp;
            p = q;
        } else if (lev == nlevel-1) { /* response */
            if (!GLMMLaplace2) {
                double mult = origRowNum - (REML ? p : 0);
                logLik +=
                    mult*(-logLikComp
                          + (log(mult)- log(2*PI)-1)/2);
            }
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
    if (stored != NULL) {
        memcpy(REAL(GET_SLOT(reStruct, install("bbetas"))),
               REAL(stored) + nrow*(ncol_levels-1), nrow*sizeof(double));
    }
    LOGICAL(GET_SLOT(reStruct, install("dirtyBbetas")))[0] = 1;
    return logLik;
}

/** 
 * Decompose the decomposed matrix and optionally populate the stored matrix
 * returning the log-likelihood
 * 
 * Only the first argument must be protected.
 *
 * @param reStruct An reStruct object
 * @param store Logical indicator of whether to populate the stored slot
 * @param pars If of positive length, use to set the parameter values
 *             of reStruct
 * 
 * @return A scalar real with the log-likelihood
 */
SEXP
nlme_logLikelihood(SEXP reStruct, const SEXPREC* pars)
{
    int newPars = length((SEXP)pars) > 0 && nlme_setParameters(&reStruct, pars);
    SEXP logLik;

    if (!newPars && !ISNA(asReal(GET_SLOT(reStruct, install("logLik")))))
        return duplicate(GET_SLOT(reStruct, install("logLik")));

    if (NAMED(reStruct) && !asLogical(GET_SLOT(reStruct, install("dontCopy"))))
        reStruct = duplicate(reStruct);
    PROTECT(reStruct);

    logLik = GET_SLOT(reStruct, install("logLik"));
    REAL(logLik)[0] = nlme_logLikelihood_internal(reStruct, 0, 0);
    UNPROTECT(1);
    return duplicate(logLik);
}

/**
 * Calculate \hat{\beta} or \hat{b_i}'s (depending on the level).
 *
 * Before the first call to this function we copy the last column of
 * the stored matrix to bbetas. We then call this function with the
 * lmeLevel objects, starting with the last but one and going back to
 * the first one. Each call puts the estimate of b_i's (beta for the
 * first call) in bbetas and modifies the bbetas values above this
 * level for use in future calls.
 *
 * Arguments need not be protected.
 *
 * @param stored The stored slot
 * @param lmeLevel The lmeLevel object for this level
 * @param bbetas The bbetas slot
 */
void
nlme_estimate_level(SEXP stored, const SEXPREC* lmeLevel, SEXP bbetas)
{
    int nlev = asInteger(GET_SLOT((SEXP) lmeLevel, install("nlev")));
    const SEXPREC* decomposedRows = GET_SLOT((SEXP) lmeLevel,
                                             install("decomposedRows"));
    const SEXPREC* storedRows = GET_SLOT((SEXP) lmeLevel,
                                         install("storedRows"));
    const SEXPREC* columns = GET_SLOT((SEXP) lmeLevel,
                                      install("columns"));
    int ldstored = INTEGER(GET_DIM(stored))[0];
    double* dBbetas = REAL(bbetas);
    int ncol = LENGTH((SEXP)columns);
    int startCol = INTEGER((SEXP)columns)[0]-1;
    double* dStored = REAL(stored) + ldstored*startCol;
    const int one = 1;
    const double one_d = 1.0;
    const double neg_one_d = -1.0;
    int j;

    for (j = 0; j < nlev; j++) {
        int startRow = INTEGER(VECTOR_ELT((SEXP)storedRows, j))[0];
        int nrow = startRow-INTEGER(VECTOR_ELT((SEXP)decomposedRows, j))[0];
        double* curr_coef = dBbetas + --startRow;
        double* R = dStored + startRow;

        /*
         * Calculate the estimate
         */
        F77_CALL(dtrsm)("L", "U", "N", "N", &ncol, &one, &one_d,
                        R, &ldstored, curr_coef, &ncol);

        /*
         * Update the bbetas values for upper levels
         */
        if (nrow > 0) {
            F77_CALL(dgemm)("N", "N", &nrow, &one, &ncol, &neg_one_d,
                            R-nrow, &ldstored, curr_coef, &ncol, &one_d,
                            curr_coef-nrow, &nrow);
        }
    }
}

/**
 * Invert the R matrix for the level
 *
 * Arguments need not be protected.
 *
 * @param stored The stored slot
 * @param lmeLevel The lmeLevel object for this level
 * @param ncol_levels Total number of columns in all the levels
 */
void
nlme_invert_level(SEXP stored, const SEXPREC* lmeLevel, int ncol_levels)
{
    int nlev = asInteger(GET_SLOT((SEXP) lmeLevel, install("nlev")));
    const SEXPREC* decomposedRows = GET_SLOT((SEXP) lmeLevel,
                                             install("decomposedRows"));
    const SEXPREC* storedRows = GET_SLOT((SEXP) lmeLevel,
                                         install("storedRows"));
    const SEXPREC* columns = GET_SLOT((SEXP) lmeLevel,
                                      install("columns"));
    int ldstored = INTEGER(GET_DIM(stored))[0];
    int ncol = LENGTH((SEXP)columns);
    int startCol = INTEGER((SEXP)columns)[0];
    int ncolRest = ncol_levels-ncol-startCol;
    double* dStored = REAL(stored) + ldstored*--startCol;
    const double one_d = 1.0;
    const double neg_one_d = -1.0;
    int j;

    for (j = 0; j < nlev; j++) {
        int startRow = INTEGER(VECTOR_ELT((SEXP)storedRows, j))[0];
        int nrow = startRow-INTEGER(VECTOR_ELT((SEXP)decomposedRows, j))[0];
        double* R = dStored + --startRow;
        int info;

        /*
         * Calculate inverse of R
         */
        F77_CALL(dtrtri)("U", "N", &ncol, R, &ldstored, &info);
        nlme_check_Lapack_error(info, "dtrtri");

        if (ncolRest > 0) {
            F77_CALL(dtrmm)("L", "U", "N", "N", &ncol, &ncolRest,
                            &neg_one_d, R, &ldstored, R+ncol*ldstored,
                            &ldstored);
            if (nrow > 0) {
                F77_CALL(dgemm)("N", "N", &nrow, &ncolRest, &ncol, &one_d,
                                R-nrow, &ldstored, R+ncol*ldstored, &ldstored,
                                &one_d, R+ncol*ldstored-nrow, &ldstored);
            }
        }
        if (nrow > 0) {
            F77_CALL(dtrmm)("R", "U", "N", "N", &nrow, &ncol,
                            &one_d, R, &ldstored, R-nrow,
                            &ldstored);
        }
    }
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
nlme_factor_level(SEXP stored, SEXP lmeLevel, SEXP bbetas,
                  double sigmainv, int fixedStartCol)
{
    int nlev = asInteger(GET_SLOT((SEXP) lmeLevel, install("nlev")));
    const SEXPREC* columns = GET_SLOT(lmeLevel, install("columns"));
    const SEXPREC* storedRows = GET_SLOT((SEXP) lmeLevel,
                                         install("storedRows"));
    int ldstored = INTEGER(GET_DIM(stored))[0];
    int ncol = LENGTH((SEXP)columns);
    int srcStartCol = INTEGER((SEXP)columns)[0] - 1;
    double* dStored = REAL(stored) + ldstored*srcStartCol - 1;
    double* dBbetas = REAL(bbetas) - 1;
    int nxcol = fixedStartCol - srcStartCol;
    int nrow = (nxcol + 1) * nlev;
    double* scratch = Calloc(nrow*ncol, double);
    double* qraux = Calloc(ncol, double);
    int* pivot = Calloc(ncol, int);
    double* work;
    double tmp;
    int lwork = -1;
    int info, i;

    for (i = 0; i < nlev; i++) {
        int srcStartRow = INTEGER(VECTOR_ELT((SEXP) storedRows, i))[0];
        double* from = dStored + srcStartRow;
        double* coefs = dBbetas + srcStartRow;
        double* to = scratch + (nxcol+1)*i;
        int j;

        for (j = 0; j < ncol; j++, from++, coefs++, to += nrow) {
            int k;

            for (k = 0; k < nxcol; k++) {
                to[k] = from[k*ldstored];
            }
            to[nxcol] = *coefs*sigmainv;
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
    nlme_copyRows(scratch, nrow, ncol,
                  REAL(GET_SLOT(lmeLevel, install("updateFactor"))),
                  ncol, ncol,
                  0, 0, ncol, 0, ncol, 0, ncol, 1, NULL, NULL);
    Free(scratch);
    Free(work);
    Free(qraux);
    Free(pivot);
}

void
nlme_hessian_level(SEXP stored, SEXP lmeLevel,
                   SEXP bbetas, double sigmainv, int fixedStartColREML,
                   int fixedStartColML, int ntotal, int REML)
{
    int nlev = asInteger(GET_SLOT((SEXP) lmeLevel, install("nlev")));
    const SEXPREC* columns = GET_SLOT(lmeLevel, install("columns"));
    const SEXPREC* storedRows = GET_SLOT((SEXP) lmeLevel,
                                         install("storedRows"));
    int ldstored = INTEGER(GET_DIM(stored))[0];
    int ncol = LENGTH((SEXP)columns);
    int srcStartCol = INTEGER((SEXP)columns)[0] - 1;

/*     int nxcolREML = fixedStartColREML - srcStartCol; */
    int nxcolML = fixedStartColML - srcStartCol;
    int p = fixedStartColREML - fixedStartColML;

    double* dStored = REAL(stored) + ldstored*srcStartCol - 1;
    double* dBbetas = REAL(bbetas) - 1;
    int ncolsq = ncol*ncol;
    double* scratch;
    double* scratch1;
    double* scratch2;
    double* scratch3;
    double* hessarray;
    double sigma2inv = sigmainv*sigmainv;
    const double one = 1.0;
    const double zero = 0.0;
    int j, k;
    int lev;

    hessarray = Calloc(ncolsq*ncolsq, double);
    scratch = Calloc(ncolsq, double);
    scratch1 = Calloc(ncolsq, double);
    scratch2 = Calloc(ncolsq, double);
    scratch3 = Calloc(p*ncolsq, double);

    for (lev = 0; lev < nlev; lev++) {
        int srcStartRow = INTEGER(VECTOR_ELT((SEXP) storedRows, lev))[0];
        double* R = dStored + srcStartRow;
        double* Rmore = R + nxcolML*ldstored;
        double* this_b = dBbetas + srcStartRow;

        for (j = 0, k = 0; j < ncol; j++) {
            double b_j = this_b[j]*sigma2inv;
            int i;
            for (i = 0; i < ncol; i++) {
                double tmp = this_b[i]*b_j;
                scratch[k] += tmp;
                scratch1[k++] = tmp;
            }
        }
        F77_CALL(dgemm)("N", "T", &ncol, &ncol, &nxcolML, &one,
                        R, &ldstored, R, &ldstored, &zero,
                        scratch2, &ncol);
        for (j = 0, k = 0; j < ncolsq; j++) {
            double RR_j = scratch2[j];
            int i;
            for (i = 0; i < ncolsq; i++) {
                hessarray[k++] += scratch2[i]*RR_j;
            }
        }

        for (j = 0, k = 0; j < ncolsq; j++) {
            double bb_j = scratch1[j];
            /* double RR_j = scratch2[j]; */
            int i;
            for (i = 0; i < ncolsq; i++) {
                hessarray[k++] += 2*scratch2[i]*bb_j;
            }
        }

        for (j = 0; j < ncol; j++) {
            double b_j = this_b[j]*sigma2inv;
            int i;
            for (i = 0; i < ncol; i++) {
                double* tmp = scratch3 + p*(ncol*j+i);
                for (k = 0; k < p; k++) {
                    tmp[k] += b_j*Rmore[ldstored*i+k];
                }
            }
        }
    }

    Free(scratch1);
    Free(scratch2);

    scratch1 = REAL(GET_SLOT(lmeLevel, install("hessianArray")));
    /* This is basically aperm(hessarray, c(4, 1:3)) */
    for (j = 0, k = 0; j < ncol; j++) {
        double* tmp = hessarray+j*ncolsq*ncol;
        int i;
        for (i = 0; i < ncolsq*ncol; i++) {
            scratch1[i*ncol+j] = tmp[i];
        }
    }
    Free(hessarray);
    hessarray = scratch1;
    for (j = 0, k = 0; j < ncolsq; j++) {
        double* tmp1 = scratch3+p*j;
        double tmp = scratch[j]/ntotal;
        int i;
        for (i = 0; i < ncolsq; i++) {
            double* tmp2 = scratch3+p*i;
            double sum = scratch[i]*tmp;
            int l;
            for (l = 0; l < p; l++) {
                sum += tmp1[k]*tmp2[k];
            }
            hessarray[k++] += sum;
        }
    }
    Free(scratch3);
    Free(scratch);
}

/**
 * Decomposition routine common to EM and gradient computation.
 *
 * Calculates the \hat{b_i}'s and R_{11(i)}^{-1}'s and A's from that.
 * 
 * Arguments must be protected.
 *
 * @param reStruct An reStruct object
 * @param pars If of positive length, use to set the parameter values
 *             of reStruct
 * 
 * @return reStruct
 */
SEXP
nlme_commonDecompose(SEXP reStruct, const SEXPREC* pars)
{
    int newPars = length((SEXP)pars) > 0 && nlme_setParameters(&reStruct, pars);
    int analyticHessian = asLogical(GET_SLOT(reStruct,
                                             install("analyticHessian")));
    int REML = asLogical(GET_SLOT(reStruct, install("REML")));
    SEXP logLik;

    if (!(newPars || ISNA(asReal(GET_SLOT(reStruct, install("logLik")))) ||
          asLogical(GET_SLOT(reStruct, install("dirtyBbetas"))) ||
          asLogical(GET_SLOT(reStruct, install("dirtyStored")))))
        return reStruct;

    if (NAMED(reStruct) && !asLogical(GET_SLOT(reStruct, install("dontCopy"))))
        reStruct = duplicate(reStruct);
    PROTECT(reStruct);
    
    logLik = GET_SLOT(reStruct, install("logLik"));

    {
        SEXP random = GET_SLOT(reStruct, install("random"));
        SEXP stored = GET_SLOT(reStruct, install("stored"));
        SEXP bbetas = GET_SLOT(reStruct, install("bbetas"));
        int* dim = INTEGER(GET_DIM(stored));
        int nlevel = LENGTH(random) - 2;
        int ncol_levels = asInteger(GET_SLOT(VECTOR_ELT(random, nlevel+1),
                                             install("columns")));
        int fixedStartColREML = 
            asInteger(GET_SLOT(VECTOR_ELT(random, nlevel+1),
                               install("columns")))-1;
        int fixedStartColML = 
            asInteger(GET_SLOT(VECTOR_ELT(random, nlevel),
                               install("columns")))-1;
        int fixedStartCol = REML?fixedStartColREML:fixedStartColML;
        int df =  asInteger(GET_DIM(GET_SLOT(reStruct,
                                             install("original"))));
        double sqrtDF, sigmainv;
        int lev;

        REAL(logLik)[0] = nlme_logLikelihood_internal(reStruct, 1, 0);

        if (REML)
            df -= LENGTH(GET_SLOT(VECTOR_ELT(random, nlevel),
                                  install("columns")));
        sqrtDF = sqrt((double) df);
        if (asLogical(GET_SLOT(reStruct, install("dirtyBbetas")))) {
            for (lev = nlevel; lev >= 0; lev--) {
                nlme_estimate_level(stored, VECTOR_ELT(random, lev), bbetas);
            }
            LOGICAL(GET_SLOT(reStruct, install("dirtyBbetas")))[0] = 0;
        }
        for (lev = nlevel; lev >= 0; lev--) {
            nlme_invert_level(stored, VECTOR_ELT(random, lev), ncol_levels);
        }
        sigmainv = REAL(stored)[dim[0]*ncol_levels-1];
        sigmainv = sqrtDF/((sigmainv < 0.0)?-sigmainv:sigmainv);
        for (lev = 0; lev < nlevel; lev++) {
            nlme_factor_level(stored, VECTOR_ELT(random, lev), bbetas,
                              sigmainv, fixedStartCol);
            if (analyticHessian) {
                nlme_hessian_level(stored,
                                   VECTOR_ELT(random, lev), bbetas,
                                   sigmainv, fixedStartColREML,
                                   fixedStartColML, df, REML);
            }
        }

        LOGICAL(GET_SLOT(reStruct, install("dirtyStored")))[0] = 0;
    }
    UNPROTECT(1);
    return reStruct;
}

/** 
 * Create an R object of storage mode integer with value start:end
 * 
 * @param start 
 * @param end 
 * 
 * @return start:end
 */
static SEXP
nlme_seq(int start, int end)
{
    int len1 = end-start;
    SEXP ans = allocVector(INTSXP, len1+1);
    int* ansp = INTEGER(ans);
    int* p = ansp+len1;

    if (start <= end) {
        while (p >= ansp) {
            *p-- = end--;
        }
    } else {
        while (p >= ansp) {
            *p-- = end++;
        }
    }
    return ans;
}

/** 
 * Fix the storedRows and decomposedRows components in each lmeLevel
 * object in the random slot of an reStruct object.
 * 
 * Arguments must be protected.
 *
 * @param reStruct an reStruct object - gets modified in the code
 * 
 * @return the reStruct object
 */
SEXP
nlme_reStructDims(SEXP reStruct)
{
    SEXP decomposedRows_sym = install("decomposedRows");
    SEXP originalRows_sym = install("originalRows");
    SEXP storedRows_sym = install("storedRows");
    SEXP columns_sym = install("columns");
    SEXP random;
    int nlevels = LENGTH(GET_SLOT(reStruct, install("random")));
    int* decStart;
    int* nextLevelEnds;
    int* nextLevelIndex;
    int* columnLengths;
    SEXP* origIndices;
    SEXP* decIndices;
    SEXP* storedIndices;
    int i, n, currow, q1;

    if (nlevels <= 2) {
        error("reStruct object not initialized correctly");
    }

    if (NAMED(reStruct) && !asLogical(GET_SLOT(reStruct, install("dontCopy")))) {
        reStruct = duplicate(reStruct);
    }
    PROTECT(reStruct);
    random = GET_SLOT(reStruct, install("random"));
    decStart = (int*) Calloc(nlevels*(4*sizeof(int)+3*sizeof(SEXP)),
                             char);
    nextLevelEnds = decStart + nlevels;
    nextLevelIndex = nextLevelEnds + nlevels;
    columnLengths = nextLevelIndex + nlevels;
    origIndices = (SEXP*) (columnLengths + nlevels);
    decIndices = origIndices + nlevels;
    storedIndices = decIndices + nlevels;

    /*
     * Go through random[i]@originalRows and figure out starting
     * indices for chunks.
     */

    for (i = 0; i < nlevels; i++) {
        SEXP lmeLevel = VECTOR_ELT(random, i);
        SEXP rows = GET_SLOT(lmeLevel, originalRows_sym);

        origIndices[i] = rows;
        nextLevelEnds[i] = LENGTH(VECTOR_ELT(rows, 0));
        columnLengths[i] = LENGTH(GET_SLOT(lmeLevel, columns_sym));

        storedIndices[i] = GET_SLOT(lmeLevel, storedRows_sym);

        decIndices[i] = GET_SLOT(lmeLevel, decomposedRows_sym);
        nextLevelIndex[i] = 0;
        decStart[i] = 1;
    }

    nextLevelEnds[0] = 0;
    n = LENGTH(origIndices[0]);
    currow = 1;
    q1 = columnLengths[0];
    for (i = 0; i < n; i++) {
        int j;
        int end = currow + q1 - 1;

        SET_VECTOR_ELT(storedIndices[0], i,
                       nlme_seq(currow, end));
        SET_VECTOR_ELT(decIndices[0], i,
                       nlme_seq(currow, end));
        currow = end+1;
        nextLevelEnds[0] += LENGTH(VECTOR_ELT(origIndices[0], i));
        for (j = 1; j < nlevels && nextLevelEnds[j] <= nextLevelEnds[0]; j++) {
            int levIndx = nextLevelIndex[j];
            int q = columnLengths[j];

            end += q;
            SET_VECTOR_ELT(storedIndices[j], levIndx,
                           nlme_seq(currow, end));
            SET_VECTOR_ELT(decIndices[j], levIndx,
                           nlme_seq(decStart[j], end));
            currow += q;
            nextLevelIndex[j] = ++levIndx;
            if (levIndx < LENGTH(origIndices[j]))
                nextLevelEnds[j] += LENGTH(VECTOR_ELT(origIndices[j],
                                                      levIndx));
        }
        while (j > 1) {
            decStart[--j] = currow;
        }
    }

    Free(decStart);
    UNPROTECT(1);
    return reStruct;
}

/** 
 * Calculate the fitted values for an reStruct object.
 *
 * Arguments must be protected.
 * 
 * @param reStruct The reStruct object
 * @param ans numeric vector of correct length for storing the fitted values
 * 
 * @return A numeric vector of fitted values
 */
SEXP
nlme_reStruct_fitted_internal(const SEXPREC* reStruct, SEXP ans,
                              const SEXPREC* level)
{
/*     int nUsedLevel = asInteger((SEXP)level); */
    const SEXPREC* offset = GET_SLOT((SEXP)reStruct, install("offset"));
    const SEXPREC* random = GET_SLOT((SEXP)reStruct, install("random"));
    const SEXPREC* original = GET_SLOT((SEXP)reStruct, install("original"));
    const double* bbetas = REAL(GET_SLOT((SEXP)reStruct, install("bbetas")))-1;
    SEXP originalRows_sym = install("originalRows");
    SEXP storedRows_sym = install("storedRows");
    SEXP columns_sym = install("columns");
    int nlevel = LENGTH((SEXP)random) - 2;
    int* dim = INTEGER(GET_DIM((SEXP)original));
    int nrow = dim[0];
/*     int ncol = dim[1]; */
    const double one_d = 1.0;
    const int one = 1;
    double* ans_d;
    int i;

/*     if (nrow == 0 || ncol == 0) { */
/*         original = GET_SLOT((SEXP)reStruct, install("original")); */
/*         dim = INTEGER(GET_DIM((SEXP)original)); */
/*         nrow = dim[0]; */
/*     } */
    ans_d = REAL(ans);
    if (LENGTH((SEXP)offset) == 0)
        memset(ans_d, 0, nrow*sizeof(double));
    else if (LENGTH((SEXP)offset) == 1) {
        double off = asReal((SEXP)offset);
        for (i = 0; i < nrow; i++)
            ans_d[i] = off;
    } else {
        double* off = REAL((SEXP)offset);
        for (i = 0; i < nrow; i++)
            ans_d[i] = off[i];
    }
    ans_d--;

    if (level == R_NilValue) {
        level = (const SEXPREC*) nlme_seq(0, nlevel);
    }
    PROTECT((SEXP)level);

    for (i = LENGTH((SEXP)level)-1; i >= 0; i--) {
        const SEXPREC* lmeLevel = VECTOR_ELT((SEXP)random,
                                             nlevel-INTEGER((SEXP)level)[i]);
        const SEXPREC* columns = GET_SLOT((SEXP)lmeLevel, columns_sym);
        const SEXPREC* originalRows = GET_SLOT((SEXP)lmeLevel,
                                               originalRows_sym);
        const SEXPREC* storedRows = GET_SLOT((SEXP)lmeLevel,
                                             storedRows_sym);
        int M = LENGTH((SEXP) originalRows);
        int q = LENGTH((SEXP)columns);
        double* X = REAL((SEXP)original)+(INTEGER((SEXP)columns)[0]-1)*nrow-1;
        int j;

        for (j = 0; j < M; j++) {
            const SEXPREC* origRowsIndx = VECTOR_ELT((SEXP)originalRows, j);
            int bstart = INTEGER(VECTOR_ELT((SEXP)storedRows, j))[0];
            int m = LENGTH((SEXP)origRowsIndx);
            int startRow = INTEGER((SEXP)origRowsIndx)[0];
            F77_CALL(dgemv)("N", &m, &q, &one_d, X+startRow, &nrow,
                            bbetas+bstart, &one, &one_d, ans_d+startRow,
                            &one);
        }
    }

    UNPROTECT(1);
    return ans;
}

/** 
 * Calculate the fitted values for an reStruct object
 * 
 * @param reStruct The reStruct object
 * 
 * @return A numeric vector of fitted values
 */
SEXP
nlme_reStruct_fitted(const SEXPREC* reStruct, const SEXPREC* level)
{
    int n = INTEGER(GET_DIM(GET_SLOT((SEXP)reStruct,
                                     install("original"))))[0];
    SEXP ans = PROTECT(allocVector(REALSXP, n));
    nlme_reStruct_fitted_internal(reStruct, ans, level);
    UNPROTECT(1);
    return ans;
}

/** 
 * Determine the penalized least squares estimates only.  This
 * function is used in the iterative steps of the second-order
 * Laplacian approximation to the log-likelihood of a GLMM.
 * 
 * Both arguments must be protected.
 *
 * @param reStruct Pointer to an reStruct object that will be modified
 * @param pars Parameter value (of type numeric if non-NULL) to
 *             set. Not used if of length zero (or is NULL).
 * 
 * @return the argument reStruct after updating the decomposed slot
 *         and the bbetas slot
 */
SEXP
nlme_solveOnly(SEXP reStruct)
{
    SEXP logLik;
    if (!(asLogical(GET_SLOT(reStruct, install("dirtyDecomposed"))) ||
          ISNA(asReal(GET_SLOT(reStruct, install("logLik")))) ||
          asLogical(GET_SLOT(reStruct, install("dirtyBbetas")))))
        return reStruct;
    if (NAMED(reStruct) && !asLogical(GET_SLOT(reStruct, install("dontCopy"))))
        reStruct = duplicate(reStruct);
    nlme_predecompose(PROTECT(reStruct));
    logLik = GET_SLOT(reStruct, install("logLik"));
    if (ISNA(REAL(logLik)[0])) {
        REAL(GET_SLOT(reStruct, install("logLik")))[0] =
            nlme_logLikelihood_internal(reStruct, 1, 0);
    }
    {
        SEXP random = GET_SLOT(reStruct, install("random"));
        SEXP stored = GET_SLOT(reStruct, install("stored"));
        SEXP bbetas = GET_SLOT(reStruct, install("bbetas"));
        int nlevel = LENGTH(random) - 2;
        int lev;

        nlme_logLikelihood_internal(reStruct, 1, 0);
        for (lev = nlevel; lev >= 0; lev--) {
            nlme_estimate_level(stored, VECTOR_ELT(random, lev), bbetas);
        }
    }
    LOGICAL(GET_SLOT(reStruct, install("dirtyBbetas")))[0] = 0;
    UNPROTECT(1);

    return reStruct;
}

SEXP
nlme_reStructEMsteps(SEXP reStruct, SEXP niter, SEXP isVerbose)
{
    int n = asInteger(niter);
    int verbose = asLogical(isVerbose);
    SEXP EMupdateCall = PROTECT(lang4(install("EMupdate<-"), R_NilValue,
                                      R_MissingArg, R_NilValue));
    SEXP coefcall = PROTECT(lang2(install("coef"), R_NilValue));
    int nlevel;
    int i;
    for (i = 0; i < n; i++) {
        int j;
        SEXP random;
        reStruct = PROTECT(nlme_commonDecompose(reStruct, R_NilValue));
        random = GET_SLOT(reStruct, install("random"));
        nlevel = LENGTH(random)-2;
        for (j = 0; j < nlevel; j++) {
            SEXP lmeLevel = VECTOR_ELT(random, j);
            SETCADR(EMupdateCall, lmeLevel);
            SETCADDDR(EMupdateCall, GET_SLOT(lmeLevel,
                                             install("updateFactor")));
            SET_VECTOR_ELT(random, j,
                           eval(EMupdateCall, R_GlobalEnv));
        }
        REAL(GET_SLOT(reStruct, install("logLik")))[0] = NA_REAL;
        if (verbose) {
            SEXP pars;
            Rprintf("\n**EM iteration %d %f\n", i,
                    asReal(nlme_logLikelihood(reStruct, R_NilValue)));
            SETCADR(coefcall, reStruct);
            pars = eval(coefcall, R_GlobalEnv);
            PrintValue(pars);
        }
        UNPROTECT(1);
    }
    UNPROTECT(2);
    return reStruct;
}

/** 
 * Check the invariance of a column of the fixed-effects design matrix
 * within the groups defined by the origRows list.
 * 
 * @param x pointer to a column of the fixed-effects design matrix
 * @param origRows a list of (1-based) row numbers that constitute groups
 * 
 * @return 0 if x varies within the groups, otherwise 1
 */

static int
nlme_check_invariance(double *x, SEXP origRows)
{
    int i, j;
    
    for(i = 0; i < LENGTH(origRows); i++) {
	SEXP rows = VECTOR_ELT(origRows, i);
	int nr = LENGTH(rows), *rr = INTEGER(rows);
	double frstval = x[rr[0] - 1];
	for (j = 1; j < nr; j++) {
	    if (x[rr[j] - 1] != frstval) return 0;
	}
    }
    return 1;
}
    
/** 
 * Calculate the denominator degrees of freedom for individual columns
 * of the X matrix and for each term in the fixed-effects specification.
 * 
 * @param reStruct Pointer to an reStruct object
 * 
 * @return A list with two components called X and terms.
 */

SEXP nlme_getFixDF(const SEXPREC* reStruct)
{
    SEXP original = GET_SLOT((SEXP) reStruct, install("original"));
    SEXP random = GET_SLOT((SEXP) reStruct, install("random"));
    SEXP assign = GET_SLOT((SEXP) reStruct, install("assign.X"));
    SEXP nlev_sym = install("nlev");
    SEXP origRows_sym = install("originalRows");
    int Q = LENGTH(random)-2;
    SEXP columns = GET_SLOT(VECTOR_ELT(random, Q), install("columns"));
    int pp = LENGTH(columns);
    SEXP valX = PROTECT(allocVector(INTSXP, pp));
/*     SEXP valXnames = PROTECT(allocVector(STRSXP, pp)); */
    SEXP valTerms;
    SEXP val = PROTECT(allocVector(VECSXP, 2));
    SEXP valnames = PROTECT(allocVector(STRSXP, 2));
    int nn = INTEGER(GET_DIM(original))[0];
    double *orig = REAL(original);
    int i, j, nterms = 0,
	*ngrps,			/* number of groups at each level */
	*dflev,			/* degrees of freedom at each level */
	*level;			/* level assigned to each coefficient */

/* We will number the levels as in the multilevel literature, but with
 * 0-based indices, not 1-based.  The observation level is level 0,
 * the level of groups with the smallest groups is level 1,  etc. */
    ngrps = Calloc(Q+1, int);
    dflev = Calloc(Q+1, int);
    level = Calloc(pp, int);
    ngrps[0] = nn;
    for(i = 1; i <= Q; i++) {
	ngrps[i] = INTEGER(GET_SLOT(VECTOR_ELT(random, i - 1), nlev_sym))[0];
	dflev[i - 1] = ngrps[i - 1] - ngrps[i];
    }
    dflev[Q] = ngrps[Q];
				/* Assign a level to each column of X */
    for(j = 0; j < pp; j++) {
	double *col = orig + (INTEGER(columns)[j] - 1) * nn;
	int asgn = INTEGER(assign)[j];

	if (nterms < asgn) nterms = asgn; /* keep track of max(asgn) */
	level[j] = 0;
	if (asgn == 0) {	/* intercept is special case */
	    dflev[Q]--;
	    continue;
	}
	for (i = 0; i < Q; i++) {
	    if (nlme_check_invariance(col,
				      GET_SLOT(VECTOR_ELT(random, i),
					       origRows_sym))) {
		level[j] = i + 1;
	    } else {
		continue;
	    }
	}
	dflev[level[j]]--;
    }	
    valTerms = PROTECT(allocVector(INTSXP, nterms + 1));
    for(j = 0; j <= nterms; j++) INTEGER(valTerms)[j] = 0;

    for(j = 0; j < pp; j++) {
	int asgn = INTEGER(assign)[j];
	int thisCoefDF = dflev[level[j]];
	int thisTermDF = INTEGER(valTerms)[asgn];

 	INTEGER(valX)[j] = thisCoefDF;
	if (thisTermDF && thisTermDF != thisCoefDF) {
	    error("Inconsistent calculation of DF within terms");
	}
	INTEGER(valTerms)[asgn] = thisCoefDF;
    }
    SET_VECTOR_ELT(val, 0, valX);
    SET_VECTOR_ELT(val, 1, valTerms);
    SET_STRING_ELT(valnames, 0, mkChar("X"));
    SET_STRING_ELT(valnames, 1, mkChar("terms"));
    UNPROTECT(4);
    Free(level); Free(dflev); Free(ngrps);
    return namesgets(val, valnames);
}

