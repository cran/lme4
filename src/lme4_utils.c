#include "lme4_utils.h"

SEXP alloc_dgeMatrix(int m, int n, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, m * n);
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dpoMatrix(int n, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, n * n);
    SET_SLOT(ans, lme4_uploSym, mkString(uplo));
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dtrMatrix(int n, char *uplo, char *diag, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dtrMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, n * n);
    SET_SLOT(ans, lme4_uploSym, mkString(uplo));
    SET_SLOT(ans, lme4_diagSym, mkString(diag));
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dsCMatrix(int n, int nz, char *uplo, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dsCMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz);
    ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz);
    ALLOC_SLOT(ans, lme4_pSym, INTSXP, n + 1);
    SET_SLOT(ans, lme4_uploSym, mkString(uplo));
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}

SEXP alloc_dgCMatrix(int m, int n, int nz, SEXP rownms, SEXP colnms)
{
    SEXP ans = PROTECT(NEW_OBJECT(MAKE_CLASS("dgCMatrix"))), dn;
    int *dims = INTEGER(ALLOC_SLOT(ans, lme4_DimSym, INTSXP, 2));

    dims[0] = m; dims[1] = n;
    ALLOC_SLOT(ans, lme4_xSym, REALSXP, nz);
    ALLOC_SLOT(ans, lme4_iSym, INTSXP, nz);
    ALLOC_SLOT(ans, lme4_pSym, INTSXP, n + 1);
    dn = ALLOC_SLOT(ans, lme4_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dn, 0, duplicate(rownms));
    SET_VECTOR_ELT(dn, 1, duplicate(colnms));
    UNPROTECT(1);
    return ans;
}
