#ifndef LME4_UTILS_H
#define LME4_UTILS_H
#include <Rdefines.h>
#include "Matrix.h"
#include "Syms.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#endif

extern cholmod_common c;

/* zero an array */
#define AZERO(x, n) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = 0;}

/**
 * Allocate an SEXP of given type and length, assign it as slot nm in
 * the object, and return the SEXP.  The validity of this function
 * depends on SET_SLOT not duplicating val when NAMED(val) == 0.  If
 * this behavior changes then ALLOC_SLOT must use SET_SLOT followed by
 * GET_SLOT to ensure that the value returned is indeed the SEXP in
 * the slot.
 *
 * @param obj object in which to assign the slot
 * @param nm name of the slot, as an R name object
 * @param type type of SEXP to allocate
 * @param length length of SEXP to allocate
 *
 * @return SEXP of given type and length assigned as slot nm in obj
 */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
    SEXP val = allocVector(type, length);

    SET_SLOT(obj, nm, val);
    return val;
}


/**
 * Check for a complete match on matrix dimensions
 *
 * @param xd dimensions of first matrix
 * @param yd dimensions of second matrix
 *
 * @return 1 if dimensions match, otherwise 0
 */
static R_INLINE
int match_mat_dims(const int xd[], const int yd[])
{
    return xd[0] == yd[0] && xd[1] == yd[1];
}

SEXP alloc_dgeMatrix(int m, int n, SEXP rownms, SEXP colnms);
SEXP alloc_dpoMatrix(int n, char *uplo, SEXP rownms, SEXP colnms);
SEXP alloc_dtrMatrix(int n, char *uplo, char *diag, SEXP rownms, SEXP colnms);
SEXP alloc_dsCMatrix(int n, int nz, char *uplo, SEXP rownms, SEXP colnms);
SEXP alloc_dgCMatrix(int m, int n, int nz, SEXP rownms, SEXP colnms);

#endif
