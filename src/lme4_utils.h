#ifndef LME4_UTILS_H
#define LME4_UTILS_H
#include <R_ext/Constants.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rversion.h>
#include "Matrix.h"
#include "Syms.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#endif

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attr_hidden __attribute__ ((visibility ("hidden")))
#else
# define attr_hidden
#endif

extern cholmod_common c;

#define flag_not_factored(x) *INTEGER(GET_SLOT(x, lme4_statusSym)) = 0

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

/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static R_INLINE double*
internal_symmetrize(double *a, int nc)
{
    for (int i = 1; i < nc; i++)
	for (int j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/**
 * Create a named vector of type TYP
 *
 * @param TYP a vector SEXP type (e.g. REALSXP)
 * @param names names of list elements with null string appended
 *
 * @return pointer to a named vector of type TYP
 */
static R_INLINE SEXP
internal_make_named(int TYP, char **names)
{
    SEXP ans, nms;
    int n;

    for (n = 0; strlen(names[n]) > 0; n++) {}
    ans = PROTECT(allocVector(TYP, n));
    nms = PROTECT(allocVector(STRSXP, n));
    for (int i = 0; i < n; i++) SET_STRING_ELT(nms, i, mkChar(names[i]));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/**
 * Return the element of a given name from a named list
 *
 * @param list
 * @param nm name of desired element
 *
 * @return element of list with name nm
 */
static R_INLINE SEXP
internal_getElement(SEXP list, char *nm) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int ll = LENGTH(list);

    for (int i = 0; i < ll; i++)
	if (!strcmp(CHAR(STRING_ELT(names, i)), nm))
	    return(VECTOR_ELT(list, i));
    return R_NilValue;
}

/**
 * Calculate the length of the parameter vector (historically called "coef"
 * even though these are not coefficients).
 *
 * @param nf number of factors
 * @param nc number of columns in the model matrices for each factor
 *
 * @return total length of the coefficient vector
 */
static R_INLINE
int coef_length(int nf, const int nc[])
{
    int i, ans = 0;
    for (i = 0; i < nf; i++) ans += (nc[i] * (nc[i] + 1))/2;
    return ans;
}

SEXP attr_hidden alloc_dgeMatrix(int m, int n, SEXP rownms, SEXP colnms);
SEXP attr_hidden alloc_dpoMatrix(int n, char *uplo, SEXP rownms, SEXP colnms);
SEXP attr_hidden alloc_dtrMatrix(int n, char *uplo, char *diag, SEXP rownms, SEXP colnms);
SEXP attr_hidden alloc_dsCMatrix(int n, int nz, char *uplo, SEXP rownms, SEXP colnms);
SEXP attr_hidden alloc_dgCMatrix(int m, int n, int nz, SEXP rownms, SEXP colnms);
SEXP attr_hidden alloc3Darray(SEXPTYPE mode, int nrow, int ncol, int nface);

/* declared here but defined in lmer.c */
SEXP mer_factor(SEXP x);
SEXP mer_secondary(SEXP x);
SEXP mer_gradComp(SEXP x);
/* declared here but defined in Wishart.c */
double attr_hidden *std_rWishart_factor(double df, int p, double ans[]);

double attr_hidden
internal_betab_update(int p, int q, double sigma, cholmod_factor *L,
		      double RZX[], double RXX[], double betahat[],
		      double bhat[], double betanew[], double bnew[]);
double attr_hidden *internal_mer_fitted(SEXP x, const double initial[], double val[]);
double attr_hidden *internal_mer_ranef(SEXP x);

int attr_hidden internal_mer_Xfactor(SEXP x);

void attr_hidden internal_ECMEsteps(SEXP x, int nEM, int verb);
void attr_hidden
internal_Omega_update(SEXP Omega, const double b[], double sigma, int nf,
		      const int nc[], const int Gp[], double *vals, int trans);
void attr_hidden internal_mer_refactor(SEXP x);
void attr_hidden internal_mer_coefGets(SEXP x, const double cc[], int ptyp);
void attr_hidden internal_mer_Zfactor(SEXP x, cholmod_factor *L);

#endif
