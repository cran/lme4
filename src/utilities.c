#include "utilities.h"

/** 
 * Check the error code returned by an Lapack routine and create an
 * appropriate error message.
 * 
 * @param info Error code as returned from the Lapack routine
 * @param laName Character string containing the name of the Lapack routine
 */
void
nlme_check_Lapack_error(int info, const char *laName)
{
    if (info != 0) {
        if (info > 0)
            error("error code %d from Lapack routine %s", info, laName);
        error("argument no. %d to Lapack routine %s is illegal",
              -info, laName);
    }
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
double *
nlme_symmetrize(double *a, const int nc)
{
    int i, j;

    for (i = 1; i < nc; i++) 
	for (j = 0; j < i; j++)
	    a[i + j*nc] = a[j + i*nc];
    return a;
}

/** 
 * Calculate the inner product of vec(nlev*D^{-1} - A'A)/2 and the
 * pdgradient array regarded as a nc*nc by plen matrix.  This
 * calculation is used in several of the LMEgradient methods.
 * 
 * @param factor The nc by nc factor of the pdMat object
 * @param A The nc by nc matrix A from the LME decomposition.
 * @param nlev The number of groups associated with the random effect
 * @param nc The number of columns in the matrix
 * @param pdgradient A pdgradient object of dimension nc by nc by plen
 * @param value An array of length plen in which the gradient will be  returned
 * 
 * @return value, with the LME gradient
 */
double *
LMEgradient(const double* factor, const double* A, const int nlev,
	     const int nc, const double* pdgradient, const int plen,
	     double* value)
{
    int info, ncsq = nc*nc, one_i = 1;
    double nlev2_d = ((double) nlev)/2., mhalf_d = -0.5, one_d = 1.0,
	zero_d = 0.0;
    double *fact = Calloc(nc * nc, double);
    
    F77_CALL(dlacpy)("U", &nc, &nc, factor, &nc, fact, &nc);
    F77_CALL(dpotri)("U", &nc, fact, &nc, &info);
    nlme_check_Lapack_error(info, "dpotri");
    F77_CALL(dsyrk)("U", "T", &nc, &nc, &mhalf_d, A, &nc, &nlev2_d,
		    fact, &nc);
    nlme_symmetrize(fact, nc);
    F77_CALL(dgemv)("T", &ncsq, &plen, &one_d, pdgradient, &ncsq,
		    fact, &one_i, &zero_d, value, &one_i);
    Free(fact);
    return value;
}

/** 
 * Replace the value of a slot or subslot of an object in place.  This
 * routine purposely does not copy the value of obj.  Use with caution.
 * 
 * @param obj object with slot to be replaced
 * @param names vector of names.  The last element is the name of the slot to replace.  The leading elements are the names of slots and subslots of obj.
 * @param value the replacement value for the slot
 * 
 * @return obj, with the named slot modified in place.
 */
SEXP
nlme_replaceSlot(SEXP obj, const SEXPREC* names, const SEXPREC* value)
{
    int lnm1 = length((SEXP) names) - 1;

    if (lnm1 >= 0) {
	SEXP comp = obj;
	int i;

	for (i = 0; i < lnm1; i++) {
	    comp = GET_SLOT(comp, install(CHAR(STRING_ELT((SEXP) names, i))));
	}
	SET_SLOT(comp, install(CHAR(STRING_ELT((SEXP) names, lnm1))),
		 (SEXP) value);
    }
    return obj;
}

